# rustar-aligner v0.1.0: two-pass sjdb root-cause and determinism check

Two short investigations on top of the PR #1855 audit. Section A finds
the root cause of the `Annotated (sjdb) = 0` / ~50% missing splices
issue documented in
[`rustar_bam_comparison.md`](rustar_bam_comparison.md) (issue 3) and
[`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md)
(section 5). Section B answers whether rustar is deterministic
across reruns with the same `--runRNGseed 0`.

## Section A: two-pass / sjdb root-cause

### Verdict: **BUG** (medium severity, single-line/coordinate-space fix at one call site, plus a second consistency fix at a parallel site)

### What the data already told us

From the existing docs (not re-derived here):

- `grep 'Annotated (sjdb)' WT_REP2.Log.final.out` -> `| 0` on every
  rustar run with `--twopassMode Basic --sjdbGTFfile genome.gtf`.
  STAR on the same inputs reports 235-1 152 annotated splices.
- Total spliced reads per sample drop by ~50% (e.g. WT_REP2: 371
  rustar vs 762 STAR).
- 234 of 340 same-position CIGAR diffs on WT_REP2 are reads where STAR
  finds an annotated GT-AG splice and rustar emits straight `M`/`S`.

### What I found in the source

The rustar source has two independent call sites that look up
annotation status in the GTF-derived `SpliceJunctionDb`:

1. **Stitch-time** (the alignment-scoring path): determines whether to
   apply the `sjdb_score` bonus to a candidate splice and whether to
   use `align_sjdb_overhang_min` (default 3) vs the stricter
   `align_sj_overhang_min` (default 5) for the overhang gate. This
   directly affects which splices survive pass 1.

2. **Stats-time** (per-aligned-transcript bookkeeping): determines the
   `annotated=0/1` column in `SJ.out.tab` and increments the
   `splices_annotated` counter that ends up in
   `Log.final.out -> Number of splices: Annotated (sjdb)`.

The two call sites pass **different coordinate spaces** to
`SpliceJunctionDb::is_annotated`, and **neither matches what the DB
stores**.

#### How the DB is keyed

`src/junction/gtf.rs` extracts junctions in **chromosome-local
1-based** STAR-convention coordinates
([`gtf.rs:154-229`](https://github.com/scverse/rustar-aligner/blob/main/src/junction/gtf.rs#L154-L229),
raw URL
[here](https://raw.githubusercontent.com/scverse/rustar-aligner/main/src/junction/gtf.rs)).
Specifically:

```rust
// gtf.rs:207-209
// Intron coordinates (1-based, STAR convention)
let intron_start = exon1.end + 1;
let intron_end = exon2.start - 1;
```

`exon1.end` is the GTF column-5 value (chromosome-local, 1-based per
GTF spec). On both index build
([`index/mod.rs:84`](https://github.com/scverse/rustar-aligner/blob/main/src/index/mod.rs#L84))
and index load
([`index/io.rs:45-56`](https://github.com/scverse/rustar-aligner/blob/main/src/index/io.rs#L45-L56))
the DB is constructed via `SpliceJunctionDb::from_raw_junctions` /
`from_gtf_configured` which stores those tuples **verbatim** as
`JunctionKey { chr_idx, intron_start, intron_end, strand }`. No
conversion. The DB is therefore keyed by **(chr_idx,
intron_start_local_1based, intron_end_local_1based, strand)**.

#### How stitch-time looks the DB up

[`src/align/stitch.rs:1305-1314`](https://github.com/scverse/rustar-aligner/blob/main/src/align/stitch.rs#L1305-L1314)
(raw
[here](https://raw.githubusercontent.com/scverse/rustar-aligner/main/src/align/stitch.rs)):

```rust
// Check annotation (needed for sjdbScore bonus and finalization check)
let is_annotated = junction_db.is_some_and(|db| {
    let junc_donor_sa = (donor_sa as i64 + jr_shift as i64) as u64;
    let donor_fwd =
        index.sa_pos_to_forward(junc_donor_sa, cluster.is_reverse, del as usize);
    let acceptor_fwd = donor_fwd + del as u64;
    db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 0)
        || db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 1)
        || db.is_annotated(cluster.chr_idx, donor_fwd, acceptor_fwd, 2)
});
```

`donor_sa` is initialised earlier in the same function:

```rust
// stitch.rs:1199-1200
// donor_sa = exclusive end of exon A = STAR's gAend+1. jr_shift = STAR's jR.
let donor_sa = last_exon.genome_end;
```

That is **genome-wide absolute 0-based** (the comment "STAR's gAend+1"

- the fact that `sa_pos_to_forward` returns a position in the padded
  genome buffer, not chr-local space; confirmed by
  [`index/mod.rs:166-175`](https://github.com/scverse/rustar-aligner/blob/main/src/index/mod.rs#L166-L175)
  which never references `chr_start`). For yeast chrI in the test
  profile, `chr_start[0] = 0`, so the genome-absolute 0-based value
  equals the chr-local 0-based value, which is **1 less** than the
  chr-local 1-based value the DB stores. For any other chromosome the
  offset compounds by `chr_start[chr_idx]`.

Result: **every stitch-time lookup misses by exactly 1 on chr 0 (the
only chromosome with annotated junctions in the test data), and by
more on higher chromosomes.** The `sjdb_score` bonus is never
applied, and the stricter `align_sj_overhang_min` gate is used for
every splice. Pass 1 therefore drops every GT-AG splice whose
overhang falls in the [3, 5) window that STAR would accept as
annotated.

#### How stats-time looks the DB up

[`src/lib.rs:1860-1894`](https://github.com/scverse/rustar-aligner/blob/main/src/lib.rs#L1860-L1894)
(raw
[here](https://raw.githubusercontent.com/scverse/rustar-aligner/main/src/lib.rs)):

```rust
// lib.rs:1858-1894 (excerpt)
let intron_start = genome_pos + 1; // 1-based, first intronic base
let intron_end = genome_pos + intron_len as u64; // 1-based, last intronic base
// ...
let strand = match motif.implied_strand() {
    Some('+') => 1u8,
    Some('-') => 2u8,
    _ => 0u8,
};
let annotated = index.junction_db.is_annotated(
    transcript.chr_idx,
    intron_start,
    intron_end,
    strand,
);
```

`genome_pos` is `transcript.genome_start + <CIGAR-accumulated>` and is
explicitly "WITHOUT n_genome offset" per the docstring at
[`stitch.rs:21`](https://github.com/scverse/rustar-aligner/blob/main/src/align/stitch.rs#L21).
So `intron_start = genome_pos + 1` is **genome-wide absolute 1-based**.
On chr 0 (chr_start=0) that equals the DB's chr-local 1-based key and
the lookup succeeds. On chr 1+ it still misses by `chr_start[chr_idx]`.

This is why the same WT_REP2 BAM produces `SJ.out.tab` with two rows
flagged `annotated=1` (chr 0 splices that the stats-time lookup found)
**while** `Log.final.out` reports `Number of splices: Annotated (sjdb)
= 0` — `Log.final.out`'s counter comes from
`transcript.junction_annotated` which is filled at stitch time
([`stats.rs:177-182`](https://github.com/scverse/rustar-aligner/blob/main/src/stats.rs#L177-L182)),
not from the stats-time lookup. The two counters disagree per file
because of the parallel-but-not-identical query bug.

Verification from the rustar test profile (yeast chrI + GFP, both
`chr_start` values readable in `chrStart.txt` -> `0\n262144`):

```
WT_REP2 rustar SJ.out.tab annotated column distribution:
   12 rows with annotated=0
    2 rows with annotated=1
WT_REP2 rustar Log.final.out:
   Number of splices: Total              371
   Number of splices: Annotated (sjdb)     0
   Number of splices: GT/AG              276
   Number of splices: Non-canonical       92
WT_REP2 STAR Log.final.out (same inputs):
   Number of splices: Total              762
   Number of splices: Annotated (sjdb)   644
   Number of splices: GT/AG              724
   Number of splices: Non-canonical       35
```

The 2 vs 0 split in rustar's own outputs is the smoking gun — both
counters should be reading the same DB, and they aren't.

### Where to fix

One of:

1. **(preferred)** Normalise the DB keys to **genome-absolute 0-based**
   at construction time. This matches the convention already used by
   `prepared_junctions` ([`index/mod.rs:92-105`](https://github.com/scverse/rustar-aligner/blob/main/src/index/mod.rs#L92-L105))
   and `SpliceJunctionStats` (which stores genome-absolute keys; the
   chr-local conversion happens at `sj_output.rs:262-264`). Then make
   the stats-time lookup pass `intron_start - 1`/`intron_end - 1`
   instead of `+1`/`+intron_len` to align with the new key scheme.
2. **(alternative)** Subtract `chr_start[chr_idx]` and add 1 at the
   stitch-time call site:
   ```rust
   let chr_off = index.genome.chr_start[cluster.chr_idx];
   let donor_local_1b = donor_fwd - chr_off + 1;
   let acceptor_local_1b = acceptor_fwd - chr_off + 1;
   db.is_annotated(cluster.chr_idx, donor_local_1b, acceptor_local_1b, ...)
   ```
   And at the stats-time call site subtract `chr_start[chr_idx]` from
   `intron_start`/`intron_end` (the `+1` offset already matches the DB
   on chr 0). This requires changes at both call sites and keeps the
   internal contract about "what coord space the DB stores"
   undocumented.

Option 1 is cleaner because the rest of the file already standardises
on genome-absolute 0-based.

### Minimal upstream reproducer (proposed)

The current rnaseq test profile is fine but heavyweight. A 50-line
yeast snippet:

- One chromosome of length ~10 kb.
- One transcript with two exons separated by a GT-AG intron of length
  ~1 kb, written into a 4-line GTF.
- A 200 nt paired-end read pair spanning the junction with 6 nt
  overhang on the downstream side (in the
  `[align_sjdb_overhang_min=3, align_sj_overhang_min=5)` gap, so
  rustar would accept it as annotated and reject it as unannotated).

Expected: STAR maps `XNYM` spliced, rustar maps `(X+Y)M` plain. After
the fix, rustar matches STAR; `Annotated (sjdb)` in Log.final.out is
1 instead of 0; the splice appears in SJ.out.tab with `annotated=1`.

### Upstream issue body (ready to paste)

<details>
<summary>scverse/rustar-aligner issue draft</summary>

**Title:** GTF-derived splice junction DB is keyed in chr-local 1-based, but is consulted in genome-absolute coords; `Annotated (sjdb) = 0` everywhere and ~50% of GT/AG splices are dropped from pass 1

**Body:**

Running rustar v0.1.0 inside nf-core/rnaseq PR #1855 with
`--twopassMode Basic --sjdbGTFfile genome.gtf`, every Log.final.out
reports `Number of splices: Annotated (sjdb) | 0`, even though
sjdbList.fromGTF.out.tab contains 100+ junctions and SJ.out.tab
contains rows that match those junctions exactly. The total spliced
read count is also ~50% lower than STAR on the same inputs
(WT_REP2: 371 rustar vs 762 STAR).

Tracing the source, `SpliceJunctionDb::is_annotated` is consulted by
two independent call sites:

- `src/align/stitch.rs:1306-1314` (alignment-time, drives the
  `sjdb_score` bonus and the `align_sjdb_overhang_min` gate), passing
  `donor_fwd, acceptor_fwd` which are **genome-absolute 0-based** SA
  positions decoded via `sa_pos_to_forward`.
- `src/lib.rs:1877-1882` (stats-time, drives the `annotated` column in
  SJ.out.tab and the `splices_annotated` counter), passing
  `intron_start = genome_pos + 1`, `intron_end = genome_pos + intron_len`,
  i.e. **genome-absolute 1-based-equivalent**.

But the DB itself is keyed by **chromosome-local 1-based** intron
coordinates extracted from the GTF
(`src/junction/gtf.rs:208-209`, `intron_start = exon1.end + 1`,
`intron_end = exon2.start - 1`; stored verbatim via
`SpliceJunctionDb::from_raw_junctions` /
`SpliceJunctionDb::from_gtf_configured` with no coordinate change).

The two call sites therefore disagree with each other **and** disagree
with the DB. On the test-profile yeast chrI (`chr_start[0]=0`):

- stitch-time misses by exactly 1 on every lookup (0-based vs 1-based).
- stats-time matches on chr 0 (0+1 == 1-based) but misses on chr 1+.

The stitch-time miss is the load-bearing one. The same WT_REP2 BAM
exhibits the inconsistency in its own outputs:

```
SJ.out.tab:  2 rows annotated=1, 12 rows annotated=0      (stats-time lookup)
Log.final.out: Number of splices: Annotated (sjdb) | 0    (stitch-time lookup)
```

Suggested fix (one of two equivalent options): either normalise
`SpliceJunctionDb` keys to **genome-absolute 0-based** at construction
(matching the rest of the file — `prepared_junctions` in
`index/mod.rs:92-105` already use genome-absolute 0-based, and
`SpliceJunctionStats` does too) and update both call sites to query in
that space; **or** subtract `chr_start[chr_idx]` (and adjust by 1) at
both call sites. Option 1 is cleaner and contains the fix to one
file.

Reproducer: branch
[`rustar-aligner`](https://github.com/nf-core/rnaseq/tree/rustar-aligner)
in nf-core/rnaseq, `-profile test,docker -params-file rustar.params.yml`
where `rustar.params.yml` is `use_rustar_star: true`. Inspect
`work/.../*.Log.final.out` and `work/.../*.SJ.out.tab` for any
paired-end sample (WT_REP1, WT_REP2, RAP1_IAA_30M_REP1).

A targeted unit test would build a `SpliceJunctionDb` from a 2-exon
toy GTF on a 10 kb single-chromosome genome, then call `is_annotated`
with both coord conventions and assert the call from the alignment
path returns `true`. The bug is in
[`stitch.rs:1306-1314`](https://github.com/scverse/rustar-aligner/blob/main/src/align/stitch.rs#L1306-L1314)
(raw URL
[here](https://raw.githubusercontent.com/scverse/rustar-aligner/main/src/align/stitch.rs))
and the inconsistency between that and
[`lib.rs:1860-1894`](https://github.com/scverse/rustar-aligner/blob/main/src/lib.rs#L1860-L1894)
(raw URL
[here](https://raw.githubusercontent.com/scverse/rustar-aligner/main/src/lib.rs)).

</details>

## Section B: determinism check

### Verdict: **deterministic except for order** (see [the bottom of this section](#verdict-deterministic-except-for-order) for the details)

### Setup

- Run 1: the original PR-1855 rustar work dirs at
  `/home/ubuntu/rnaseq-rustar-aligner/work/{e7,d7,cc,e1,ed}/...`
  (catalogued in `rustar_bam_comparison.md`). Initial `mtime` 14:14.
- Run 2: new rustar run launched in this investigation with
  `nextflow run . -profile test,docker -params-file rustar-det.params.yml`
  (same params file content, only `outdir` changed to
  `results-rustar-det`, so all upstream work was forced to re-execute
  rather than resume). Work dirs:
  - WT_REP1: `work/db/1067f7af22e0060f5dd51e3e3567da`
  - WT_REP2: `work/16/d9d802566f49cb84bfccd08ddef610`
  - RAP1_IAA_30M_REP1: `work/a7/3cb7403e663404029828bf259b75a2`
  - RAP1_UNINDUCED_REP1: `work/df/f3cec00c49545fa21e27ddaca2abe2`
  - RAP1_UNINDUCED_REP2: `work/ae/dd8f27fcf749ee0da2bd1d5e268a5c`

Both runs invoke the same rustar image at the same `--runRNGseed 0`
with all other CLI args identical (verified via `.command.sh`).

### Plain-file comparison

`md5sum` per file, per sample:

| Sample              | Log.final.out | SJ.out.tab | SJ.pass1.out.tab |
| ------------------- | ------------- | ---------- | ---------------- |
| WT_REP1             | DIFFER\*      | IDENTICAL  | IDENTICAL        |
| WT_REP2             | DIFFER\*      | IDENTICAL  | IDENTICAL        |
| RAP1_IAA_30M_REP1   | DIFFER\*      | IDENTICAL  | IDENTICAL        |
| RAP1_UNINDUCED_REP1 | DIFFER\*      | IDENTICAL  | IDENTICAL        |
| RAP1_UNINDUCED_REP2 | DIFFER\*      | IDENTICAL  | IDENTICAL        |

\* `Log.final.out` differs in **timestamps only**:

```
< Started job on |  May 12 14:14:23
< Started mapping on |  May 12 14:14:23
< Finished on |  May 12 14:14:58
< Mapping speed, Million of reads per hour |  5.14
---
> Started job on |  May 12 15:44:09
> Started mapping on |  May 12 15:44:09
> Finished on |  May 12 15:44:45
> Mapping speed, Million of reads per hour |  4.86
```

All numerical mapping statistics (including the `Number of splices:
Annotated (sjdb) | 0` line) are byte-identical. The wall-clock-derived
"Mapping speed" varies because the second run shared the VM with five
other concurrent pipelines, not because of any aligner behaviour.

`ReadsPerGene.out.tab` was not produced (pipeline doesn't pass
`--quantMode GeneCounts`).

### BAM comparison

Raw BAM bytes differ on every file. That's expected at the gzip-block
level even without record-level differences. The interesting tests
are SAM-text equality (header stripped) at same record order and at
name-sorted order.

#### Record counts (record-count parity)

Same per file between the two runs, on every BAM:

| Sample              | Aligned.out.bam (run1/run2) | Aligned.toTranscriptome.out.bam (run1/run2) |
| ------------------- | --------------------------: | ------------------------------------------: |
| WT_REP1             |           188 322 / 188 322 |                           161 724 / 161 724 |
| WT_REP2             |             94 458 / 94 458 |                             82 744 / 82 744 |
| RAP1_IAA_30M_REP1   |             94 950 / 94 950 |                             85 152 / 85 152 |
| RAP1_UNINDUCED_REP1 |             48 952 / 48 952 |                             45 609 / 45 609 |
| RAP1_UNINDUCED_REP2 |             98 201 / 98 201 |                             92 675 / 92 675 |

No record loss or gain across reruns.

#### Same-order SAM body

`docker run ... samtools view <bam> | md5sum` — all 10 BAMs (5 samples

- 2 BAM types) produce DIFFERENT md5s between the two runs. The
  records are emitted in different order across reruns.

#### Name-sorted SAM body

`docker run ... samtools sort -n -O sam <bam> | grep -v '^@' | md5sum`:

| Sample              | Aligned.out.bam | Aligned.toTranscriptome.out.bam |
| ------------------- | --------------- | ------------------------------- |
| WT_REP1             | IDENTICAL       | IDENTICAL                       |
| WT_REP2             | IDENTICAL       | IDENTICAL                       |
| RAP1_IAA_30M_REP1   | IDENTICAL       | IDENTICAL                       |
| RAP1_UNINDUCED_REP1 | IDENTICAL       | IDENTICAL                       |
| RAP1_UNINDUCED_REP2 | IDENTICAL       | IDENTICAL                       |

After name-sorting, every record-level field (flag, MAPQ, CIGAR,
sequence, quality, every tag including `AS`, `nM`, `MD`, `NH`, `HI`,
`RG`) matches byte-for-byte across the two runs on every BAM.

### Verdict: **deterministic except for order**

rustar v0.1.0 with `--runRNGseed 0` produces the same records, same
counts, and same per-record fields and tags on a rerun. What changes
between runs:

- the **order** records are written to `Aligned.out.bam` /
  `Aligned.toTranscriptome.out.bam` (the BAM body is not in input-read
  order, and that order varies between runs).
- the **wall-clock timestamps** in `Log.final.out` (and the derived
  `Mapping speed` field). Every other line is byte-identical.

`SJ.out.tab` and `SJ.pass1.out.tab` are byte-identical out-of-the-box
(they're sorted before write — see `sj_output.rs:248-253`).

This is good news for downstream tools that care about content but not
emission order — the pipeline name-sorts via `SAMTOOLS_SORT`
immediately downstream of the aligner anyway, so the order divergence
is invisible to everything after that step.

**Verdict per `rustar_bam_comparison.md` vocabulary**: BAM
record-order divergence between same-seed reruns is **FLOATING-POINT
/ ORDER**. There is no record-level non-determinism. Not worth filing
upstream unless the rustar maintainers care about strict byte-level
reproducibility for forensic / hash-based testing (which is a
reasonable thing to want but isn't on the critical path for nf-core
consumers).

Root cause is almost certainly the parallel reader/writer threading in
the aligner: reads finish in a different order each run as the OS
schedules workers, and the records are flushed in completion order
rather than input order. STAR has the same issue without
`--outSAMmultOutputOrder Old_2.4` and equivalent.

No upstream issue body needed; this isn't a bug.
