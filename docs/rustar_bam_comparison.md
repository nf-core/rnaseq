# rustar-aligner v0.1.0 vs STAR 2.7.11b: categorical BAM comparison

Date: 2026-05-12. Pipeline: nf-core/rnaseq branch `rustar-aligner` (PR #1855). Aligner under test: rustar-aligner v0.1.0 from https://github.com/scverse/rustar-aligner. Reference: STAR 2.7.11b. Both aligners invoked with identical CLI args except the binary name (see `rustar_investigation_wt_rep2.md` for the full args; same nf-core test profile, paired-end yeast subset + GFP).

## TL;DR (net-new findings beyond the transcriptome-BAM mate-pair bug)

Six new, file-able issues surfaced on top of the previously-characterised paired-end transcriptome BAM mate-pair bug (already documented in [`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md) and not re-derived here). In rough severity order:

1. **`NM` tag is replaced with `nM`, with the wrong semantics.** rustar emits `nM` everywhere STAR emits `NM`. STAR's `NM` follows the SAM spec (edit distance = mismatches + indel bases); rustar's `nM` only counts substitutions. ~2% of identical-CIGAR records disagree (1 663 / 89 107 in WT_REP2; 850 / 46 869 in RAP1_UNINDUCED_REP1). The records with disagreement all contain indels, and the rustar value is always lower. Any downstream tool that reads `NM:i:` (samtools stats, Picard, MultiQC's BAM-stats parsers, mappability filters keyed on edit distance) will silently get nothing back. **BUG (high)**.
2. **`XS` tag (strand) is never emitted.** Pipeline ran with `--outSAMstrandField intronMotif`; STAR adds `XS:A:+` / `XS:A:-` on spliced reads (1 300 records in WT_REP2 genome BAM, 356 in RAP1_UNINDUCED_REP1). rustar emits zero `XS` tags. StringTie, Cufflinks, and rseqc's `infer_experiment.py` all need `XS` for stranded transcript assembly / strand inference. **BUG (high)**.
3. **GTF-supplied splice junctions are not seeded into pass 1, dropping ~50% of splices.** Total spliced reads drops by ~half on every paired-end sample and `Number of splices: Annotated (sjdb) = 0` on every sample (STAR: 235-1 152 annotated splices). On WT_REP2, 234 of 340 same-position CIGAR diffs are reads where STAR finds a splice and rustar doesn't. This is the same root issue flagged briefly in the prior investigation; we now show it changes the per-read CIGAR, not just the log header. **BUG (medium)** (mapping rate impact in test profile is < 0.5 pp, but downstream junction-aware tools see different evidence).
4. **`@PG` header line is content-free.** rustar emits `@PG\tID:rustar-aligner` with no `PN`, no `VN`, no `CL`. STAR emits a fully populated PG record. Provenance is lost for any downstream tool that reads PG (e.g. MultiQC's program-version table). **BUG (low)**.
5. **Per-record `RG:Z:` tag is missing from the transcriptome BAM** (despite the `@RG` header being present). STAR writes `RG:Z:WT_REP2` on every transcriptome record; rustar writes none. Any tool that splits multi-sample transcriptome BAMs by RG (Salmon's `--validateMappings` flow with bundled BAMs, custom QC) won't be able to. **BUG (low)**. Genome BAM is unaffected; per-record `RG:Z:` is present there.
6. **rustar emits ~17% more secondary alignments than STAR**, with NH values up to 20 vs STAR's max of 7 on the same data. Symmetric: 368 reads move from `NH=1` -> `NH=2` and 318 from `NH=2` -> `NH=1` in WT_REP2; the rest is rustar finding extra noisy hits at the long tail. No `outFilterMultimapScoreRange` analogue appears to be applied. **BEHAVIOURAL / possibly BUG (medium)** depending on intended scoring threshold.

Issues 1-2 are likely the most impactful for downstream nf-core/rnaseq consumers (`samtools stats` rolls up `NM` into the MultiQC report; `infer_experiment.py` needs `XS` for stranded protocol detection). Issue 6 is a behavioural choice that could be intentional but should at minimum be documented.

## Methodology

### Samples and what was characterised

| Sample              | Layout | Genome BAM categories 1-8 | Transcriptome BAM                   | Spot-check |
| ------------------- | ------ | ------------------------- | ----------------------------------- | ---------- |
| WT_REP2             | PE     | full                      | known-bug confirmed (see prior doc) | -          |
| RAP1_UNINDUCED_REP1 | SE     | full                      | not characterised                   | -          |
| WT_REP1             | PE     | flagstat + tag inventory  | known-bug confirmed                 | yes        |
| RAP1_IAA_30M_REP1   | PE     | flagstat + tag inventory  | known-bug confirmed                 | yes        |
| RAP1_UNINDUCED_REP2 | SE     | flagstat + tag inventory  | not paired                          | yes        |

### Commands and helper scripts

All BAM inspection was done with samtools 1.21 from `community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4` (the STAR Wave container the pipeline uses). All artefacts are on the `nf-dev-rnaseq` VM at `/home/ubuntu/rustar_bam_inv/`:

- `flagstats/{star,rus}_<sample>_{Aligned.out,Aligned.toTranscriptome.out}.txt` - flagstat output for every BAM.
- `sorted/{star,rus}_<sample>_genome_all.sam` - name-sorted SAM (header stripped) for the two fully-characterised samples.
- `readnames/{star,rus}_<sample>_genome.names` - sorted, unique primary-read-name sets.
- `analyze.py` / `analyze2.py` / `analyze3.py` - the three Python scripts that produced the per-category counts; copies in `/tmp/rustar_bam_inv/`.

Top-level invocation:

```bash
# On the VM, primary read name set comparison:
docker run --rm -v /home/ubuntu:/home/ubuntu community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4 \
    samtools view -F 0x100 <BAM> | cut -f1 | sort -u > <names>

# Name-sort for joint comparison:
samtools sort -n -@ 4 -O sam <BAM> | grep -v "^@" > sorted/<key>.sam

# Then the analyze*.py scripts joined on (read_name, mate_bit).
```

The work directories used:

| Sample              | STAR work dir                            | rustar work dir                          |
| ------------------- | ---------------------------------------- | ---------------------------------------- |
| WT_REP1             | `work/ee/09a9d7045e7a6cbb606372f97420eb` | `work/e7/91f0997ccfaa8ad4e623a304820459` |
| WT_REP2             | `work/e0/b8c327bfa0a964eecd894bcb05b569` | `work/d7/43755befdffb99383bedb820e900f9` |
| RAP1_IAA_30M_REP1   | `work/e2/49de29a5e06418584554f94d388406` | `work/cc/81563a4b47309a642a01cfbbcc5be1` |
| RAP1_UNINDUCED_REP1 | `work/89/bfc76336d60abc43e6f1a6135778a4` | `work/e1/28e07ef7b8fcf635042af437f24c57` |
| RAP1_UNINDUCED_REP2 | `work/d9/76ac6573acf0d54637352628612807` | `work/ed/529e7f6245e1e01fb8037310826ba9` |

(Paths under `/home/ubuntu/rnaseq-rustar-aligner/`.)

## Category-by-category results

### 1. Reads present in one BAM but not the other (primary records, by name + mate)

Both BAMs cover ~99.9% of the same fragments. Per sample:

| Sample              | STAR primary keys | rustar primary keys |   both | STAR-only | rustar-only |
| ------------------- | ----------------: | ------------------: | -----: | --------: | ----------: |
| WT_REP2             |            90 709 |              90 716 | 90 659 |        50 |          57 |
| RAP1_UNINDUCED_REP1 |            47 713 |              47 693 | 47 688 |        25 |           5 |

Most of the 50-100 "X-only" reads on WT_REP2 are reads where STAR found a spliced alignment that rustar reports as unmapped, or vice versa. The numbers are too small to chase individually for a 5-sample test profile, but a >0.1% read-set divergence rate at this depth means real production runs will see hundreds of differing reads.

Spot check on WT_REP1: STAR `180 597` primary, rustar `180 596` (1 fewer); same single-fragment-different pattern.

**Verdict: BEHAVIOURAL** for the bulk (low-mapability reads where the splice issue in category 4 is the proximate cause). The known sjdb bug (issue 3 in TL;DR) explains most of these; once that's filed, this category is downstream of it.

### 2. SAM flag-bit divergence (reads present in both)

In WT_REP2 (PE): only 2 reads out of 90 659 have any flag-bit difference, both flipping `0x10` + `0x20` (strand of self and mate). RAP1_UNINDUCED_REP1 (SE) has 13 reads flipping `0x10` alone.

| Bit                |   WT_REP2 differs | RAP1_UN1 differs |
| ------------------ | ----------------: | ---------------: |
| 0x1 paired         |                 0 |                0 |
| 0x2 proper-pair    | 0 (genome BAM!)\* |              n/a |
| 0x4 unmapped       |                 0 |                0 |
| 0x8 mate-unmap     |                 0 |              n/a |
| 0x10 reverse       |                 2 |               13 |
| 0x20 mate-rev      |                 2 |              n/a |
| 0x40 / 0x80        |                 0 |              n/a |
| 0x100 secondary    |                 0 |                0 |
| 0x200 QC fail      |                 0 |                0 |
| 0x400 duplicate    |                 0 |                0 |
| 0x800 supplemental |                 0 |                0 |

\* Genome BAM `0x2` is correctly populated (90 714 records in WT*REP2). The proper-pair bug is \_transcriptome-BAM specific*. This is the key data point for "is the prior bug also in the genome BAM?" - it isn't.

The strand-flip cases are reads with ambiguous orientation (e.g. palindromic-ish soft-clipped ends); the per-mate position usually swaps as well. Looks like RNG-driven tie-breaking on which mate-orientation to call primary.

**Verdict: RNG / TIE-BREAKING** (15 reads total across the two analysed samples).

### 3. MAPQ differences

Both aligners use the same `{0, 1, 3, 255}` STAR-style scheme keyed on multi-mapper count (NH=1 -> 255; NH=2 -> 3; NH=3 -> 1; NH>=4 -> 0). MAPQ deltas are driven entirely by NH transitions:

| Sample              | MAPQ identical | MAPQ differs |
| ------------------- | -------------: | -----------: |
| WT_REP2             |         89 805 |          854 |
| RAP1_UNINDUCED_REP1 |         47 496 |          192 |

The MAPQ-delta distribution on WT_REP2 is symmetric around zero (`{-255: 14, -254: 96, -252: 368, -3: 10, -2: 28, -1: 14, +2: 4, +252: 318, +255: 2}`), consistent with reads moving symmetrically between NH classes (368 went NH 1->2, 318 went NH 2->1 - see category 7). The +/- 252 entries are reads going NH=1 (mapq 255) <-> NH=2 (mapq 3).

**Verdict: BEHAVIOURAL / RNG** (driven by category 4 - missing splice causes some unique mappers to drop to NH>1 and vice versa - and category 7 - multi-mapper sampling).

### 4. CIGAR differences (same name + flag + position)

WT_REP2: 340 reads have identical position + orientation but different CIGAR.

| CIGAR shape category                              | WT_REP2 | RAP1_UN1 |
| ------------------------------------------------- | ------: | -------: |
| STAR has splice (N), rustar has only M/S          |     234 |       45 |
| rustar has splice (N), STAR has only M/S          |      80 |        4 |
| Same op set, lengths differ (mostly M/S boundary) |      13 |        0 |
| STAR has soft-clip (S), rustar doesn't            |       9 |        0 |
| rustar has soft-clip (S), STAR doesn't            |       2 |        1 |
| Same op set with indels, lengths differ           |       2 |        1 |

The dominant pattern (>70%) is **STAR finds a splice through a known annotated junction, rustar emits a straight `101M` (or `99M2S`)**. Example `SRR6357072.32572100`: STAR maps `96M14674N5M` across a yeast intron (one of the annotated junctions in `genome_gfp.gtf`); rustar maps `101M` at the same exonic start - it just doesn't try the splice. NH for this read drops from 2 (STAR, where the spliced and unspliced both score) to 1 (rustar, only the unspliced version found).

This is the same root cause as the `Annotated (sjdb) = 0` finding in the prior doc, but now visible at the per-read CIGAR level.

**Verdict: BUG (medium)** - same root cause as the sjdb-not-seeded issue. File a separate issue: rustar isn't using `--sjdbGTFfile` to seed pass-1 alignment, which costs ~50% of GT/AG splices.

### 5. Position / RNAME differences

Same name + mate-bit, both primary:

| Sample              | Same rname, diff pos | Different rname |
| ------------------- | -------------------: | --------------: |
| WT_REP2             |                1 212 |               0 |
| RAP1_UNINDUCED_REP1 |                  768 |               0 |

No read maps to a different reference between the two aligners. The 1 212 same-chr position differences are mostly NH>=2 reads where rustar picks a different alignment as primary (see category 7).

**Verdict: BEHAVIOURAL / RNG** - driven by primary-locus selection on multi-mappers (category 7).

### 6. Optional SAM tags

Genome BAM tags emitted by each aligner (per `--outSAMattributes NH HI AS NM MD` + the `--outSAMstrandField intronMotif` directive):

| Tag                         | STAR present | rustar present | Notes                                                                                                                                                                                                              |
| --------------------------- | :----------: | :------------: | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `NH`                        |     yes      |      yes       | values differ on 870 reads (WT_REP2); see category 7                                                                                                                                                               |
| `HI`                        |     yes      |      yes       | rare disagreement (8 reads WT_REP2) - tie-breaking on which hit is index 1                                                                                                                                         |
| `AS`                        |     yes      |      yes       | **864 same-CIGAR records disagree** (1 % of all comparable reads); rustar AS is always lower by 2-5. Worth filing - byte-level identity on AS for uniquely-mapped reads is explicitly promised in rustar's README. |
| `NM`                        |   **yes**    |     **no**     | rustar emits `nM` instead - see issue 1 in TL;DR                                                                                                                                                                   |
| `nM`                        |      no      |      yes       | rustar's `nM` counts substitutions only, not indels. 1 663 records have `NM != nM` even on identical CIGAR (WT_REP2).                                                                                              |
| `MD`                        |     yes      |      yes       | identical where CIGAR is identical (89 107 / 89 107 in WT_REP2). Differs only when CIGAR differs.                                                                                                                  |
| `XS` (strand)               |     yes      |     **no**     | 1 300 records in WT_REP2 lose strand info. **BUG (high)**                                                                                                                                                          |
| `RG:Z:` (genome BAM)        |     yes      |      yes       | both populated on every record                                                                                                                                                                                     |
| `RG:Z:` (transcriptome BAM) |     yes      |     **no**     | rustar @RG header exists but per-record tag is missing - issue 5                                                                                                                                                   |

STAR-specific `jM` (junction motif) and `jI` (intron coordinates) are not requested by the pipeline (`--outSAMattributes` excludes them) so neither aligner emits them; **NOT TESTED**.

**Verdict per tag**:

- `NM`/`nM`: **BUG (high)** - file separately.
- `XS`: **BUG (high)** - file separately.
- `AS`: **BUG (low)** - small per-read divergence (<=5 score units), needs an issue against rustar's promised byte-equivalence claim.
- `RG:Z:` transcriptome: **BUG (low)** - file separately.

### 7. Primary/secondary assignment for multi-mappers

For reads where STAR reports NH>=2 _and_ rustar reports the same NH (i.e. comparable cardinality):

| Sample              | Multi-mapper reads | Same primary locus | Different primary locus | Same full locus set |
| ------------------- | -----------------: | -----------------: | ----------------------: | ------------------: |
| WT_REP2             |              1 572 |          671 (43%) |               901 (57%) |         1 560 (99%) |
| RAP1_UNINDUCED_REP1 |              1 013 |          337 (33%) |               676 (67%) |         1 010 (99%) |

**99% of multi-mappers have the same locus set in both BAMs** - the aligners agree on which positions are tied. They disagree only on which tied alignment is reported with `HI:i:1`. This is consistent with rustar using a different RNG (ChaCha) than STAR (Mersenne Twister) for tie-breaking; the call sites are deterministic with `--runRNGseed 0` but the bit-level streams are not aligned.

**Verdict: RNG / TIE-BREAKING.** Not a bug per se; the rustar README documents this explicitly. The 1% locus-set difference (12 / 1 572 in WT_REP2) is where rustar genuinely found a different alignment - small enough to fold into category 4 (mostly splice differences).

For multi-mappers where STAR has NH=1 and rustar has NH=2+ (368 reads in WT_REP2), see issue 6 in TL;DR: rustar finds extra noisy hits. The extra hits show NH values up to 20 (vs STAR max 7) on the same data.

### 8. Sort order

Both genome BAMs have `@HD VN:1.4` (STAR) / `@HD VN:1.6` (rustar), both without `SO:` tag (= unsorted). The aligners emit records in input-read order; pair-internal ordering occasionally differs (STAR sometimes emits mate2 before mate1; rustar always emits mate1 first), but this has no downstream impact because the pipeline immediately name-sorts via SAMTOOLS_SORT.

**Verdict: BEHAVIOURAL / equivalent.** rustar's choice of always emitting mate1 first is arguably nicer than STAR's; no downstream consequence.

### 9. Unmapped-read FASTQ contents

**NOT TESTED.** STAR was not invoked with `--outReadsUnmapped Fastx`; both aligners produced an empty unmapped-FASTQ stream and the pipeline correctly handles this. To test this category one would need to add `--outReadsUnmapped Fastx` to the ALIGN_STAR module CLI (and verify rustar honours the flag at all).

Reproducer if needed:

```bash
# Add --outReadsUnmapped Fastx to the args in both modules/nf-core/star/align/main.nf
#  and modules/local/rustar/align/main.nf, then compare .Unmapped.out.mate{1,2}.
diff <(zcat <STAR>.unmapped_1.fastq.gz | head -1000) \
     <(zcat <RUS>.unmapped_1.fastq.gz | head -1000)
```

### 10. Transcriptome BAM-specific behaviour

Beyond the known paired-end mate-fields bug, the transcriptome BAM also:

- **Omits per-record `RG:Z:` tags** despite the `@RG` header line being present (issue 5 in TL;DR).
- **Adds `AS:i:` and `nM:i:`** that STAR's transcriptome BAM doesn't emit. STAR's transcriptome BAM tags are exactly `{NH, HI, RG}`; rustar's are `{NH, HI, AS, nM}`. Whether this is a feature or a divergence depends on downstream tool expectations. Salmon ignores both `AS` and `NM` in alignment-mode quant, so this is downstream-neutral for the current pipeline. **BEHAVIOURAL** (potentially nice-to-have for other consumers).
- Sort order matches STAR's (input-read-order, unsorted).

Primary / secondary record counts are within 2-3% of STAR's on every sample - the 1 600-record gap on WT_REP2 (74 772 vs 73 172) is real reads that rustar didn't pass to the transcriptome filter; almost certainly tied to the splice-miss issue (category 4).

## Consolidated issues to file at https://github.com/scverse/rustar-aligner

| #   | Severity   | Title                                                                                                              | One-line repro                                                                                                                                |
| --- | ---------- | ------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | **high**   | `NM` tag is renamed to `nM` and re-defined to exclude indels                                                       | `samtools view <rustar>.Aligned.out.bam \| head \| grep -oE '[nN]M:i:[0-9]+'` returns `nM:i:N` only, even when the CIGAR contains `D`/`I` ops |
| 2   | **high**   | `XS` tag (intron-motif strand) never emitted with `--outSAMstrandField intronMotif`                                | `samtools view <rustar>.Aligned.out.bam \| grep -c 'XS:A:'` returns 0; STAR-equivalent returns ~1 % of reads (1 300 / 90 716 on WT_REP2)      |
| 3   | **medium** | `--sjdbGTFfile` not seeded into pass-1 alignment; `Annotated (sjdb)` always 0, ~50 % of splices missed             | `grep 'Annotated (sjdb)' <rustar>.Log.final.out` -> `\| 0` on every sample with `--twopassMode Basic --sjdbGTFfile genome.gtf`                |
| 4   | **medium** | More secondary alignments emitted than STAR with same `--outFilterMultimapNmax`; NH tail extends to 20 vs STAR's 7 | `samtools view -f 0x100 -c <rustar>.Aligned.out.bam` returns ~17 % more than STAR on identical input                                          |
| 5   | **low**    | Transcriptome BAM omits per-record `RG:Z:` tag despite `@RG` header being present                                  | `samtools view <rustar>.Aligned.toTranscriptome.out.bam \| head -1 \| grep -c 'RG:Z:'` returns 0; STAR returns 1                              |
| 6   | **low**    | `@PG` header line lacks `PN`, `VN`, `CL` fields                                                                    | `samtools view -H <rustar>.Aligned.out.bam \| grep '^@PG'` shows `@PG ID:rustar-aligner` only; provenance is lost                             |
| 7   | **low**    | `AS` tag value disagrees with STAR by 2-5 on identical-CIGAR records                                               | 864 records in WT_REP2; rustar's README claims byte-equivalence on uniquely-mapped AS, so this contradicts the spec                           |

Plus the existing high-severity issue (paired-end transcriptome BAM mate-pair fields), already drafted in `rustar_investigation_wt_rep2.md` and ready to file.

**Filing strategy**: 1 and 2 are independent and both affect the nf-core/rnaseq downstream stack directly (samtools stats, infer_experiment.py); file as separate issues. 3-5 are tied to per-feature scope in the rustar source so are each a discrete issue. 4 is borderline-behavioural and may want a clarifying issue before a fix request. 6-7 can be folded into a single "BAM provenance / fidelity" issue if the maintainer prefers consolidated reports.

## What I couldn't measure and why

- **Unmapped FASTQ contents (category 9)**: pipeline doesn't pass `--outReadsUnmapped Fastx`. Needs a one-line patch to both ALIGN_STAR modules to generate output; cmds in the category-9 section.
- **`jM` / `jI` (junction-motif and intron-coordinate tags)**: not requested by `--outSAMattributes NH HI AS NM MD`. Add them to the CLI on both aligners and re-compare; expectation is that rustar omits both (consistent with the strand-tag gap), but unverified.
- **Strand-orientation correctness of `XS` once it lands**: cannot be measured until `XS` is emitted by rustar at all. Once category-2 issue is fixed, follow up with a comparison of `XS:A:+/-` values vs STAR's on reads with annotated splices.
- **Coordinate-sorted BAM byte-equivalence**: not measured because rnaseq's downstream `SAMTOOLS_SORT` step rewrites the BAM anyway, so any byte-level divergence in the aligner output is moot. The interesting comparison is at the post-sort BAM stage and would need the published `results-*/star_salmon/<sample>.markdup.bam` rather than the work-dir output.
- **MD tag correctness on the rustar-only splice case** (where rustar emits unspliced `101M` for a read STAR splices): would need to align both manually against the reference to score correctness. Visually, both look internally consistent; the rustar `MD:Z:101` (no mismatches in the unspliced alignment) is just less informative than STAR's splice-aware alignment.

## Appendix: raw analysis output

Per-category counts, full distributions, examples per category for the two fully-analysed samples are on disk at `/home/ubuntu/rustar_bam_inv/reports/wt_rep2_genome.txt` and `/home/ubuntu/rustar_bam_inv/reports/rap1un1_genome.txt`. The driver scripts are in `/home/ubuntu/rustar_bam_inv/analyze.py` (top-level per-category counts), `analyze2.py` (NH transition matrix + CIGAR shape categorisation), `analyze3.py` (multi-mapper primary-locus + NM/nM/MD/AS equivalence on identical-CIGAR records). All three are deterministic and re-runnable from the per-sample name-sorted SAMs in `/home/ubuntu/rustar_bam_inv/sorted/`.
