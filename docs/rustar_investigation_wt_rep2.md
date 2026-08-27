# rustar-aligner v0.1.0: transcriptome BAM is missing mate-pair fields, breaking Salmon TPM (WT_REP2 case study)

Date of investigation: 2026-05-12. Pipeline branch: `rustar-aligner` (nf-core/rnaseq PR #1855). Aligner under test: `rustar-aligner` v0.1.0 from https://github.com/scverse/rustar-aligner. Reference aligner: STAR 2.7.11b.

> **Status:** filed upstream as [scverse/rustar-aligner#22](https://github.com/scverse/rustar-aligner/issues/22).

## TL;DR

In nf-core/rnaseq `-profile test,docker`, gene-level Salmon TPMs correlate >= 0.995 between STAR and rustar for four of five samples, but drop to **Pearson r = 0.985040 on WT_REP2** while gene-level NumReads stays at r = 0.999842. The split between very-high read agreement and noticeably-lower TPM agreement points at Salmon's effective-length / fragment-length-distribution stage, not at the alignments themselves.

Root cause, confirmed at three levels (Salmon meta_info, raw BAM records, rustar source):

- **rustar's `Aligned.toTranscriptome.out.bam` does not populate mate-pair fields** (`RNEXT`, `PNEXT`, `TLEN`) and **does not set the proper-pair flag (0x2)**. Every paired-end record looks single-end to Salmon.
- Salmon falls back to its default fragment-length prior (mean = 250 bp, SD = 25 bp) instead of inferring the actual distribution (~168 bp, SD ~71 bp on these data).
- The wrong fragment-length distribution shrinks Salmon's `EffectiveLength` for short transcripts disproportionately (e.g. snR18 transcript: STAR 13.01 -> rustar 3.54, a 73 % drop), so per-transcript TPMs balloon while NumReads stay nearly identical.
- The pattern reproduces on **every paired-end sample** in the test profile (WT_REP1, WT_REP2, RAP1_IAA_30M_REP1 all show `frag_length_mean = 250.0, sd = 25.0` from rustar but inferred values from STAR). WT_REP2 just happens to have the largest TPM Pearson hit because of how its multi-mapper / short-feature mass distributes.
- Bug location in rustar source: `src/lib.rs::build_transcriptome_records_pe` lines 762-768 only OR-stamps `SEGMENTED | FIRST_SEGMENT` / `LAST_SEGMENT`; the underlying `src/io/sam.rs::SamWriter::build_transcriptome_records` (lines 566-660) never sets `PROPERLY_SEGMENTED`, `mate_reference_sequence_id_mut`, `mate_alignment_start_mut`, or `template_length_mut`.

## Reproduction from a clean clone

```bash
# 1. Clone the pipeline branch
git clone --branch rustar-aligner https://github.com/nf-core/rnaseq.git
cd rnaseq

# 2. params file (rustar.params.yml)
cat > rustar.params.yml <<'YAML'
use_rustar_star: true
save_align_intermeds: true   # keep BAMs for inspection; off in normal CI
outdir: results-rustar
YAML

# 3. Run rustar
nextflow run . -profile test,docker -params-file rustar.params.yml -resume

# 4. Reference STAR run (same inputs, save_align_intermeds=true)
nextflow run . -profile test,docker --outdir results-star --save_align_intermeds -resume
```

Both runs use:

- Samplesheet: `https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a/samplesheet/v3.10/samplesheet_test.csv` (yeast subset + GFP transgene, 5 samples, 3 paired-end, 2 single-end).
- STAR / rustar CLI args (from `.command.sh` in both work dirs): identical apart from the binary name -

  ```
  --genomeDir star --readFilesIn <m1> <m2> --runThreadN 4 --outFileNamePrefix WT_REP2.
  --sjdbGTFfile genome_gfp.gtf
  --quantMode TranscriptomeSAM --outSAMtype BAM Unsorted
  --outSAMattributes NH HI AS NM MD --readFilesCommand zcat
  --twopassMode Basic --runRNGseed 0 --outFilterMultimapNmax 20
  --alignSJDBoverhangMin 1 --outSAMstrandField intronMotif
  --quantTranscriptomeSAMoutput BanSingleEnd
  --outSAMattrRGline 'ID:WT_REP2' 'SM:WT_REP2'
  ```

## Methodology

All artefacts live on the `nf-dev-rnaseq` VM at:

- `/home/ubuntu/rnaseq-rustar-aligner/results-star/`
- `/home/ubuntu/rnaseq-rustar-aligner/results-rustar/`

Work directories used in this investigation (BAMs only exist in `work/` because `save_align_intermeds=false` was the default on the run):

- STAR transcriptome BAM: `work/e0/b8c327bfa0a964eecd894bcb05b569/WT_REP2.Aligned.toTranscriptome.out.bam`
- rustar transcriptome BAM: `work/d7/43755befdffb99383bedb820e900f9/WT_REP2.Aligned.toTranscriptome.out.bam`

What was measured:

1. **Salmon `meta_info.json` per sample** (`results-*/star_salmon/<sample>/aux_info/meta_info.json`): `frag_length_mean`, `frag_length_sd`, `num_eq_classes`, `num_mapped`.
2. **Per-transcript / per-gene divergence** from `quant.sf` (script at `/tmp/rustar_inv/analyze.py`, embedded in the Appendix).
3. **STAR-style alignment stats** from `WT_REP2.Log.final.out` in the corresponding rustar/STAR work dirs.
4. **BAM-level inspection** via samtools 1.23 in the `community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5` container:
   - `samtools flagstat` on both transcriptome BAMs.
   - `samtools view ... | head` to inspect raw FLAG / RNEXT / PNEXT / TLEN fields.
   - NH-tag distributions via `samtools view | grep -oE 'NH:i:[0-9]+' | sort | uniq -c`.
5. **Source inspection** of:
   - rustar: `src/io/sam.rs`, `src/lib.rs`, `src/quant/transcriptome.rs` at https://github.com/scverse/rustar-aligner @ `main`.
   - Salmon: `include/salmon/internal/config/SalmonDefaults.hpp`, `src/alignment/SalmonQuantifyAlignments.cpp`, `include/salmon/internal/alignment/BAMQueue.tpp` at https://github.com/COMBINE-lab/salmon @ `master`.

## Evidence

### 1. Salmon's frag-length inference fails on every rustar paired-end sample

From `aux_info/meta_info.json`:

| Sample              | Layout | STAR frag mean | STAR frag SD | rustar frag mean | rustar frag SD | STAR eq classes | rustar eq classes |
| ------------------- | ------ | -------------- | ------------ | ---------------- | -------------- | --------------- | ----------------- |
| RAP1_IAA_30M_REP1   | PE     | 168.84         | 66.09        | **250.00**       | **25.00**      | 113             | **429**           |
| RAP1_UNINDUCED_REP1 | SE     | 250.00         | 25.00        | 250.00           | 25.00          | 94              | 95                |
| RAP1_UNINDUCED_REP2 | SE     | 250.00         | 25.00        | 250.00           | 25.00          | 107             | 104               |
| WT_REP1             | PE     | 168.18         | 68.07        | **250.00**       | **25.00**      | 115             | **476**           |
| WT_REP2             | PE     | 167.85         | 70.64        | **250.00**       | **25.00**      | 100             | **392**           |

`fragLenPriorMean = 250.0` and `fragLenPriorSD = 25.0` are Salmon's hard-coded defaults (`SalmonDefaults.hpp:58-59`). For paired-end samples, STAR's BAM lets Salmon infer the real distribution; rustar's BAM does not. Single-end samples land on the prior for both (Salmon cannot infer from single-end). The blow-up in equivalence-class count for rustar paired-end samples (3-4x) is a downstream consequence of the wrong frag length distribution changing which (transcript, position) combinations look plausible.

### 2. flagstat: rustar's transcriptome BAM has **zero properly-paired** reads

```text
=== STAR WT_REP2.Aligned.toTranscriptome.out.bam ===
84406 total, 74772 primary, 9634 secondary
74772 paired in sequencing  (37386 read1 + 37386 read2)
74772 properly paired (100.00%)
0 with mate mapped to a different chr

=== rustar WT_REP2.Aligned.toTranscriptome.out.bam ===
82744 total, 73172 primary, 9572 secondary
73172 paired in sequencing  (36586 read1 + 36586 read2)
0 properly paired (0.00%)
73172 with mate mapped to a different chr   <-- 100% of mates
63114 with mate mapped to a different chr (mapQ>=5)
```

Verbatim record comparison for fragment `SRR6357072.1900203`:

| Field     | STAR (r1 / r2)        | rustar (r1 / r2)      |
| --------- | --------------------- | --------------------- |
| FLAG      | 163 / 83              | 81 / 129              |
| RNAME     | `YAR010C` / `YAR010C` | `YAR010C` / `YAR010C` |
| POS       | 1045 / 1103           | 1103 / 1045           |
| MAPQ      | 255 / 255             | 255 / 255             |
| CIGAR     | `101M` / `101M`       | `101M` / `101M`       |
| **RNEXT** | `=` / `=`             | **`*`** / **`*`**     |
| **PNEXT** | 1103 / 1045           | **0** / **0**         |
| **TLEN**  | 159 / -159            | **0** / **0**         |

The alignment positions are equivalent (within mate-strand swap) - rustar genuinely finds the same fragment. It just refuses to tell Salmon that the two records are mates. STAR's FLAG 163 contains `PROPERLY_SEGMENTED` (0x2); rustar's FLAG 81 does not.

### 3. NH/HI multi-mapper structure is essentially the same

| Tag    | STAR records | rustar records |
| ------ | ------------ | -------------- |
| NH:i:1 | 65 200       | 63 646         |
| NH:i:2 | 19 084       | 19 004         |
| NH:i:3 | 18           | 18             |
| NH:i:4 | 80           | 64             |
| NH:i:6 | 24           | 12             |

So this is not a tie-breaking / multi-mapper assignment problem. The total multi-mapper population is nearly identical; the issue is the SAM bookkeeping for paired records.

### 4. Per-gene divergence in WT_REP2 maps cleanly onto EffectiveLength change

Distribution of `EffectiveLength_rustar / EffectiveLength_star` across 125 transcripts in WT_REP2 (script in Appendix):

```
count    125.000000
mean       0.743381
std        0.255925
min        0.139372
5%         0.277408
25%        0.541054
50%        0.865689
75%        0.947846
95%        0.977230
max        0.981523
```

Top 6 genes by absolute TPM delta in WT_REP2 (NumReads agreement is excellent; EffectiveLength change drives the TPM swing):

| gene_id | TPM_star | TPM_rustar | TPM_delta | NumReads_star | NumReads_rustar | Length | EffLen_star | EffLen_rustar | EffLen_ratio |
| ------- | -------: | ---------: | --------: | ------------: | --------------: | -----: | ----------: | ------------: | -----------: |
| snR18   |   28 641 | **93 080** |   +64 439 |          8.00 |            8.00 |    102 |       13.01 |      **3.54** |         0.27 |
| YAR009C |  229 302 |    201 000 |   -28 302 |        16 852 |          16 290 |  3 591 |    3 423.15 |      3 341.00 |         0.98 |
| YAL038W |  220 550 |    207 930 |   -12 620 |         6 322 |           6 320 |  1 503 |    1 335.15 |      1 253.00 |         0.94 |
| YAR010C |  245 119 |    235 742 |    -9 376 |         6 079 |           6 136 |  1 323 |    1 155.15 |      1 073.00 |         0.93 |
| YAL005C |  109 599 |    102 066 |    -7 533 |         4 144 |           4 157 |  1 929 |    1 761.15 |      1 679.00 |         0.95 |
| YAL003W |   72 309 |     65 003 |    -7 306 |           704 |             585 |    621 |      453.49 |        371.00 |         0.82 |

Gene NumReads agreement summary:

- 100 / 125 genes agree to < 1 read.
- 121 / 125 to < 5 reads.
- Max |NumReads_delta| over all 125 genes = 562 (YAR009C, 16 852 vs 16 290 - 3.3 %).

Top by relative TPM delta (filtered to NumReads_star >= 5): snR18 +225 %, YAL046C +61 %, YAL026C-A +27 %, YAL039C +22 %. All of these are short transcripts whose `EffectiveLength` is dominated by the assumed fragment-length distribution: the prior (mean 250) chops a fixed-ish absolute number of bases off each transcript regardless of length, but a 102 bp transcript gets reduced from 13 bp -> 3.54 bp (73 % shrink) while a 3 591 bp transcript only shrinks from 3 423 -> 3 341 (2 % shrink).

Because TPMs are normalised so that `sum(TPM) = 1e6`, a small over-shrinking of effective length on short transcripts steals TPM mass from long ones. That's exactly the per-gene pattern we see in the table above (snR18 +225 %, while every long transcript is mildly negative).

### 5. rustar `Log.final.out` reports `Annotated (sjdb) = 0` on WT_REP2

| Metric                                  | STAR           | rustar         |
| --------------------------------------- | -------------- | -------------- |
| Number of splices: Total                | 762            | 371            |
| **Number of splices: Annotated (sjdb)** | **644**        | **0**          |
| Number of splices: GT/AG                | 724            | 276            |
| Number of splices: Non-canonical        | 35             | 92             |
| Reads unmapped: too short               | 1 540 (3.11 %) | 4 193 (8.46 %) |
| Reads unmapped: other                   | 2 656 (5.36 %) | 0              |

This is consistent with rustar honouring `--sjdbGTFfile` for index lookup but never crediting annotation-derived junctions to the `sjdb` bucket in the log. The total-splice gap (371 vs 762) and the higher non-canonical rate (92 vs 35) are also consistent with rustar's two-pass step not seeding pass-1 from the GTF SJ database, but this is a **second, lower-severity issue** and is _not_ the cause of the TPM divergence; the BAM bookkeeping is.

## Root cause in rustar source

Tracing the BAM output path:

- Entry point for paired transcriptome BAM emission: `src/lib.rs::build_transcriptome_records_pe` ( https://github.com/scverse/rustar-aligner/blob/main/src/lib.rs ), lines 677-777.
- Per-mate record construction: calls `SamWriter::build_transcriptome_records` twice, once for mate1 (line 743) and once for mate2 (line 752), each as if it were a separate single-end alignment list.
- After both calls return, the wrapper stamps paired flags at lines 762-768:

```rust
use noodles::sam::alignment::record::Flags;
for r in rec1s.iter_mut() {
    *r.flags_mut() |= Flags::SEGMENTED | Flags::FIRST_SEGMENT;
}
for r in rec2s.iter_mut() {
    *r.flags_mut() |= Flags::SEGMENTED | Flags::LAST_SEGMENT;
}
```

That is the only post-processing applied to paired transcriptome records. Notably **absent**:

- `Flags::PROPERLY_SEGMENTED` (0x2) - never set on any transcriptome record.
- `Flags::MATE_REVERSE_COMPLEMENTED` (0x20) - never set; STAR sets this on the mate that's on the opposite strand.
- `record.mate_reference_sequence_id_mut()` - never set, so RNEXT remains `*` (the default).
- `record.mate_alignment_start_mut()` - never set, so PNEXT remains 0.
- `record.template_length_mut()` - never set, so TLEN remains 0.

Compare with the _genome_ BAM emission, which does set these correctly: `src/io/sam.rs::build_paired_records` line 276 onward, and `build_paired_mate_record` lines 1163-1260 which explicitly sets PROPERLY_SEGMENTED at line 1187-1188:

```rust
if is_proper_pair {
    flags |= sam::alignment::record::Flags::PROPERLY_SEGMENTED; // 0x2
}
```

and mate fields at lines 1239-1251:

```rust
*record.mate_reference_sequence_id_mut() = Some(mate_transcript.chr_idx);
*record.mate_alignment_start_mut()       = Some(mate_pos.try_into()...);
*record.template_length_mut()            = insert_size;
```

So the bookkeeping logic exists; it just isn't applied to the transcriptome-space path.

The underlying record builder, `SamWriter::build_transcriptome_records` ( `src/io/sam.rs:566-660` ), constructs flags as:

```rust
let mut flags = sam::alignment::record::Flags::empty();
if t.is_reverse {
    flags |= sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
}
if hit_idx != primary_hit_idx {
    flags |= sam::alignment::record::Flags::SECONDARY;
}
*record.flags_mut() = flags;
```

i.e. it has no concept of "mate".

## Why this breaks Salmon specifically

From Salmon master:

- `include/salmon/internal/config/SalmonDefaults.hpp:58-59`:

  ```cpp
  constexpr const double fragLenPriorMean{250.0};
  constexpr const double fragLenPriorSD{25.0};
  ```

- `src/alignment/SalmonQuantifyAlignments.cpp` `processMiniBatch`, lines 1218-1223:

  ```cpp
  if (aln->isPaired() and !salmonOpts.noFragLengthDist) {
      double fragLength = aln->fragLengthPedantic(transcript.RefLength);
      if (fragLength > 0) {
          fragLengthDist.addVal(fragLength, logForgettingMass);
      }
  }
  ```

  `fragLengthPedantic` derives the length from TLEN / mate position. With TLEN = 0 and PNEXT = 0, it returns 0, and the update is skipped.

- `src/alignment/SalmonQuantifyAlignments.cpp` line ~1107: pair-conditional fragment-probability uses `aln->isInward()`, which also requires populated mate-position info.

- `include/salmon/internal/alignment/BAMQueue.tpp` line 360 onward: pair classification uses `BAM_FPROPER_PAIR` to decide `AlignmentType::MappedConcordantPair`. Without 0x2, records get classified differently and the equivalence-class assignment shifts (consistent with the 100 -> 392 eq-class blow-up we see).

Net effect: Salmon never updates the empirical fragment-length distribution, falls back to its prior (mean 250, SD 25), and computes `EffectiveLength` from the prior. STAR's TLEN-bearing transcriptome BAM avoids this entirely.

## Hypothesis

Two issues, only the first matters for this PR:

1. **(Blocker)** `build_transcriptome_records_pe` in rustar v0.1.0 emits transcriptome BAM records without the proper-pair flag and without RNEXT/PNEXT/TLEN populated. This makes the BAM unusable as input to Salmon's alignment-mode quant for paired-end libraries (TPMs systematically wrong; NumReads close-but-not-identical). This is not a Salmon-specific concern: any downstream tool that derives fragment-length statistics from the transcriptome BAM (e.g. RSEM, kallisto's BAM mode, custom fragment-size QC) will be affected.

2. **(Lower severity; not yet investigated end-to-end)** rustar's `Log.final.out` reports `Number of splices: Annotated (sjdb) = 0` even when `--sjdbGTFfile` is provided with `--twopassMode Basic`. The total-splice count is also ~half of STAR's. This is consistent with the GTF-derived SJ database not being used to seed pass-1 alignment in the way STAR does. The mapping-rate impact in our test (-0.15 to -0.21 pp) is tiny, but on larger genomes with more annotated alternative splicing this could matter more.

## What would falsify hypothesis 1

The cleanest one-shot falsification: pick any paired-end BAM produced by rustar v0.1.0 with `--quantMode TranscriptomeSAM`, count records with the 0x2 flag set:

```bash
samtools view -c -f 2 WT_REP2.Aligned.toTranscriptome.out.bam   # rustar
samtools view -c -f 2 WT_REP2.Aligned.toTranscriptome.out.bam   # STAR (reference)
```

If rustar's count is > 0, the hypothesis is wrong as stated. With our run it is exactly 0.

A secondary, fragment-length-distribution check (also one-line):

```bash
jq '.frag_length_mean, .frag_length_sd' results-*/star_salmon/WT_REP2/aux_info/meta_info.json
```

If rustar reports anything other than 250.0 / 25.0 _and_ WT_REP2 Pearson r still drops below 0.99 on gene TPM, the proximate cause is elsewhere.

A third check, downstream of the fix: re-run Salmon on rustar's existing transcriptome BAM with `--fldMean 168 --fldSD 71` (forcing the right prior) and verify that gene TPM Pearson rises to >= 0.995. If it does, that confirms the EffectiveLength path is the sole driver; if it doesn't, there's a second problem hiding behind this one.

## Suggested upstream issue

**Title**: `transcriptome BAM (--quantMode TranscriptomeSAM) omits mate-pair fields for paired-end reads, breaks Salmon TPM`

**Body**:

````markdown
## Summary

For paired-end input, rustar-aligner v0.1.0's `Aligned.toTranscriptome.out.bam`
emits each mate as if it were a separate single-end alignment record. The
proper-pair flag (0x2), `RNEXT`, `PNEXT`, and `TLEN` are not set, even though
both mates are present in the file and aligned to the same transcript at
sensible positions. This makes the BAM unsuitable as input to Salmon's
alignment-mode quantification (and any other paired-aware downstream tool).

## Repro

Standard nf-core/rnaseq `-profile test,docker` with `--use_rustar_star`
(PR nf-core/rnaseq#1855), or equivalently any paired-end FASTQ + GTF +
`--quantMode TranscriptomeSAM --twopassMode Basic --runRNGseed 0
--quantTranscriptomeSAMoutput BanSingleEnd`.

## Evidence

`samtools flagstat WT_REP2.Aligned.toTranscriptome.out.bam`:

| Field                               |   STAR 2.7.11b | rustar 0.1.0 |
| ----------------------------------- | -------------: | -----------: |
| properly paired                     | 74 772 (100 %) |  **0 (0 %)** |
| with mate mapped to a different chr |              0 |       73 172 |

First record pair (`SRR6357072.1900203`):

| Field | STAR        | rustar    |
| ----- | ----------- | --------- |
| FLAG  | 163 / 83    | 81 / 129  |
| RNAME | `YAR010C`   | `YAR010C` |
| RNEXT | `=`         | **`*`**   |
| PNEXT | 1103 / 1045 | **0 / 0** |
| TLEN  | 159 / -159  | **0 / 0** |

Salmon `aux_info/meta_info.json` `frag_length_mean` / `frag_length_sd`:

- STAR: inferred 167.85 / 70.64
- rustar: default prior 250.00 / 25.00

The Salmon prior (`SalmonDefaults.hpp:58-59`) is `fragLenPriorMean=250.0,
fragLenPriorSD=25.0` - exactly the values rustar's output causes Salmon to
fall back on. The `EffectiveLength` shift inflates TPMs of short transcripts
disproportionately, e.g. snR18 (102 bp): EffLen 13.01 -> 3.54, TPM 28 641 ->
93 080 with NumReads unchanged at 8.

## Root cause

`src/lib.rs::build_transcriptome_records_pe`, lines 762-768:

```rust
for r in rec1s.iter_mut() {
    *r.flags_mut() |= Flags::SEGMENTED | Flags::FIRST_SEGMENT;
}
for r in rec2s.iter_mut() {
    *r.flags_mut() |= Flags::SEGMENTED | Flags::LAST_SEGMENT;
}
```
````

This is the only post-processing applied to paired transcriptome records.
`PROPERLY_SEGMENTED`, `MATE_REVERSE_COMPLEMENTED`,
`mate_reference_sequence_id_mut`, `mate_alignment_start_mut`, and
`template_length_mut` are never set. The genome-space paired path in
`src/io/sam.rs::build_paired_mate_record` (lines 1163-1260) shows the
intended logic; it just isn't wired into the transcriptome path.

## Suggested fix

In `build_transcriptome_records_pe`, after the two `build_transcriptome_records`
calls and the existing flag-stamping loops, iterate over the interleaved
(rec1, rec2) pairs and set:

- `PROPERLY_SEGMENTED` on both (both mates landed on the same transcript by
  construction at line 725 of `lib.rs`, so concordance is guaranteed).
- `MATE_REVERSE_COMPLEMENTED` on r1 if r2 is reverse, and vice versa
  (`REVERSE_COMPLEMENTED` is already set per-record by
  `build_transcriptome_records`).
- `mate_reference_sequence_id_mut(Some(<paired transcript chr_idx>))`.
- `mate_alignment_start_mut(Some(<paired transcript POS>))`.
- `template_length_mut(<signed TLEN>)` - same sign convention as STAR
  (leftmost mate gets +TLEN, rightmost gets -TLEN; TLEN = mate_end -
  this_start where applicable).

Alternatively, extract a shared helper from the genome-space
`build_paired_mate_record` so the transcriptome path goes through the same
mate-bookkeeping code rather than duplicating it.

## Why this matters

Any downstream tool that interprets the transcriptome BAM as paired (Salmon
alignment-mode, RSEM, tximeta wrappers, custom fragment-size QC) will either
crash, warn, or silently degrade. Salmon silently degrades, which is
arguably the worst outcome - users get plausible-looking TPMs that are
systematically biased toward short transcripts.

## Test plan

Add an integration test that:

1. Runs rustar with `--quantMode TranscriptomeSAM` on a paired-end FASTQ.
2. Asserts `samtools view -c -f 2 <bam>` == count of primary paired records.
3. Asserts no records have `RNEXT = '*'` for both-mate-mapped reads.
4. (Optional) Runs Salmon downstream and asserts `meta_info.json
.frag_length_mean` is not equal to 250.0 (i.e. inference succeeded).

````

## What a follow-up investigator should pick up

If the upstream issue confirms the fix scope, the nf-core/rnaseq side needs:

- A version bump or a containers pin to a rustar release that includes the
  fix.
- A focused regression test that runs the test profile with
  `--use_rustar_star` and checks gene-TPM Pearson against the saved snapshot
  / STAR run (threshold > 0.995 across all samples). The current PR's
  `bin/compare_aligner_runs.py` is sufficient for the check; just add an
  assertion at the top of its main on the gene-TPM column for paired-end
  samples.

If the upstream maintainer asks for a minimal failing input, the WT_REP2
FASTQs are at:

- `https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_1.fastq.gz`
- `https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_2.fastq.gz`

The yeast + GFP reference and GTF are reachable from `conf/test.config` in
this branch.

The SJ-annotation issue (rustar's `Number of splices: Annotated (sjdb) = 0`
with `--sjdbGTFfile`) is separate and lower-severity and should be filed
as its own issue once the BAM fix lands. Mapping-rate impact in our test
profile is < 0.25 pp.

## Appendix: helper script

`/tmp/rustar_inv/analyze.py` (run with the `nf-core` micromamba env active):

```python
#!/usr/bin/env python3
"""Per-transcript / per-gene divergence analysis for WT_REP2 STAR vs rustar."""
import pandas as pd
import numpy as np

star   = pd.read_csv("/tmp/rustar_inv/quant.star.sf",   sep="\t")
rustar = pd.read_csv("/tmp/rustar_inv/quant.rustar.sf", sep="\t")
tx2g   = pd.read_csv("/tmp/rustar_inv/tx2gene.tsv", sep="\t", header=None,
                     names=["transcript_id", "gene_id", "gene_name"])

m = star.merge(rustar, on="Name", suffixes=("_star", "_rustar"))
m = m.merge(tx2g, left_on="Name", right_on="transcript_id", how="left")
m["EffLen_ratio"]  = m["EffectiveLength_rustar"] / m["EffectiveLength_star"]
m["TPM_delta"]     = m["TPM_rustar"] - m["TPM_star"]
m["TPM_abs_delta"] = m["TPM_delta"].abs()
m["NumReads_delta"]= m["NumReads_rustar"] - m["NumReads_star"]

g = m.groupby("gene_id", dropna=False).agg(
    TPM_star=("TPM_star", "sum"),
    TPM_rustar=("TPM_rustar", "sum"),
    NumReads_star=("NumReads_star", "sum"),
    NumReads_rustar=("NumReads_rustar", "sum"),
    Length_star=("Length_star", "max"),
    EffLen_star=("EffectiveLength_star", "max"),
    EffLen_rustar=("EffectiveLength_rustar", "max"),
)
g["TPM_delta"]      = g["TPM_rustar"] - g["TPM_star"]
g["TPM_abs_delta"]  = g["TPM_delta"].abs()
g["EffLen_ratio"]   = g["EffLen_rustar"] / g["EffLen_star"]

print(g.sort_values("TPM_abs_delta", ascending=False).head(20).round(2))
````

Files needed (copied from the VM):

```bash
scp nf-dev-rnaseq:/home/ubuntu/rnaseq-rustar-aligner/results-star/star_salmon/WT_REP2/quant.sf   /tmp/rustar_inv/quant.star.sf
scp nf-dev-rnaseq:/home/ubuntu/rnaseq-rustar-aligner/results-rustar/star_salmon/WT_REP2/quant.sf /tmp/rustar_inv/quant.rustar.sf
scp nf-dev-rnaseq:/home/ubuntu/rnaseq-rustar-aligner/results-star/star_salmon/salmon.merged.tx2gene.tsv /tmp/rustar_inv/tx2gene.tsv
```

BAM inspection (samtools is not on the VM; use the htslib_samtools wave container):

```bash
ssh nf-dev-rnaseq
IMG=community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5
STAR_BAM=/home/ubuntu/rnaseq-rustar-aligner/work/e0/b8c327bfa0a964eecd894bcb05b569/WT_REP2.Aligned.toTranscriptome.out.bam
RUSTAR_BAM=/home/ubuntu/rnaseq-rustar-aligner/work/d7/43755befdffb99383bedb820e900f9/WT_REP2.Aligned.toTranscriptome.out.bam

docker run --rm -v /home/ubuntu/rnaseq-rustar-aligner/work:/work $IMG \
    samtools flagstat /work/e0/b8c327bfa0a964eecd894bcb05b569/WT_REP2.Aligned.toTranscriptome.out.bam
docker run --rm -v /home/ubuntu/rnaseq-rustar-aligner/work:/work $IMG \
    samtools flagstat /work/d7/43755befdffb99383bedb820e900f9/WT_REP2.Aligned.toTranscriptome.out.bam
```
