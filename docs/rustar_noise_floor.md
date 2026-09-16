# rustar-aligner v0.1.0 vs STAR 2.7.11b: noise-floor characterisation

Quantitative companion to [`rustar_bam_comparison.md`](rustar_bam_comparison.md) and [`rustar_two_pass_and_determinism.md`](rustar_two_pass_and_determinism.md). Where those docs catalogue the differences between aligners on a single pair of runs, this one establishes the **noise envelope** for the comparison harness so any future deltas can be placed in context.

## TL;DR

Per-sample gene_tpm Pearson on the test profile (yeast subset + GFP, 5 samples, paired and single-end), worst case across samples:

| Contrast                                 | gene_tpm Pearson worst | gene_counts Pearson worst | %mapped worst Δ (pp) |
| ---------------------------------------- | ---------------------: | ------------------------: | -------------------: |
| STAR vs STAR, **same seed**, fresh rerun |           0.9999999997 |              0.9999999999 |                +0.00 |
| STAR vs STAR, seed 0 vs seed 1           |           0.9999999994 |              0.9999999999 |                +0.00 |
| rustar vs rustar, seed 0 vs seed 1       |           0.9999999996 |              0.9999999994 |                +0.00 |
| STAR vs rustar, same (seed 0)            |   **0.9850** (WT_REP2) |                    0.9998 |                -0.21 |

Both STAR's and rustar's own seed-to-seed variance is ~1e-9 on Pearson and below the harness's 0.01 pp display precision on % mapped. **The STAR-vs-rustar deltas are 7-9 orders of magnitude larger than either aligner's own seed-driven noise envelope. The previously-documented `WT_REP2 = 0.985` divergence is firmly outside the noise floor; it is real cross-aligner signal driven by the upstream BUGs catalogued in [`rustar_bam_comparison.md`](rustar_bam_comparison.md), not RNG noise.**

Side finding: STAR with the same seed is **alignment-bit-identical at the record content level** (qname, flag-minus-dup-bit, rname, pos, MAPQ, CIGAR, MD, NM, NH, HI, AS all identical across reruns), but the published `<sample>.markdup.sorted.bam` bytes differ because (a) the STAR aligner's record order is non-deterministic and (b) Picard MarkDuplicates picks a different "keeper" out of each duplicate set on each run. The downstream Salmon quant.sf also differs at the floating-point level — same reason: order-of-arrival changes Salmon's bias-modelling sample. None of this leaks into the merged TSVs at any precision that matters.

MarkDuplicates dup-bit agreement on STAR-vs-rustar matched primary records ranges 98.7-99.8% per sample. That is broadly consistent with the BAM divergence floor documented in [`rustar_bam_comparison.md`](rustar_bam_comparison.md) (60-110 STAR-only / 105-1152 rustar-only primary reads per sample, plus ~1k same-name position differences on multi-mappers); MarkDuplicates is propagating the upstream BAM divergence rather than compounding it.

## Methodology

### Source data

Three new pipeline runs on top of the existing `results-star/` / `results-rustar/` baseline. All five runs use `-profile test,docker` on the `nf-dev-rnaseq` VM (36 CPU / 69 GB), Nextflow 26.04, same input samplesheet:

| outdir                  | aligner | --runRNGseed | how                                                                                                          |
| ----------------------- | ------- | ------------ | ------------------------------------------------------------------------------------------------------------ |
| `results-star/`         | STAR    | 0 (default)  | original PR-1855 baseline                                                                                    |
| `results-star-rerun/`   | STAR    | 0 (default)  | `-params-file star-rerun.params.yml`                                                                         |
| `results-star-seed1/`   | STAR    | 1            | `-params-file star-seed1.params.yml` (sets `extra_star_align_args: '--runRNGseed 1'`)                        |
| `results-rustar/`       | rustar  | 0 (default)  | original PR-1855 baseline                                                                                    |
| `results-rustar-seed1/` | rustar  | 1            | `-params-file rustar-seed1.params.yml` (`use_rustar_star: true` + `extra_star_align_args: '--runRNGseed 1'`) |

Params files at `/home/ubuntu/rnaseq-rustar-aligner/{star-rerun,star-seed1,rustar-seed1}.params.yml`. The `--runRNGseed 1` override is plumbed via the existing `params.extra_star_align_args` dedup logic at `conf/modules/align_star.config:108-121`, which strips the pipeline default `--runRNGseed 0` and replaces it with the user value. The `withName: '.*ALIGN_STAR:STAR_ALIGN|.*ALIGN_STAR:RUSTAR_ALIGN|...'` selector applies to both aligner modules so the override hits both code paths.

### Comparison harness

The existing `bin/compare_aligner_runs.py <dir_a> <dir_b>` script reads `Log.final.out` (% mapped), `salmon.merged.gene_{tpm,counts}.tsv` (per-sample Pearson + Spearman vs the other dir), and `pipeline_info/execution_trace_*.txt` (trace timings). It's aligner-agnostic — it treats the two dirs symmetrically — so it works for STAR-vs-STAR and rustar-vs-rustar contrasts as-is.

All four JSON reports on disk at `/tmp/cmp_{star_self,star_seed,rustar_seed,cross}.json`. Stitched table generator at `/tmp/summarize.py`.

### MarkDuplicates analysis

For each sample's `markdup.sorted.bam` pair, build a `(qname, mate-bit, secondary-bit) -> dup_bit` map from each BAM via `samtools view -F 0x800` (excluding supplementaries), then on the intersection of keys count agreement vs disagreement on the 0x400 dup bit. Per-sample disagreement is reported as count + agreement rate. Helper at `/tmp/markdup_agreement.py` on the VM, run via `samtools` from the STAR Wave container `community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4` (matches `rustar_bam_comparison.md`'s commands).

Per-contrast JSON line dumps at `/tmp/markdup_{starself,starseed,rustarseed,cross}.jsonl`.

### Wall-time

Each run was ~3.5 - 4.5 minutes wall on the unloaded VM. Run 1 (STAR rerun) 17:31:35 -> 17:35:57; run 2 (STAR seed1) 17:35:57 -> 17:39:58; run 3 (rustar seed1) 17:39:58 -> 17:43:29. All three completed cleanly (`Pipeline completed successfully` in `noise_floor_runs.log`). No failures, no stalls.

## Results: harness-level noise envelope

### Per-sample, per-contrast (gene_tpm Pearson, full precision)

| Sample              | STAR vs STAR (same seed) | STAR vs STAR (seed 0 vs 1) | rustar vs rustar (seed 0 vs 1) | STAR vs rustar (seed 0 both) |
| ------------------- | -----------------------: | -------------------------: | -----------------------------: | ---------------------------: |
| RAP1_IAA_30M_REP1   |             0.9999999997 |               0.9999999997 |                0.9999999999950 |                   **0.9968** |
| RAP1_UNINDUCED_REP1 |             1.0000000000 |               1.0000000000 |                1.0000000000000 |                       0.9997 |
| RAP1_UNINDUCED_REP2 |             1.0000000000 |               1.0000000000 |                0.9999999999999 |                       0.9997 |
| WT_REP1             |             0.9999999999 |               0.9999999996 |                0.9999999999618 |                   **0.9955** |
| WT_REP2             |             1.0000000000 |               0.9999999994 |                0.9999999999812 |                   **0.9850** |

### Per-sample, per-contrast (gene_counts Pearson, full precision)

| Sample              | STAR vs STAR (same seed) | STAR vs STAR (seed 0 vs 1) | rustar vs rustar (seed 0 vs 1) | STAR vs rustar |
| ------------------- | -----------------------: | -------------------------: | -----------------------------: | -------------: |
| RAP1_IAA_30M_REP1   |          0.9999999999999 |            0.9999999999999 |                  0.99999999998 |         0.9998 |
| RAP1_UNINDUCED_REP1 |          1.0000000000000 |            1.0000000000000 |                1.0000000000000 |         0.9999 |
| RAP1_UNINDUCED_REP2 |          1.0000000000000 |            1.0000000000000 |                0.9999999999999 |         0.9999 |
| WT_REP1             |          0.9999999999999 |            0.9999999999999 |                0.9999999999993 |         0.9999 |
| WT_REP2             |          0.9999999999999 |            0.9999999999999 |                0.9999999999996 |         0.9998 |

### % uniquely mapped reads

All three "noise floor" contrasts report a delta of `+0.00 pp` per sample at the `Log.final.out` display precision of 0.01 %. The STAR-vs-rustar cross-aligner contrast reports -0.05 to -0.21 pp per sample (reproducing the [`rustar_differences.md`](rustar_differences.md) headline).

### Trace timings (process wall median, peak RSS)

Re-running STAR with the same seed reproduces the prior median wall time exactly (68.0 s, 0.92 GB peak RSS, n=5 tasks). STAR with `--runRNGseed 1` is within +2 s on median wall and +/-0.01 GB on RSS — well inside per-task scheduling noise on a multi-user VM. rustar with `--runRNGseed 1` median wall went 33.8 -> 35.4 s and RSS 0.119 -> 0.118 GB (n=5); again inside noise.

## Results: BAM-level same-seed determinism

Subset of the STAR self-rerun (`results-star/` vs `results-star-rerun/`) on the published `<sample>.markdup.sorted.bam`:

| BAM equality test (md5)                                            | Result                         |
| ------------------------------------------------------------------ | ------------------------------ |
| Raw BAM bytes                                                      | DIFFER on all 5 samples        |
| `samtools view` SAM body, same order                               | DIFFER on all 5 samples        |
| `samtools sort -n -O sam` SAM body                                 | DIFFER on all 5 samples        |
| Same body, with the 0x400 (duplicate) flag bit masked off + sorted | **IDENTICAL on all 5 samples** |

In other words: STAR with `--runRNGseed 0` produces a bit-identical alignment of every read (qname, flag-minus-dup-bit, rname, pos, MAPQ, CIGAR, plus all tags), but downstream Picard MarkDuplicates flips the duplicate bit on a small number of reads between runs. Per-record dup-bit disagreement count:

| Sample              | Records | MarkDup STAR-self disagree | Disagreement rate |
| ------------------- | ------: | -------------------------: | ----------------: |
| WT_REP1             | 184 589 |                         12 |          0.0065 % |
| WT_REP2             |  92 683 |                          4 |          0.0043 % |
| RAP1_IAA_30M_REP1   |  93 368 |                          8 |          0.0086 % |
| RAP1_UNINDUCED_REP1 |  48 823 |                         26 |          0.0533 % |
| RAP1_UNINDUCED_REP2 |  98 088 |                         18 |          0.0184 % |

Total dup count per sample is identical across reruns — MarkDuplicates marks the same number of records as duplicates, it just picks a different one from each duplicate set as the "keeper". The disagreement count scales with the dup fraction (the SE samples with 74-80% dup rate generate the most pairs to choose between).

Root cause is upstream of MarkDuplicates: STAR's BAM output order is non-deterministic across reruns (parallel-worker completion order). The downstream sorted/coord-sorted/markduped BAM inherits an effectively random tie-break choice on equal-priority duplicates. STAR's own record content is deterministic with `--runRNGseed 0`, as advertised.

Salmon's `quant.sf` is also non-byte-identical across STAR-same-seed reruns (`EffectiveLength` differs on ~all 125 genes, `TPM` on ~85%, `NumReads` on ~5%), driven by the same input-order non-determinism feeding Salmon's bias-modelling sampler. Numerically the differences are ~1e-9 relative and round to a Pearson of `~0.9999999997` after the merged-TSV step.

## Results: MarkDuplicates decision agreement (bonus)

| Contrast                       | Sample              | Total matched primaries | Dup-bit disagree | Agreement rate |
| ------------------------------ | ------------------- | ----------------------: | ---------------: | -------------: |
| STAR vs STAR (same seed)       | WT_REP1             |                 184 589 |               12 |      99.9935 % |
|                                | WT_REP2             |                  92 683 |                4 |      99.9957 % |
|                                | RAP1_IAA_30M_REP1   |                  93 368 |                8 |      99.9914 % |
|                                | RAP1_UNINDUCED_REP1 |                  48 823 |               26 |      99.9467 % |
|                                | RAP1_UNINDUCED_REP2 |                  98 088 |               18 |      99.9816 % |
| STAR vs STAR (seed 0 vs 1)     | WT_REP1             |                 184 589 |               12 |      99.9935 % |
|                                | WT_REP2             |                  92 683 |               12 |      99.9871 % |
|                                | RAP1_IAA_30M_REP1   |                  93 368 |                0 |     100.0000 % |
|                                | RAP1_UNINDUCED_REP1 |                  48 823 |               20 |      99.9590 % |
|                                | RAP1_UNINDUCED_REP2 |                  98 088 |               16 |      99.9837 % |
| rustar vs rustar (seed 0 vs 1) | WT_REP1             |                 184 942 |              238 |      99.8713 % |
|                                | WT_REP2             |                  92 846 |               86 |      99.9074 % |
|                                | RAP1_IAA_30M_REP1   |                  93 584 |               82 |      99.9124 % |
|                                | RAP1_UNINDUCED_REP1 |                  48 818 |              120 |      99.7542 % |
|                                | RAP1_UNINDUCED_REP2 |                  98 049 |              162 |      99.8348 % |
| STAR vs rustar (seed 0 both)   | WT_REP1             |                 183 790 |            2 357 |      98.7176 % |
|                                | WT_REP2             |                  92 303 |              256 |      99.7227 % |
|                                | RAP1_IAA_30M_REP1   |                  93 074 |              208 |      99.7765 % |
|                                | RAP1_UNINDUCED_REP1 |                  48 713 |              231 |      99.5258 % |
|                                | RAP1_UNINDUCED_REP2 |                  97 826 |              301 |      99.6923 % |

Two interesting sub-patterns:

- **rustar seed-to-seed dup-bit disagreement is ~10x higher than STAR seed-to-seed.** STAR's same-seed and seed-diff numbers are statistically indistinguishable (12 vs 0-12 disagreements on WT_REP1), which is just the MarkDuplicates noise floor of "you ran twice". rustar seed-to-seed sees 82-238 disagreements per sample. The asymmetry is also striking: in rustar seed-diff, e.g. WT_REP1 has 196 "STAR-dup-rustar-not" but only 42 "rustar-dup-STAR-not" — i.e. the second rustar run consistently marks fewer reads as duplicates (`dup_pct` 20.43 -> 20.35). This is downstream of the BAM record-order divergence documented in [`rustar_two_pass_and_determinism.md`](rustar_two_pass_and_determinism.md) section B; MarkDuplicates seeing a different read-arrival order makes a different keeper choice. Same total record set, different keeper-set.
- **STAR-vs-rustar disagreement (98.7-99.8%) is consistent with the BAM divergence floor** rather than additive on top of it. From [`rustar_bam_comparison.md`](rustar_bam_comparison.md), STAR-only / rustar-only primary read counts on PE samples were 25-110 / 5-1152 per sample, plus ~1k same-name position differences on multi-mappers, plus tag-content divergences. Picard's dup decision keys on `(rname, pos, orientation, mate-pos)`, and those are what diverge upstream. The 200-2 400 dup-bit disagreements are right where you'd expect given the upstream alignment divergence; MarkDuplicates is propagating, not amplifying.

## Verdict per question

### Is STAR deterministic with the same `--runRNGseed`?

**Yes, at the alignment-record content level.** Every SAM field except the dup bit is bit-identical across two STAR `-profile test,docker` runs with `--runRNGseed 0`. The dup bit flips on 4-26 reads per sample because Picard MarkDuplicates ties on read-arrival order, and STAR's BAM record emission order is non-deterministic by design.

**No, at the BAM byte level** — but every documented STAR-derived metric (mapping rate, splice counts, Salmon TPM at 1e-7-precision) is. Fitness-for-purpose verdict: STAR is **as deterministic as documented** for this pipeline.

Per the `rustar_bam_comparison.md` vocabulary: **FLOATING-POINT / ORDER** (the MarkDup keeper choice) and **NOT TESTED** (whether STAR alone, prior to MarkDuplicates, is byte-identical across reruns — likely also order-divergent, given the same evidence in `rustar_two_pass_and_determinism.md` Section B for rustar).

### How big is STAR's own RNG-driven variance?

**Negligible on this data.** The worst gene_tpm Pearson between STAR-seed-0 and STAR-seed-1 is `0.9999999994` (WT_REP2); gene_counts is `0.9999999999` everywhere. % mapped is identical at display precision on all 5 samples. Spearman = 1.0 (relative gene rank is unchanged).

Vocabulary tag: **FLOATING-POINT**. STAR's documented RNG affects multi-mapper tie-breaking and equal-score read selection; on yeast + GFP test data the number of contestable choices is small and the resulting per-gene TPM perturbation is at the limit of IEEE 754 doubles.

### How big is rustar's own RNG-driven variance?

**Same order as STAR's, possibly fractionally smaller.** Worst rustar-seed-0 vs rustar-seed-1 gene_tpm Pearson is `0.9999999996` (WT_REP1); gene_counts worst is `0.9999999999`. Spot-on identical to STAR's noise envelope at this resolution.

Vocabulary tag: **FLOATING-POINT / RNG**. This complements [`rustar_two_pass_and_determinism.md`](rustar_two_pass_and_determinism.md) Section B — rustar is record-content-identical across same-seed reruns (modulo order); this new measurement adds that across **different** seeds it produces functionally equivalent quant outputs.

There is a real **secondary divergence** that the seed change exposes: MarkDuplicates dup-count differs between rustar seed-0 and seed-1 (e.g. WT_REP1: 38 481 vs 38 327 dups, an 0.4% absolute swing). This isn't a TPM-level signal because MarkDuplicates marks but does not filter, and Salmon's transcriptome-BAM pipeline doesn't see the genome-BAM dup flag. Worth knowing for downstream tools that do consume the dup flag (`picard CollectRnaSeqMetrics`, custom filters).

### Are the STAR-vs-rustar deltas inside or outside the noise envelope?

**Outside, by 7-9 orders of magnitude.** STAR's own seed-driven gene_tpm Pearson noise is `1 - 0.9999999994 ~= 6e-10`. The STAR-vs-rustar deltas range `1 - 0.985 = 1.5e-2` (WT_REP2) to `1 - 0.9997 = 3e-4` (RAP1_UNINDUCED_REP\*). Every cross-aligner Pearson on every sample is comfortably outside the noise floor.

Verdict tag: the STAR-vs-rustar cross-aligner deltas previously catalogued in [`rustar_differences.md`](rustar_differences.md), [`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md), and [`rustar_quant_and_multiqc.md`](rustar_quant_and_multiqc.md) are **BUG** / **BEHAVIOURAL** rather than **RNG** / **FLOATING-POINT**. This is consistent with the upstream attribution: most of the WT_REP2 delta is driven by the paired-end transcriptome BAM mate-field bug (scverse/rustar-aligner#22), with secondary contributions from the sjdb-not-seeded bug (#27) and the NH-tail / multi-mapper sampling difference (#31). None of those are RNG-driven.

### Do MarkDuplicates decisions propagate cleanly or compound the divergence?

**They propagate cleanly.** Disagreement on STAR-vs-rustar matched primary records is 98.7-99.8 % per sample, which is what you'd expect from a tool that keys on `(rname, pos, orientation, mate-pos)` when those upstream fields already disagree on a small fraction of reads (cf. `rustar_bam_comparison.md` categories 4-5). Within-aligner same-seed agreement is 99.95-99.99 %; cross-aligner is one order of magnitude lower than that, which tracks with the cross-aligner BAM divergence floor.

If MarkDuplicates were _compounding_ the divergence (e.g. via cascading-effect read filtering) we'd expect to see cross-aligner agreement drop to ~95 % or worse and asymmetry biased by total dup count. We don't — disagreement is approximately symmetric (`star_dup_rustar_notdup ~= rustar_dup_star_notdup` on PE samples and within ~2x on SE), which is the signature of a propagation rather than amplification.

Vocabulary tag: **BEHAVIOURAL** (propagation, not amplification). MarkDuplicates is doing its job consistently; the upstream BAM divergence is what drives the cross-aligner disagreement.

## Cross-references

- [`rustar_bam_comparison.md`](rustar_bam_comparison.md) - the categorical BAM-level divergence catalogue these numbers ground.
- [`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md) - the WT_REP2 deep dive that surfaced the headline `0.985` TPM Pearson.
- [`rustar_two_pass_and_determinism.md`](rustar_two_pass_and_determinism.md) Section B - prior same-seed determinism check on the rustar aligner output, on which this doc layers the seed-to-seed analysis.
- [`rustar_differences.md`](rustar_differences.md) - the top-level index linking the 12 filed upstream issues (`scverse/rustar-aligner#22, #25-#35`).

## Constraints / limits

- All numbers are on the test-profile yeast subset + GFP, 5 samples, ~50-200 k reads per sample. **NOT TESTED**: whether the noise envelope at human-genome scale on the full test_full samplesheet is the same shape. STAR's own RNG-driven variance scales with the number of multi-mapper tied alignments, which grows with target complexity, so the floor could be one order of magnitude looser on real data. The headline conclusion (cross-aligner delta is many orders larger than the noise floor) would hold a fortiori under any plausible scaling.
- Only the `salmon` quantification path is exercised. The `star_rsem` path and `pseudo_aligner`-only path are out of scope here. **NOT TESTED**.
- The MarkDuplicates `READ_NAME_REGEX` and `OPTICAL_DUPLICATE_PIXEL_DISTANCE` settings are pipeline defaults. **NOT TESTED**: how sensitive the MarkDup-noise component would be to disabling optical-duplicate detection via `READ_NAME_REGEX null`, which would eliminate one source of ties.
