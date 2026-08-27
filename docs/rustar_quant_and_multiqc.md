# rustar-aligner v0.1.0 vs STAR: per-transcript quant and MultiQC report

Date: 2026-05-12. Pipeline: nf-core/rnaseq PR #1855 (`rustar-aligner`). Both
sections operate on the existing `results-{star,rustar}/` artefacts on
`nf-dev-rnaseq`; no re-run. Vocabulary follows
[`rustar_bam_comparison.md`](rustar_bam_comparison.md): **BUG / RNG /
FLOATING-POINT / BEHAVIOURAL / NOT TESTED**.

## A. Per-transcript Salmon `quant.sf` divergence on single-end samples

### Setup

Both `RAP1_UNINDUCED_REP1` and `RAP1_UNINDUCED_REP2` are single-end and
quantified against the same 125-transcript yeast+GFP transcriptome (one
isoform per gene; gene-level rollup hides nothing in this profile).
Source files:

| Sample              | STAR `quant.sf`                                         | rustar `quant.sf`                                         |
| ------------------- | ------------------------------------------------------- | --------------------------------------------------------- |
| RAP1_UNINDUCED_REP1 | `results-star/star_salmon/RAP1_UNINDUCED_REP1/quant.sf` | `results-rustar/star_salmon/RAP1_UNINDUCED_REP1/quant.sf` |
| RAP1_UNINDUCED_REP2 | `results-star/star_salmon/RAP1_UNINDUCED_REP2/quant.sf` | `results-rustar/star_salmon/RAP1_UNINDUCED_REP2/quant.sf` |

Driver script: `/tmp/rustar_taskA/analyze.py` (not committed). Per-transcript
deltas computed for `Length`, `EffectiveLength`, `TPM`, `NumReads`.

### Headline numbers

| Metric                                         | REP1                             | REP2                               |
| ---------------------------------------------- | -------------------------------- | ---------------------------------- |
| Transcripts (both runs)                        | 125                              | 125                                |
| Pearson r, `Length`                            | identical                        | identical                          |
| Pearson r, `EffectiveLength`                   | **1.000000**                     | **1.000000**                       |
| Pearson r, `NumReads`                          | 0.999904                         | 0.999906                           |
| Pearson r, `TPM`                               | 0.999673                         | 0.999746                           |
| Sum `NumReads` STAR -> rustar (delta)          | 40 158 -> 39 537 (-621, -1.55 %) | 80 630 -> 79 471 (-1 159, -1.44 %) |
| Transcripts with non-zero \|ΔEffectiveLength\| | **0 / 125**                      | **0 / 125**                        |
| Transcripts with non-zero ΔLength              | 0 / 125                          | 0 / 125                            |

**`EffectiveLength` matches to 6 decimal places on every transcript.** The
PE mate-fields bug
([`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md))
does not bleed into SE samples and produces no new EL shrinkage here.
This is the key negative finding: **no new BUG** in the SE quant path.

### Top 10 by |ΔTPM|, RAP1_UNINDUCED_REP1

| Transcript | Gene name |  Len |   TPM STAR | TPM rustar |      ΔTPM |  rel% | NR STAR | NR rustar |  ΔNR |
| ---------- | --------- | ---: | ---------: | ---------: | --------: | ----: | ------: | --------: | ---: |
| YAL003W    | EFB1      |  621 |  79 842.26 |  72 347.32 | -7 494.94 |  -9.4 |     809 |       721 |  -88 |
| YAR010C    | YAR010C   | 1323 | 206 210.65 | 210 318.41 | +4 107.77 |  +2.0 |    6043 |      6062 |  +19 |
| YAL038W    | CDC19     | 1503 | 186 084.50 | 188 572.57 | +2 488.07 |  +1.3 |    6368 |      6347 |  -21 |
| YAR009C    | YAR009C   | 3591 | 185 145.95 | 182 715.54 | -2 430.41 |  -1.3 |   16894 |     16398 | -496 |
| YAL005C    | SSA1      | 1929 | 132 742.75 | 134 585.75 | +1 843.00 |  +1.4 |    6087 |      6070 |  -17 |
| snR18      | snR18     |  102 |  82 673.20 |  84 055.78 | +1 382.59 |  +1.7 |       8 |         8 |    0 |
| YAL030W    | SNC1      |  354 |   4 576.80 |   3 579.49 |   -997.31 | -21.8 |      13 |        10 |   -3 |
| YAL069W    | YAL069W   |  315 |     560.51 |       0.00 |   -560.51 |  -100 |       1 |         0 |   -1 |
| YAL012W    | CYS3      | 1185 |   9 046.04 |   9 197.32 |   +151.28 |  +1.7 |     231 |       231 |    0 |
| YAR002C-A  | ERP1      |  660 |   6 876.46 |   6 991.46 |   +115.00 |  +1.7 |      77 |        77 |    0 |

(All single-isoform genes - no multi-isoform parent in test profile.)

### Top relative ΔTPM (require STAR TPM > 0.5 to drop noise)

REP1 outliers above the +1.67 % baseline (see "Two-band structure"):

| Transcript | Gene    |  Len |   TPM STAR | TPM rustar |        rel% | NR STAR | NR rustar |
| ---------- | ------- | ---: | ---------: | ---------: | ----------: | ------: | --------: |
| YAL069W    | YAL069W |  315 |     560.51 |       0.00 | **-100.00** |       1 |         0 |
| YAL063C    | FLO9    | 3969 |      96.00 |     119.27 |      +24.23 |     9.8 |      11.9 |
| YAL030W    | SNC1    |  354 |    4576.80 |    3579.49 |      -21.79 |      13 |        10 |
| YAL062W    | GDH3    | 1374 |     586.36 |     529.93 |       -9.62 |      18 |        16 |
| YAL003W    | EFB1    |  621 |  79 842.26 |  72 347.32 |       -9.39 |     809 |       721 |
| YAL011W    | SWC3    | 1878 |     539.78 |     525.94 |       -2.56 |      24 |        23 |
| YAR010C    | YAR010C | 1323 | 206 210.65 | 210 318.41 |       +1.99 |    6043 |      6062 |

REP2 has the same shape: `EFB1 -144 reads / -7.1 % TPM`, `SNC1 -11
reads / -26.4 % TPM`, `YAL003W EFB1` again the largest absolute. Zero
multi-isoform genes either side - gene-level rollup hides no isoform
artefact.

### Two-band structure and length banding (REP1)

Almost every transcript sits in one of two bands: a "global scale-shift"
band at +1.67 % (90+ of 125 transcripts), driven by rustar's ~1.5 %
total-read deficit dropping the TPM denominator uniformly; and a
sparse set of off-band outliers where reads were lost faster than the
global rate. Every outlier traces to known issues 3 (sjdb-not-seeded:
`EFB1`, `YAR009C`) or 6 (multi-mapper sampling: `SNC1`, `GDH3`,
`YAL069W`) in `rustar_bam_comparison.md`.

| Length band |   n | median \|ΔTPM\| | max \|ΔTPM\| | median rel TPM % | median \|ΔNR\| |
| ----------- | --: | --------------: | -----------: | ---------------: | -------------: |
| 0-200       |   9 |            0.00 |      1 382.6 |             1.67 |              0 |
| 200-500     |  33 |            0.00 |        997.3 |             1.67 |              0 |
| 500-1 000   |  31 |           14.83 |      7 494.9 |             1.67 |              0 |
| 1 000-2 000 |  28 |           13.84 |      4 107.8 |             1.67 |              0 |
| 2 000-5 000 |  24 |           19.38 |      2 430.4 |             1.67 |              0 |

No length-band concentration: the +1.67 % constant lives in every band
and the outliers are scattered. No EL shrinkage analogue of the PE
mate-fields effect.

### Verdict for Task A

`Length`, `EffectiveLength`: **clean** (identical to 6 dp; the PE
fragment-length-prior bug has no SE analogue). `NumReads`: ~1.5 %
systematic under-count, distributed non-uniformly - downstream of
known issues 3 and 6 in `rustar_bam_comparison.md`. `TPM`: mostly the
renormalisation echo of NumReads, with `EFB1` -9 % and `YAL069W`
-100 % as worst outliers. **No new BUG file-able from the SE quant
inspection.**

## B. MultiQC report content comparison

### Inventory

Both runs emit identical module sets in `multiqc_report_data/`:
`bbmap`, `Cutadapt`, `FastQC (raw|trimmed|filtered)`, `Picard`,
`QualiMap`, `RSeQC`, `STAR`, `Salmon`, `Samtools`. No module absent on
either side. Driver: `/tmp/rustar_taskB_compare.py` (not committed).

Per-module deltas roll up as follows (sample shown is WT_REP2 unless
noted; all paired-end samples show the same pattern, magnitudes vary).

### Per-module summary

| MultiQC module                | Status                          | Worst-case user-visible misread                                                                                                                                  |
| ----------------------------- | ------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| STAR Alignment Plot           | **MATERIALLY-DIFFERENT**        | Annotated splices = 0 (vs 235-1 152); non-canonical splice rate inflated 1.6-4x; mapping-rate ±0.2 pp; unmapped reclassified entirely from "other" -> "tooshort" |
| Samtools Stats                | **MATERIALLY-DIFFERENT**        | `error_rate` reads as 0.000000 (vs ~0.01) and `mismatches` reads as 0; `average_quality` reads ~2x (68.3 vs 35.3)                                                |
| Samtools Flagstat             | numerically-close-but-different | Total record count inflated +0.5-+0.7 % from extra secondary alignments (16-19 % more)                                                                           |
| Samtools Idxstats             | numerically-close-but-different | Per-chromosome read totals shift with secondary-alignment count                                                                                                  |
| Salmon (MultiQC's parse)      | numerically-close-but-different | frag_length_mean from flenDist matches STAR; **but Salmon's internal effective-length used the 250-bp default - hidden from MultiQC**                            |
| Picard MarkDuplicates         | numerically-close-but-different | %duplication shifts -0.02 to +0.13 pp; `SECONDARY_OR_SUPPLEMENTARY_RDS` +17 %                                                                                    |
| RSeQC infer_experiment        | same / numerically-close        | strand fractions stable within ±0.005; XS-tag absence does NOT break this tool                                                                                   |
| RSeQC bam_stat                | **MATERIALLY-DIFFERENT**        | `splice_reads` halved (~ -50 to -53 %); `proper_pairs_percent` drops by ~0.8 pp; `non_primary_hits` inflated 16-19 %                                             |
| RSeQC read_distribution       | numerically-close               | Intron/intergenic tag distribution shifts up to ±11 %, mostly tiny numerators (6-17 introns total)                                                               |
| RSeQC junction_annotation     | **MATERIALLY-DIFFERENT**        | `known_splicing_events_pct` 46-55 % -> 15-27 %; novel splice rate is inversely inflated                                                                          |
| RSeQC junction_saturation     | numerically-different           | Subsample curves swing ±50 % - artefact of tiny denominators (3-4 known junctions in test data)                                                                  |
| RSeQC inner_distance          | same / RNG                      | Per-bin inner-distance histogram shifts by ≤1 %; reads-paired shifts ±0.07 %                                                                                     |
| Qualimap RNA-seq              | numerically-close-but-different | `non_unique_alignments` +14-15 % from extra multi-mappers; `5_3_bias` -1.7-1.9 %                                                                                 |
| dupRadar                      | numerically-close               | Per-bin model shifts; bin x-coords move by 1-3 % from upstream gene-quant changes                                                                                |
| featureCounts biotype         | same / numerically-close        | Counts shift by < 0.2 %                                                                                                                                          |
| FastQC (raw/trimmed/filtered) | **identical**                   | (parses input FASTQs only, untouched by aligner choice)                                                                                                          |
| cutadapt                      | **identical**                   | Same input FASTQs                                                                                                                                                |
| BBMap (bbsplit)               | **identical**                   | Pre-alignment                                                                                                                                                    |
| Salmon DESeq2 PCA (gene)      | numerically-close               | PC1/PC2 shift by 0.06-1.6 % - same sample order, same cluster shape                                                                                              |
| STAR_SALMON DESeq2 PCA        | **materially-different**        | PC2 swing on RAP1_UNINDUCED_REP2 78 % and on WT_REP1 ~6 800 % (sign flip near zero) - cosmetic, axes are noisy near origin                                       |

### Quantitative drill-down on the modules that misread

#### STAR Alignment Plot (parses `Log.final.out`)

| Sample              | Metric                   |  STAR | rustar | Delta   |
| ------------------- | ------------------------ | ----: | -----: | ------- |
| WT_REP1             | num_annotated_splices    | 1 152 |      0 | -100 %  |
| WT_REP1             | num_splices (total)      | 1 406 |    636 | -54.8 % |
| WT_REP1             | num_noncanonical_splices |    68 |    149 | +119 %  |
| WT_REP1             | unmapped_other           | 5 421 |      0 | -100 %  |
| WT_REP1             | unmapped_tooshort        | 3 504 |  8 927 | +154 %  |
| RAP1_IAA_30M_REP1   | num_annotated_splices    |   528 |      0 | -100 %  |
| RAP1_UNINDUCED_REP1 | num_annotated_splices    |   235 |      0 | -100 %  |

`num_annotated_splices = 0` is the loudest signal in the MultiQC report
and is the direct downstream of **issue 3 (sjdb not seeded)**. The
~5 000-read reclassification `unmapped_other -> unmapped_tooshort` is
rustar's internal categorisation; totals match but the stacked bar
chart shifts visually. **BEHAVIOURAL + downstream of BUG**.

#### Samtools Stats (parses BAM `NM:i` and quality strings)

| Sample  | Metric          |             STAR |          rustar | Note                                                                                                                                                                                                                                                 |
| ------- | --------------- | ---------------: | --------------: | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| All     | mismatches      | 54 926 - 173 449 |           **0** | samtools stats reads `NM:i:`; rustar emits `nM`                                                                                                                                                                                                      |
| All     | error_rate      |  0.0091 - 0.0117 |         **0.0** | derived from `mismatches / bases_mapped_(cigar)`                                                                                                                                                                                                     |
| All     | average_quality |      35.3 - 35.5 | **68.3 - 68.5** | exactly +33 offset; root cause confirmed by the verification session as [scverse/rustar-aligner#34](https://github.com/scverse/rustar-aligner/issues/34) — rustar writes Phred+33 ASCII bytes into the BAM `QUAL` column instead of raw Phred values |
| WT_REP2 | reads_MQ0       |               18 |              54 | +200 %; downstream of category-7 multi-mapper redistribution                                                                                                                                                                                         |

`error_rate: 0` is the most actively misleading number in the report.
**A user looking at MultiQC will conclude rustar BAMs are error-free.**
This is the user-visible projection of **issue 1 (NM->nM)** in
`rustar_bam_comparison.md` (filed upstream as [#29](https://github.com/scverse/rustar-aligner/issues/29)).
The `average_quality` doubling — flagged here as "new, separate symptom"
when this doc was written — was diagnosed later as the BAM `QUAL` column
being written as Phred+33 ASCII rather than raw Phred values, filed as
[#34](https://github.com/scverse/rustar-aligner/issues/34). 68 − 33 = 35,
which is exactly the offset.

#### Salmon (MultiQC parses `flenDist.txt`, not `meta_info.json`)

| Sample               | frag_length_mean STAR | frag_length_mean rustar | Delta     |
| -------------------- | --------------------: | ----------------------: | --------- |
| WT_REP2              |                165.86 |                  165.88 | +0.014 %  |
| WT_REP1              |                165.65 |                  165.65 | 0.000 %   |
| RAP1_IAA_30M_REP1    |                167.41 |                  167.46 | +0.029 %  |
| RAP1_UNINDUCED_REP\* |     250 (default, SE) |       250 (default, SE) | identical |

**Most actively misleading entry in the report.** MultiQC's
`frag_length_mean` comes from the empirical distribution Salmon writes
to `flenDist.txt` after quantification (Salmon can infer it
post-hoc by re-pairing in the transcriptome BAM at QC time, regardless
of what it used internally). The PE rustar flenDist agrees with STAR
within 0.03 %, **but the per-transcript `EffectiveLength` in
`quant.sf` was computed using the 250 bp default** (see
`rustar_investigation_wt_rep2.md`). A user reading MultiQC sees
"rustar fragment length looks identical to STAR" and concludes TPMs
are reliable; they are not. DOCUMENTATION GAP - MultiQC has no
surface for `meta_info.json::frag_length_mean` (the prior Salmon
actually used).

#### RSeQC bam_stat / junction_annotation (parses BAM)

| Sample  | Metric                    |  STAR | rustar | Delta       |
| ------- | ------------------------- | ----: | -----: | ----------- |
| WT_REP1 | splice_reads              | 1 103 |    522 | -52.7 %     |
| WT_REP1 | proper_pairs_percent      | 74.19 |  73.59 | -0.82 pp    |
| WT_REP1 | non_primary_hits          | 6 508 |  7 726 | +18.7 %     |
| WT_REP1 | known_splicing_events_pct |  45.6 |   20.2 | -55.8 % rel |
| WT_REP2 | known_splicing_events_pct |  52.6 |   21.3 | -59.5 % rel |

All downstream of **issue 3 (sjdb)** and **issue 6 (extra secondary
alignments)** in `rustar_bam_comparison.md`. No new BUG; existing
upstream bugs surface here as user-readable headline numbers.

#### Picard MarkDuplicates, Qualimap RNA-seq

Picard %duplication shifts -0.02 to +0.13 pp;
`SECONDARY_OR_SUPPLEMENTARY_RDS` +17 % mirrors flagstat. Qualimap
`non_unique_alignments +14 %`, `5_3_bias` -1.7 to -4 %. All
BEHAVIOURAL, all downstream of issue 6 (extra multi-mappers) - no new
bug, but the +14 % Qualimap multi-map rate is the kind of number a
user might flag as "unusual sample quality".

#### RSeQC infer_experiment

**Correctness preserved despite the missing `XS` tag.** infer_experiment
uses the BAM strand bit (`0x10`) plus the BED/GTF reference; it does
not consume `XS:A:`. Strand fractions in MultiQC match STAR within
±0.005 on every sample. Tools that **do** need `XS` (StringTie,
Cufflinks) remain at risk per issue 2 in `rustar_bam_comparison.md`,
but that risk is not reflected in MultiQC.

### Headline answer

**Running `--use_rustar_star` today, the misread MultiQC numbers, by
visibility:**

1. **`error_rate: 0` and `mismatches: 0` (Samtools Stats)** -
   conclude rustar BAMs are pristine. They are not; samtools-stats
   indexes off `NM:i:`. (Issue 1, NM->nM.)
2. **`average_quality = 68` vs STAR's 35 (Samtools Stats)** - a
   consistent doubling across all five samples. New symptom on top of
   the NM bug; root cause not yet localised.
3. **`Annotated (sjdb) = 0` in STAR Alignment Plot** - looks like no
   GTF was provided. It was. (Issue 3.)
4. **RSeQC `known_splicing_events_pct` halved**, `splice_reads`
   halved. Same root cause as (3).
5. **`non_primary_hits` and `non_unique_alignments` inflated 14-19 %**
   in flagstat / bam_stat / Qualimap - would prompt a user to flag
   the sample for unusual multi-mapping. Real, behavioural (issue 6).
6. **`frag_length_mean` (Salmon) looks normal but is a lie about what
   Salmon used internally.** MultiQC has no surface for the real
   prior; user has no signal that PE TPMs are distorted.

Every misread is the user-visible projection of a bug already
documented in `rustar_bam_comparison.md` /
`rustar_investigation_wt_rep2.md`. **No new BUG file-able from
MultiQC alone.** The report does give materially weaker QC for
rustar than for STAR until issues 1 and 3 are fixed upstream.

### What I couldn't measure

- **Unmapped category re-routing** in `Log.final.out` -
  `unmapped_other` -> 0 vs `unmapped_tooshort` +154 %: totals match,
  internal categorisation differs. Joint-walking unmapped FASTQs is
  needed to characterise; same blocker as
  `rustar_bam_comparison.md` category 9. **NOT TESTED.**
- **`average_quality` doubling root cause**: worth a one-line
  reproducer against samtools 1.21 outside the pipeline once rustar
  emits `NM`.
