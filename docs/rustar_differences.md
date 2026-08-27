# rustar-aligner integration: observed differences vs STAR

Running notes captured while wiring up the experimental `--use_rustar_star`
path in PR #1855. The intent here is to track every divergence we observe so
that nothing surprises us at review and so we can file targeted upstream
issues at https://github.com/scverse/rustar-aligner (or fix things on our
side) rather than discovering them in production. This will be cleaned up
before merge.

The verification setup: standard `-profile test,docker` on the
`nf-dev-rnaseq` VM (36 CPU / 69 GB), back-to-back STAR and rustar runs,
identical inputs.

**2026-08-27 re-verification**: every issue filed below has since been
closed upstream (rustar-aligner moved 0.1.0 → 0.2.0 over the same period).
[`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md)
re-checks each one directly against the current image rather than trusting
the "closed" label - most are genuinely fixed (including the headline #22
TPM bug), two are improved but short of full parity, and one (#31, NH-tail
depth) shows no measurable change despite being closed with no linked PR.
The tables below are left as the historical record of what the original
investigation found; the "Tracked upstream" table at the bottom is updated
with current status and links to the re-verification evidence.

## Verified

### Wall-time and RAM (test profile, one tile of yeast + GFP)

From `pipeline_info/execution_trace_*.txt`, comparing the per-task medians:

| Process              | n | Wall (s) STAR → rustar | Peak RSS (GB) STAR → rustar |
|----------------------|---|------------------------|------------------------------|
| `STAR_GENOMEGENERATE` / `RUSTAR_GENOMEGENERATE` | 1 | 0.3 → 0.3 | 0.01 → 0.02 |
| `STAR_ALIGN` / `RUSTAR_ALIGN`                   | 5 | 68.0 → 33.8 | 0.92 → 0.12 |

Caveat: this is on the tiny test genome (a yeast subset plus GFP transgene)
with ≤10 k reads per sample, run inside Docker. The absolute numbers say
nothing about human-scale performance. Re-running on the `test_full`
samplesheet on AWS is a follow-up.

### Mapping rate (per `Log.final.out`)

Original (v0.1.0, 2026-05-12):

| Sample              | STAR  | rustar | Δ (pp) |
|---------------------|-------|--------|--------|
| RAP1_IAA_30M_REP1   | 90.44 | 90.23  | -0.21  |
| RAP1_UNINDUCED_REP1 | 95.96 | 95.88  | -0.08  |
| RAP1_UNINDUCED_REP2 | 95.85 | 95.80  | -0.05  |
| WT_REP1             | 88.99 | 88.81  | -0.18  |
| WT_REP2             | 89.54 | 89.39  | -0.15  |

Re-verified (v0.2.0, 2026-08-27):

| Sample              | STAR  | rustar | Δ (pp) |
|---------------------|-------|--------|--------|
| RAP1_IAA_30M_REP1   | 90.44 | 90.24  | -0.20  |
| RAP1_UNINDUCED_REP1 | 95.96 | 95.91  | -0.05  |
| RAP1_UNINDUCED_REP2 | 95.85 | 95.83  | -0.02  |
| WT_REP1             | 88.99 | 88.82  | -0.17  |
| WT_REP2             | 89.54 | 89.40  | -0.14  |

Unchanged, within ±0.25 pp of STAR both times.

### Quantification concordance (per-sample Pearson on merged Salmon matrices)

Original (v0.1.0, 2026-05-12):

| Sample              | gene_tpm | gene_counts |
|---------------------|----------|-------------|
| RAP1_IAA_30M_REP1   | 0.996808 | 0.999848    |
| RAP1_UNINDUCED_REP1 | 0.999673 | 0.999904    |
| RAP1_UNINDUCED_REP2 | 0.999746 | 0.999906    |
| WT_REP1             | 0.995496 | 0.999890    |
| WT_REP2             | **0.985040** | 0.999842 |

`gene_counts` (raw `NumReads`) was essentially identical across both
runs even then. `gene_tpm` diverged materially on `WT_REP2`, with
`RAP1_IAA_30M_REP1` and `WT_REP1` showing the same effect at smaller
magnitude - all three are paired-end. The two single-end samples
(`RAP1_UNINDUCED_REP1/2`) were already clean, which was the tell: see
[`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md) for
the root-cause deep dive (rustar v0.1.0's `Aligned.toTranscriptome.out.bam`
didn't populate mate-pair fields or the proper-pair flag on paired-end
records, so Salmon fell back to its default fragment-length prior,
distorting `EffectiveLength`/TPM for short transcripts). Filed as
[scverse/rustar-aligner#22](https://github.com/scverse/rustar-aligner/issues/22).

**Fixed as of v0.2.0** (re-verified 2026-08-27, see
[`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md)):

| Sample              | gene_tpm | gene_counts |
|---------------------|----------|-------------|
| RAP1_IAA_30M_REP1   | 0.999841 | 0.999852    |
| RAP1_UNINDUCED_REP1 | 0.999914 | 0.999912    |
| RAP1_UNINDUCED_REP2 | 0.999916 | 0.999911    |
| WT_REP1             | 0.999853 | 0.999890    |
| WT_REP2             | 0.999839 | 0.999858    |

All five samples now clear 0.9998 on both matrices. The paired-end
transcriptome BAM records carry proper `RNEXT`/`PNEXT`/`TLEN` and the
proper-pair flag; confirmed directly with `samtools view`, not just
inferred from the TPM improvement.

## Module-level workarounds

**Update 2026-08-27: all three workarounds below have been removed.**
rustar-aligner 0.2.0 fixed the upstream behaviour each one compensated
for (see [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md)
for the direct verification of each). Removing them was confirmed
behaviour-neutral: `bin/compare_aligner_runs.py` produced identical
pass/fail results and Pearson correlations before and after removal.
Kept here as a historical record of what v0.1.0 needed and why.

### `--limitGenomeGenerateRAM` was not accepted (v0.1.0)

STAR exposes `--limitGenomeGenerateRAM`; the upstream `STAR_GENOMEGENERATE`
module derives a value from `task.memory` and passes it. rustar v0.1.0
rejected this flag at startup (`error: unexpected argument
'--limitGenomeGenerateRAM' found`), so `modules/local/rustar_align/genomegenerate/main.nf`
omitted it. Filed as [scverse/rustar-aligner#25](https://github.com/scverse/rustar-aligner/issues/25).
**Fixed in 0.2.0** - the flag is now accepted (`--help` lists a
`31G` default) and is restored in the module, computed identically to
`STAR_GENOMEGENERATE` (`task.memory.toBytes() - 100000000`).

### `--outFileNamePrefix` ending in `.` was treated as a directory (v0.1.0)

STAR treats `--outFileNamePrefix SAMPLE.` as a literal string prefix and
writes `SAMPLE.Aligned.out.bam`, `SAMPLE.Log.final.out`, etc. side by
side in the work directory. rustar v0.1.0 instead interpreted the same
value as a directory name and wrote bare-named files inside it, so
`modules/local/rustar_align/align/main.nf` post-processed by flattening
that directory back into STAR-style filenames. Filed as
[scverse/rustar-aligner#26](https://github.com/scverse/rustar-aligner/issues/26).
**Fixed in 0.2.0** (via upstream PR #46) - verified directly with a
standalone CLI invocation; the flattening step is removed.

### `Log.out` and `Log.progress.out` were not written (v0.1.0)

STAR emits three log files: `Log.final.out` (summary stats, MultiQC
input), `Log.out` (verbose run log) and `Log.progress.out` (per-chunk
progress). rustar v0.1.0 only wrote `Log.final.out`, so both were
marked `optional: true` in `RUSTAR_ALIGN`. Filed as
[scverse/rustar-aligner#28](https://github.com/scverse/rustar-aligner/issues/28)
(folded into [#55](https://github.com/scverse/rustar-aligner/issues/55),
fixed via upstream PR #82). **Fixed in 0.2.0** - both files now contain
real STAR-equivalent content (parameter dump, per-phase progress,
`ALL DONE!`), not stubs; `optional: true` removed from the module.

### Extra `SJ.pass1.out.tab` at the top level (v0.1.0)

rustar wrote both `SJ.out.tab` and `SJ.pass1.out.tab` (the two-pass
intermediate) at the top level. STAR keeps the intermediate inside
`<prefix>_STARpass1/`. Covered by the same #28/#55 fix above -
**fixed in 0.2.0**: the pass-1 table now lives in `<prefix>_STARpass1/`
like STAR's, and the top-level duplicate is gone (confirmed by the
nf-test snapshot diff losing those entries entirely on re-generation).

### Version reporting

The rustar container (`ghcr.io/scverse/rustar-aligner:dev` on debian-slim)
does not bundle `samtools` or `gawk`, which are present in the STAR Wave
container. STAR_GENOMEGENERATE uses `samtools faidx` + `gawk` to
auto-compute `--genomeSAindexNbases`.

To avoid adding a `samtools`/`gawk` dependency to the rustar image,
`RUSTAR_GENOMEGENERATE` does the same heuristic in Groovy from the
on-disk FASTA size. The approximation is well inside the floor() of
`log2(len)/2 - 1` so the chosen index size matches.

`RUSTAR_ALIGN` emits only the `rustar-aligner` version through the
topic-based versions channel - no `samtools` / `gawk` emissions.

## Nextflow-side, not rustar's fault, but bites us anyway

### Boolean CLI flags get coerced to the string `"true"`

`--use_rustar_star`, `--use_rustar_star=true`, and
`--use_rustar_star true` all fail nf-schema validation with `Value is
[string] but should be [boolean]` on Nextflow 26.04 + nf-schema 2.6.1.
This is not rustar-specific; the same error occurs for
`--use_parabricks_star`. A YAML params file works:

```yaml
use_rustar_star: true
outdir: results-rustar
```

then `nextflow run ... -params-file rustar.params.yml`. Worth raising
upstream (Nextflow / nf-schema), separately from rustar.

## Still to verify

- Full-size run on the `test_full` samplesheet (GRCh37, larger reads) to
  produce performance and concordance numbers that map to user
  expectations. The test-profile numbers above are not load-bearing.
- Whether rustar's `--quantTranscriptomeSAMoutput BanSingleEnd` matches
  STAR's interpretation byte-for-byte.

## Tracked upstream

All filed against [scverse/rustar-aligner](https://github.com/scverse/rustar-aligner/issues).
Every row below is now closed upstream; the **Status** column and
[`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) record what re-checking
each one directly on rustar-aligner 0.2.0 actually found - "closed" isn't taken at face value.

| # | Severity | Summary | Status (2026-08-27) | Evidence |
|---|---|---|---|---|
| [#22](https://github.com/scverse/rustar-aligner/issues/22) | high | Paired-end transcriptome BAM omits mate fields (`RNEXT`/`PNEXT`/`TLEN`) + proper-pair flag, Salmon falls back to its default fragment-length prior and distorts TPM. | ✅ Fixed, verified | [`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md), [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#25](https://github.com/scverse/rustar-aligner/issues/25) | medium | `--limitGenomeGenerateRAM` rejected by the CLI parser. | ✅ Fixed, verified | [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#26](https://github.com/scverse/rustar-aligner/issues/26) | medium | `--outFileNamePrefix SAMPLE.` treated as a directory rather than a string prefix. | ✅ Fixed (upstream PR #46), verified | [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#27](https://github.com/scverse/rustar-aligner/issues/27) | medium | `Log.final.out` always reports `Annotated (sjdb) = 0` despite `--sjdbGTFfile`; ~50% of splices missing. Root cause: `is_annotated()` coord-space bug at `src/align/stitch.rs:1306-1314`. | 🟡 Improved, not at parity | [`rustar_two_pass_and_determinism.md`](rustar_two_pass_and_determinism.md), [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#28](https://github.com/scverse/rustar-aligner/issues/28) | low | Output-shape gaps: `Log.out` / `Log.progress.out` not written, `SJ.pass1.out.tab` lives at the top level instead of under `<prefix>_STARpass1/`. | ✅ Fixed (via #55 → upstream PR #82), verified | [`rustar_bam_comparison.md`](rustar_bam_comparison.md), [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#29](https://github.com/scverse/rustar-aligner/issues/29) | high | `--outSAMattributes NM` emits `nM:i:` instead of `NM:i:`, with different semantics (substitutions only, no indels). Breaks samtools stats, Picard, MultiQC. | ✅ Fixed, verified (indel semantics checked directly) | [`rustar_bam_comparison.md`](rustar_bam_comparison.md), [`rustar_quant_and_multiqc.md`](rustar_quant_and_multiqc.md), [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#30](https://github.com/scverse/rustar-aligner/issues/30) | high | `--outSAMstrandField intronMotif` accepted but no `XS:A:` tags ever emitted. Breaks StringTie, Cufflinks. (RSeQC `infer_experiment` uses the BAM strand bit instead so is unaffected.) | 🟡 Improved, not at parity (64% vs STAR's 100% coverage on spliced reads) | [`rustar_bam_comparison.md`](rustar_bam_comparison.md), [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#31](https://github.com/scverse/rustar-aligner/issues/31) | medium | Multi-mapper NH cap extends to 20 vs STAR's 7; ~17% more secondaries on identical input. Possibly missing an `--outFilterMultimapScoreRange`-equivalent threshold. | 🔴 Closed with no linked PR; **no measurable change** on re-test | [`rustar_bam_comparison.md`](rustar_bam_comparison.md), [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#32](https://github.com/scverse/rustar-aligner/issues/32) | low | Transcriptome BAM lacks per-record `RG:Z:` despite the `@RG` header being present. Genome BAM is fine. | ✅ Fixed, verified | [`rustar_bam_comparison.md`](rustar_bam_comparison.md), [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#33](https://github.com/scverse/rustar-aligner/issues/33) | low | `@PG` header is content-free (just `ID:rustar-aligner`, no `PN`/`VN`/`CL`); `AS:i:` values disagree by 2-5 units on 864 records with identical CIGAR. | ✅ Fixed (header), verified | [`rustar_bam_comparison.md`](rustar_bam_comparison.md), [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#34](https://github.com/scverse/rustar-aligner/issues/34) | high | BAM `QUAL` field is offset by +33 (Phred+33 ASCII bytes written instead of raw Phred values). Explains the "average_quality = 68 vs STAR's 35" symptom in MultiQC; spotted by the verification session, not our own audits. Highest-impact BAM-correctness issue after #22 because every downstream tool that reads QUAL is wrong. | ✅ Fixed, verified (`average_quality` now 35.5 for both) | [`rustar_quant_and_multiqc.md`](rustar_quant_and_multiqc.md), [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#35](https://github.com/scverse/rustar-aligner/issues/35) | medium | `--chimSegmentMin > 0` + `--twopassMode Basic` aborts the run when `--outFileNamePrefix` doesn't end in `/`. Silent run-killer: no `Aligned.out.bam`, no `Log.final.out`. | ✅ Fixed, verified | [`rustar_cli_compat.md`](rustar_cli_compat.md), [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#47](https://github.com/scverse/rustar-aligner/issues/47) / [#50](https://github.com/scverse/rustar-aligner/issues/50) | medium | Follow-ups on #27: pass-1 doesn't seed candidates from `--sjdbGTFfile`; `sjdb_score` bonus added instead of replacing motif score. | 🟡 Improved, not at parity (same gap as #27) | [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#48](https://github.com/scverse/rustar-aligner/issues/48) | medium | `Log.final.out` folds all unmapped reads into `too short`; `other` bucket always 0. | 🔴 Closed with no observed fix; **no measurable change** on re-test | [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#53](https://github.com/scverse/rustar-aligner/issues/53) | low | `Number of reads mapped to too many loci` always 0; reads exceeding `--outFilterMultimapNmax` dropped from accounting. | ⚪ Untestable on the yeast test set (nothing exceeds the threshold either way) | [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |
| [#55](https://github.com/scverse/rustar-aligner/issues/55) | low | `Log.out`/`Log.progress.out` need real STAR-equivalent content, not stubs (a rejected earlier PR attempted a stub fix). | ✅ Fixed (upstream PR #82), verified | [`rustar_reverification_2026-08.md`](rustar_reverification_2026-08.md) |

## Fixed in this PR (was originally suspected upstream)

- **Prokaryotic mode + rustar produced an empty transcriptome BAM**. `conf/modules/prepare_genome.config`'s `withName:` selector for `--sjdbGTFfeatureExon CDS` listed STAR + Parabricks but not `RUSTAR_GENOMEGENERATE`, so the flag was silently dropped from rustar's index build. Adding `RUSTAR_GENOMEGENERATE` to the selector makes rustar byte-equivalent to STAR on the same inputs (13 `@SQ`, 8 082 records). Originally diagnosed as a rustar transcriptome-projection bug; reclassified after the verification session showed rustar honours the flag fine when it's plumbed through. See [`rustar_mode_smoke_tests.md`](rustar_mode_smoke_tests.md).
