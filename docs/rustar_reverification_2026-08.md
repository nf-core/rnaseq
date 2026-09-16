# Re-verification against rustar-aligner 0.2.0 (2026-08-27)

All 17 issues filed against [scverse/rustar-aligner](https://github.com/scverse/rustar-aligner)
during the original PR #1855 investigation (`rustar_differences.md` et al., May 2026) are
closed. This doc re-verifies each one directly, on identical `-profile test,docker` inputs,
rather than trusting the tracker's "closed" label - per the standing rule for this integration,
the bar is methodological parity (same computation, same categorisation), not "issue says
fixed". Several issues turned out to be genuinely fixed; a few are improved but not at parity;
two (#31, #48) show no measurable change - traced to premature closes: the PRs that
closed them (#77, #141) fix a different, adjacent bug, not the one reported. See below.

Setup: `ghcr.io/scverse/rustar-aligner:dev`, image built 2026-08-23T22:29Z, `rustar-aligner
0.2.0` (git `207773b`), on `nf-dev-rnaseq` (36 CPU / 69 GB). STAR and rustar run back-to-back,
`-profile test,docker`, `-params-file {star,rustar}.params.yml` (`save_align_intermeds: true`,
boolean params must go through a params file - see the existing nf-schema note below). Compared
with `bin/compare_aligner_runs.py` plus direct `samtools view`/`stats` inspection using
`community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4`.

## Fully resolved (verified, not just closed)

| #                                                                                                                       | Original finding                                                                                            | Verification                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| ----------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [#22](https://github.com/scverse/rustar-aligner/issues/22)                                                              | Paired-end transcriptome BAM omitted `RNEXT`/`PNEXT`/`TLEN` and the proper-pair flag, distorting Salmon TPM | `WT_REP2.Aligned.toTranscriptome.out.bam` records now carry `RNEXT`/`PNEXT`/`TLEN` and flags `83`/`163` (proper-pair bit set). `gene_tpm` Pearson on WT_REP2: **0.985040 → 0.999839**. All five samples now ≥0.9998 on both `gene_tpm` and `gene_counts`.                                                                                                                                                                                             |
| [#25](https://github.com/scverse/rustar-aligner/issues/25)                                                              | `--limitGenomeGenerateRAM` rejected at startup                                                              | `--help` now lists it (`[default: 31G]`). Restored to `RUSTAR_GENOMEGENERATE` mirroring `STAR_GENOMEGENERATE`'s exact computation (`task.memory.toBytes() - 100000000`); ran clean in the full pipeline with `--limitGenomeGenerateRAM 16006127360` on the task's `.command.sh`.                                                                                                                                                                      |
| [#26](https://github.com/scverse/rustar-aligner/issues/26)                                                              | `--outFileNamePrefix SAMPLE.` treated as a directory                                                        | Direct CLI test: `--outFileNamePrefix TESTSAMPLE.` now produces flat `TESTSAMPLE.Aligned.out.bam` etc., no directory. Directory-flattening workaround removed from `RUSTAR_ALIGN`; a full pipeline run with the workaround gone produced no stray `SAMPLE.` directory anywhere in `work/`.                                                                                                                                                            |
| [#28](https://github.com/scverse/rustar-aligner/issues/28) / [#55](https://github.com/scverse/rustar-aligner/issues/55) | `Log.out`/`Log.progress.out` not written; `SJ.pass1.out.tab` at top level instead of `<prefix>_STARpass1/`  | Both files now written with real STAR-equivalent content (parameter dump, per-phase timing, `ALL DONE!` progress line) - not stubs. `optional: true` removed from both outputs in `RUSTAR_ALIGN`. `WT_REP2._STARpass1/` now exists as a real directory holding the pass-1 SJ table, matching STAR's layout; the top-level duplicate is gone (confirmed via the nf-test snapshot diff: the `SJ.pass1.out.tab` top-level entries disappeared entirely). |
| [#29](https://github.com/scverse/rustar-aligner/issues/29)                                                              | `--outSAMattributes NM` emitted `nM:i:` (substitutions only) instead of spec `NM:i:`                        | Genome BAM records now carry `NM:i:`. Verified semantics, not just presence: a record with CIGAR `69M2D32M` (a 2-base deletion, 0 mismatches, confirmed via `MD:Z:69^AT32`) reports `NM:i:2` - correct SAM-spec edit distance (mismatches + indel bases), matching STAR's definition.                                                                                                                                                                 |
| [#32](https://github.com/scverse/rustar-aligner/issues/32)                                                              | Transcriptome BAM missing per-record `RG:Z:`                                                                | Present on every checked transcriptome-BAM record.                                                                                                                                                                                                                                                                                                                                                                                                    |
| [#33](https://github.com/scverse/rustar-aligner/issues/33)                                                              | `@PG` header content-free                                                                                   | Now `ID:rustar-aligner PN:rustar-aligner VN:0.2.0 CL:<full command line>`, matching STAR's `@PG` shape field-for-field.                                                                                                                                                                                                                                                                                                                               |
| [#34](https://github.com/scverse/rustar-aligner/issues/34)                                                              | BAM `QUAL` bytes offset by +33 (Phred+33 ASCII instead of raw Phred)                                        | `samtools stats average_quality`: **35.5 for both aligners** (was 68.3 rustar vs 35.3 STAR). `error_rate` now 0.90% rustar vs 0.93% STAR (was reading 0 for rustar due to the compounding QUAL bug).                                                                                                                                                                                                                                                  |
| [#35](https://github.com/scverse/rustar-aligner/issues/35)                                                              | `--chimSegmentMin > 0` + `--twopassMode Basic` aborted when `--outFileNamePrefix` doesn't end in `/`        | Direct CLI test with `--outFileNamePrefix CHIMTEST.` (no trailing slash) + `--twopassMode Basic --chimSegmentMin 12`: completed cleanly, exit 0.                                                                                                                                                                                                                                                                                                      |

## Improved but not at parity - flagging honestly rather than closing the file

| #                                                                                                                                                                                    | Original finding                                                         | Where it stands now                                                                                                                                                                                                                                                                                                                                                                                                     |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [#27](https://github.com/scverse/rustar-aligner/issues/27) / [#47](https://github.com/scverse/rustar-aligner/issues/47) / [#50](https://github.com/scverse/rustar-aligner/issues/50) | `Annotated (sjdb)` always 0; ~50% of splices dropped from pass-1 seeding | Massively better - no longer 0 - but still short of STAR. WT_REP2: rustar reports 674 total splices / 388 annotated (57.6%); STAR reports 762 total / 644 annotated (84.5%). Total-splice recovery is 88% of STAR's count; annotated-fraction is 27pp below STAR's. Not a regression from "closed", but not parity either.                                                                                              |
| [#30](https://github.com/scverse/rustar-aligner/issues/30)                                                                                                                           | `--outSAMstrandField intronMotif` never emitted `XS:A:`                  | XS is now emitted, but not on every spliced record: 896/1393 (64%) of N-CIGAR reads in WT_REP2's genome BAM carry `XS:A:`, vs 1100/1100 (100%) for STAR on the same reads. The shortfall tracks the #27 gap above - reads whose splice rustar can't credit as annotated/canonical likely fall into `Non-canonical` (92 rustar vs 35 STAR) and lose the strand call along with it, rather than being an independent bug. |

## No measurable change - and now we know why: the closing PRs fixed a different bug

Both re-confirmed live against `ghcr.io/scverse/rustar-aligner:dev` pulled fresh on
2026-08-27T15:01Z - git `3713dfb`, built minutes before this check, the literal tip
of `main` at the time (top commit: "perf(build): enable geometric LCP memoization for
genome indexing (#228)", same day). This rules out image staleness entirely: every fix
PR below has been in every `:dev` build for weeks, and the symptom is unchanged.

| #                                                          | Original finding                                                                   | Why the closing PR doesn't fix it                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| ---------------------------------------------------------- | ---------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [#31](https://github.com/scverse/rustar-aligner/issues/31) | NH tail extends deeper than STAR's (20 vs 6-7)                                     | Closed citing PR #77 ("Integration/pr batch 1", a 14-PR batch merge). Read directly: #77 bundles PR #54, which fixes `too_many_loci` **accounting** (routing reads over `--outFilterMultimapNmax` to the right counter) - it does not touch the NH-depth cap itself. Re-ran on the fresh image: rustar still reaches `NH:20` on WT_REP2 (40 primary alignments at NH:20, entries at NH:11/15/20 that STAR never produces), STAR's max on the same input is `NH:7`. The underlying divergence #31 reported was never addressed by any PR in the batch - the issue was closed on the strength of adjacent fixes in the same integration PR, not this one.                                                                                                                                                                                    |
| [#48](https://github.com/scverse/rustar-aligner/issues/48) | `Log.final.out` folds all unmapped reads into `too short`; `other` bucket always 0 | Closed citing PR #141 ("check too many loci after the quality filter, matching STAR's"). Read directly: #141 fixes reads being mislabeled `too many loci` instead of `too short` when they fail the quality filter first (its own repro: 269 vs 26,285 in `too many loci` before the fix, on a chr21 human test) - a `too_many_loci` vs `too_short` bug. It says nothing about, and doesn't touch, the `too_short` vs `other` split that #48 actually reported. Re-ran on the fresh image: WT_REP2 still reports `too short: 4189 (8.45%)`, `other: 0 (0.00%)`; STAR on the identical input reports `too short: 1540 (3.11%)`, `other: 2656 (5.36%)`. Total unmapped-non-mismatch count matches closely (4189 vs 4196) so it's a mis-categorisation, not a counting bug - and it's a different mis-categorisation than the one #141 fixed. |

Both look like premature closes: the maintainer (or an automated linker) associated the
issue with a PR that fixed a neighbouring symptom in the same subsystem, not the one
reported. Recommend filing fresh, narrowly-scoped issues rather than reopening the old
ones, since the old repros (0.1.0-era) no longer match the current codebase's internals -
these fresh numbers are the right repro to attach.

## Untestable on this dataset

- [#53](https://github.com/scverse/rustar-aligner/issues/53) ("too many loci" always 0, reads exceeding `--outFilterMultimapNmax` dropped from accounting): both STAR and rustar report 0 reads mapped to too many loci on this sample. Nothing in the yeast test set exceeds `--outFilterMultimapNmax 20`, so this can't be exercised here either way. Not verified fixed or broken; would need a larger/repeat-rich reference to trigger.

## What this means for the module workarounds

All three module-level workarounds documented in the original `rustar_differences.md` are now
removed - upstream fixed the behaviour they were compensating for (see #25/#26/#28 above), and
removing them was verified behaviour-neutral end to end: `bin/compare_aligner_runs.py` produced
byte-identical pass/fail results and the same Pearson correlations with the workarounds present
vs removed. See `modules/local/rustar_align/{align,genomegenerate}/main.nf`.

## Fresh comparison numbers (test profile, docker, rustar-aligner 0.2.0)

### % Uniquely mapped reads

| Sample              | STAR  | rustar | Δ (pp) |
| ------------------- | ----- | ------ | ------ |
| RAP1_IAA_30M_REP1   | 90.44 | 90.24  | -0.20  |
| RAP1_UNINDUCED_REP1 | 95.96 | 95.91  | -0.05  |
| RAP1_UNINDUCED_REP2 | 95.85 | 95.83  | -0.02  |
| WT_REP1             | 88.99 | 88.82  | -0.17  |
| WT_REP2             | 89.54 | 89.40  | -0.14  |

### Salmon merged matrices (per-sample Pearson)

| Sample              | gene_tpm | gene_counts |
| ------------------- | -------- | ----------- |
| RAP1_IAA_30M_REP1   | 0.999841 | 0.999852    |
| RAP1_UNINDUCED_REP1 | 0.999914 | 0.999912    |
| RAP1_UNINDUCED_REP2 | 0.999916 | 0.999911    |
| WT_REP1             | 0.999853 | 0.999890    |
| WT_REP2             | 0.999839 | 0.999858    |

Every sample now clears the harness's 0.999 pass threshold on both matrices; the WT_REP2 outlier
that drove the original `#22` investigation is gone. `nf-test tests/rustar_default.nf.test
--profile=+test,docker` passes 4/4 against a regenerated snapshot (the snapshot changed:
`rustar-aligner` version string 0.1.0 → 0.2.0, `SJ.out.tab`/`SJ.pass1.out.tab`/featureCounts/
StringTie content hashes shifted downstream of the real alignment and splice-detection changes
above, and the duplicate top-level `SJ.pass1.out.tab` entries are gone per the #28/#55 fix).
