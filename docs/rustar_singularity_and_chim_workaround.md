# rustar-aligner: Singularity/Apptainer parity and `--chimSegmentMin` workaround

Two verifications for PR nf-core/rnaseq#1855 at HEAD `9f87cfa9a` on `nf-dev-rnaseq` (Ubuntu 22.04, Apptainer 1.5.0). Index: [`rustar_differences.md`](rustar_differences.md). The BAM divergences and the chim-path finding referenced below sit in [`rustar_bam_comparison.md`](rustar_bam_comparison.md) and [`rustar_cli_compat.md`](rustar_cli_compat.md) ([scverse/rustar-aligner#35](https://github.com/scverse/rustar-aligner/issues/35) upstream).

## A. Singularity / Apptainer container parity

VM had neither `singularity` nor `apptainer`; installed Apptainer 1.5.0 from the official PPA (`add-apt-repository ppa:apptainer/ppa && apt-get install apptainer`); the apptainer package symlinks `singularity` -> `apptainer`, so the pipeline's `singularity` profile works unchanged. Setup: worktree `/home/ubuntu/rnaseq-rustar-sing`, `-profile test,singularity`, fresh `NXF_SINGULARITY_CACHEDIR`, same `rustar.params.yml` as the Docker baseline at `results-rustar/`. Pass criteria: TPM Pearson >= 0.999, % mapped within +-0.5 pp. `compare_aligner_runs.py results-rustar results-rustar-sing`:

| Check                                              | Verdict              | Detail                                                                                                                                                                                       |
| -------------------------------------------------- | -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Image pull via `singularity_pull_docker_container` | IDENTICAL            | 33 MB SIF, single `singularity pull docker://...` invocation, no manual conversion or auth                                                                                                   |
| `rustar-aligner --version` inside SIF vs Docker    | IDENTICAL            | Both `0.1.0` / `70be24d` / built `2026-05-12T15:14:41Z`                                                                                                                                      |
| % uniquely mapped, all 5 samples                   | IDENTICAL            | 0.00 pp delta everywhere (RAP1_IAA_30M_REP1: 90.23%, RAP1_UNINDUCED_REP1: 95.88%, RAP1_UNINDUCED_REP2: 95.80%, WT_REP1: 88.81%, WT_REP2: 89.39%)                                             |
| `gene_tpm` Pearson                                 | IDENTICAL            | Min 0.99999999968 (WT_REP1), 1.0 on three samples, Spearman 1.0 on all                                                                                                                       |
| `gene_counts` Pearson                              | IDENTICAL            | 1.000000 on all five                                                                                                                                                                         |
| `samtools flagstat` on markdup BAMs                | IDENTICAL            | Matching primary/secondary/mapped/duplicate on WT_REP1 (180596 / 7726 / 100.00% / 38481) and WT_REP2 (90716 / 3742 / 100.00% / 12183)                                                        |
| BAM byte equality                                  | FLOATING-POINT/ORDER | 50-260 B size delta per sample, consistent with the rustar threading non-determinism in [`rustar_two_pass_and_determinism.md`](rustar_two_pass_and_determinism.md); not Singularity-specific |

Output JSON `/tmp/rustar_sing_vs_docker.json`. Overall: **IDENTICAL**.

### Gotchas for HPC users

- **No SIF is published for rustar.** `singularity pull docker://ghcr.io/scverse/rustar-aligner:dev` works, but HPC sites that block `docker://` URLs on shared nodes will need to pre-pull on a head node and stage the SIF into `NXF_SINGULARITY_CACHEDIR`. Worth a line in `docs/usage.md`.
- **The `:dev` tag floats.** Mid-test, the Singularity pull fetched `70be24d` while a stale Docker image on the same host still resolved to `5f8ad08`; a `docker pull` brought them into sync. Users should pin to a digest once rustar cuts a tagged release - otherwise the same `:dev` reference produces different outputs on consecutive runs.
- **No writable-filesystem or `--nv` issues.** `singularity.autoMounts = true` is sufficient; no extra `runOptions`. rustar is CPU-only so GPU container flags are not in scope.

## B. `--chimSegmentMin > 0` workaround verification ([scverse/rustar-aligner#35](https://github.com/scverse/rustar-aligner/issues/35))

The chim path-builder appends `/Chimeric.out.junction` to `--outFileNamePrefix` unconditionally; pipeline-style `SAMPLE.` crashes with `No such file or directory`. Documented workaround in [`rustar_cli_compat.md`](rustar_cli_compat.md) Detail B: `--outFileNamePrefix dir/` with parent pre-created. Confirms the workaround end-to-end on PE (cli_compat probed SE) and re-checks the dot-prefix regression. Image `70be24d`, sample `WT_REP2`, existing index, `genome_gfp.gtf`.

### Test 1: directory-style prefix (the workaround), paired-end

`--outFileNamePrefix WT_REP2/ --chimSegmentMin 12 --chimOutType Junctions --twopassMode Basic --outSAMtype BAM Unsorted --sjdbGTFfile genome_gfp.gtf`

| Check                                     | Result                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ----------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Exit status                               | 0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `WT_REP2/Aligned.out.bam`                 | 6 682 971 B, non-empty                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `WT_REP2/Chimeric.out.junction`           | exists; 0 records on this PE yeast subset (Log.final.out: `Number of chimeric reads = 0`). Re-ran SE-only on the same `WT_REP2_primary_1.fastq.gz` for chim-data sanity: 5393 records, 14 columns, header sample `I 160892 + I 108830 - 0 0 0 SRR6357072.6691718 160799 7S94M 108810 70S21M10S` - matches STAR's 14-column format (chr_donor, brkpoint_donor, strand_donor, chr_acceptor, brkpoint_acceptor, strand_acceptor, junction_type, repeat_left, repeat_right, read_name, segment1_start, segment1_cigar, segment2_start, segment2_cigar) |
| `WT_REP2/Log.final.out`                   | 1984 B, 89.39% uniquely mapped (matches the non-chim runs)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `WT_REP2/SJ.out.tab` + `SJ.pass1.out.tab` | both present                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |

The workaround produces every standard rustar output; the chim file is well-formed when chimeric reads exist (cli_compat used SE to surface them since the PE yeast subset has none).

### Test 2: trailing-dot prefix (regression check)

Same invocation, `--outFileNamePrefix WT_REP2.` instead of `WT_REP2/`. As predicted:

```
[INFO  rustar_aligner] Chimeric detection enabled (chimSegmentMin=12)
Error: I/O error: No such file or directory (os error 2) (WT_REP2./Chimeric.out.junction)

Caused by:
    No such file or directory (os error 2)
```

Exit 1, no `WT_REP2.*` files produced. Regression for [#35](https://github.com/scverse/rustar-aligner/issues/35) is still present at `70be24d`.

### Verdict: workaround **works**, but is not pipeline-compatible as-is

The directory-prefix shape is the only rustar invocation that survives `--chimSegmentMin > 0`. A user adding `--chimSegmentMin 12` to `extra_star_align_args` today hits the regression: `modules/local/rustar_align/align/main.nf` builds `--outFileNamePrefix ${prefix}.` to match STAR, with a post-step `mv` flattening the resulting directory. Switching to `${prefix}/` unblocks chim but diverges further from STAR on a flag whose upstream fix is one line (`create_dir_all` before opening the chim writer; see [#35](https://github.com/scverse/rustar-aligner/issues/35)).

Recommended: wait for [#35](https://github.com/scverse/rustar-aligner/issues/35) upstream; no chim-specific module hack in this PR. If a user reports the failure first, suggest `extra_star_align_args = '--chimSegmentMin 12 --chimOutType WithinBAM'` (chim records in the genome BAM, sidestepping the path-builder) - `WithinBAM` semantics in rustar are not yet verified, flag as a follow-up probe.

## Files of interest

- Singularity worktree on VM: `/home/ubuntu/rnaseq-rustar-sing/` (log `rustar_sing.log`, results `results-rustar-sing/`).
- `compare_aligner_runs.py` output: `/tmp/rustar_sing_vs_docker.{json,md}` on the VM.
- Chim probe artefacts: `/home/ubuntu/rustar_chim_workaround/{WT_REP2/, WT_REP2_SE/, chim_dir.log, chim_dir_se.log, chim_dot.log}`.
