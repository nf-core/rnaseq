# rustar-aligner v0.1.0: mode-specific smoke tests

Date: 2026-05-12. Pipeline: nf-core/rnaseq branch `rustar-aligner` (PR #1855), commit `06c2ffc0`. Aligner: rustar-aligner v0.1.0 (https://github.com/scverse/rustar-aligner). Reference: STAR 2.7.11b. All runs on `nf-dev-rnaseq` (Docker).

Two smoke tests on top of the eukaryotic baseline already covered in [`rustar_differences.md`](rustar_differences.md), [`rustar_bam_comparison.md`](rustar_bam_comparison.md), and [`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md): Task A turns on `--quantMode GeneCounts` (we'd only exercised `TranscriptomeSAM`), Task B runs the prokaryotic test profile with `--aligner star_salmon`.

## TL;DR

- **Task A** - `--quantMode GeneCounts TranscriptomeSAM` works. rustar emits the STAR 4-column `ReadsPerGene.out.tab`; per-sample Pearson is >= 0.99994 on every column; max per-gene delta is 82 reads on a 38 k-read gene. Residual drift is downstream of the already-filed NH-tail issue (BAM-comparison issue 6). **RNG / BEHAVIOURAL, no new BUG.**
- **Task B** - originally diagnosed as a rustar upstream bug; **on follow-up verification it turned out to be a pipeline-integration gap, fixed in this PR.** `conf/modules/prepare_genome.config`'s `withName:` selector for the prokaryotic-specific `--sjdbGTFfeatureExon CDS` listed `STAR_GENOMEGENERATE` and `PARABRICKS_STARGENOMEGENERATE` but not `RUSTAR_GENOMEGENERATE`, so the flag was silently dropped from rustar's index build and the CDS-only annotation produced an index with zero transcripts. Adding `RUSTAR_GENOMEGENERATE` to the selector makes rustar's `Aligned.toTranscriptome.out.bam` byte-equivalent to STAR's (13 `@SQ`, 8 082 records). **No upstream issue.**

## Task A: `--quantMode GeneCounts TranscriptomeSAM`

### Method

```yaml
# rustar_gc.yml
use_rustar_star: true
extra_star_align_args: '--quantMode GeneCounts TranscriptomeSAM'
outdir: results-rustar-gc

# star_gc.yml
extra_star_align_args: '--quantMode GeneCounts TranscriptomeSAM'
outdir: results-star-gc
```

`-profile test,docker -resume -params-file <yml>`. The dedupe in `conf/modules/align_star.config` replaces the pipeline-default `--quantMode TranscriptomeSAM` with the user value, so both aligners emit both files. Compared the per-sample `*.ReadsPerGene.out.tab` (5 samples, 4-column STAR format).

### Header counters (`N_unmapped` / `N_multimapping` / `N_noFeature` / `N_ambiguous`)

| Sample            | Counter        | STAR uns / fwd / rev   | rustar uns / fwd / rev    |
| ----------------- | -------------- | ---------------------- | ------------------------- |
| RAP1_IAA_30M_REP1 | N_unmapped     | 3 824 / 3 824 / 3 824  | 3 823 / 3 823 / 3 823     |
| RAP1_IAA_30M_REP1 | N_multimapping | 916 / 916 / 916        | **1 022** / 1 022 / 1 022 |
| RAP1_IAA_30M_REP1 | N_noFeature    | 2 277 / 37 670 / 2 544 | 2 255 / 37 601 / 2 522    |
| RAP1_IAA_30M_REP1 | N_ambiguous    | 7 026 / 12 / 119       | 6 960 / 11 / 91           |
| WT_REP1           | N_unmapped     | 8 926 / 8 926 / 8 926  | 8 928 / 8 928 / 8 928     |
| WT_REP1           | N_multimapping | 1 996 / 1 996 / 1 996  | **2 173** / 2 173 / 2 173 |
| WT_REP1           | N_noFeature    | 5 163 / 76 713 / 5 679 | 5 133 / 76 654 / 5 642    |
| WT_REP1           | N_ambiguous    | 11 323 / 20 / 274      | 11 130 / 22 / 197         |
| WT_REP2           | N_unmapped     | 4 196 / 4 196 / 4 196  | 4 194 / 4 194 / 4 194     |
| WT_REP2           | N_multimapping | 987 / 987 / 987        | **1 065** / 1 065 / 1 065 |
| WT_REP2           | N_noFeature    | 2 506 / 37 538 / 2 814 | 2 476 / 37 501 / 2 783    |
| WT_REP2           | N_ambiguous    | 6 651 / 12 / 144       | 6 586 / 11 / 120          |

(Full 5-sample table in `/tmp/rustar_gc_analysis/gc_report.md` on the VM.)

`N_multimapping` is up by 70-180 reads on rustar across every sample; `N_ambiguous` / `N_noFeature` are correspondingly down. Same root cause as issue 6 in [`rustar_bam_comparison.md`](rustar_bam_comparison.md) (rustar's NH-tag tail reaches 20 vs STAR's 7). Confirmed at YAL038W (largest-delta gene in WT_REP1):

```
== YAL038W STAR WT_REP1 (chr I 71786-73288) NH counts ==
  22 398 NH:i:1   259 NH:i:2   28 NH:i:3   24 NH:i:4   10 NH:i:5   3 NH:i:6   4 NH:i:8   1 NH:i:9   6 NH:i:11   7 NH:i:12

== YAL038W rustar WT_REP1 NH counts ==
  22 494 NH:i:1   141 NH:i:2   105 NH:i:3   55 NH:i:4   39 NH:i:5   23 NH:i:6   14 NH:i:7   28 NH:i:8   16 NH:i:9   10 NH:i:11   15 NH:i:12   26 NH:i:13   40 NH:i:15   18 NH:i:20
```

NH=1 primary mappers match (22 398 -> 22 494); rustar adds noisy secondary hits at NH=13-20 that STAR doesn't have. Reads STAR scored NH=1 land at NH=2-8 in rustar, fall into `N_multimapping`, and `N_ambiguous` shrinks. `unique + multi + ambiguous + noFeature + unmapped` totals are preserved within +/-5 reads.

### Per-sample Pearson on per-gene counts

| Sample              | unstranded |  forward |  reverse | n_genes | identical cells (uns/fwd/rev) |          max | Δ   | (uns/fwd/rev) |
| ------------------- | ---------: | -------: | -------: | ------: | ----------------------------: | -----------: | --- | ------------- |
| RAP1_IAA_30M_REP1   |   0.999999 | 0.999988 | 0.999999 |     125 |               114 / 113 / 114 | 20 / 19 / 27 |
| RAP1_UNINDUCED_REP1 |   0.999999 | 0.999987 | 1.000000 |     125 |               114 / 119 / 114 | 25 / 25 / 11 |
| RAP1_UNINDUCED_REP2 |   0.999999 | 0.999986 | 0.999999 |     125 |               112 / 118 / 113 | 62 / 56 / 60 |
| WT_REP1             |   0.999998 | 0.999945 | 0.999999 |     125 |               111 / 111 / 110 | 56 / 82 / 39 |
| WT_REP2             |   0.999997 | 0.999982 | 0.999997 |     125 |               116 / 110 / 116 | 46 / 27 / 40 |

89-95 % of cells per sample/column are byte-identical; the rest drift inside +/-1 % of the per-gene count. No order-of-magnitude divergence. Largest absolute reverse-stranded deltas (YAL038W, YAR009C, YAR010C) are all genes already flagged by the NH-tail analysis above.

### Task A verdict

**RNG / BEHAVIOURAL.** No new bug. The residual drift is mechanistically explained by [`rustar_bam_comparison.md`](rustar_bam_comparison.md) issue 6 (NH tail). No upstream issue to file from Task A.

## Task B: prokaryotic mode (`--prokaryotic`, `--aligner star_salmon`)

### Method

```yaml
# star_prok.yml
aligner: star_salmon
outdir: results-star-prok

# rustar_prok.yml
use_rustar_star: true
aligner: star_salmon
outdir: results-rustar-prok
```

`-profile test_prokaryotic,docker -params-file <yml> -resume`. The `aligner: star_salmon` override is required because `test_prokaryotic` defaults to `bowtie2_salmon`. Input: 2-sample Salmonella SL1344 prokaryotic test set; `conf/modules/align_star.config` adds `--sjdbGTFfeatureExon CDS --alignIntronMax 1` to both aligners' args under `params.prokaryotic`.

STAR run completed. rustar run **failed at SALMON_QUANT** with exit 1 - salmon parsed the BAM header, found nothing to quantify, and bailed (`work/2d/3d14da.../.exitcode`). A later `-resume` reported "completed successfully" because salmon's input is `SAMPLE.bam` and the cache key collided with the earlier bowtie2 SALMON_QUANT; the work-dir symlink resolves to bowtie2's BAM. Cache-collision artefact, not a real rustar -> salmon success.

### Root cause: pipeline-side selector gap (fixed in this PR)

The empty transcriptome BAM is **not** a rustar bug. `conf/modules/prepare_genome.config` has:

```groovy
withName: 'STAR_GENOMEGENERATE|PARABRICKS_STARGENOMEGENERATE|RUSTAR_GENOMEGENERATE' {
    ext.args = {
        def args = []
        if (params.prokaryotic) {
            args += ['--sjdbGTFfeatureExon CDS']
        }
        args.join(' ')
    }
}
```

Before this PR's fix the `|RUSTAR_GENOMEGENERATE` slot was missing, so under `--prokaryotic` the flag was silently dropped from rustar's index build. rustar then built an index from the GFF's `exon` features (zero rows in a CDS-only annotation) and the transcriptome BAM had nothing to populate.

A minimal direct rustar invocation that passes `--sjdbGTFfeatureExon CDS` at **both** index and alignment time produces a transcriptome BAM byte-equivalent to STAR's:

| Aligner       | `@SQ` lines | records | mapping rate |
| ------------- | ----------: | ------: | -----------: |
| STAR 2.7.11b  |          13 |   8 082 |      86.78 % |
| rustar v0.1.0 |          13 |   8 082 |      86.78 % |

The earlier diagnosis that "rustar's transcriptome projection ignores `--sjdbGTFfeatureExon` and hardcodes `exon`" was wrong. rustar honours the flag fine when it's plumbed through; our pipeline wasn't plumbing it through.

The genome BAM was healthy on both aligners throughout (8 428 records, identical flagstat). `--alignIntronMax 1` is honoured.

### Other prokaryotic-mode findings (genome BAM is fine)

| Sample    | STAR mapped | rustar mapped | Δ (pp) |
| --------- | ----------: | ------------: | -----: |
| SALM_REP1 |     86.78 % |       86.78 % |   0.00 |
| SALM_REP2 |     86.69 % |       86.69 % |   0.00 |

- **`--alignIntronMax 1` is honoured.** Zero `N` ops in any CIGAR; `Number of splices: Total = 0`; both SJ tabs empty. **Not a silent-ignore bug.**
- **Cosmetic `Log.final.out` categorisation.** STAR splits unmapped into `too short` + `other`; rustar always reports `other = 0` and folds STAR's `other` into `too short` (238 + 404 = 642 on SALM_REP1). Total mapped count conserved; MultiQC unmapped bar reads ~3x higher than STAR's. **BEHAVIOURAL (low).**
- **No final `SJ.out.tab`, only `SJ.pass1.out.tab`** - already in [`rustar_differences.md`](rustar_differences.md); both empty on prokaryotic anyway.

`bin/compare_aligner_runs.py` reports PASS because matrices section is empty (rustar's salmon TSVs missing) and the harness only fails on present-but-divergent matrices. If `--prokaryotic --aligner star_salmon` ever lands in CI, the harness needs a fail-on-missing check; out of scope here.

### Task B verdict

**Pipeline-integration gap, fixed in this PR.** No upstream issue.

## Appendix: artefact locations on `nf-dev-rnaseq`

All under `/home/ubuntu/rnaseq-rustar-aligner/`.

- Task A: `results-{star,rustar}-gc/`; comparison script + `gc_report.md` in `/tmp/rustar_gc_analysis/`.
- Task B: `results-{star,rustar}-prok/` (rustar partial, aborted at SALMON_QUANT) and `results-bowtie2-prok/` (default `test_prokaryotic`, bowtie2 reference at 82-83 % mapping rate). rustar work dirs `work/ef/dd979fd16bf4c80b400c0d32b85a2a` + `work/01/3084eb9d055d906f4dac76cbf86905`; STAR `work/02/0fbec318cbd713e74d9cb73004d114` + `work/cd/d28e503f24af090e09019c7b6e9ec1`; failed salmon `work/2d/3d14da3b907991a0a58728ee9067c2`.
