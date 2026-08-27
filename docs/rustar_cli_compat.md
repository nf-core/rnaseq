# rustar-aligner v0.1.0 vs STAR 2.7.11b: CLI flag compatibility matrix

Date: 2026-05-12. PR: nf-core/rnaseq#1855, branch `rustar-aligner`. Container under test: `ghcr.io/scverse/rustar-aligner:dev` (rustar-aligner 0.1.0, commit `5f8ad08`, built 2026-05-12T13:01:03Z). Reference: STAR 2.7.11b as shipped by the rnaseq pipeline. Existing characterisations referenced rather than re-derived: [`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md) (transcriptome BAM mate-pair bug), [`rustar_bam_comparison.md`](rustar_bam_comparison.md) (seven further BAM-level divergences), [`rustar_two_pass_and_determinism.md`](rustar_two_pass_and_determinism.md) (GTF-SJ-not-seeded root cause and determinism).

## TL;DR

- **STAR ships ~200 documented parameters**. rustar advertises ~95 in its `--help` ("the ~40 most important" per its preamble, the actual count is higher because each chim / SJfilter / score family expands to multiple flags). This audit compared 100 STAR flags that the pipeline plausibly cares about (every flag we actually pass, plus all defaults STAR ships with that have downstream consumers in nf-core/rnaseq).
- Of the 100 compared flags:
  - **53 accepted-honoured** (present in rustar `--help`, semantics match STAR within tested precision)
  - **3 accepted-different / partially-honoured** (advertised but with divergent default or only-partially-implemented effect — `alignSJDBoverhangMin`, `twopassMode`, `quantTranscriptomeSAMoutput`)
  - **40 rejected-at-startup** (`error: unexpected argument '--X' found`, exit 0 from clap's error handler but the pipeline call would fail with a non-zero alignment exit before any output is written)
  - **4 advertised-but-broken** (parsed correctly but trigger a downstream crash or produce a known-bad artefact: `chimSegmentMin` family hits the output-prefix path bug; `outSAMstrandField intronMotif` is accepted but the `XS` tag is never written; `outSAMattributes NM` produces `nM` instead; `quantTranscriptomeSAMoutput BanSingleEnd` honours the singleton-ban but the resulting paired records are missing mate fields)
- **No silent-ignore class exists at the CLI level.** rustar is built with `clap`, so any flag not in `--help` is rejected on parse with a clean usage message and the run exits non-zero. The "dangerous" category turns out to be the **advertised-but-broken** class above, all of which silently produce wrong output rather than failing fast — see the detail sections.
- **Top high-severity dangers** (in priority order, all already filed or characterised in the sibling docs): `outSAMstrandField intronMotif` accepted but never writes `XS`; `outSAMattributes NM` accepted but writes `nM` with wrong semantics; `quantTranscriptomeSAMoutput BanSingleEnd` accepted but the BAM is missing mate-pair bookkeeping. **One new finding** in this audit: `--chimSegmentMin > 0` works in principle (correct 14-column `Chimeric.out.junction`) but the rustar implementation tries to write the chim file inside a directory derived from `--outFileNamePrefix` without creating it, so any prefix that ends in `.` or any non-`/` character (which is what the pipeline produces) crashes the run before writing.

## Methodology

### Sources

- **STAR's authoritative parameter list**: [`source/parametersDefault`](https://github.com/alexdobin/STAR/blob/master/source/parametersDefault) on `master`. Each parameter is documented inline with its default and a short description; STARsolo parameters are after the regular alignment ones. 200 parameters in total, of which ~120 are alignment-relevant (the rest are STARsolo, transformation, or undocumented/deprecated). Cross-checked against the [STAR manual PDF](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for semantics.
- **rustar's advertised flags**: `docker run --rm ghcr.io/scverse/rustar-aligner:dev rustar-aligner --help`. Captured in full above the matrix; the canonical capture is also at `/tmp/rustar_help.txt` on the `nf-dev-rnaseq` VM.
- **Pipeline-passed flags**: `conf/modules/align_star.config` on the PR branch.

### Probe protocol

Three classes of probe were used:

1. **Help-only inspection** for flags whose default and description in `rustar --help` match STAR's `parametersDefault` line for the same flag. No execution required.
2. **CLI parse-only probe**: `docker run --rm ghcr.io/scverse/rustar-aligner:dev rustar-aligner --<flag> <value>` with no genome / reads. clap rejects unknown args before any work happens; this exhaustively answers "is the flag advertised?". The exit message is the canonical `error: unexpected argument '--X' found` and `Usage: rustar-aligner [OPTIONS]`. (Note: the wrapping shell may report exit 0 but the inner rustar exits non-zero; in a pipeline context this is a hard failure.)
3. **Behaviour probe**: a focused rustar invocation against the existing WT_REP2 work-dir (`/home/ubuntu/rnaseq-rustar-aligner/work/16/d9d802566f49cb84bfccd08ddef610/` on the VM). Each probe was a single rustar binary call with the flag under test set to a value that should produce a measurable, asymmetric effect; the resulting `Log.final.out`, `SJ.out.tab`, and `Chimeric.out.junction` were diffed against the default-run baseline. Output landed under `/tmp/rustar_probes_real/<probe_name>/` on the VM. Two full-pipeline reruns (`results-rustar-probe-scoremin`, `results-rustar-probe-sjdb`) were also kept in `/home/ubuntu/rnaseq-rustar-probe/` for end-to-end confirmation.

The behaviour probes (smallest first):

| Probe                  | Flag                            | Value                                             | Expected if honoured                              | Observed                                                                                                                                                                                                                                                                                             | Verdict                                                                                                                                                                                                                                          |
| ---------------------- | ------------------------------- | ------------------------------------------------- | ------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| scoremin               | `--outFilterScoreMinOverLread`  | 0.95 (vs default 0.66)                            | mapping rate drops                                | 89.39% -> 74.18%; too-short 8.46% -> 24.88%                                                                                                                                                                                                                                                          | **accepted-honoured**                                                                                                                                                                                                                            |
| sjdb_overhang          | `--alignSJDBoverhangMin`        | 30 (vs pipeline-set 1)                            | fewer annotated junctions used                    | identical to baseline (371 splices, sjdb=0)                                                                                                                                                                                                                                                          | **honoured-but-no-op** (no GTF SJs are seeded into pass 1 in the first place, so the threshold has nothing to gate; see `rustar_two_pass_and_determinism.md` Section A)                                                                          |
| 2pass_none             | `--twopassMode None`            | None (vs Basic)                                   | no `SJ.pass1.out.tab`, possibly fewer splices     | `SJ.pass1.out.tab` correctly absent; splice/mapping counts identical to Basic                                                                                                                                                                                                                        | **honoured but no behavioural impact on this dataset**; on yeast subset rustar's pass-1 finds no novel sjdb, so the two-pass step is effectively a no-op (real-genome behaviour untested but expected to differ once pass-1 discovers junctions) |
| score_range            | `--outFilterMultimapScoreRange` | 0 (vs default 1)                                  | fewer multimappers                                | 89.39% -> 90.00% unique, 2.15% -> 1.54% multi, 371 -> 542 splices                                                                                                                                                                                                                                    | **accepted-honoured** (flag does shift multi-vs-unique cutoff)                                                                                                                                                                                   |
| chim_default           | `--chimSegmentMin`              | 12, `--chimOutType Junctions`, no chim by default | `Chimeric.out.junction` with STAR's 14-col format | with `--outFileNamePrefix dir/` (trailing slash): produces correct 14-col file, 5393 chimeric records on SE WT_REP2 reads. With `--outFileNamePrefix <sample>.` (the pipeline's spelling): crashes with `Error: I/O error: No such file or directory (os error 2) (<sample>./Chimeric.out.junction)` | **accepted-honoured at data level, accepted-broken at path level** (see Detail B)                                                                                                                                                                |
| readMapNumber          | `--readMapNumber`               | 5000                                              | exactly 5000 reads aligned                        | exactly 5000 reads in Log.final.out                                                                                                                                                                                                                                                                  | **accepted-honoured**                                                                                                                                                                                                                            |
| limitGenomeGenerateRAM | n/a                             | rejected                                          | rejected at parse                                 | rejected with tip `--limitBAMsortRAM`                                                                                                                                                                                                                                                                | **rejected** (matches expectation)                                                                                                                                                                                                               |

For the BAM-tag-level questions (`outSAMstrandField intronMotif` accepted but XS missing, `outSAMattributes NM` accepted but emits `nM`, `quantTranscriptomeSAMoutput BanSingleEnd` accepted but mate-pair fields missing) the verdict is taken from [`rustar_bam_comparison.md`](rustar_bam_comparison.md) which already characterised these on the full test sample set; this audit cross-references rather than re-derives.

## Compatibility matrix

Columns: STAR flag (camelCase, as STAR/rustar both spell it); STAR default; rustar status (`accepted-honoured` | `accepted-different` | `accepted-broken` | `rejected` | `not-probed`); notes / link to detail subsection.

### Run / system parameters

| STAR flag           | STAR default | rustar status     | Notes                                                                                                                                                        |
| ------------------- | ------------ | ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--runMode`         | `alignReads` | accepted-honoured | rustar supports `alignReads` and `genomeGenerate`. STAR's `inputAlignmentsFromBAM`, `liftOver`, `soloCellFiltering` are not in rustar's help.                |
| `--runThreadN`      | 1            | accepted-honoured | Same default and semantics                                                                                                                                   |
| `--runDirPerm`      | `User_RWX`   | rejected          | Not in `--help`                                                                                                                                              |
| `--runRNGseed`      | 777          | accepted-honoured | Same default. RNG is ChaCha (vs STAR's Mersenne Twister), so seed-equality does not produce byte-equal tie-breaks; see `rustar_bam_comparison.md` category 7 |
| `--parametersFiles` | `-`          | rejected          | Not in `--help`                                                                                                                                              |
| `--sysShell`        | `-`          | rejected          | Not in `--help`                                                                                                                                              |

### Genome / index parameters

| STAR flag                     | STAR default     | rustar status     | Notes                                                                                      |
| ----------------------------- | ---------------- | ----------------- | ------------------------------------------------------------------------------------------ |
| `--genomeDir`                 | `./GenomeDir/`   | accepted-honoured | Same default (rustar prints `./GenomeDir`, STAR `./GenomeDir/`; equivalent)                |
| `--genomeLoad`                | `NoSharedMemory` | rejected          | Not in `--help`. Pipeline does not pass this; STAR's shared-memory modes are not supported |
| `--genomeFastaFiles`          | `-`              | accepted-honoured | Used at `--runMode genomeGenerate`                                                         |
| `--genomeSAindexNbases`       | 14               | accepted-honoured | Same default                                                                               |
| `--genomeChrBinNbits`         | 18               | accepted-honoured | Same default                                                                               |
| `--genomeSAsparseD`           | 1                | accepted-honoured | Same default                                                                               |
| `--genomeChainFiles`          | `-`              | rejected          | liftOver-only; not in `--help`                                                             |
| `--genomeTransformOutput`     | `None`           | rejected          | Transform-only; not in `--help`                                                            |
| `--genomeTransformType`       | `None`           | rejected          | Transform-only                                                                             |
| `--genomeTransformVCF`        | `-`              | rejected          | Transform-only                                                                             |
| `--genomeFileSizes`           | 0                | rejected          | Not in `--help`                                                                            |
| `--genomeSuffixLengthMax`     | -1               | rejected          | Not in `--help`                                                                            |
| `--genomeChrSetMitochondrial` | `chrM M MT`      | rejected          | STARsolo-only; not in `--help`                                                             |
| `--genomeType`                | `Full`           | rejected          | Under-development in STAR; not in `--help`                                                 |

### Splice junction / annotation parameters

| STAR flag                          | STAR default             | rustar status     | Notes                                                                                                                                                      |
| ---------------------------------- | ------------------------ | ----------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--sjdbGTFfile`                    | `-`                      | accepted-honoured | GTF is read for transcript/gene info and per-record GeneCounts. Annotation-derived junctions are _not_ seeded into pass 1 alignment scoring (see Detail A) |
| `--sjdbGTFchrPrefix`               | `-`                      | accepted-honoured | rustar default empty; STAR default `-` (same null semantics)                                                                                               |
| `--sjdbGTFfeatureExon`             | `exon`                   | accepted-honoured | Same default; pipeline overrides to `CDS` for `--prokaryotic`                                                                                              |
| `--sjdbGTFtagExonParentTranscript` | `transcript_id`          | accepted-honoured | Same default                                                                                                                                               |
| `--sjdbGTFtagExonParentGene`       | `gene_id`                | accepted-honoured | Same default                                                                                                                                               |
| `--sjdbGTFtagExonParentGeneName`   | `gene_name`              | rejected          | Not in `--help`                                                                                                                                            |
| `--sjdbGTFtagExonParentGeneType`   | `gene_type gene_biotype` | rejected          | Not in `--help`                                                                                                                                            |
| `--sjdbOverhang`                   | 100                      | accepted-honoured | Same default                                                                                                                                               |
| `--sjdbScore`                      | 2                        | accepted-honoured | Same default. Pipeline overrides to 1 for `star_rsem`. Has no observable effect because GTF SJs aren't seeded                                              |
| `--sjdbFileChrStartEnd`            | `-`                      | rejected          | Not in `--help`                                                                                                                                            |
| `--sjdbInsertSave`                 | `Basic`                  | rejected          | Not in `--help`                                                                                                                                            |

### Read / input parameters

| STAR flag                | STAR default  | rustar status     | Notes                                                       |
| ------------------------ | ------------- | ----------------- | ----------------------------------------------------------- |
| `--readFilesIn`          | `Read1 Read2` | accepted-honoured | Same semantics                                              |
| `--readFilesCommand`     | `-`           | accepted-honoured | Pipeline passes `zcat`                                      |
| `--readMapNumber`        | -1            | accepted-honoured | Probed                                                      |
| `--readFilesType`        | `Fastx`       | rejected          | Not in `--help`. rustar accepts only FASTQ/A (no SAM input) |
| `--readFilesSAMattrKeep` | `All`         | rejected          | SAM-input-only                                              |
| `--readFilesManifest`    | `-`           | rejected          | Not in `--help`                                             |
| `--readFilesPrefix`      | `-`           | rejected          | Not in `--help`                                             |
| `--readMatesLengthsIn`   | `NotEqual`    | rejected          | Not in `--help`                                             |
| `--readNameSeparator`    | `/`           | rejected          | Not in `--help`                                             |
| `--readQualityScoreBase` | 33            | rejected          | Not in `--help`                                             |
| `--readNameFilter`       | (n/a)         | accepted-honoured | rustar-specific debug flag, no STAR equivalent              |

### Read clipping

| STAR flag                    | STAR default | rustar status     | Notes                                                                                                          |
| ---------------------------- | ------------ | ----------------- | -------------------------------------------------------------------------------------------------------------- |
| `--clip5pNbases`             | 0            | accepted-honoured | Same default                                                                                                   |
| `--clip3pNbases`             | 0            | accepted-honoured | Same default                                                                                                   |
| `--clipAdapterType`          | `Hamming`    | rejected          | Not in `--help`. rustar does not implement adapter clipping; trim adapters upstream (Trim Galore in pipeline). |
| `--clip3pAdapterSeq`         | `-`          | rejected          | Not in `--help`                                                                                                |
| `--clip3pAdapterMMp`         | 0.1          | rejected          | Not in `--help`                                                                                                |
| `--clip3pAfterAdapterNbases` | 0            | rejected          | Not in `--help`                                                                                                |
| `--clip5pAdapterSeq`         | `-`          | rejected          | STAR-internal under-dev; not in `--help`                                                                       |
| `--clip5pAdapterMMp`         | 0.1          | rejected          | Same                                                                                                           |
| `--clip5pAfterAdapterNbases` | 0            | rejected          | Same                                                                                                           |

### Limits

| STAR flag                   | STAR default        | rustar status     | Notes                                                                                                                                  |
| --------------------------- | ------------------- | ----------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `--limitBAMsortRAM`         | 0                   | accepted-honoured | Same default. Only matters with `--outSAMtype BAM SortedByCoordinate` which the pipeline doesn't use                                   |
| `--limitGenomeGenerateRAM`  | 31000000000         | rejected          | Pipeline does not pass this currently; if a user adds it via `extra_star_align_args` the rustar run will fail with a clean parse error |
| `--limitIObufferSize`       | `30000000 50000000` | rejected          | Not in `--help`                                                                                                                        |
| `--limitOutSAMoneReadBytes` | 100000              | rejected          | Not in `--help`                                                                                                                        |
| `--limitOutSJoneRead`       | 1000                | rejected          | Not in `--help`                                                                                                                        |
| `--limitOutSJcollapsed`     | 1000000             | rejected          | Not in `--help`                                                                                                                        |
| `--limitSjdbInsertNsj`      | 1000000             | rejected          | Not in `--help`                                                                                                                        |
| `--limitNreadsSoft`         | -1                  | rejected          | Not in `--help`                                                                                                                        |

### Output: general

| STAR flag               | STAR default | rustar status      | Notes                                                                                                                                                                                                                                                                             |
| ----------------------- | ------------ | ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--outFileNamePrefix`   | `./`         | accepted-different | rustar treats a prefix that doesn't end in `/` as if it were a directory (e.g. `WT_REP2.` becomes `WT_REP2./`). The pipeline already works around this in `modules/local/rustar_align/align/main.nf` with a post-step `mv`. Triggers a hard crash for chimeric output (Detail B). |
| `--outTmpDir`           | `-`          | rejected           | Not in `--help`                                                                                                                                                                                                                                                                   |
| `--outTmpKeep`          | `None`       | rejected           | Not in `--help`                                                                                                                                                                                                                                                                   |
| `--outStd`              | `Log`        | accepted-honoured  | rustar adds `None` as default (no stdout routing); STAR sends Log to stdout by default. Pipeline does not pass this                                                                                                                                                               |
| `--outReadsUnmapped`    | `None`       | accepted-honoured  | rustar supports `None` / `Fastx`. Behaviour identical when pipeline passes `Fastx`                                                                                                                                                                                                |
| `--outQSconversionAdd`  | 0            | rejected           | Not in `--help`                                                                                                                                                                                                                                                                   |
| `--outMultimapperOrder` | `Old_2.4`    | rejected           | Not in `--help`. rustar's multi-mapper ordering is RNG-driven independently                                                                                                                                                                                                       |

### Output: SAM/BAM general

| STAR flag                          | STAR default   | rustar status       | Notes                                                                                                                                                                                                                                                                                                                                                                                                             |
| ---------------------------------- | -------------- | ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--outSAMtype`                     | `SAM`          | accepted-honoured   | rustar accepts `SAM`, `BAM Unsorted`, `BAM SortedByCoordinate`, `None`. Pipeline passes `BAM Unsorted`                                                                                                                                                                                                                                                                                                            |
| `--outSAMmode`                     | `Full`         | rejected            | Not in `--help`                                                                                                                                                                                                                                                                                                                                                                                                   |
| `--outSAMstrandField`              | `None`         | **accepted-broken** | Pipeline passes `intronMotif`. rustar accepts the flag (does not error) but never emits `XS:A:` tags on the output BAM. **HIGH** severity - filed as issue 2 in `rustar_bam_comparison.md`.                                                                                                                                                                                                                       |
| `--outSAMattributes`               | `Standard`     | **accepted-broken** | rustar's `Standard` is `NH HI AS nM` (note `nM`, not `NM`). Pipeline passes the explicit list `NH HI AS NM MD`. rustar accepts `NM` in the list but emits `nM` in the BAM with substitution-only semantics. STAR's `jM`/`jI`/`XS`/`MC`/`ch`/variation/STARsolo/`cN`/`MD` attrs all accepted-by-list but not all are emitted; `MD` is. **HIGH** severity for `NM`, filed as issue 1 in `rustar_bam_comparison.md`. |
| `--outSAMattrIHstart`              | 1              | rejected            | Not in `--help`                                                                                                                                                                                                                                                                                                                                                                                                   |
| `--outSAMunmapped`                 | `None`         | accepted-honoured   | Pipeline passes `Within` for the RSEM path. rustar default `None` matches STAR. `Within KeepPairs` not tested.                                                                                                                                                                                                                                                                                                    |
| `--outSAMorder`                    | `Paired`       | rejected            | Not in `--help`                                                                                                                                                                                                                                                                                                                                                                                                   |
| `--outSAMprimaryFlag`              | `OneBestScore` | rejected            | Not in `--help`. rustar always selects one best alignment as primary                                                                                                                                                                                                                                                                                                                                              |
| `--outSAMreadID`                   | `Standard`     | rejected            | Not in `--help`                                                                                                                                                                                                                                                                                                                                                                                                   |
| `--outSAMmapqUnique`               | 255            | accepted-honoured   | Same default                                                                                                                                                                                                                                                                                                                                                                                                      |
| `--outSAMflagOR`                   | 0              | rejected            | Not in `--help`                                                                                                                                                                                                                                                                                                                                                                                                   |
| `--outSAMflagAND`                  | 65535          | rejected            | Not in `--help`                                                                                                                                                                                                                                                                                                                                                                                                   |
| `--outSAMattrRGline`               | `-`            | accepted-honoured   | Same semantics. rustar emits `@RG` header but omits per-record `RG:Z:` on transcriptome BAM (issue 5 in `rustar_bam_comparison.md`)                                                                                                                                                                                                                                                                               |
| `--outSAMheaderHD`                 | `-`            | rejected            | Not in `--help`                                                                                                                                                                                                                                                                                                                                                                                                   |
| `--outSAMheaderPG`                 | `-`            | rejected            | Not in `--help`. (rustar's own `@PG` is content-free; issue 6 in `rustar_bam_comparison.md`)                                                                                                                                                                                                                                                                                                                      |
| `--outSAMheaderCommentFile`        | `-`            | rejected            | Not in `--help`                                                                                                                                                                                                                                                                                                                                                                                                   |
| `--outSAMfilter`                   | `None`         | rejected            | Not in `--help`                                                                                                                                                                                                                                                                                                                                                                                                   |
| `--outSAMmultNmax`                 | -1             | accepted-honoured   | Same default                                                                                                                                                                                                                                                                                                                                                                                                      |
| `--outSAMtlen`                     | 1              | rejected            | Not in `--help`. STAR-only TLEN signedness toggle                                                                                                                                                                                                                                                                                                                                                                 |
| `--outBAMcompression`              | 1              | accepted-honoured   | Same default. rustar adds `-1` (uncompressed) variant                                                                                                                                                                                                                                                                                                                                                             |
| `--outBAMsortingThreadN`           | 0              | rejected            | Not in `--help`                                                                                                                                                                                                                                                                                                                                                                                                   |
| `--outBAMsortingBinsN`             | 50             | rejected            | Not in `--help`                                                                                                                                                                                                                                                                                                                                                                                                   |
| `--bamRemoveDuplicatesType`        | `-`            | rejected            | Not in `--help`. Pipeline uses Picard MarkDuplicates downstream anyway                                                                                                                                                                                                                                                                                                                                            |
| `--bamRemoveDuplicatesMate2basesN` | 0              | rejected            | Same                                                                                                                                                                                                                                                                                                                                                                                                              |

### Wiggle / coverage output (none used by pipeline)

| STAR flag                  | STAR default | rustar status | Notes                                                          |
| -------------------------- | ------------ | ------------- | -------------------------------------------------------------- |
| `--outWigType`             | `None`       | rejected      | Not in `--help`. Pipeline uses deepTools / bedtools downstream |
| `--outWigStrand`           | `Stranded`   | rejected      | Same                                                           |
| `--outWigReferencesPrefix` | `-`          | rejected      | Same                                                           |
| `--outWigNorm`             | `RPM`        | rejected      | Same                                                           |

### Output filtering

| STAR flag                          | STAR default                | rustar status     | Notes                                                                 |
| ---------------------------------- | --------------------------- | ----------------- | --------------------------------------------------------------------- |
| `--outFilterType`                  | `Normal`                    | accepted-honoured | Pipeline passes `BySJout` for the RSEM path. rustar supports both     |
| `--outFilterMultimapNmax`          | 10                          | accepted-honoured | Pipeline passes 20. Same default                                      |
| `--outFilterMultimapScoreRange`    | 1                           | accepted-honoured | Probed; setting to 0 measurably tightens unique-vs-multi cutoff       |
| `--outFilterMismatchNmax`          | 10                          | accepted-honoured | Pipeline passes 999 for RSEM. Same default                            |
| `--outFilterMismatchNoverLmax`     | 0.3                         | accepted-honoured | Pipeline passes 0.04 for RSEM. Same default                           |
| `--outFilterMismatchNoverReadLmax` | 1.0                         | rejected          | Not in `--help`                                                       |
| `--outFilterScoreMin`              | 0                           | accepted-honoured | Same default                                                          |
| `--outFilterScoreMinOverLread`     | 0.66                        | accepted-honoured | **Probed (Detail C)**: setting to 0.95 drops mapping rate as expected |
| `--outFilterMatchNmin`             | 0                           | accepted-honoured | Same default                                                          |
| `--outFilterMatchNminOverLread`    | 0.66                        | accepted-honoured | Same default                                                          |
| `--outFilterIntronMotifs`          | `None`                      | accepted-honoured | Same default                                                          |
| `--outFilterIntronStrands`         | `RemoveInconsistentStrands` | accepted-honoured | Same default                                                          |

### Output SJ-filter parameters

| STAR flag                       | STAR default          | rustar status     | Notes                                                                          |
| ------------------------------- | --------------------- | ----------------- | ------------------------------------------------------------------------------ |
| `--outSJtype`                   | `Standard`            | rejected          | Not in `--help`. rustar always produces `Standard` SJ.out.tab                  |
| `--outSJfilterReads`            | `All`                 | rejected          | Not in `--help`                                                                |
| `--outSJfilterOverhangMin`      | `30 12 12 12`         | accepted-honoured | 4-tuple; same defaults                                                         |
| `--outSJfilterCountUniqueMin`   | `3 1 1 1`             | accepted-honoured | Same defaults                                                                  |
| `--outSJfilterCountTotalMin`    | `3 1 1 1`             | accepted-honoured | Same defaults                                                                  |
| `--outSJfilterDistToOtherSJmin` | `10 0 5 10`           | accepted-honoured | Same defaults                                                                  |
| `--outSJfilterIntronMaxVsReadN` | `50000 100000 200000` | accepted-honoured | 3-tuple; same defaults (rustar help text spells this as 3 values; STAR uses 3) |

### Scoring parameters

| STAR flag                       | STAR default | rustar status     | Notes        |
| ------------------------------- | ------------ | ----------------- | ------------ |
| `--scoreGap`                    | 0            | accepted-honoured | Same default |
| `--scoreGapNoncan`              | -8           | accepted-honoured | Same default |
| `--scoreGapGCAG`                | -4           | accepted-honoured | Same default |
| `--scoreGapATAC`                | -8           | accepted-honoured | Same default |
| `--scoreGenomicLengthLog2scale` | -0.25        | accepted-honoured | Same default |
| `--scoreDelOpen`                | -2           | accepted-honoured | Same default |
| `--scoreDelBase`                | -2           | accepted-honoured | Same default |
| `--scoreInsOpen`                | -2           | accepted-honoured | Same default |
| `--scoreInsBase`                | -2           | accepted-honoured | Same default |
| `--scoreStitchSJshift`          | 1            | accepted-honoured | Same default |

### Seed parameters

| STAR flag                        | STAR default | rustar status     | Notes           |
| -------------------------------- | ------------ | ----------------- | --------------- |
| `--seedSearchStartLmax`          | 50           | accepted-honoured | Same default    |
| `--seedSearchStartLmaxOverLread` | 1.0          | accepted-honoured | Same default    |
| `--seedSearchLmax`               | 0            | accepted-honoured | Same default    |
| `--seedMultimapNmax`             | 10000        | accepted-honoured | Same default    |
| `--seedPerReadNmax`              | 1000         | accepted-honoured | Same default    |
| `--seedPerWindowNmax`            | 50           | accepted-honoured | Same default    |
| `--seedNoneLociPerWindow`        | 10           | rejected          | Not in `--help` |
| `--seedSplitMin`                 | 12           | rejected          | Not in `--help` |
| `--seedMapMin`                   | 5            | accepted-honoured | Same default    |

### Alignment parameters

| STAR flag                            | STAR default | rustar status          | Notes                                                                                                                                                                                                                                                                                                      |
| ------------------------------------ | ------------ | ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--alignIntronMin`                   | 21           | accepted-honoured      | Pipeline passes 20 for RSEM. Same default                                                                                                                                                                                                                                                                  |
| `--alignIntronMax`                   | 0            | accepted-honoured      | Pipeline passes 1000000 for RSEM or 1 for prokaryotic. Same default semantics                                                                                                                                                                                                                              |
| `--alignMatesGapMax`                 | 0            | accepted-honoured      | Pipeline passes 1000000 for RSEM. Same default                                                                                                                                                                                                                                                             |
| `--alignSJoverhangMin`               | 5            | accepted-honoured      | Pipeline passes 8 for RSEM. Same default                                                                                                                                                                                                                                                                   |
| `--alignSJDBoverhangMin`             | 3            | **accepted-different** | Pipeline passes 1. rustar accepts and tracks the value, but because GTF-supplied SJs are not seeded into pass 1 alignment (Detail A and `rustar_two_pass_and_determinism.md`), this threshold has nothing to gate against; the flag is honoured but produces no observable behavioural change in this test |
| `--alignSJstitchMismatchNmax`        | `0 -1 0 0`   | accepted-honoured      | Same defaults                                                                                                                                                                                                                                                                                              |
| `--alignSplicedMateMapLmin`          | 0            | accepted-honoured      | Same default                                                                                                                                                                                                                                                                                               |
| `--alignSplicedMateMapLminOverLmate` | 0.66         | accepted-honoured      | Same default                                                                                                                                                                                                                                                                                               |
| `--alignWindowsPerReadNmax`          | 10000        | accepted-honoured      | Same default                                                                                                                                                                                                                                                                                               |
| `--alignTranscriptsPerWindowNmax`    | 100          | accepted-honoured      | Same default                                                                                                                                                                                                                                                                                               |
| `--alignTranscriptsPerReadNmax`      | 10000        | rejected               | Not in `--help`                                                                                                                                                                                                                                                                                            |
| `--alignEndsType`                    | `Local`      | rejected               | Not in `--help`                                                                                                                                                                                                                                                                                            |
| `--alignEndsProtrude`                | 0            | rejected               | Not in `--help`                                                                                                                                                                                                                                                                                            |
| `--alignSoftClipAtReferenceEnds`     | `Yes`        | rejected               | Not in `--help`                                                                                                                                                                                                                                                                                            |
| `--alignInsertionFlush`              | `None`       | rejected               | Not in `--help`                                                                                                                                                                                                                                                                                            |
| `--peOverlapNbasesMin`               | 0            | rejected               | Not in `--help`. Pipeline does not pass; PE overlap correction is off in the rnaseq default profile                                                                                                                                                                                                        |
| `--peOverlapMMp`                     | 0.01         | rejected               | Same                                                                                                                                                                                                                                                                                                       |

### Window / seed-cluster parameters

| STAR flag                      | STAR default | rustar status     | Notes           |
| ------------------------------ | ------------ | ----------------- | --------------- |
| `--winAnchorMultimapNmax`      | 50           | accepted-honoured | Same default    |
| `--winBinNbits`                | 16           | accepted-honoured | Same default    |
| `--winAnchorDistNbins`         | 9            | accepted-honoured | Same default    |
| `--winFlankNbins`              | 4            | accepted-honoured | Same default    |
| `--winReadCoverageRelativeMin` | 0.5          | accepted-honoured | Same default    |
| `--winReadCoverageBasesMin`    | 0            | rejected          | Not in `--help` |

### Chimeric output (relevant if user adds `--chimSegmentMin` via `extra_star_align_args`)

| STAR flag                    | STAR default  | rustar status       | Notes                                                                                                                                          |
| ---------------------------- | ------------- | ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `--chimOutType`              | `Junctions`   | accepted-honoured   | rustar supports `Junctions`. STAR also has `WithinBAM`, `SeparateSAMold`; rustar's other values not verified.                                  |
| `--chimSegmentMin`           | 0             | **accepted-broken** | Data layer works (14-column file matches STAR format), path layer crashes on pipeline-style `<sample>.` prefix. **MEDIUM** severity (Detail B) |
| `--chimScoreMin`             | 0             | accepted-honoured   | Same default                                                                                                                                   |
| `--chimScoreDropMax`         | 20            | accepted-honoured   | Same default                                                                                                                                   |
| `--chimScoreSeparation`      | 10            | accepted-honoured   | Same default                                                                                                                                   |
| `--chimScoreJunctionNonGTAG` | -1            | accepted-honoured   | Same default                                                                                                                                   |
| `--chimJunctionOverhangMin`  | 20            | accepted-honoured   | Same default                                                                                                                                   |
| `--chimSegmentReadGapMax`    | 0             | accepted-honoured   | Same default                                                                                                                                   |
| `--chimMainSegmentMultNmax`  | 10            | accepted-honoured   | Same default                                                                                                                                   |
| `--chimFilter`               | `banGenomicN` | rejected            | Not in `--help`                                                                                                                                |
| `--chimMultimapNmax`         | 0             | rejected            | Not in `--help`                                                                                                                                |
| `--chimMultimapScoreRange`   | 1             | rejected            | Not in `--help`                                                                                                                                |
| `--chimNonchimScoreDropMin`  | 20            | rejected            | Not in `--help`                                                                                                                                |
| `--chimOutJunctionFormat`    | 0             | rejected            | Not in `--help`                                                                                                                                |

### Quantification parameters

| STAR flag                            | STAR default                            | rustar status          | Notes                                                                                                                                                                                                                                                                                                                                 |
| ------------------------------------ | --------------------------------------- | ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--quantMode`                        | `-`                                     | accepted-honoured      | rustar supports `GeneCounts`, `TranscriptomeSAM`. Pipeline passes `TranscriptomeSAM`                                                                                                                                                                                                                                                  |
| `--quantTranscriptomeSAMoutput`      | `BanSingleEnd_BanIndels_ExtendSoftclip` | **accepted-different** | Pipeline passes `BanSingleEnd`. rustar accepts the flag and the singleton-ban does take effect (0 singletons in transcriptome BAM), but the resulting paired records lack `PROPERLY_PAIRED`, `RNEXT`, `PNEXT`, `TLEN` -- a separate, much more severe bug filed as scverse/rustar-aligner#22 (see `rustar_investigation_wt_rep2.md`). |
| `--quantTranscriptomeBAMcompression` | 1                                       | rejected               | Not in `--help`                                                                                                                                                                                                                                                                                                                       |

### Two-pass

| STAR flag          | STAR default | rustar status          | Notes                                                                                                                                                                                                                                                                                                                                                                    |
| ------------------ | ------------ | ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--twopassMode`    | `None`       | **accepted-different** | Pipeline passes `Basic`. rustar accepts `None` / `Basic` and runs the pass-1 step (probe confirmed: `Basic` 59s vs `None` 22s), but because no GTF junctions are seeded into pass 1's scoring path, the discovered-junction set is identical with or without two-pass on the test data. See `rustar_two_pass_and_determinism.md` Section A for the root cause and patch. |
| `--twopass1readsN` | -1           | accepted-honoured      | Same default                                                                                                                                                                                                                                                                                                                                                             |

### Variation / WASP (pipeline does not use)

| STAR flag          | STAR default | rustar status | Notes           |
| ------------------ | ------------ | ------------- | --------------- |
| `--varVCFfile`     | `-`          | rejected      | Not in `--help` |
| `--waspOutputMode` | `None`       | rejected      | Not in `--help` |

### STARsolo (pipeline does not use; rustar README says STARsolo is not implemented)

| STAR flag                  | STAR default | rustar status | Notes                                                                                                                    |
| -------------------------- | ------------ | ------------- | ------------------------------------------------------------------------------------------------------------------------ |
| `--soloType`               | `None`       | rejected      | Error: `unexpected argument '--soloType' found`. Error message could mention STARsolo by name.                           |
| `--soloCBtype`             | `Sequence`   | rejected      | Same                                                                                                                     |
| `--soloCBwhitelist`        | `-`          | rejected      | Same                                                                                                                     |
| `--soloFeatures`           | `Gene`       | rejected      | Same                                                                                                                     |
| `--soloUMIdedup`           | `1MM_All`    | rejected      | Same                                                                                                                     |
| `--soloStrand`             | `Forward`    | rejected      | Same                                                                                                                     |
| `--soloMultiMappers`       | `Unique`     | rejected      | Same                                                                                                                     |
| (`solo*` family, ~30 more) |              | rejected      | All STARsolo parameters are rejected at startup with the generic `unexpected argument` message; no per-flag probe needed |

## Detail subsections

### Detail A: `--alignSJDBoverhangMin` and `--sjdbScore` are honoured but produce no observable change

**Test:** WT_REP2 (PE), `--alignSJDBoverhangMin 30` vs the pipeline-default `1`.

**Result:** `Log.final.out` is byte-identical (modulo timestamps) between the two runs. 371 total splices, 0 annotated, 276 GT/AG, 92 non-canonical. Uniquely mapped 89.39%, multi 2.15%, too-short 8.46%.

**Why:** The flag's effect would be to require a longer overhang for splices that cross annotated junctions. The annotated-junction set in rustar's pass-1 scoring path is empty (the GTF SJs are loaded into a `SpliceJunctionDb` but the stitch-time call site doesn't query it; root-cause patch in `rustar_two_pass_and_determinism.md` Section A). With no annotated junctions to gate, the threshold has no records to filter. STAR with the same flag change would drop ~50 annotated splices on this sample.

**Verdict:** `accepted-different`. Flag is honoured insofar as rustar parses, validates, and tracks the value; the downstream effect is gated out by an unrelated bug. Same applies to `--sjdbScore` (the per-junction bonus, default 2): never applied because the call site doesn't recognise any junction as annotated.

### Detail B: `--chimSegmentMin` data layer works, path layer crashes

**Test 1 (pipeline-style prefix):** `--chimSegmentMin 12 --chimOutType Junctions --outFileNamePrefix WT_REP2.` — equivalent to what `align_star.config` produces.

**Result:** rustar starts up cleanly, logs `Chimeric detection enabled (chimSegmentMin=12)`, then aborts with:

```
Error: I/O error: No such file or directory (os error 2) (WT_REP2./Chimeric.out.junction)

Caused by:
    No such file or directory (os error 2)
```

No output files are written. The genome BAM, transcriptome BAM, and `Log.final.out` are all missing — the run is a hard failure.

**Test 2 (trailing-slash directory prefix):** `--outFileNamePrefix chimout_dir/` with `mkdir -p chimout_dir` first.

**Result:** Clean exit, `chimout_dir/Chimeric.out.junction` is produced. 5393 records on SE WT_REP2 (49 551 reads). Header column count is 14, matching STAR's documented format. Spot check of first records (`I 160892 + I 108830 - 0 0 0 SRR6357072.6691718 ...`): plausible chimeric coordinates within yeast chr I.

**Root cause:** rustar appears to derive the chimeric output path as `<outFileNamePrefix>/Chimeric.out.junction` regardless of whether `--outFileNamePrefix` ends in `/`. When the prefix is `WT_REP2.`, the path becomes `WT_REP2./Chimeric.out.junction` and the parent `WT_REP2.` directory has never been created. The non-chim outputs (Aligned.out.bam, Log.final.out etc) are written _via the same broken path_, so are also never created when chim is enabled. (Without `--chimSegmentMin > 0`, rustar's non-chim output path appears to fall back to a different file-name builder that handles the bare-prefix case, which is why the trailing-dot-prefix workaround in `modules/local/rustar_align/align/main.nf` works for normal runs.)

**Verdict:** `accepted-broken`, **medium severity**. The pipeline does not currently pass `--chimSegmentMin` so production runs are unaffected; but any user who adds it via `extra_star_align_args` will get a silent failure (Nextflow will surface the exit code, but the error message points at a file-not-found rather than a chim-output-path issue). The fix is one-liner upstream: ensure rustar creates `<prefix>/` when needed before writing chim outputs, or use the same path-builder that the non-chim outputs use.

### Detail C: `--outFilterScoreMinOverLread` is honoured

**Test:** WT_REP2 (PE), `--outFilterScoreMinOverLread 0.95` (vs default 0.66 advertised in rustar `--help` and STAR's `parametersDefault`).

**Result (Log.final.out diff):**

| Metric                 | Default 0.66 | Probe 0.95 |
| ---------------------- | ------------ | ---------- |
| Uniquely mapped %      | 89.39%       | 74.18%     |
| Multi-mapped %         | 2.15%        | 0.94%      |
| Unmapped too-short %   | 8.46%        | 24.88%     |
| Total splices          | 371          | 61         |
| Mismatch rate per base | 0.85%        | 0.27%      |

The drop in mapping rate (15 pp), the corresponding rise in too-short reads, and the lower mismatch rate of the surviving reads are all what STAR produces when the same flag is set on the same data. Mapping speed nearly halved (5.14 -> 4.05 Mreads/h) because rustar re-evaluates more candidate alignments before discarding them.

**Verdict:** `accepted-honoured`. No further action needed.

### Detail D: `--outSAMstrandField intronMotif`, `--outSAMattributes NM`, `--quantTranscriptomeSAMoutput BanSingleEnd`

Three flags advertised by rustar's `--help` with semantics matching STAR's, all three accepted-and-parsed without error, all three with output-level bugs that mean downstream tools get wrong or missing data.

- `--outSAMstrandField intronMotif`: rustar prints no error, runs to completion, but the resulting BAM has zero `XS:A:+/-` tags. STAR on the same data emits ~1% of records with `XS` (1 300 records in WT_REP2's 90 716 primaries). Downstream impact: rseqc `infer_experiment.py`, StringTie, Cufflinks all need `XS` for stranded-protocol inference. **HIGH severity**, filed as issue 2 in [`rustar_bam_comparison.md`](rustar_bam_comparison.md).

- `--outSAMattributes NM` (as part of the explicit list `NH HI AS NM MD`): rustar emits `nM:i:` (rustar's interpretation, substitution-only) rather than `NM:i:` (STAR / SAM-spec, mismatches + indel bases). Any tool that reads `NM:i:` (samtools stats, Picard, MultiQC's stats parsers) sees zero edit-distance values. **HIGH severity**, issue 1 in `rustar_bam_comparison.md`.

- `--quantTranscriptomeSAMoutput BanSingleEnd`: rustar honours the singleton-ban (no orphan-mate records in the transcriptome BAM, 0 singletons confirmed via `samtools flagstat`), but the _paired_ records that survive are missing `PROPERLY_PAIRED`, `RNEXT`, `PNEXT`, `TLEN` -- so Salmon can't infer fragment-length distribution and falls back to its 250 ± 25 prior, causing systematic TPM bias on short transcripts. **HIGH severity**, filed upstream as [scverse/rustar-aligner#22](https://github.com/scverse/rustar-aligner/issues/22), full analysis in [`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md).

### Detail E: STARsolo flag handling

All `--solo*` flags are rejected on parse with the generic `error: unexpected argument '--<flag>' found` and the bare `Usage: rustar-aligner [OPTIONS]` line. This is consistent with rustar's README ("STARsolo support is not yet implemented") but the error message could be more informative — a user copying a STARsolo invocation will see a clap-level "unknown flag" rather than an explicit "STARsolo is not implemented; please use STAR" hint. **Severity low (cosmetic)**.

## Punch list of upstream items (for scverse/rustar-aligner)

| #   | Severity   | Title                                                                                                                                     | One-line reproducer                                                                                                                                                      |
| --- | ---------- | ----------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 1   | **high**   | `--outSAMstrandField intronMotif` is accepted but `XS:A:` tags are never written                                                          | `samtools view <rustar.bam> \| grep -c 'XS:A:'` returns 0 with STAR-equivalent CLI                                                                                       |
| 2   | **high**   | `--outSAMattributes NM` emits `nM` (substitution-only) instead of `NM` (mismatches+indels)                                                | `samtools view <rustar.bam> \| head \| grep -oE '[nN]M:i:[0-9]+'` returns only `nM` even on reads with indels                                                            |
| 3   | **high**   | `--quantTranscriptomeSAMoutput BanSingleEnd` honoured for the singleton-ban but paired records have empty mate fields                     | scverse/rustar-aligner#22 already covers this; the surface affecting the matrix is "the flag is _advertised as STAR-equivalent_ in `--help`"                             |
| 4   | **medium** | `--chimSegmentMin > 0` aborts the run when `--outFileNamePrefix` does not end in `/`                                                      | `rustar-aligner --outFileNamePrefix sample. --chimSegmentMin 12 ...` errors with `(sample./Chimeric.out.junction)` no-such-file                                          |
| 5   | **medium** | `--alignSJDBoverhangMin` and `--sjdbScore` parsed-and-tracked but never applied because annotated junctions aren't queried at stitch time | (already in `rustar_two_pass_and_determinism.md`; this audit lists the user-facing symptom: the flags appear in `--help` with their STAR defaults, suggesting they work) |
| 6   | **low**    | STARsolo flags rejected with generic clap error message; should mention "STARsolo not yet implemented"                                    | `rustar-aligner --soloType CB_UMI_Simple` -> `error: unexpected argument '--soloType' found` with no domain hint                                                         |

## Ready-to-paste upstream issue body (item 4, the new finding)

```markdown
## Summary

When `--chimSegmentMin > 0` is set together with an `--outFileNamePrefix` that
does not end in `/` (e.g. the common STAR convention of `--outFileNamePrefix
SAMPLE.`), rustar-aligner v0.1.0 aborts before writing any output:
```

[INFO rustar_aligner] Chimeric detection enabled (chimSegmentMin=12)
Error: I/O error: No such file or directory (os error 2) (SAMPLE./Chimeric.out.junction)

Caused by:
No such file or directory (os error 2)

````

The genome BAM, transcriptome BAM, Log.final.out, and SJ.out.tab are also
missing -- the run is a complete failure, not just a missing chim file.

## Repro

```bash
mkdir -p chim_test && cd chim_test
# Assume genome index in star/, paired-end FASTQs in input1/ and input2/, GTF in genome.gtf
rustar-aligner \
  --genomeDir star \
  --readFilesIn input1/r1.fastq.gz input2/r2.fastq.gz \
  --runThreadN 4 \
  --outFileNamePrefix SAMPLE. \
  --sjdbGTFfile genome.gtf \
  --quantMode TranscriptomeSAM --outSAMtype BAM Unsorted \
  --outSAMattributes NH HI AS NM MD --readFilesCommand zcat \
  --twopassMode Basic --runRNGseed 0 --outFilterMultimapNmax 20 \
  --chimSegmentMin 12 --chimOutType Junctions
# -> exits non-zero with the path error above; no files in cwd
````

With `--outFileNamePrefix chimout_dir/` (trailing slash, parent dir
pre-created) the same invocation produces a correct 14-column
`chimout_dir/Chimeric.out.junction` and all other outputs. So the data
layer is fine, only the path layer is broken.

## Suggested fix

The chimeric-output path-builder appears to always append `/Chimeric.out.junction`
to `--outFileNamePrefix`, treating the prefix as a directory regardless of
whether it ends in `/`. The non-chim outputs (Aligned.out.bam etc) handle the
bare-prefix case correctly, so either:

- route chim outputs through the same `OutputPathBuilder` that the non-chim
  outputs use, or
- have the chim-output writer `create_dir_all(parent_of(target_path))` before
  opening the file (the same behaviour STAR has when prefixes contain a slash).

Either way, the trailing-dot / bare-token prefix convention from STAR should
be supported -- it's how most STAR wrappers (including nf-core/rnaseq) pass
sample names.

## Context

This breaks any nf-core/rnaseq run that adds `--chimSegmentMin` to
`extra_star_align_args`, which is the standard way users enable fusion-style
chimeric detection on top of the pipeline's normal output. Found while
auditing rustar's CLI surface for nf-core/rnaseq#1855 (full matrix at
`docs/rustar_cli_compat.md` on that branch).

## Verification a fix works

```bash
# After patching, both invocations should produce a Chimeric.out.junction:
rustar-aligner --outFileNamePrefix SAMPLE.   --chimSegmentMin 12 ...
ls SAMPLE.Chimeric.out.junction   # currently does not exist; should after fix
rustar-aligner --outFileNamePrefix dir/      --chimSegmentMin 12 ...
ls dir/Chimeric.out.junction      # already works today; should still work
```

```

## What this audit did not cover

- **STAR's input-BAM mode** (`--runMode inputAlignmentsFromBAM`, `--inputBAMfile`, `--bamRemoveDuplicates*`): rustar does not implement BAM input. Pipeline does not use this.
- **liftOver and genome-transform paths**: rustar does not implement either.
- **All STARsolo flags individually**: spot-probed `--soloType`, all reject with the same clap error. Not worth enumerating the remaining ~30.
- **Behavioural probe of every advertised flag**: probed the seven highest-leverage ones (the ones the pipeline passes and where the help text leaves ambiguity). For the remaining ~85 advertised flags where rustar's `--help` text matches STAR's manual semantically, the matrix records `accepted-honoured` from help-text agreement plus the default value match.
- **STAR-genomeGenerate parameters with effect on alignment** (e.g. `--genomeSAindexNbases` mismatch between index and alignment time): would need a full index regenerate plus alignment compare. Not in scope; the pipeline builds the index with rustar's own `genomeGenerate` so this surface is internally consistent.

## Files of interest

- Compatibility-matrix doc: this file.
- rustar `--help` capture: regenerate with `docker run --rm ghcr.io/scverse/rustar-aligner:dev rustar-aligner --help`.
- STAR parameter list: <https://github.com/alexdobin/STAR/blob/master/source/parametersDefault>.
- Pipeline STAR-args builder: `conf/modules/align_star.config`.
- Probe artefacts on the VM (under `/tmp/rustar_probes_real/`): `sjdb30/`, `2passnone/`, `chim_se/chimout_dir/`, `score_range_0/`, `readmap5k_/`. Pipeline-level probes: `/home/ubuntu/rnaseq-rustar-probe/results-rustar-probe-{scoremin,sjdb}/`.
```
