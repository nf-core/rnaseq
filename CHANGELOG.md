# nf-core/rnaseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[3.15.1](https://github.com/nf-core/rnaseq/releases/tag/3.15.1)] - 2024-09-16

### Enhancements & fixes

- [PR #1374](https://github.com/nf-core/rnaseq/pull/1374) - Bump pipeline version to 3.16.0dev
- [PR #1380](https://github.com/nf-core/rnaseq/pull/1380) - Fix issues with R modules changing sample names
- [PR #1381](https://github.com/nf-core/rnaseq/pull/1381) - Update all modules following massive conda usage update in nf-core modules
- [PR #1382](https://github.com/nf-core/rnaseq/pull/1381) - Slight fixes for rnaseq preprocessing
- [P$ #1383](https://github.com/nf-core/rnaseq/pull/1383) - Prerelease 3.15.1 version bumps

## [[3.15.0](https://github.com/nf-core/rnaseq/releases/tag/3.15.0)] - 2024-09-04

### Credits

Special thanks to the following for their contributions to the release:

- [Adam Talbot](https://github.com/adamrtalbot)
- [David Carlson](https://github.com/davidecarlson)
- [Edmund Miller](https://github.com/edmundmiller)
- [Jonathan Manning](https://github.com/pinin4fjords)
- [Laramie Lindsey](https://github.com/laramiellindsey)
- [Luke Zappia](https://github.com/lazappi)
- [Matthias Zepper](https://github.com/MatthiasZepper)
- [Maxime Garcia](https://github.com/maxulysse)
- [Pieter Moris](https://github.com/pmoris)
- [Rob Syme](https://github.com/robsyme)
- [Thomas Danhorn](https://github.com/tdanhorn)

Thank you to everyone else that has contributed by reporting bugs, enhancements or in any other way, shape or form.

### Enhancements & fixes

- [PR #1180](https://github.com/nf-core/rnaseq/pull/1180) - Bump pipeline version to 3.15.0dev
- [PR #1186](https://github.com/nf-core/rnaseq/pull/1186) - Properly update qualimap/rnaseq module (ie not patch)
- [PR #1197](https://github.com/nf-core/rnaseq/pull/1197) - Delete lib directory and replace with utils\_\* subworkflows
- [PR #1199](https://github.com/nf-core/rnaseq/pull/1199) - Replace modules.config with more modular config files per module/subworkflow/workflow
- [PR #1201](https://github.com/nf-core/rnaseq/pull/1201) - Template update for nf-core/tools v2.12
- [PR #1206](https://github.com/nf-core/rnaseq/pull/1206) - Remove `lib/` directory and `modules.config`
- [PR #1210](https://github.com/nf-core/rnaseq/pull/1210) - Use pseudoalignment subworkflow components from nf-core/modules
- [PR #1212](https://github.com/nf-core/rnaseq/pull/1212) - Update all modules and subworkflows to latest versions
- [PR #1213](https://github.com/nf-core/rnaseq/pull/1213) - Pass transcriptome fasta through to samtools stats
- [PR #1214](https://github.com/nf-core/rnaseq/pull/1214) - Bump umitools + delocalise prepareforrsem ([#831](https://github.com/nf-core/rnaseq/issues/831))
- [PR #1216](https://github.com/nf-core/rnaseq/pull/1216) - Delocalise catadditionalfasta ([#1162](https://github.com/nf-core/rnaseq/issues/1162))
- [PR #1217](https://github.com/nf-core/rnaseq/pull/1217) - Update Emiller88 => edmundmiller in README
- [PR #1218](https://github.com/nf-core/rnaseq/pull/1218) - Template update for nf-core/tools v2.13
- [PR #1220](https://github.com/nf-core/rnaseq/pull/1220) - Initialise nf-test and add pipeline level test
- [PR #1221](https://github.com/nf-core/rnaseq/pull/1221) - Use nf-test test for all nf-core components
- [PR #1226](https://github.com/nf-core/rnaseq/pull/1226) - Reuse bbsplit index and don't keep overwriting ([#1225](https://github.com/nf-core/rnaseq/issues/1225))
- [PR #1228](https://github.com/nf-core/rnaseq/pull/1228) - Make README usage consistent with docs/usage.md
- [PR #1229](https://github.com/nf-core/rnaseq/pull/1229) - Template update for nf-core/tools v2.13.1
- [PR #1231](https://github.com/nf-core/rnaseq/pull/1231) - Add sortmerna index possibilities
- [PR #1232](https://github.com/nf-core/rnaseq/pull/1232) - Add nf-test tests to star_genomegenerate_igenomes
- [PR #1233](https://github.com/nf-core/rnaseq/pull/1233) - Add nf-test tests to star_align_igenomes
- [PR #1234](https://github.com/nf-core/rnaseq/pull/1234) - Use genomecov from nf-core/modules
- [PR #1235](https://github.com/nf-core/rnaseq/pull/1235) - Add nf-test tests to utils_nfcore_rnaseq_pipeline tests
- [PR #1236](https://github.com/nf-core/rnaseq/pull/1236) - Add nf-test tests to gtf_filter
- [PR #1237](https://github.com/nf-core/rnaseq/pull/1237) - Fix concurrency error in Github CI workflow
- [PR #1238](https://github.com/nf-core/rnaseq/pull/1238) - Add nf-test tests to preprocess_transcripts_fasta_gencode
- [PR #1239](https://github.com/nf-core/rnaseq/pull/1239) - Add nf-test tests to align_star
- [PR #1240](https://github.com/nf-core/rnaseq/pull/1240) - Fix reference files params usage
- [PR #1241](https://github.com/nf-core/rnaseq/pull/1241) - Add nf-test tests to deseq2_qc
- [PR #1242](https://github.com/nf-core/rnaseq/pull/1242) - Use dupradar from nf-core/modules
- [PR #1243](https://github.com/nf-core/rnaseq/pull/1243) - Add nf-test for module MULTIQC_CUSTOM_BIOTYPE
- [PR #1244](https://github.com/nf-core/rnaseq/pull/1244) - Add gtf2bed tests
- [PR #1245](https://github.com/nf-core/rnaseq/pull/1245) - nf test quantify rsem
- [PR #1246](https://github.com/nf-core/rnaseq/pull/1246) - nf-test quantify pseudoalignment
- [PR #1247](https://github.com/nf-core/rnaseq/pull/1247) - nf-test prepare_genome
- [PR #1248](https://github.com/nf-core/rnaseq/pull/1248) - Improved ext.args consolidation for STAR and TRIMGALORE extra_args parameters
- [PR #1249](https://github.com/nf-core/rnaseq/pull/1249) - Include nf-tests for rsem_merge_counts module
- [PR #1250](https://github.com/nf-core/rnaseq/pull/1250) - Remove all tags.yml files because the testing system has changed
- [PR #1251](https://github.com/nf-core/rnaseq/pull/1251) - Replace deseq2qc paths
- [PR #1252](https://github.com/nf-core/rnaseq/pull/1252) - Fix genomeAttribute usage
- [PR #1253](https://github.com/nf-core/rnaseq/pull/1253) - Use nf-test files as matrix to test over in CI/CD for efficiency
- [PR #1260](https://github.com/nf-core/rnaseq/pull/1260) - Update CHANGELOG
- [PR #1261](https://github.com/nf-core/rnaseq/pull/1261) - Add more tests for PREPARE_GENOME
- [PR #1262](https://github.com/nf-core/rnaseq/pull/1262) - Fix CI pipeline
- [PR #1264](https://github.com/nf-core/rnaseq/pull/1264) - Add GHA files back into include statement of Github workflow change detection
- [PR #1265](https://github.com/nf-core/rnaseq/pull/1265) - Small updates noticed during code review
- [PR #1266](https://github.com/nf-core/rnaseq/pull/1266) - Delete unecessary tags from nf.test files for modules and subworkflows
- [PR #1271](https://github.com/nf-core/rnaseq/pull/1271) - Update trimming subworkflow to include more tests
- [PR #1272](https://github.com/nf-core/rnaseq/pull/1272) - Simple pipeline level nf-tests
- [PR #1274](https://github.com/nf-core/rnaseq/pull/1274) - Update bam_markduplicates_picard subworkflow
- [PR #1278](https://github.com/nf-core/rnaseq/pull/1278) - Delocalise pseudo quant workflow
- [PR #1279](https://github.com/nf-core/rnaseq/pull/1279) - Add pseudoaligner pipeline level tests to test suite
- [PR #1280](https://github.com/nf-core/rnaseq/pull/1280) - Reorganise pipeline level tests into flat directory structure
- [PR #1282](https://github.com/nf-core/rnaseq/pull/1282) - Fix CHANGELOG error
- [PR #1283](https://github.com/nf-core/rnaseq/pull/1283) - Add output files to nf-test snapshot
- [PR #1293](https://github.com/nf-core/rnaseq/pull/1293) - Update subworkflow utils_nfcore_pipeline
- [PR #1297](https://github.com/nf-core/rnaseq/pull/1297) - Important! Template update for nf-core/tools v2.14.1
- [PR #1302](https://github.com/nf-core/rnaseq/pull/1302) - Add missing files from Tximport processing
- [PR #1304](https://github.com/nf-core/rnaseq/pull/1304) - Remove redundant gene TPM outputs
- [PR #1306](https://github.com/nf-core/rnaseq/pull/1306) - Overhaul strandedness detection / comparison
- [PR #1307](https://github.com/nf-core/rnaseq/pull/1307) - Clarify infer strandedness step in subway map and text
- [PR #1308](https://github.com/nf-core/rnaseq/pull/1308) - Various MultiQC issues: FastQC sections for raw and trimmed reads // umi-tools dedup and extraction plots, custom content styling.
- [PR #1309](https://github.com/nf-core/rnaseq/pull/1309) - Document FASTP sampling
- [PR #1310](https://github.com/nf-core/rnaseq/pull/1310) - Reinstate pseudoalignment subworkflow config
- [PR #1312](https://github.com/nf-core/rnaseq/pull/1312) - Fix issues with unzipping of GTF/ GFF files without absolute paths
- [PR #1314](https://github.com/nf-core/rnaseq/pull/1314) - Add reference genome recommendations to usage docs
- [PR #1317](https://github.com/nf-core/rnaseq/pull/1317) - Strip problematic ifEmpty()
- [PR #1319](https://github.com/nf-core/rnaseq/pull/1319) - Reinstate oncomplete error messages
- [PR #1321](https://github.com/nf-core/rnaseq/pull/1321) - Remove push and release triggers from CI
- [PR #1322](https://github.com/nf-core/rnaseq/pull/1322) - Use pre-built Github Action to detect nf-test changes
- [PR #1323](https://github.com/nf-core/rnaseq/pull/1323) - Update actions/checkout to v4
- [PR #1324](https://github.com/nf-core/rnaseq/pull/1324) - Fix tags entries and rename pipeline level tests
- [PR #1325](https://github.com/nf-core/rnaseq/pull/1325) - Minor fixes to strandedness settings and messaging
- [PR #1326](https://github.com/nf-core/rnaseq/pull/1326) - Move Conda dependencies for local modules to individual environment file
- [PR #1328](https://github.com/nf-core/rnaseq/pull/1328) - Add pipeline level test for STAR-RSEM and HISAT2
- [PR #1329](https://github.com/nf-core/rnaseq/pull/1329) - Remove tags from all nf-test files
- [PR #1330](https://github.com/nf-core/rnaseq/pull/1330) - Update all nf-core/modules and subworkflows
- [PR #1331](https://github.com/nf-core/rnaseq/pull/1331) - Adding stubs for local modules
- [PR #1334](https://github.com/nf-core/rnaseq/pull/1334) - Update all nf-core/modules and subworkflows with stubs
- [PR #1335](https://github.com/nf-core/rnaseq/pull/1335) - Adding stubs at all levels
- [PR #1336](https://github.com/nf-core/rnaseq/pull/1336) - Use nf-core/setup-nf-test to install nf-test from cache during CI/CD
- [PR #1340](https://github.com/nf-core/rnaseq/pull/1340) - Remove out-of-date Azure specific guidance
- [PR #1341](https://github.com/nf-core/rnaseq/pull/1341) - Add rename in the MultiQC report for samples without techreps
- [PR #1342](https://github.com/nf-core/rnaseq/pull/1342) - Factor out preprocessing
- [PR #1345](https://github.com/nf-core/rnaseq/pull/1345) - Fix preprocessing call
- [PR #1350](https://github.com/nf-core/rnaseq/pull/1350) - Reduce resource usage for sort process in bedtools/genomecov
- [PR #1352](https://github.com/nf-core/rnaseq/pull/1352) - Assorted fixes to MultiQC usage
- [PR #1353](https://github.com/nf-core/rnaseq/pull/1353) - Correct conditional for salmon indexing in preprocessing workflow
- [PR #1355](https://github.com/nf-core/rnaseq/pull/1355) - Make all curves on subway map better looking, and all lines now have the same width
- [PR #1357](https://github.com/nf-core/rnaseq/pull/1357) - Fix anchor issue in multiqc
- [PR #1358](https://github.com/nf-core/rnaseq/pull/1358) - Update test profiles to restore a static URI for megatests
- [PR #1359](https://github.com/nf-core/rnaseq/pull/1359) - Update MultiQC and revert unnecessary workaround
- [PR #1360](https://github.com/nf-core/rnaseq/pull/1360) - More complete snapshots for pipeline level tests
- [PR #1361](https://github.com/nf-core/rnaseq/pull/1361) - Animate subway map
- [PR #1362](https://github.com/nf-core/rnaseq/pull/1362) - Move multiqc module prefix for nf-test to module
- [PR #1363](https://github.com/nf-core/rnaseq/pull/1363) - Minor updates of nf-core modules and subworkflows
- [PR #1363](https://github.com/nf-core/rnaseq/pull/1363) - Update dupradar script
- [PR #1366](https://github.com/nf-core/rnaseq/pull/1366) - Clarify docs on different tximport count files
- [PR #1367](https://github.com/nf-core/rnaseq/pull/1367) - Clarify design formula and blind dispersion estimation
- [PR #1370](https://github.com/nf-core/rnaseq/pull/1370) - Bump versions for 3.15.0
- [PR #1371](https://github.com/nf-core/rnaseq/pull/1371) - Apply Maxime's CHANGELOG edits
- [PR #1372](https://github.com/nf-core/rnaseq/pull/1372) - Bump tximeta/tximport for gene table row names fix

### Parameters

| Old parameter      | New parameter                    |
| ------------------ | -------------------------------- |
|                    | `--pipelines_testdata_base_path` |
|                    | `--sortmerna_index`              |
|                    | `--stranded_threshold`           |
|                    | `--unstranded_threshold`         |
| `--test_data_base` |                                  |

### Software dependencies

| Dependency                          | Old version | New version |
| ----------------------------------- | ----------- | ----------- |
| `bedtools`                          | 2.30.0      | 2.31.1      |
| `bioconductor-dupradar`             | 1.28.0      | 1.32.0      |
| `bioconductor-summarizedexperiment` | 1.24.0      | 1.32.0      |
| `bioconductor-tximeta`              | 1.12.0      | 1.20.1      |
| `gffread`                           | 0.12.1      | 0.12.7      |
| `multiqc`                           | 1.20        | 1.22.3      |
| `picard`                            | 3.0.0       | 3.1.1       |
| `samtools`                          | 1.17        | 1.20        |
| `sortmerna`                         | 4.3.4       | 4.3.6       |
| `umi_tools`                         | 1.14        | 1.15        |
| `untar`                             | 1.3         | 1.34        |

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

## [[3.14.0](https://github.com/nf-core/rnaseq/releases/tag/3.14.0)] - 2024-01-08

### Credits

Special thanks to the following for their contributions to the release:

- [Adam Talbot](https://github.com/adamrtalbot)
- [Jonathan Manning](https://github.com/pinin4fjords)
- [Mahesh Binzer-Panchal](https://github.com/mahesh-panchal)
- [Matthias Zepper](https://github.com/MatthiasZepper)
- [Maxime Garcia](https://github.com/maxulysse)
- [Phil Ewels](https://github.com/ewels)
- [Vlad Savelyev](https://github.com/vladsavelyev)

Thank you to everyone else that has contributed by reporting bugs, enhancements or in any other way, shape or form.

### Enhancements & fixes

- [PR #1135](https://github.com/nf-core/rnaseq/pull/1135) - Update [action-tower-launch](https://github.com/marketplace/actions/action-tower-launch) to v2 which supports more variable handling
- [PR #1141](https://github.com/nf-core/rnaseq/pull/1141) - Important! Template update for nf-core/tools v2.11
- [PR #1143](https://github.com/nf-core/rnaseq/pull/1143) - Move fasta check back to Groovy ([#1142](https://github.com/nf-core/rnaseq/issues/1142))
- [PR #1144](https://github.com/nf-core/rnaseq/pull/1144) - Interface to kmer size for pseudoaligners ([#1111](https://github.com/nf-core/rnaseq/issues/1111))
- [PR #1149](https://github.com/nf-core/rnaseq/pull/1149) - Fix and patch version commands for Fastp, FastQC and UMI-tools modules ([#1103](https://github.com/nf-core/rnaseq/issues/1103))
- [PR #1150](https://github.com/nf-core/rnaseq/pull/1150) - Be more flexible on attribute values in GTFs ([#1132](https://github.com/nf-core/rnaseq/issues/1132))
- [PR #1151](https://github.com/nf-core/rnaseq/pull/1151) - fix to #1150: reinstate conditional
- [PR #1152](https://github.com/nf-core/rnaseq/pull/1152) - Bump container versions for tools using Docker V1 manifest ([#1140](https://github.com/nf-core/rnaseq/issues/1140))
- [PR #1154](https://github.com/nf-core/rnaseq/pull/1154) - Prerelease 3.14.0 fixes ([#1111](https://github.com/nf-core/rnaseq/issues/1111), [#1153](https://github.com/nf-core/rnaseq/issues/1153))
- [PR #1157](https://github.com/nf-core/rnaseq/pull/1157) - Add slash to `--outdir` for cloud tests to fix Azure validation issue
- [PR #1159](https://github.com/nf-core/rnaseq/pull/1159) - Issues loading MultiQC report ([#1158](https://github.com/nf-core/rnaseq/issues/1158))

### Parameters

| Old parameter | New parameter                |
| ------------- | ---------------------------- |
|               | `--pseudo_aligner_kmer_size` |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `multiqc`  | 1.17        | 1.19        |
| `qualimap` | 2.2.2d      | 2.3         |
| `rseqc`    | 3.0.1       | 5.0.3       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

## [[3.13.2](https://github.com/nf-core/rnaseq/releases/tag/3.13.2)] - 2023-11-21

### Credits

Special thanks to the following for their contributions to the release:

- [Jonathan Manning](https://github.com/pinin4fjords)
- [Regina Hertfelder Reynolds](https://github.com/RHReynolds)
- [Matthias Zepper](https://github.com/MatthiasZepper)

### Enhancements & fixes

- [PR #1123](https://github.com/nf-core/rnaseq/pull/1123) - Overhaul tximport.r, output length tables
- [PR #1124](https://github.com/nf-core/rnaseq/pull/1124) - Ensure pseudoaligner is set if pseudoalignment is not skipped
- [PR #1126](https://github.com/nf-core/rnaseq/pull/1126) - Pipeline fails if transcript_fasta not provided and `skip_gtf_filter = true`.
- [PR #1127](https://github.com/nf-core/rnaseq/pull/1127) - Enlarge sampling to determine the number of columns in `filter_gtf.py` script.

## [[3.13.1](https://github.com/nf-core/rnaseq/releases/tag/3.13.1)] - 2023-11-17

### Enhancements and fixes

- [PR #1121](https://github.com/nf-core/rnaseq/pull/1121) - Changes for 3.13.1 patch release incl. igenomes star fix

## [[3.13.0](https://github.com/nf-core/rnaseq/releases/tag/3.13.0)] - 2023-11-17

### Credits

Special thanks to the following for their contributions to the release:

- [Adam Talbot](https://github.com/adamrtalbot)
- [hmehlan](https://github.com/hmehlan)
- [Jonathan Manning](https://github.com/pinin4fjords)
- [Júlia Mir Pedrol](https://github.com/mirpedrol)
- [Matthias Zepper](https://github.com/MatthiasZepper)
- [Maxime Garcia](https://github.com/maxulysse)
- [Steffen Möller](https://github.com/smoe)

Thank you to everyone else that has contributed by reporting bugs, enhancements or in any other way, shape or form.

### Enhancements & fixes

- [PR #1049](https://github.com/nf-core/rnaseq/pull/1049) - Display a warning when `--extra_star_align_args` are used with `--aligner star_rsem`
- [PR #1051](https://github.com/nf-core/rnaseq/pull/1051) - Remove `public_aws_ecr` profile
- [PR #1054](https://github.com/nf-core/rnaseq/pull/1054) - Template update to nf-core/tools v2.9
- [PR #1058](https://github.com/nf-core/rnaseq/pull/1058) - Use `nf-validation` plugin for parameter and samplesheet validation
- [PR #1068](https://github.com/nf-core/rnaseq/pull/1068) - Update `grep` version for `untar` module
- [PR #1073](https://github.com/nf-core/rnaseq/pull/1073) - Update documentation to discourage use of `--genome`
- [PR #1078](https://github.com/nf-core/rnaseq/pull/1078) - Updated pipeline template to [nf-core/tools 2.10](https://github.com/nf-core/tools/releases/tag/2.10)
- [PR #1083](https://github.com/nf-core/rnaseq/pull/1083) - Move local modules and subworkflows to subfolders
- [PR #1088](https://github.com/nf-core/rnaseq/pull/1088) - Updates contributing and code of conduct documents with nf-core template 2.10
- [PR #1091](https://github.com/nf-core/rnaseq/pull/1091) - Reorganise parameters in schema for better usability
- [PR #1106](https://github.com/nf-core/rnaseq/pull/1106) - Kallisto quantification
- [PR #1107](https://github.com/nf-core/rnaseq/pull/1107) - Expand GTF filtering to remove rows with empty transcript ID when required, fix STAR GTF usage
- [#976](https://github.com/nf-core/rnaseq/issues/976) - Add author and licenses for all custom scripts
- [#1050](https://github.com/nf-core/rnaseq/issues/1050) - Provide custom prefix/suffix for summary files to avoid overwriting
- [#1074](https://github.com/nf-core/rnaseq/issues/1074) - Enable quantification using StringTie AND a custom
- [#1082](https://github.com/nf-core/rnaseq/issues/1082) - More informative error message for `filter_gtf_for_genes_in_genome.py`
- [#1102](https://github.com/nf-core/rnaseq/issues/1102) - gene entries with empty transcript_id fields

### Software dependencies

| Dependency              | Old version | New version |
| ----------------------- | ----------- | ----------- |
| `fastqc`                | 0.11.9      | 0.12.1      |
| `multiqc`               | 1.14        | 1.17        |
| `ucsc-bedgraphtobigwig` | 377         | 445         |

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

### Modules / Subworkflows

| Script                   | Old name          | New name                    |
| ------------------------ | ----------------- | --------------------------- |
| `local/gtf_filter`       | `GTF_GENE_FILTER` | `GTF_FILTER`                |
| `local/tx2gene`          | `SALMON_TX2GENE`  | `TX2GENE`                   |
| `local/tximport`         | `SALMON_TXIMPORT` | `TXIMPORT`                  |
| `local/quantify_salmon`  | `QUANTIFY_SALMON` | `QUANTIFY_PSEUDO_ALIGNMENT` |
| `nf-core/kallisto_index` |                   | `KALLISTO_INDEX`            |
| `nf-core/kallisto_quant` |                   | `KALLISTO_QUANT`            |

## [[3.12.0](https://github.com/nf-core/rnaseq/releases/tag/3.12.0)] - 2023-06-02

### Credits

Special thanks to the following for their contributions to the release:

- [Adam Talbot](https://github.com/adamrtalbot)
- [Esha Joshi](https://github.com/ejseqera)
- [Ghepardo](https://github.com/Ghepardo)
- [Matthias Zepper](https://github.com/MatthiasZepper)
- [Maxime Garcia](https://github.com/maxulysse)
- [Rob Syme](https://github.com/robsyme)

Thank you to everyone else that has contributed by reporting bugs, enhancements or in any other way, shape or form.

### Enhancements & fixes

- [[#1011](https://github.com/nf-core/rnaseq/issues/1011)] - FastQ files from UMI-tools not being passed to fastp
- [[#1018](https://github.com/nf-core/rnaseq/issues/1018)] - Ability to skip both alignment and pseudoalignment to only run pre-processing QC steps.
- [PR #1016](https://github.com/nf-core/rnaseq/pull/1016) - Updated pipeline template to [nf-core/tools 2.8](https://github.com/nf-core/tools/releases/tag/2.8)
- [PR #1025](https://github.com/nf-core/fetchngs/pull/1025) - Add `public_aws_ecr.config` to source mulled containers when using `public.ecr.aws` Docker Biocontainer registry
- [PR #1038](https://github.com/nf-core/rnaseq/pull/1038) - Updated error log for count values when supplying `--additional_fasta`
- [PR #1042](https://github.com/nf-core/rnaseq/pull/1042) - revert samtools_sort modules to no memory assignement

### Parameters

| Old parameter | New parameter             |
| ------------- | ------------------------- |
|               | `--skip_pseudo_alignment` |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `fastp`    | 0.23.2      | 0.23.4      |
| `samtools` | 1.16.1      | 1.17        |

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

## [[3.11.2](https://github.com/nf-core/rnaseq/releases/tag/3.11.2)] - 2023-04-25

### Credits

Special thanks to the following for their contributions to the release:

- [Jonathan Manning](https://github.com/pinin4fjords)
- [Maxime Garcia](https://github.com/maxulysse)
- [Rob Syme](https://github.com/robsyme)
- [W. Lee Pang](https://github.com/wleepang)

Thank you to everyone else that has contributed by reporting bugs, enhancements or in any other way, shape or form.

### Enhancements & fixes

- [[#1003](https://github.com/nf-core/rnaseq/pull/1003)] - `FASTQ_SUBSAMPLE_FQ_SALMON:SALMON_INDEX` is launched multiple times and fails

## [[3.11.1](https://github.com/nf-core/rnaseq/releases/tag/3.11.1)] - 2023-03-31

### Credits

Special thanks to the following for their code contributions to the release:

- [Adam Talbot](https://github.com/adamrtalbot)
- [Rob Syme](https://github.com/robsyme)
- [suhrig](https://github.com/suhrig)

### Enhancements & fixes

- [[#987](https://github.com/nf-core/rnaseq/pull/987)] - Fix issue with incorrect cacheing of test datasets during CI/CD
- [[#988](https://github.com/nf-core/rnaseq/issues/988)] - `DESEQ2_QC_STAR_SALMON` fails when sample names have many components
- Remove `wait: false` option from Tower Actions which is the default
- Fix release trigger for full-sized multi-cloud tests
- Adding `[ci fast]` to commit message now skips all tests except for standard `-profile test` pipeline run

## [[3.11.0](https://github.com/nf-core/rnaseq/releases/tag/3.11.0)] - 2023-03-30

### Credits

Special thanks to the following for their code contributions to the release:

- [J Lorent](https://github.com/jlorent)
- [Luca Beltrame](https://github.com/lbeltrame)
- [Matthias Zepper](https://github.com/MatthiasZepper)
- [Maxime Garcia](https://github.com/maxulysse)
- [Ryan Yordanoff](https://github.com/ryanyord)
- [Thomas Sandmann](https://github.com/tomsing1)

Thank you to everyone else that has contributed by reporting bugs, enhancements or in any other way, shape or form.

### Enhancements & fixes

- Add infrastructure and CI for multi-cloud full-sized tests run via Nextflow Tower (see [#981](https://github.com/nf-core/rnaseq/pull/981))
- Added fastp support.
  - Users can now select between `--trimmer trimgalore` (default) and `--trimmer fastp`.
  - Trim Galore! specific pipeline parameters have been deprecated: `--clip_r1`, `--clip_r2`, `--three_prime_clip_r1`, `--three_prime_clip_r2` and `--trim_nextseq`
  - Any additional options can now be specified via the `--extra_trimgalore_args` and `--extra_fastp_args` parameters, respectively.
- [[#663](https://github.com/nf-core/rnaseq/issues/663)] - Alternative trimming step for polyA/T removal
- [[#781](https://github.com/nf-core/rnaseq/issues/781)] - Add Warning for poly(A) libraries
- [[#878](https://github.com/nf-core/rnaseq/issues/878)] - Allow tabs in fasta header when creating decoys for salmon index
- [[#931](https://github.com/nf-core/rnaseq/issues/931)] - Save transcriptome BAM files when using `--save_umi_intermeds` / `--save_align_intermeds`
- [[#934](https://github.com/nf-core/rnaseq/pull/934)] - Union of `ext.args` and `params.extra_star_align_args` prevents parameter clashes in the STAR module
- [[#940](https://github.com/nf-core/rnaseq/issues/940)] - Bugfix in `salmon_summarizedexperiment.r` to ensure `rbind` doesn't fail when `rowdata` has no `tx` column.
- [[#944](https://github.com/nf-core/rnaseq/issues/944)] - Read clipping using clip_r1, clip_r2, three_prime_clip_r1, three_prime_clip_r2 disabled in 3.10
- [[#956](https://github.com/nf-core/rnaseq/pull/956)] - Implement 'auto' as default strandedness argument in `fastq_dir_to_samplesheet.py` script
- [[#960](https://github.com/nf-core/rnaseq/issues/960)] - Failure with awsbatch when running processes that are using `executor: local`
- [[#961](https://github.com/nf-core/rnaseq/issues/961)] - Add warnings to STDOUT for all skipped and failed strandedness check samples
- [[#975](https://github.com/nf-core/rnaseq/issues/975)] - `SALMON_INDEX` runs when using `--aligner star_rsem` even if samples have explicit strandedness
- Remove HISAT2 from automated AWS full-sized tests

### Parameters

| Old parameter           | New parameter             |
| ----------------------- | ------------------------- |
|                         | `--trimmer`               |
|                         | `--extra_trimgalore_args` |
| `--clip_r1`             |                           |
| `--clip_r2`             |                           |
| `--three_prime_clip_r1` |                           |
| `--three_prime_clip_r2` |                           |
| `--tracedir`            |                           |
| `--trim_nextseq`        |                           |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| `fastp`     |             | 0.23.2      |
| `multiqc`   | 1.13        | 1.14        |
| `picard`    | 2.27.4      | 3.0.0       |
| `salmon`    | 1.9.0       | 1.10.1      |
| `umi_tools` | 1.1.2       | 1.1.4       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

## [[3.10.1](https://github.com/nf-core/rnaseq/releases/tag/3.10.1)] - 2023-01-05

### Enhancements & fixes

- [[#919](https://github.com/nf-core/rnaseq/issues/919)] - Salmon quant not run after FastQ subsampling if index not provided
- [[#922](https://github.com/nf-core/rnaseq/issues/922)] - Passing TrimGalore `--hardtrim3` / `--hardtrim5` via custom config raises missing output filename error

## [[3.10](https://github.com/nf-core/rnaseq/releases/tag/3.10)] - 2022-12-21

### Enhancements & fixes

- Bump minimum Nextflow version from `21.10.3` -> `22.10.1`
- Updated pipeline template to [nf-core/tools 2.7.2](https://github.com/nf-core/tools/releases/tag/2.7.2)
- [[#729](https://github.com/nf-core/rnaseq/issues/729)] - Add 'auto' option to samplesheet to automatically detect strandedness for samples
- [[#889](https://github.com/nf-core/rnaseq/issues/889)] - Document valid options for `--genome` parameter
- [[#891](https://github.com/nf-core/rnaseq/issues/891)] - Skip MarkDuplicates when UMIs are used
- [[#896](https://github.com/nf-core/rnaseq/issues/896)] - Remove `copyTo` call for iGenomes README
- [[#897](https://github.com/nf-core/rnaseq/issues/897)] - Use `--skip_preseq` by default
- [[#898](https://github.com/nf-core/rnaseq/issues/898)] - Documentation on salmon decoy-aware index creation, gcbias and seqbias
- [[#900](https://github.com/nf-core/rnaseq/issues/900)] - Add `--recursive` option to `fastq_dir_to_samplesheet.py` script
- [[#902](https://github.com/nf-core/rnaseq/issues/902)] - `check_samplesheet.py` script doesn't output optional columns in samplesheet
- [[#907](https://github.com/nf-core/rnaseq/issues/907)] - Add `--extra_star_align_args` and `--extra_salmon_quant_args` parameter
- [[#912](https://github.com/nf-core/rnaseq/issues/912)] - Add UMI deduplication before quantification in tube map

### Parameters

| Old parameter    | New parameter               |
| ---------------- | --------------------------- |
| `--enable_conda` |                             |
|                  | `--extra_star_align_args`   |
|                  | `--extra_salmon_quant_args` |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency                          | Old version | New version |
| ----------------------------------- | ----------- | ----------- |
| `bbmap`                             | 38.93       | 39.01       |
| `bioconductor-dupradar`             | 1.18.0      | 1.28.0      |
| `bioconductor-summarizedexperiment` | 1.20.0      | 1.24.0      |
| `bioconductor-tximeta`              | 1.8.0       | 1.12.0      |
| `fq`                                |             | 0.9.1       |
| `salmon`                            | 1.5.2       | 1.9.0       |
| `samtools`                          | 1.15.1      | 1.16.1      |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [[3.9](https://github.com/nf-core/rnaseq/releases/tag/3.9)] - 2022-09-30

### Enhancements & fixes

- [[#746](https://github.com/nf-core/rnaseq/issues/746)] - Add `tin.py` output to MultiQC report
- [[#841](https://github.com/nf-core/rnaseq/issues/841)] - Turn `--deseq2_vst` on by default
- [[#853](https://github.com/nf-core/rnaseq/issues/853)] - Pipeline fails at email step: Failed to invoke `workflow.onComplete` event handler
- [[#857](https://github.com/nf-core/rnaseq/issues/857)] - Missing parameter required by StringTie if using STAR as aligner
- [[#862](https://github.com/nf-core/rnaseq/issues/862)] - Filter samples that have no reads after trimming
- [[#864](https://github.com/nf-core/rnaseq/issues/864)] - Pre-process transcripts fasta when using `--gencode`
- Expose additional arguments to UMI-tools as pipeline params: `--umitools_bc_pattern2` is required if the UMI is located on read 2. `--umitools_umi_separator` will often be needed in conjunction with `--skip_umi_extract` as most other tools such as Illumina's `BCL Convert` use a colon instead of an underscore to separate the UMIs. The `--umitools_grouping_method` allows to fine-tune handling of similar but non-identical UMIs.
- Updated pipeline template to [nf-core/tools 2.5.1](https://github.com/nf-core/tools/releases/tag/2.5.1)

### Parameters

| Old parameter | New parameter                |
| ------------- | ---------------------------- |
|               | `--umitools_bc_pattern2`     |
|               | `--umitools_umi_separator`   |
|               | `--umitools_grouping_method` |

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `hisat2`   | 2.2.0       | 2.2.1       |
| `multiqc`  | 1.11        | 1.13        |
| `picard`   | 2.26.10     | 2.27.4      |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [[3.8.1](https://github.com/nf-core/rnaseq/releases/tag/3.8.1)] - 2022-05-27

- [[#834](https://github.com/nf-core/rnaseq/issues/834)] - `nf-core download` fails with version 3.8 of the pipeline

## [[3.8](https://github.com/nf-core/rnaseq/releases/tag/3.8)] - 2022-05-25

### :warning: Major enhancements

Fixed quite a well hidden bug in the UMI processing mode of the pipeline when using `--with_umi --aligner star_salmon` as reported by [Lars Roed Ingerslev](https://github.com/lars-work-sund). Paired-end BAM files were not appropriately name sorted after `umi_tools dedup` which ultimately resulted in incorrect reading and quantification with Salmon. If you have used previous versions of the pipeline to analyse paired-end UMI data it will need to be reprocessed using this version of the pipeline. See [#828](https://github.com/nf-core/rnaseq/issues/828) for more context.

### Enhancements & fixes

- [[#824](https://github.com/nf-core/rnaseq/issues/824)] - Add explicit docs for usage of featureCounts in the pipeline
- [[#825](https://github.com/nf-core/rnaseq/issues/825)] - Pipeline fails due to trimming related removal of all reads from a sample
- [[#827](https://github.com/nf-core/rnaseq/issues/827)] - Control generation of --output-stats when running umi-tools dedup
- [[#828](https://github.com/nf-core/rnaseq/issues/828)] - Filter BAM output of UMI-tools dedup before passing to Salmon quant
- Updated pipeline template to [nf-core/tools 2.4.1](https://github.com/nf-core/tools/releases/tag/2.4.1)

### Parameters

| Old parameter | New parameter            |
| ------------- | ------------------------ |
|               | `--min_trimmed_reads`    |
|               | `--umitools_dedup_stats` |

## [[3.7](https://github.com/nf-core/rnaseq/releases/tag/3.7)] - 2022-05-03

### :warning: Major enhancements

- Updated default STAR version to latest available (`2.7.10a`; see [#808](https://github.com/nf-core/rnaseq/issues/808]))
- Vanilla Linux Docker container changed from `biocontainers/biocontainers:v1.2.0_cv1` to `ubuntu:20.04` to fix issues observed on GCP (see [#764](https://github.com/nf-core/rnaseq/issues/764]))

### Enhancements & fixes

- [[#762](https://github.com/nf-core/rnaseq/issues/762)] - Explicitly set `--skip_bbsplit false` with `--bbsplit_fasta_list` to use BBSplit
- [[#764](https://github.com/nf-core/rnaseq/issues/764)] - Test fails when using GCP due to missing tools in the basic biocontainer
- [[#765](https://github.com/nf-core/rnaseq/issues/765)] - Add docs for the usage of nf-core/rnaseq with prokaryotic data
- [[#775](https://github.com/nf-core/rnaseq/issues/775)] - Incorrect columns in Salmon transcript files
- [[#791](https://github.com/nf-core/rnaseq/issues/791)] - Add outputs for umitools dedup summary stats
- [[#797](https://github.com/nf-core/rnaseq/issues/797)] - Add `--skip_umi_extract` to account for pre-existing UMIs header embeddings.
- [[#798](https://github.com/nf-core/rnaseq/issues/798)] - Decompress transcript fasta error
- [[#799](https://github.com/nf-core/rnaseq/issues/799)] - Issue with using `--retain_unpaired` with the `FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE` module
- [[#802](https://github.com/nf-core/rnaseq/issues/802)] - `--bam_csi_index` error generated if `--skip_alignment` specified
- [[#808](https://github.com/nf-core/rnaseq/issues/808)] - Auto-detect usage of Illumina iGenomes reference
- [[#809](https://github.com/nf-core/rnaseq/issues/809)] - Add metro map for pipeline
- [[#814](https://github.com/nf-core/rnaseq/issues/814)] - Use decimal values for `--min_mapped_reads`
- Updated pipeline template to [nf-core/tools 2.3.2](https://github.com/nf-core/tools/releases/tag/2.3.2)

### Parameters

| Old parameter | New parameter        |
| ------------- | -------------------- |
|               | `--skip_umi_extract` |

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| `samtools`  | 1.14        | 1.15.1      |
| `star`      | 2.6.1d      | 2.7.10a     |
| `stringtie` | 2.1.7       | 2.2.1       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [[3.6](https://github.com/nf-core/rnaseq/releases/tag/3.6)] - 2022-03-04

### Enhancements & fixes

- [nf-core/tools#1415](https://github.com/nf-core/tools/issues/1415) - Make `--outdir` a mandatory parameter
- [[#734](https://github.com/nf-core/rnaseq/issues/734)] - Is a vulnerable picard still used ? log4j vulnerability
- [[#744](https://github.com/nf-core/rnaseq/issues/744)] - Auto-detect and raise error if CSI is required for BAM indexing
- [[#750](https://github.com/nf-core/rnaseq/issues/750)] - Optionally ignore R1 / R2 after UMI extraction process
- [[#752](https://github.com/nf-core/rnaseq/issues/752)] - How to set publishing mode for all processes?
- [[#753](https://github.com/nf-core/rnaseq/issues/753)] - Add warning when user provides `--transcript_fasta`
- [[#754](https://github.com/nf-core/rnaseq/issues/754)] - DESeq2 QC issue linked to `--count_col` parameter
- [[#755](https://github.com/nf-core/rnaseq/issues/755)] - Rename RSEM_PREPAREREFERENCE_TRANSCRIPTS process
- [[#759](https://github.com/nf-core/rnaseq/issues/759)] - Empty lines in samplesheet.csv cause a crash
- [[#769](https://github.com/nf-core/rnaseq/issues/769)] - Do not run RSeQC tin.py by default

### Parameters

| Old parameter | New parameter        |
| ------------- | -------------------- |
|               | `--publish_dir_mode` |
|               | `--umi_discard_read` |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
>
> **NB:** Parameter has been **added** if just the new parameter information is present.
>
> **NB:** Parameter has been **removed** if new parameter information isn't present.

## [[3.5](https://github.com/nf-core/rnaseq/releases/tag/3.5)] - 2021-12-17

### Enhancements & fixes

- Port pipeline to the updated Nextflow DSL2 syntax adopted on nf-core/modules
  - Removed `--publish_dir_mode` as it is no longer required for the new syntax
- Bump minimum Nextflow version from `21.04.0` -> `21.10.3`
- Updated pipeline template to [nf-core/tools 2.2](https://github.com/nf-core/tools/releases/tag/2.2)
- [[#664](https://github.com/nf-core/rnaseq/issues/664)] - Conflict of library names for technical replicates
- [[#720](https://github.com/nf-core/rnaseq/issues/720)] - KeyError 'gene_id' in salmon_tx2gene.py
- [[#724](https://github.com/nf-core/rnaseq/issues/724)] - Deal with warnings generated when native NF processes are used
- [[#725](https://github.com/nf-core/rnaseq/issues/725)] - Untar needs `--no-same-owner` on DNAnexus
- [[#727](https://github.com/nf-core/rnaseq/issues/727)] - Fix transcriptome staging issues on DNAnexus for rsem/prepareference
- [[#728](https://github.com/nf-core/rnaseq/issues/728)] - Add RSeQC TIN.py as a quality metric for the pipeline

## [[3.4](https://github.com/nf-core/rnaseq/releases/tag/3.4)] - 2021-10-05

### Enhancements & fixes

- Software version(s) will now be reported for every module imported during a given pipeline execution
- Added `python3` shebang to appropriate scripts in `bin/` directory
- [[#407](https://github.com/nf-core/rnaseq/issues/407)] - Filter mouse reads from PDX samples
- [[#570](https://github.com/nf-core/rnaseq/issues/570)] - Update SortMeRNA to use SilvaDB 138 (for commercial use)
- [[#690](https://github.com/nf-core/rnaseq/issues/690)] - Error with post-trimmed read 2 sample names from FastQC in MultiQC
- [[#693](https://github.com/nf-core/rnaseq/issues/693)] - Cutadapt version missing from MultiQC report
- [[#697](https://github.com/nf-core/rnaseq/issues/697)] - pipeline_report.{txt,html} missing from pipeline_info directory
- [[#705](https://github.com/nf-core/rnaseq/issues/705)] - Sample sheet error check false positive

### Parameters

| Old parameter | New parameter          |
| ------------- | ---------------------- |
|               | `--bbsplit_fasta_list` |
|               | `--bbsplit_index`      |
|               | `--save_bbsplit_reads` |
|               | `--skip_bbsplit`       |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency    | Old version | New version |
| ------------- | ----------- | ----------- |
| `bbmap`       |             | 38.93       |
| `hisat2`      | 2.2.0       | 2.2.1       |
| `picard`      | 2.23.9      | 2.25.7      |
| `salmon`      | 1.4.0       | 1.5.2       |
| `samtools`    | 1.12        | 1.13        |
| `sortmerna`   | 4.2.0       | 4.3.4       |
| `trim-galore` | 0.6.6       | 0.6.7       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [[3.3](https://github.com/nf-core/rnaseq/releases/tag/3.3)] - 2021-07-29

### Enhancements & fixes

- Updated pipeline template to [nf-core/tools 2.1](https://github.com/nf-core/tools/releases/tag/2.1)
- [[#556](https://github.com/nf-core/rnaseq/issues/556)] - Genome index is not recreated with --additional_fasta unless --star_index false
- [[#668](https://github.com/nf-core/rnaseq/issues/668)] - Salmon quant with UMI-tools does not work
- [[#674](https://github.com/nf-core/rnaseq/issues/674)] - Launch pipeline regex fails

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| `samtools`  | 1.10        | 1.12        |
| `stringtie` | 2.1.4       | 2.1.7       |
| `umi_tools` | 1.1.1       | 1.1.2       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [[3.2](https://github.com/nf-core/rnaseq/releases/tag/3.2)] - 2021-06-18

### Enhancements & fixes

- Removed workflow to download data from public databases in favour of using [nf-core/fetchngs](https://nf-co.re/fetchngs)
- Added a stand-alone Python script [`bin/fastq_dir_to_samplesheet.py`](https://github.com/nf-core/rnaseq/blob/master/bin/fastq_dir_to_samplesheet.py) to auto-create samplesheet from a directory of FastQ files
- Added docs about overwriting default container definitions to use latest versions e.g. Pangolin
- [[#637](https://github.com/nf-core/rnaseq/issues/637)] - Add `--salmon_quant_libtype` parameter to provide the `--libType` option to salmon quantification
- [[#645](https://github.com/nf-core/rnaseq/issues/645)] - Remove trailing slash from `params.igenomes_base`
- [[#649](https://github.com/nf-core/rnaseq/issues/649)] - DESeq2 fails with only one sample
- [[#652](https://github.com/nf-core/rnaseq/issues/652)] - Results files have incorrect file names
- [[nf-core/viralrecon#201](https://github.com/nf-core/viralrecon/issues/201)] - Conditional include are not expected to work

### Parameters

| Old parameter               | New parameter            |
| --------------------------- | ------------------------ |
| `--public_data_ids`         |                          |
| `--skip_sra_fastq_download` |                          |
|                             | `--salmon_quant_libtype` |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if parameter information isn't present.

## [[3.1](https://github.com/nf-core/rnaseq/releases/tag/3.1)] - 2021-05-13

### :warning: Major enhancements

- Samplesheet format has changed from `group,replicate,fastq_1,fastq_2,strandedness` to `sample,fastq_1,fastq_2,strandedness`
  - This gives users the flexibility to name their samples however they wish (see [#550](https://github.com/nf-core/rnaseq/issues/550]))
  - PCA generated by DESeq2 will now be monochrome and will not be grouped by using the replicate id
- Updated Nextflow version to `v21.04.0` (see [nextflow#572](https://github.com/nextflow-io/nextflow/issues/1964))
- Restructure pipeline scripts into `modules/`, `subworkflows/` and `workflows/` directories

### Enhancements & fixes

- Updated pipeline template to nf-core/tools `1.14`
- Initial implementation of a standardised samplesheet JSON schema to use with user interfaces and for validation
- Only FastQ files that require to be concatenated will be passed to `CAT_FASTQ` process
- [[#449](https://github.com/nf-core/modules/pull/449)] - `--genomeSAindexNbases` will now be auto-calculated before building STAR indices
- [[#460](https://github.com/nf-core/rnaseq/issues/460)] - Auto-detect and bypass featureCounts execution if biotype doesn't exist in GTF
- [[#544](https://github.com/nf-core/rnaseq/issues/544)] - Update test-dataset for pipeline
- [[#553](https://github.com/nf-core/rnaseq/issues/553)] - Make tximport output files using all the samples; identified by @j-andrews7
- [[#561](https://github.com/nf-core/rnaseq/issues/561)] - Add gene symbols to merged output; identified by @grst
- [[#563](https://github.com/nf-core/rnaseq/issues/563)] - samplesheet.csv merge error
- [[#567](https://github.com/nf-core/rnaseq/issues/567)] - Update docs to mention trimgalore core usage nuances
- [[#568](https://github.com/nf-core/rnaseq/issues/568)] - `--star_index` argument is ignored with `--aligner star_rsem` option
- [[#569](https://github.com/nf-core/rnaseq/issues/569)] - nextflow edge release documentation for running 3.0
- [[#575](https://github.com/nf-core/rnaseq/issues/575)] - Remove duplicated salmon output files
- [[#576](https://github.com/nf-core/rnaseq/issues/576)] - umi_tools dedup : Run before salmon to dedup counts
- [[#582](https://github.com/nf-core/rnaseq/issues/582)] - Generate a separate bigwig tracks for each strand
- [[#583](https://github.com/nf-core/rnaseq/issues/583)] - Samtools error during run requires use of BAM CSI index
- [[#585](https://github.com/nf-core/rnaseq/issues/585)] - Clarify salmon uncertainty for some transcripts
- [[#604](https://github.com/nf-core/rnaseq/issues/604)] - Additional fasta with GENCODE annotation results in biotype error
- [[#610](https://github.com/nf-core/rnaseq/issues/610)] - save R objects as RDS
- [[#619](https://github.com/nf-core/rnaseq/issues/619)] - implicit declaration of the workflow in main
- [[#629](https://github.com/nf-core/rnaseq/pull/629)] - Add and fix EditorConfig linting in entire pipeline
- [[nf-core/modules#423](https://github.com/nf-core/modules/pull/423)] - Replace `publish_by_id` module option to `publish_by_meta`
- [[nextflow#2060](https://github.com/nextflow-io/nextflow/issues/2060)] - Pipeline execution hang when native task fail to be submitted

### Parameters

| Old parameter               | New parameter                  |
| --------------------------- | ------------------------------ |
| `--hisat_build_memory`      | `--hisat2_build_memory`        |
| `--gtf_count_type`          | `--featurecounts_feature_type` |
| `--gtf_group_features_type` | `--featurecounts_group_type`   |
|                             | `--bam_csi_index`              |
|                             | `--schema_ignore_params`       |
|                             | `--show_hidden_params`         |
|                             | `--validate_params`            |
| `--clusterOptions`          |                                |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `bedtools` | 2.29.2      | 2.30.0      |
| `multiqc`  | 1.9         | 1.10.1      |
| `preseq`   | 2.0.3       | 3.1.2       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [[3.0](https://github.com/nf-core/rnaseq/releases/tag/3.0)] - 2020-12-15

### :warning: Major enhancements

- You will need to install Nextflow `>=20.11.0-edge` to run the pipeline. If you are using Singularity, then features introduced in that release now enable the pipeline to directly download Singularity images hosted by Biocontainers as opposed to performing a conversion from Docker images (see [#496](https://github.com/nf-core/rnaseq/issues/496)).
- The previous default of aligning BAM files using STAR and quantifying using featureCounts (`--aligner star`) has been removed. The new default is to align with STAR and quantify using Salmon (`--aligner star_salmon`).
  - This decision was made primarily because of the limitations of featureCounts to appropriately quantify gene expression data. Please see [Zhao et al., 2015](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0141910#pone-0141910-t001) and [Soneson et al., 2015](https://f1000research.com/articles/4-1521/v1)).
- For similar reasons, **quantification will not be performed** if using `--aligner hisat2` due to the lack of an appropriate option to calculate accurate expression estimates from HISAT2 derived genomic alignments.
  - This pipeline option is still available for those who have a preference for the alignment, QC and other types of downstream analysis compatible with the output of HISAT2. No gene-level quantification results will be generated.
  - In a future release we hope to add back quantitation for HISAT2 using different tools.

### Enhancements & fixes

- Updated pipeline template to nf-core/tools `1.12.1`
- Bumped Nextflow version `20.07.1` -> `20.11.0-edge`
- Added UCSC `bedClip` module to restrict bedGraph file coordinates to chromosome boundaries
- Check if Bioconda and conda-forge channels are set-up correctly when running with `-profile conda`
- Use `rsem-prepare-reference` and not `gffread` to create transcriptome fasta file
- [[#494](https://github.com/nf-core/rnaseq/issues/494)] - Issue running rnaseq v2.0 (DSL2) with test profile
- [[#496](https://github.com/nf-core/rnaseq/issues/496)] - Direct download of Singularity images via HTTPS
- [[#498](https://github.com/nf-core/rnaseq/issues/498)] - Significantly different versions of STAR in star_rsem (2.7.6a) and star (2.6.1d)
- [[#499](https://github.com/nf-core/rnaseq/issues/499)] - Use of salmon counts for DESeq2
- [[#500](https://github.com/nf-core/rnaseq/issues/500), [#509](https://github.com/nf-core/rnaseq/issues/509)] - Error with AWS batch params
- [[#511](https://github.com/nf-core/rnaseq/issues/511)] - rsem/star index fails with large genome
- [[#515](https://github.com/nf-core/rnaseq/issues/515)] - Add decoy-aware indexing for salmon
- [[#516](https://github.com/nf-core/rnaseq/issues/516)] - Unexpected error [InvocationTargetException]
- [[#525](https://github.com/nf-core/rnaseq/issues/525)] - sra_ids_to_runinfo.py UnicodeEncodeError
- [[#550](https://github.com/nf-core/rnaseq/issues/525)] - handle samplesheets with replicate=0

### Parameters

| Old parameter               | New parameter                          |
| --------------------------- | -------------------------------------- |
| `--fc_extra_attributes`     | `--gtf_extra_attributes`               |
|  `--fc_group_features`      |  `--gtf_group_features`                |
|  `--fc_count_type`          |  `--gtf_count_type`                    |
|  `--fc_group_features_type` |  `--gtf_group_features_type`           |
|                             |  `--singularity_pull_docker_container` |
|  `--skip_featurecounts`     |                                        |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency                          | Old version | New version |
| ----------------------------------- | ----------- | ----------- |
| `bioconductor-summarizedexperiment` | 1.18.1      | 1.20.0      |
| `bioconductor-tximeta`              | 1.6.3       | 1.8.0       |
| `picard`                            | 2.23.8      | 2.23.9      |
| `requests`                          |             | 2.24.0      |
| `salmon`                            | 1.3.0       | 1.4.0       |
| `ucsc-bedclip`                      |             | 377         |
| `umi_tools`                         | 1.0.1       | 1.1.1       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [[2.0](https://github.com/nf-core/rnaseq/releases/tag/2.0)] - 2020-11-12

### Major enhancements

- Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
- All software containers are now exclusively obtained from [Biocontainers](https://biocontainers.pro/#/registry)
- Added a separate workflow to download FastQ files via SRA, ENA or GEO ids and to auto-create the input samplesheet ([`ENA FTP`](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html); see [`--public_data_ids`](https://nf-co.re/rnaseq/parameters#public_data_ids) parameter)
- Added and refined a Groovy `lib/` of functions that include the automatic rendering of parameters defined in the JSON schema for the help and summary log information
- Replace [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) for the generation of PCA and heatmaps (also included in the MultiQC report)
- Creation of bigWig coverage files using [BEDTools](https://github.com/arq5x/bedtools2/) and [bedGraphToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/)
- [[#70](https://github.com/nf-core/rnaseq/issues/70)] - Added new genome mapping and quantification route with [RSEM](https://github.com/deweylab/RSEM) via the `--aligner star_rsem` parameter
- [[#72](https://github.com/nf-core/rnaseq/issues/72)] - Samples skipped due to low alignment reported in the MultiQC report
- [[#73](https://github.com/nf-core/rnaseq/issues/73), [#435](https://github.com/nf-core/rnaseq/pull/435)] - UMI barcode support
- [[#91](https://github.com/nf-core/rnaseq/issues/91)] - Ability to concatenate multiple runs of the same samples via the input samplesheet
- [[#123](https://github.com/nf-core/rnaseq/issues/123)] - The primary input for the pipeline has changed from `--reads` glob to samplesheet `--input`. See [usage docs](https://nf-co.re/rnaseq/docs/usage#introduction).
- [[#197](https://github.com/nf-core/rnaseq/issues/197)] - Samples failing strand-specificity checks reported in the MultiQC report
- [[#227](https://github.com/nf-core/rnaseq/issues/227)] - Removal of ribosomal RNA via [SortMeRNA](https://github.com/biocore/sortmerna)
- [[#419](https://github.com/nf-core/rnaseq/pull/419)] - Add `--additional_fasta` parameter to provide ERCC spike-ins, transgenes such as GFP or CAR-T as additional sequences to align to

### Other enhancements & fixes

- Updated pipeline template to nf-core/tools `1.11`
- Optimise MultiQC configuration for faster run-time on huge sample numbers
- Add information about SILVA licensing when removing rRNA to `usage.md`
- Fixed ansi colours for pipeline summary, added summary logs of alignment results
- [[#281](https://github.com/nf-core/rnaseq/issues/281)] - Add nag to cite the pipeline in summary
- [[#302](https://github.com/nf-core/rnaseq/issues/302)] - Fixed MDS plot axis labels
- [[#338](https://github.com/nf-core/rnaseq/issues/338)] - Add option for turning on/off STAR command line option (--sjdbGTFfile)
- [[#344](https://github.com/nf-core/rnaseq/issues/344)] - Added multi-core TrimGalore support
- [[#351](https://github.com/nf-core/rnaseq/issues/351)] - Fixes missing Qualimap parameter `-p`
- [[#353](https://github.com/nf-core/rnaseq/issues/353)] - Fixes an issue where MultiQC fails to run with `--skip_biotype_qc` option
- [[#357](https://github.com/nf-core/rnaseq/issues/357)] - Fixes broken links
- [[#362](https://github.com/nf-core/rnaseq/issues/362)] - Fix error with gzipped annotation file
- [[#384](https://github.com/nf-core/rnaseq/issues/384)] - Changed SortMeRNA reference dbs path to use stable URLs (v4.2.0)
- [[#396](https://github.com/nf-core/rnaseq/issues/396)] - Deterministic mapping for STAR aligner
- [[#412](https://github.com/nf-core/rnaseq/issues/412)] - Fix Qualimap not being passed on correct strand-specificity parameter
- [[#413](https://github.com/nf-core/rnaseq/issues/413)] - Fix STAR unmapped reads not output
- [[#434](https://github.com/nf-core/rnaseq/issues/434)] - Fix typo reported for work-dir
- [[#437](https://github.com/nf-core/rnaseq/issues/434)] - FastQC uses correct number of threads now
- [[#440](https://github.com/nf-core/rnaseq/issues/440)] - Fixed issue where featureCounts process fails when setting `--fc_count_type` to gene
- [[#452](https://github.com/nf-core/rnaseq/issues/452)] - Fix `--gff` input bug
- [[#345](https://github.com/nf-core/rnaseq/pull/345)] - Fixes label name in FastQC process
- [[#391](https://github.com/nf-core/rnaseq/pull/391)] - Make publishDir mode configurable
- [[#431](https://github.com/nf-core/rnaseq/pull/431)] - Update AWS GitHub actions workflow with organization level secrets
- [[#435](https://github.com/nf-core/rnaseq/pull/435)] - Fix a bug where gzipped references were not extracted when `--additional_fasta` was not specified
- [[#435](https://github.com/nf-core/rnaseq/pull/435)] - Fix a bug where merging of RSEM output would fail if only one fastq provided as input
- [[#435](https://github.com/nf-core/rnaseq/pull/435)] - Correct RSEM output name (was saving counts but calling them TPMs; now saving both properly labelled)
- [[#436](https://github.com/nf-core/rnaseq/pull/436)] - Fix a bug where the RSEM reference could not be built
- [[#458](https://github.com/nf-core/rnaseq/pull/458)] - Fix `TMP_DIR` for process MarkDuplicates and Qualimap

### Parameters

#### Updated

| Old parameter                 | New parameter               |
| ----------------------------- | --------------------------- |
| `--reads`                     | `--input`                   |
|  `--igenomesIgnore`           |  `--igenomes_ignore`        |
|  `--removeRiboRNA`            |  `--remove_ribo_rna`        |
|  `--rRNA_database_manifest`   |  `--ribo_database_manifest` |
|  `--save_nonrRNA_reads`       |  `--save_non_ribo_reads`    |
|  `--saveAlignedIntermediates` |  `--save_align_intermeds`   |
|  `--saveReference`            |  `--save_reference`         |
|  `--saveTrimmed`              |  `--save_trimmed`           |
|  `--saveUnaligned`            |  `--save_unaligned`         |
|  `--skipAlignment`            |  `--skip_alignment`         |
|  `--skipBiotypeQC`            |  `--skip_biotype_qc`        |
|  `--skipDupRadar`             |  `--skip_dupradar`          |
|  `--skipFastQC`               |  `--skip_fastqc`            |
|  `--skipMultiQC`              |  `--skip_multiqc`           |
|  `--skipPreseq`               |  `--skip_preseq`            |
|  `--skipQC`                   |  `--skip_qc`                |
|  `--skipQualimap`             |  `--skip_qualimap`          |
|  `--skipRseQC`                |  `--skip_rseqc`             |
|  `--skipTrimming`             |  `--skip_trimming`          |
|  `--stringTieIgnoreGTF`       |  `--stringtie_ignore_gtf`   |

#### Added

- `--additional_fasta` - FASTA file to concatenate to genome FASTA file e.g. containing spike-in sequences
- `--deseq2_vst` - Use vst transformation instead of rlog with DESeq2
- `--enable_conda` - Run this workflow with Conda. You can also use '-profile conda' instead of providing this parameter
- `--min_mapped_reads` - Minimum percentage of uniquely mapped reads below which samples are removed from further processing
- `--multiqc_title` - MultiQC report title. Printed as page header, used for filename if not otherwise specified
- `--public_data_ids` - File containing SRA/ENA/GEO identifiers one per line in order to download their associated FastQ files
- `--publish_dir_mode` - Method used to save pipeline results to output directory
- `--rsem_index` - Path to directory or tar.gz archive for pre-built RSEM index
- `--rseqc_modules` - Specify the RSeQC modules to run
- `--save_merged_fastq` - Save FastQ files after merging re-sequenced libraries in the results directory
- `--save_umi_intermeds` - If this option is specified, intermediate FastQ and BAM files produced by UMI-tools are also saved in the results directory
- `--skip_bigwig` - Skip bigWig file creation
- `--skip_deseq2_qc` - Skip DESeq2 PCA and heatmap plotting
- `--skip_featurecounts` - Skip featureCounts
- `--skip_markduplicates` - Skip picard MarkDuplicates step
- `--skip_sra_fastq_download` - Only download metadata for public data database ids and don't download the FastQ files
- `--skip_stringtie` - Skip StringTie
- `--star_ignore_sjdbgtf` - See [#338](https://github.com/nf-core/rnaseq/issues/338)
- `--umitools_bc_pattern` - The UMI barcode pattern to use e.g. 'NNNNNN' indicates that the first 6 nucleotides of the read are from the UMI
- `--umitools_extract_method` - UMI pattern to use. Can be either 'string' (default) or 'regex'
- `--with_umi` - Enable UMI-based read deduplication

#### Removed

- `--awsqueue` can now be provided via nf-core/configs if using AWS
- `--awsregion` can now be provided via nf-core/configs if using AWS
- `--compressedReference` now auto-detected
- `--markdup_java_options` in favour of updating centrally on nf-core/modules
- `--project` parameter from old NGI template
- `--readPaths` is not required since these are provided from the input samplesheet
- `--sampleLevel` not required
- `--singleEnd` is now auto-detected from the input samplesheet
- `--skipEdgeR` qc not performed by DESeq2 instead
- `--star_memory` in favour of updating centrally on nf-core/modules if required
- Strandedness is now specified at the sample-level via the input samplesheet
  - `--forwardStranded`
  - `--reverseStranded`
  - `--unStranded`
  - `--pico`

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency                          | Old version | New version |
| ----------------------------------- | ----------- | ----------- |
| `bioconductor-dupradar`             | 1.14.0      | 1.18.0      |
| `bioconductor-summarizedexperiment` | 1.14.0      | 1.18.1      |
| `bioconductor-tximeta`              | 1.2.2       | 1.6.3       |
| `fastqc`                            | 0.11.8      | 0.11.9      |
| `gffread`                           | 0.11.4      | 0.12.1      |
| `hisat2`                            | 2.1.0       | 2.2.0       |
| `multiqc`                           | 1.7         | 1.9         |
| `picard`                            | 2.21.1      | 2.23.8      |
| `qualimap`                          | 2.2.2c      | 2.2.2d      |
| `r-base`                            | 3.6.1       | 4.0.3       |
| `salmon`                            | 0.14.2      | 1.3.0       |
| `samtools`                          | 1.9         | 1.10        |
| `sortmerna`                         | 2.1b        | 4.2.0       |
| `stringtie`                         | 2.0         | 2.1.4       |
| `subread`                           | 1.6.4       | 2.0.1       |
| `trim-galore`                       | 0.6.4       | 0.6.6       |
| `bedtools`                          | -           | 2.29.2      |
| `bioconductor-biocparallel`         | -           | 1.22.0      |
| `bioconductor-complexheatmap`       | -           | 2.4.2       |
| `bioconductor-deseq2`               | -           | 1.28.0      |
| `bioconductor-tximport`             | -           | 1.16.0      |
| `perl`                              | -           | 5.26.2      |
| `python`                            | -           | 3.8.3       |
| `r-ggplot2`                         | -           | 3.3.2       |
| `r-optparse`                        | -           | 1.6.6       |
| `r-pheatmap`                        | -           | 1.0.12      |
| `r-rcolorbrewer`                    | -           | 1.1_2       |
| `rsem`                              | -           | 1.3.3       |
| `ucsc-bedgraphtobigwig`             | -           | 377         |
| `umi_tools`                         | -           | 1.0.1       |
| `bioconductor-edger`                | -           | -           |
| `deeptools`                         | -           | -           |
| `matplotlib`                        | -           | -           |
| `r-data.table`                      | -           | -           |
| `r-gplots`                          | -           | -           |
| `r-markdown`                        | -           | -           |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [[1.4.2](https://github.com/nf-core/rnaseq/releases/tag/1.4.2)] - 2019-10-18

- Minor version release for keeping Git History in sync
- No changes with respect to 1.4.1 on pipeline level

## [[1.4.1](https://github.com/nf-core/rnaseq/releases/tag/1.4.1)] - 2019-10-17

Major novel changes include:

- Update `igenomes.config` with NCBI `GRCh38` and most recent UCSC genomes
- Set `autoMounts = true` by default for `singularity` profile

### Pipeline enhancements & fixes

- Fixed parameter warnings [#316](https://github.com/nf-core/rnaseq/issues/316) and [318](https://github.com/nf-core/rnaseq/issues/318)
- Fixed [#307](https://github.com/nf-core/rnaseq/issues/307) - Confusing Info Printout about GFF and GTF

## [[1.4](https://github.com/nf-core/rnaseq/releases/tag/1.4)] - 2019-10-15

Major novel changes include:

- Support for Salmon as an alternative method to STAR and HISAT2
- Several improvements in `featureCounts` handling of types other than `exon`. It is possible now to handle nuclearRNAseq data. Nuclear RNA has un-spliced RNA, and the whole transcript, including the introns, needs to be counted, e.g. by specifying `--fc_count_type transcript`.
- Support for [outputting unaligned data](https://github.com/nf-core/rnaseq/issues/277) to results folders.
- Added options to skip several steps
  - Skip trimming using `--skipTrimming`
  - Skip BiotypeQC using `--skipBiotypeQC`
  - Skip Alignment using `--skipAlignment` to only use pseudoalignment using Salmon

### Documentation updates

- Adjust wording of skipped samples [in pipeline output](https://github.com/nf-core/rnaseq/issues/290)
- Fixed link to guidelines [#203](https://github.com/nf-core/rnaseq/issues/203)
- Add `Citation` and `Quick Start` section to `README.md`
- Add in documentation of the `--gff` parameter

### Reporting Updates

- Generate MultiQC plots in the results directory [#200](https://github.com/nf-core/rnaseq/issues/200)
- Get MultiQC to save plots as [standalone files](https://github.com/nf-core/rnaseq/issues/183)
- Get MultiQC to write out the software versions in a `.csv` file [#185](https://github.com/nf-core/rnaseq/issues/185)
- Use `file` instead of `new File` to create `pipeline_report.{html,txt}` files, and properly create subfolders

### Pipeline enhancements & fixes

- Restore `SummarizedExperimment` object creation in the salmon_merge process avoiding increasing memory with sample size.
- Fix sample names in feature counts and dupRadar to remove suffixes added in other processes
- Removed `genebody_coverage` process [#195](https://github.com/nf-core/rnaseq/issues/195)
- Implemented Pearsons correlation instead of Euclidean distance [#146](https://github.com/nf-core/rnaseq/issues/146)
- Add `--stringTieIgnoreGTF` parameter [#206](https://github.com/nf-core/rnaseq/issues/206)
- Removed unused `stringtie` channels for `MultiQC`
- Integrate changes in `nf-core/tools v1.6` template which resolved [#90](https://github.com/nf-core/rnaseq/issues/90)
- Moved process `convertGFFtoGTF` before `makeSTARindex` [#215](https://github.com/nf-core/rnaseq/issues/215)
- Change all boolean parameters from `snake_case` to `camelCase` and vice versa for value parameters
- Add SM ReadGroup info for QualiMap compatibility[#238](https://github.com/nf-core/rnaseq/issues/238)
- Obtain edgeR + dupRadar version information [#198](https://github.com/nf-core/rnaseq/issues/198) and [#112](https://github.com/nf-core/rnaseq/issues/112)
- Add `--gencode` option for compatibility of Salmon and featureCounts biotypes with GENCODE gene annotations
- Added functionality to accept compressed reference data in the pipeline
- Check that gtf features are on chromosomes that exist in the genome fasta file [#274](https://github.com/nf-core/rnaseq/pull/274)
- Maintain all gff features upon gtf conversion (keeps `gene_biotype` or `gene_type` to make `featureCounts` happy)
- Add SortMeRNA as an optional step to allow rRNA removal [#280](https://github.com/nf-core/rnaseq/issues/280)
- Minimal adjustment of memory and CPU constraints for clusters with locked memory / CPU relation
- Cleaned up usage, `parameters.settings.json` and the `nextflow.config`

### Dependency Updates

- Dependency list is now sorted appropriately
- Force matplotlib=3.0.3

#### Updated Packages

- Picard 2.20.0 -> 2.21.1
- bioconductor-dupradar 1.12.1 -> 1.14.0
- bioconductor-edger 3.24.3 -> 3.26.5
- gffread 0.9.12 -> 0.11.4
- trim-galore 0.6.1 -> 0.6.4
- gffread 0.9.12 -> 0.11.4
- rseqc 3.0.0 -> 3.0.1
- R-Base 3.5 -> 3.6.1

#### Added / Removed Packages

- Dropped CSVtk in favor of Unix's simple `cut` and `paste` utilities
- Added Salmon 0.14.2
- Added TXIMeta 1.2.2
- Added SummarizedExperiment 1.14.0
- Added SortMeRNA 2.1b
- Add tximport and summarizedexperiment dependency [#171](https://github.com/nf-core/rnaseq/issues/171)
- Add Qualimap dependency [#202](https://github.com/nf-core/rnaseq/issues/202)

## [[1.3](https://github.com/nf-core/rnaseq/releases/tag/1.3)] - 2019-03-26

### Pipeline Updates

- Added configurable options to specify group attributes for featureCounts [#144](https://github.com/nf-core/rnaseq/issues/144)
- Added support for RSeqC 3.0 [#148](https://github.com/nf-core/rnaseq/issues/148)
- Added a `parameters.settings.json` file for use with the new `nf-core launch` helper tool.
- Centralized all configuration profiles using [nf-core/configs](https://github.com/nf-core/configs)
- Fixed all centralized configs [for offline usage](https://github.com/nf-core/rnaseq/issues/163)
- Hide %dup in [multiqc report](https://github.com/nf-core/rnaseq/issues/150)
- Add option for Trimming NextSeq data properly ([@jburos work](https://github.com/jburos))

### Bug fixes

- Fixing HISAT2 Index Building for large reference genomes [#153](https://github.com/nf-core/rnaseq/issues/153)
- Fixing HISAT2 BAM sorting using more memory than available on the system
- Fixing MarkDuplicates memory consumption issues following [#179](https://github.com/nf-core/rnaseq/pull/179)
- Use `file` instead of `new File` to create the `pipeline_report.{html,txt}` files to avoid creating local directories when outputting to AWS S3 folders
- Fix SortMeRNA default rRNA db paths specified in assets/rrna-db-defaults.txt

### Dependency Updates

- RSeQC 2.6.4 -> 3.0.0
- Picard 2.18.15 -> 2.20.0
- r-data.table 1.11.4 -> 1.12.2
- bioconductor-edger 3.24.1 -> 3.24.3
- r-markdown 0.8 -> 0.9
- csvtk 0.15.0 -> 0.17.0
- stringtie 1.3.4 -> 1.3.6
- subread 1.6.2 -> 1.6.4
- gffread 0.9.9 -> 0.9.12
- multiqc 1.6 -> 1.7
- deeptools 3.2.0 -> 3.2.1
- trim-galore 0.5.0 -> 0.6.1
- qualimap 2.2.2b
- matplotlib 3.0.3
- r-base 3.5.1

## [[1.2](https://github.com/nf-core/rnaseq/releases/tag/1.2)] - 2018-12-12

### Pipeline updates

- Removed some outdated documentation about non-existent features
- Config refactoring and code cleaning
- Added a `--fcExtraAttributes` option to specify more than ENSEMBL gene names in `featureCounts`
- Remove legacy rseqc `strandRule` config code. [#119](https://github.com/nf-core/rnaseq/issues/119)
- Added STRINGTIE ballgown output to results folder [#125](https://github.com/nf-core/rnaseq/issues/125)
- HiSAT index build now requests `200GB` memory, enough to use the exons / splice junction option for building.
  - Added documentation about the `--hisatBuildMemory` option.
- BAM indices are stored and re-used between processes [#71](https://github.com/nf-core/rnaseq/issues/71)

### Bug Fixes

- Fixed conda bug which caused problems with environment resolution due to changes in bioconda [#113](https://github.com/nf-core/rnaseq/issues/113)
- Fixed wrong gffread command line [#117](https://github.com/nf-core/rnaseq/issues/117)
- Added `cpus = 1` to `workflow summary process` [#130](https://github.com/nf-core/rnaseq/issues/130)

## [[1.1](https://github.com/nf-core/rnaseq/releases/tag/1.1)] - 2018-10-05

### Pipeline updates

- Wrote docs and made minor tweaks to the `--skip_qc` and associated options
- Removed the depreciated `uppmax-modules` config profile
- Updated the `hebbe` config profile to use the new `withName` syntax too
- Use new `workflow.manifest` variables in the pipeline script
- Updated minimum nextflow version to `0.32.0`

### Bug Fixes

- [#77](https://github.com/nf-core/rnaseq/issues/77): Added back `executor = 'local'` for the `workflow_summary_mqc`
- [#95](https://github.com/nf-core/rnaseq/issues/95): Check if task.memory is false instead of null
- [#97](https://github.com/nf-core/rnaseq/issues/97): Resolved edge-case where numeric sample IDs are parsed as numbers causing some samples to be incorrectly overwritten.

## [[1.0](https://github.com/nf-core/rnaseq/releases/tag/1.0)] - 2018-08-20

This release marks the point where the pipeline was moved from [SciLifeLab/NGI-RNAseq](https://github.com/SciLifeLab/NGI-RNAseq)
over to the new [nf-core](http://nf-co.re/) community, at [nf-core/rnaseq](https://github.com/nf-core/rnaseq).

View the previous changelog at [SciLifeLab/NGI-RNAseq/CHANGELOG.md](https://github.com/SciLifeLab/NGI-RNAseq/blob/master/CHANGELOG.md)

In addition to porting to the new nf-core community, the pipeline has had a number of major changes in this version.
There have been 157 commits by 16 different contributors covering 70 different files in the pipeline: 7,357 additions and 8,236 deletions!

In summary, the main changes are:

- Rebranding and renaming throughout the pipeline to nf-core
- Updating many parts of the pipeline config and style to meet nf-core standards
- Support for GFF files in addition to GTF files
  - Just use `--gff` instead of `--gtf` when specifying a file path
- New command line options to skip various quality control steps
- More safety checks when launching a pipeline
  - Several new sanity checks - for example, that the specified reference genome exists
- Improved performance with memory usage (especially STAR and Picard)
- New BigWig file outputs for plotting coverage across the genome
- Refactored gene body coverage calculation, now much faster and using much less memory
- Bugfixes in the MultiQC process to avoid edge cases where it wouldn't run
- MultiQC report now automatically attached to the email sent when the pipeline completes
- New testing method, with data on GitHub
  - Now run pipeline with `-profile test` instead of using bash scripts
- Rewritten continuous integration tests with Travis CI
- New explicit support for Singularity containers
- Improved MultiQC support for DupRadar and featureCounts
  - Now works for all users instead of just NGI Stockholm
- New configuration for use on AWS batch
- Updated config syntax to support latest versions of Nextflow
- Built-in support for a number of new local HPC systems
  - CCGA, GIS, UCT HEX, updates to UPPMAX, CFC, BINAC, Hebbe, c3se
- Slightly improved documentation (more updates to come)
- Updated software packages

...and many more minor tweaks.

Thanks to everyone who has worked on this release!
