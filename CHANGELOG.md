# nf-core/rnaseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.1dev - [date]

### Enhancements & fixes

* Updated pipeline template to nf-core/tools `1.12`
* [[#500](https://github.com/nf-core/rnaseq/issues/500), [#509](https://github.com/nf-core/rnaseq/issues/509)] - Error with AWS batch params

## [[2.0](https://github.com/nf-core/rnaseq/releases/tag/2.0)] - 2020-11-12

### Major enhancements

* Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
* All software containers are now exclusively obtained from [Biocontainers](https://biocontainers.pro/#/registry)
* Added a separate workflow to download FastQ files via SRA, ENA or GEO ids and to auto-create the input samplesheet ([`ENA FTP`](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html); see [`--public_data_ids`](https://nf-co.re/rnaseq/parameters#public_data_ids) parameter)
* Added and refined a Groovy `lib/` of functions that include the automatic rendering of parameters defined in the JSON schema for the help and summary log information
* Replace [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) for the generation of PCA and heatmaps (also included in the MultiQC report)
* Creation of bigWig coverage files using [BEDTools](https://github.com/arq5x/bedtools2/) and [bedGraphToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/)
* [[#70](https://github.com/nf-core/rnaseq/issues/70)] - Added new genome mapping and quantification route with [RSEM](https://github.com/deweylab/RSEM) via the `--aligner star_rsem` parameter
* [[#72](https://github.com/nf-core/rnaseq/issues/72)] - Samples skipped due to low alignment reported in the MultiQC report
* [[#73](https://github.com/nf-core/rnaseq/issues/73), [#435](https://github.com/nf-core/rnaseq/pull/435)] - UMI barcode support
* [[#91](https://github.com/nf-core/rnaseq/issues/91)] - Ability to concatenate multiple runs of the same samples via the input samplesheet
* [[#123](https://github.com/nf-core/rnaseq/issues/123)] - The primary input for the pipeline has changed from `--reads` glob to samplesheet `--input`. See [usage docs](https://nf-co.re/rnaseq/docs/usage#introduction).
* [[#197](https://github.com/nf-core/rnaseq/issues/197)] - Samples failing strand-specificity checks reported in the MultiQC report
* [[#227](https://github.com/nf-core/rnaseq/issues/227)] - Removal of ribosomal RNA via [SortMeRNA](https://github.com/biocore/sortmerna)
* [[#419](https://github.com/nf-core/rnaseq/pull/419)] - Add `--additional_fasta` parameter to provide ERCC spike-ins, transgenes such as GFP or CAR-T as additional sequences to align to

### Other enhancements & fixes

* Updated pipeline template to nf-core/tools `1.11`
* Optimise MultiQC configuration for faster run-time on huge sample numbers
* Add information about SILVA licensing when removing rRNA to `usage.md`
* Fixed ansi colours for pipeline summary, added summary logs of alignment results
* [[#281](https://github.com/nf-core/rnaseq/issues/281)] - Add nag to cite the pipeline in summary
* [[#302](https://github.com/nf-core/rnaseq/issues/302)] - Fixed MDS plot axis labels
* [[#338](https://github.com/nf-core/rnaseq/issues/338)] - Add option for turning on/off STAR command line option (--sjdbGTFfile)
* [[#344](https://github.com/nf-core/rnaseq/issues/344)] - Added multi-core TrimGalore support
* [[#351](https://github.com/nf-core/rnaseq/issues/351)] - Fixes missing Qualimap parameter `-p`
* [[#353](https://github.com/nf-core/rnaseq/issues/353)] - Fixes an issue where MultiQC fails to run with `--skip_biotype_qc` option
* [[#357](https://github.com/nf-core/rnaseq/issues/357)] - Fixes broken links
* [[#362](https://github.com/nf-core/rnaseq/issues/362)] - Fix error with gzipped annotation file
* [[#384](https://github.com/nf-core/rnaseq/issues/384)] - Changed SortMeRNA reference dbs path to use stable URLs (v4.2.0)
* [[#396](https://github.com/nf-core/rnaseq/issues/396)] - Deterministic mapping for STAR aligner
* [[#412](https://github.com/nf-core/rnaseq/issues/412)] - Fix Qualimap not being passed on correct strand-specificity parameter
* [[#413](https://github.com/nf-core/rnaseq/issues/413)] - Fix STAR unmapped reads not output
* [[#434](https://github.com/nf-core/rnaseq/issues/434)] - Fix typo reported for work-dir
* [[#437](https://github.com/nf-core/rnaseq/issues/434)] - FastQC uses correct number of threads now
* [[#440](https://github.com/nf-core/rnaseq/issues/440)] - Fixed issue where featureCounts process fails when setting `--fc_count_type` to gene
* [[#452](https://github.com/nf-core/rnaseq/issues/452)] - Fix `--gff` input bug
* [[#345](https://github.com/nf-core/rnaseq/pull/345)] - Fixes label name in FastQC process
* [[#391](https://github.com/nf-core/rnaseq/pull/391)] - Make publishDir mode configurable
* [[#431](https://github.com/nf-core/rnaseq/pull/431)] - Update AWS GitHub actions workflow with organization level secrets
* [[#435](https://github.com/nf-core/rnaseq/pull/435)] - Fix a bug where gzipped references were not extracted when `--additional_fasta` was not specified
* [[#435](https://github.com/nf-core/rnaseq/pull/435)] - Fix a bug where merging of RSEM output would fail if only one fastq provided as input
* [[#435](https://github.com/nf-core/rnaseq/pull/435)] - Correct RSEM output name (was saving counts but calling them TPMs; now saving both properly labelled)
* [[#436](https://github.com/nf-core/rnaseq/pull/436)] - Fix a bug where the RSEM reference could not be built
* [[#458](https://github.com/nf-core/rnaseq/pull/458)] - Fix `TMP_DIR` for process MarkDuplicates and Qualimap

### Parameters

#### Updated

| Old parameter                | New parameter              |
|------------------------------|----------------------------|
| `--reads`                    | `--input`                  |
| `--igenomesIgnore`           | `--igenomes_ignore`        |
| `--removeRiboRNA`            | `--remove_ribo_rna`        |
| `--rRNA_database_manifest`   | `--ribo_database_manifest` |
| `--save_nonrRNA_reads`       | `--save_non_ribo_reads`    |
| `--saveAlignedIntermediates` | `--save_align_intermeds`   |
| `--saveReference`            | `--save_reference`         |
| `--saveTrimmed`              | `--save_trimmed`           |
| `--saveUnaligned`            | `--save_unaligned`         |
| `--skipAlignment`            | `--skip_alignment`         |
| `--skipBiotypeQC`            | `--skip_biotype_qc`        |
| `--skipDupRadar`             | `--skip_dupradar`          |
| `--skipFastQC`               | `--skip_fastqc`            |
| `--skipMultiQC`              | `--skip_multiqc`           |
| `--skipPreseq`               | `--skip_preseq`            |
| `--skipQC`                   | `--skip_qc`                |
| `--skipQualimap`             | `--skip_qualimap`          |
| `--skipRseQC`                | `--skip_rseqc`             |
| `--skipTrimming`             | `--skip_trimming`          |
| `--stringTieIgnoreGTF`       | `--stringtie_ignore_gtf`   |

#### Added

* `--additional_fasta` - FASTA file to concatenate to genome FASTA file e.g. containing spike-in sequences
* `--deseq2_vst` - Use vst transformation instead of rlog with DESeq2
* `--enable_conda` - Run this workflow with Conda. You can also use '-profile conda' instead of providing this parameter
* `--min_mapped_reads` - Minimum percentage of uniquely mapped reads below which samples are removed from further processing
* `--multiqc_title` - MultiQC report title. Printed as page header, used for filename if not otherwise specified
* `--public_data_ids` - File containing SRA/ENA/GEO identifiers one per line in order to download their associated FastQ files
* `--publish_dir_mode` - Method used to save pipeline results to output directory
* `--rsem_index` - Path to directory or tar.gz archive for pre-built RSEM index
* `--rseqc_modules` - Specify the RSeQC modules to run
* `--save_merged_fastq` - Save FastQ files after merging re-sequenced libraries in the results directory
* `--save_umi_intermeds` - If this option is specified, intermediate FastQ and BAM files produced by UMI-tools are also saved in the results directory
* `--skip_bigwig` - Skip bigWig file creation
* `--skip_deseq2_qc` - Skip DESeq2 PCA and heatmap plotting
* `--skip_featurecounts` - Skip featureCounts
* `--skip_markduplicates` - Skip picard MarkDuplicates step
* `--skip_sra_fastq_download` - Only download metadata for public data database ids and don't download the FastQ files
* `--skip_stringtie` - Skip StringTie
* `--star_ignore_sjdbgtf` - See [#338](https://github.com/nf-core/rnaseq/issues/338)
* `--umitools_bc_pattern` - The UMI barcode pattern to use e.g. 'NNNNNN' indicates that the first 6 nucleotides of the read are from the UMI
* `--umitools_extract_method` - UMI pattern to use. Can be either 'string' (default) or 'regex'
* `--with_umi` - Enable UMI-based read deduplication

#### Removed

* `--awsqueue` can now be provided via nf-core/configs if using AWS
* `--awsregion` can now be provided via nf-core/configs if using AWS
* `--compressedReference` now auto-detected
* `--markdup_java_options` in favour of updating centrally on nf-core/modules
* `--project` parameter from old NGI template
* `--readPaths` is not required since these are provided from the input samplesheet
* `--sampleLevel` not required
* `--singleEnd` is now auto-detected from the input samplesheet
* `--skipEdgeR` qc not performed by DESeq2 instead
* `--star_memory` in favour of updating centrally on nf-core/modules if required
* Strandedness is now specified at the sample-level via the input samplesheet
  * `--forwardStranded`
  * `--reverseStranded`
  * `--unStranded`
  * `--pico`

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency                          | Old version | New version |
|-------------------------------------|-------------|-------------|
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

> **NB:** Dependency has been __updated__ if both old and new version information is present.  
> **NB:** Dependency has been __added__ if just the new version information is present.  
> **NB:** Dependency has been __removed__ if version information isn't present.  

## [[1.4.2](https://github.com/nf-core/rnaseq/releases/tag/1.4.2)] - 2019-10-18

* Minor version release for keeping Git History in sync
* No changes with respect to 1.4.1 on pipeline level

## [[1.4.1](https://github.com/nf-core/rnaseq/releases/tag/1.4.1)] - 2019-10-17

Major novel changes include:

* Update `igenomes.config` with NCBI `GRCh38` and most recent UCSC genomes
* Set `autoMounts = true` by default for `singularity` profile

### Pipeline enhancements & fixes

* Fixed parameter warnings [#316](https://github.com/nf-core/rnaseq/issues/316) and [318](https://github.com/nf-core/rnaseq/issues/318)
* Fixed [#307](https://github.com/nf-core/rnaseq/issues/307) - Confusing Info Printout about GFF and GTF

## [[1.4](https://github.com/nf-core/rnaseq/releases/tag/1.4)] - 2019-10-15

Major novel changes include:

* Support for Salmon as an alternative method to STAR and HISAT2
* Several improvements in `featureCounts` handling of types other than `exon`. It is possible now to handle nuclearRNAseq data. Nuclear RNA has un-spliced RNA, and the whole transcript, including the introns, needs to be counted, e.g. by specifying `--fc_count_type transcript`.
* Support for [outputting unaligned data](https://github.com/nf-core/rnaseq/issues/277) to results folders.
* Added options to skip several steps
  * Skip trimming using `--skipTrimming`
  * Skip BiotypeQC using `--skipBiotypeQC`
  * Skip Alignment using `--skipAlignment` to only use pseudo-alignment using Salmon

### Documentation updates

* Adjust wording of skipped samples [in pipeline output](https://github.com/nf-core/rnaseq/issues/290)
* Fixed link to guidelines [#203](https://github.com/nf-core/rnaseq/issues/203)
* Add `Citation` and `Quick Start` section to `README.md`
* Add in documentation of the `--gff` parameter

### Reporting Updates

* Generate MultiQC plots in the results directory [#200](https://github.com/nf-core/rnaseq/issues/200)
* Get MultiQC to save plots as [standalone files](https://github.com/nf-core/rnaseq/issues/183)
* Get MultiQC to write out the software versions in a `.csv` file [#185](https://github.com/nf-core/rnaseq/issues/185)
* Use `file` instead of `new File` to create `pipeline_report.{html,txt}` files, and properly create subfolders

### Pipeline enhancements & fixes

* Restore `SummarizedExperimment` object creation in the salmon_merge process avoiding increasing memory with sample size.
* Fix sample names in feature counts and dupRadar to remove suffixes added in other processes
* Removed `genebody_coverage` process [#195](https://github.com/nf-core/rnaseq/issues/195)
* Implemented Pearsons correlation instead of Euclidean distance [#146](https://github.com/nf-core/rnaseq/issues/146)
* Add `--stringTieIgnoreGTF` parameter [#206](https://github.com/nf-core/rnaseq/issues/206)
* Removed unused `stringtie` channels for `MultiQC`
* Integrate changes in `nf-core/tools v1.6` template which resolved [#90](https://github.com/nf-core/rnaseq/issues/90)
* Moved process `convertGFFtoGTF` before `makeSTARindex` [#215](https://github.com/nf-core/rnaseq/issues/215)
* Change all boolean parameters from `snake_case` to `camelCase` and vice versa for value parameters
* Add SM ReadGroup info for QualiMap compatibility[#238](https://github.com/nf-core/rnaseq/issues/238)
* Obtain edgeR + dupRadar version information [#198](https://github.com/nf-core/rnaseq/issues/198) and [#112](https://github.com/nf-core/rnaseq/issues/112)
* Add `--gencode` option for compatibility of Salmon and featureCounts biotypes with GENCODE gene annotations
* Added functionality to accept compressed reference data in the pipeline
* Check that gtf features are on chromosomes that exist in the genome fasta file [#274](https://github.com/nf-core/rnaseq/pull/274)
* Maintain all gff features upon gtf conversion (keeps `gene_biotype` or `gene_type` to make `featureCounts` happy)
* Add SortMeRNA as an optional step to allow rRNA removal [#280](https://github.com/nf-core/rnaseq/issues/280)
* Minimal adjustment of memory and CPU constraints for clusters with locked memory / CPU relation
* Cleaned up usage, `parameters.settings.json` and the `nextflow.config`

### Dependency Updates

* Dependency list is now sorted appropriately
* Force matplotlib=3.0.3

#### Updated Packages

* Picard 2.20.0 -> 2.21.1
* bioconductor-dupradar 1.12.1 -> 1.14.0
* bioconductor-edger 3.24.3 -> 3.26.5
* gffread 0.9.12 -> 0.11.4
* trim-galore 0.6.1 -> 0.6.4
* gffread 0.9.12 -> 0.11.4
* rseqc 3.0.0 -> 3.0.1
* R-Base 3.5 -> 3.6.1

#### Added / Removed Packages

* Dropped CSVtk in favor of Unix's simple `cut` and `paste` utilities
* Added Salmon 0.14.2
* Added TXIMeta 1.2.2
* Added SummarizedExperiment 1.14.0
* Added SortMeRNA 2.1b
* Add tximport and summarizedexperiment dependency [#171](https://github.com/nf-core/rnaseq/issues/171)
* Add Qualimap dependency [#202](https://github.com/nf-core/rnaseq/issues/202)

## [[1.3](https://github.com/nf-core/rnaseq/releases/tag/1.3)] - 2019-03-26

### Pipeline Updates

* Added configurable options to specify group attributes for featureCounts [#144](https://github.com/nf-core/rnaseq/issues/144)
* Added support for RSeqC 3.0 [#148](https://github.com/nf-core/rnaseq/issues/148)
* Added a `parameters.settings.json` file for use with the new `nf-core launch` helper tool.
* Centralized all configuration profiles using [nf-core/configs](https://github.com/nf-core/configs)
* Fixed all centralized configs [for offline usage](https://github.com/nf-core/rnaseq/issues/163)
* Hide %dup in [multiqc report](https://github.com/nf-core/rnaseq/issues/150)
* Add option for Trimming NextSeq data properly ([@jburos work](https://github.com/jburos))

### Bug fixes

* Fixing HISAT2 Index Building for large reference genomes [#153](https://github.com/nf-core/rnaseq/issues/153)
* Fixing HISAT2 BAM sorting using more memory than available on the system
* Fixing MarkDuplicates memory consumption issues following [#179](https://github.com/nf-core/rnaseq/pull/179)
* Use `file` instead of `new File` to create the `pipeline_report.{html,txt}` files to avoid creating local directories when outputting to AWS S3 folders
* Fix SortMeRNA default rRNA db paths specified in assets/rrna-db-defaults.txt

### Dependency Updates

* RSeQC 2.6.4 -> 3.0.0
* Picard 2.18.15 -> 2.20.0
* r-data.table 1.11.4 -> 1.12.2
* bioconductor-edger 3.24.1 -> 3.24.3
* r-markdown 0.8 -> 0.9
* csvtk 0.15.0 -> 0.17.0
* stringtie 1.3.4 -> 1.3.6
* subread 1.6.2 -> 1.6.4
* gffread 0.9.9 -> 0.9.12
* multiqc 1.6 -> 1.7
* deeptools 3.2.0 -> 3.2.1
* trim-galore 0.5.0 -> 0.6.1
* qualimap 2.2.2b
* matplotlib 3.0.3
* r-base 3.5.1

## [[1.2](https://github.com/nf-core/rnaseq/releases/tag/1.2)] - 2018-12-12

### Pipeline updates

* Removed some outdated documentation about non-existent features
* Config refactoring and code cleaning
* Added a `--fcExtraAttributes` option to specify more than ENSEMBL gene names in `featureCounts`
* Remove legacy rseqc `strandRule` config code. [#119](https://github.com/nf-core/rnaseq/issues/119)
* Added STRINGTIE ballgown output to results folder [#125](https://github.com/nf-core/rnaseq/issues/125)
* HiSAT index build now requests `200GB` memory, enough to use the exons / splice junction option for building.
  * Added documentation about the `--hisatBuildMemory` option.
* BAM indices are stored and re-used between processes [#71](https://github.com/nf-core/rnaseq/issues/71)

### Bug Fixes

* Fixed conda bug which caused problems with environment resolution due to changes in bioconda [#113](https://github.com/nf-core/rnaseq/issues/113)
* Fixed wrong gffread command line [#117](https://github.com/nf-core/rnaseq/issues/117)
* Added `cpus = 1` to `workflow summary process` [#130](https://github.com/nf-core/rnaseq/issues/130)

## [[1.1](https://github.com/nf-core/rnaseq/releases/tag/1.1)] - 2018-10-05

### Pipeline updates

* Wrote docs and made minor tweaks to the `--skip_qc` and associated options
* Removed the depreciated `uppmax-modules` config profile
* Updated the `hebbe` config profile to use the new `withName` syntax too
* Use new `workflow.manifest` variables in the pipeline script
* Updated minimum nextflow version to `0.32.0`

### Bug Fixes

* [#77](https://github.com/nf-core/rnaseq/issues/77): Added back `executor = 'local'` for the `workflow_summary_mqc`
* [#95](https://github.com/nf-core/rnaseq/issues/95): Check if task.memory is false instead of null
* [#97](https://github.com/nf-core/rnaseq/issues/97): Resolved edge-case where numeric sample IDs are parsed as numbers causing some samples to be incorrectly overwritten.

## [[1.0](https://github.com/nf-core/rnaseq/releases/tag/1.0)] - 2018-08-20

This release marks the point where the pipeline was moved from [SciLifeLab/NGI-RNAseq](https://github.com/SciLifeLab/NGI-RNAseq)
over to the new [nf-core](http://nf-co.re/) community, at [nf-core/rnaseq](https://github.com/nf-core/rnaseq).

View the previous changelog at [SciLifeLab/NGI-RNAseq/CHANGELOG.md](https://github.com/SciLifeLab/NGI-RNAseq/blob/master/CHANGELOG.md)

In addition to porting to the new nf-core community, the pipeline has had a number of major changes in this version.
There have been 157 commits by 16 different contributors covering 70 different files in the pipeline: 7,357 additions and 8,236 deletions!

In summary, the main changes are:

* Rebranding and renaming throughout the pipeline to nf-core
* Updating many parts of the pipeline config and style to meet nf-core standards
* Support for GFF files in addition to GTF files
  * Just use `--gff` instead of `--gtf` when specifying a file path
* New command line options to skip various quality control steps
* More safety checks when launching a pipeline
  * Several new sanity checks - for example, that the specified reference genome exists
* Improved performance with memory usage (especially STAR and Picard)
* New BigWig file outputs for plotting coverage across the genome
* Refactored gene body coverage calculation, now much faster and using much less memory
* Bugfixes in the MultiQC process to avoid edge cases where it wouldn't run
* MultiQC report now automatically attached to the email sent when the pipeline completes
* New testing method, with data on GitHub
  * Now run pipeline with `-profile test` instead of using bash scripts
* Rewritten continuous integration tests with Travis CI
* New explicit support for Singularity containers
* Improved MultiQC support for DupRadar and featureCounts
  * Now works for all users instead of just NGI Stockholm
* New configuration for use on AWS batch
* Updated config syntax to support latest versions of Nextflow
* Built-in support for a number of new local HPC systems
  * CCGA, GIS, UCT HEX, updates to UPPMAX, CFC, BINAC, Hebbe, c3se
* Slightly improved documentation (more updates to come)
* Updated software packages

...and many more minor tweaks.

Thanks to everyone who has worked on this release!
