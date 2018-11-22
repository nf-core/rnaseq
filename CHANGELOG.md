# nf-core/rnaseq

## Version 1.2dev

#### Bug Fixes
* Fixed conda bug which caused problems with environment resolution due to changes in bioconda [#113](https://github.com/nf-core/rnaseq/issues/113)
* Remove legacy rseqc `strandRule` config code. [#119](https://github.com/nf-core/rnaseq/issues/119)
* Removed some outdated documentation about non-existent features
* Config refactoring and code cleaning
* Added a `--fcExtraAttributes` option to specify more than ENSEMBL gene names in `featureCounts`

## [Version 1.1](https://github.com/nf-core/rnaseq/releases/tag/1.1) - 2018-10-05

#### Pipeline updates
* Wrote docs and made minor tweaks to the `--skip_qc` and associated options
* Removed the depreciated `uppmax-modules` config profile
* Updated the `hebbe` config profile to use the new `withName` syntax too
* Use new `workflow.manifest` variables in the pipeline script
* Updated minimum nextflow version to `0.32.0`

#### Bug Fixes
* [#77](https://github.com/nf-core/rnaseq/issues/77): Added back `executor = 'local'` for the `workflow_summary_mqc`
* [#95](https://github.com/nf-core/rnaseq/issues/95): Check if task.memory is false instead of null
* [#97](https://github.com/nf-core/rnaseq/issues/97): Resolved edge-case where numeric sample IDs are parsed as numbers causing some samples to be incorrectly overwritten.


## [Version 1.0](https://github.com/nf-core/rnaseq/releases/tag/1.0) - 2018-08-20

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
