# NGI-RNAseq

## 1.3dev

* Updated HISAT2 from v2.0.5 to v2.1.0
  * Uses `--new-summary` and `--summary-file` to give output that will work with MultiQC,
  * UPPMAX environment module load and Docker image
* Moved `pipefail` statement into config, applies to all processes
* Rewrote summary e-mail commands. Now uses `sendmail`
  * Proper multipart text/html e-mail with embedded images
  * If sendmail fails, falls back to sending plaintext using `mail`
* MultiQC process now runs using `local` executor for internet access on some UPPMAX clusters.
* UPPMAX featureCounts process now loads `python/2.7.11` environment module
* New test script for uppmax with HISAT2
* Timeline and trace now always generated for every run
* Script now checks that the version of Nextflow is recent enough and warns if not

## [1.2](https://github.com/SciLifeLab/NGI-RNAseq/releases/tag/1.2) - 2017-06-13

* MultiQC now runs using `local` instead of `slurm` with default config
  * Means that it stands a better chance of getting information from a remote database with plugins
* MultiQC now uses the workflow name for report title and filename if specified
* Fixed error where trimming parameters weren't set properly for `--pico`
* Made pipeline run in forward-stranded mode when using `--pico`


## [1.1](https://github.com/SciLifeLab/NGI-RNAseq/releases/tag/1.1) - 2017-05-23
New release of the pipeline with a few additions:

* Summary HTML report now created upon pipeline completion
  * This is saved in the results folder and also e-mailed if `--email` is set.
  * Nice way to alert the person running the pipeline that it completed and shows whether any errors were encountered.
  * Good for reproducibility as it logs all pipeline parameters to a easy to interpret file.
* Pipeline assumes that it's running with Paired-end files now
  * An error is now raised if the file glob doesn't give pairs of files
  * If running with single-end data, use the `--singleEnd` command-line option.
* Timelines and traces now created by default for the testing configs
* New configs and documentation about running the pipeline on AWS
* Made sure that the `.bam` files ended up in the main STAR directory when `--saveAlignedIntermediates` is used, instead of `STAR/logs`
* Made MultiQC load its config file through a channel instead of directly copying from `baseDir`


## [1.0.4](https://github.com/SciLifeLab/NGI-RNAseq/releases/tag/1.0.4) - 2017-04-21
RseQC hotfix, input file was not supplied properly to one of the scripts

## [1.0.3](https://github.com/SciLifeLab/NGI-RNAseq/releases/tag/1.0.3) - 2017-04-19
Hotfix to fix minor bug affecting strandedness for StringTie run.

## [1.0.2](https://github.com/SciLifeLab/NGI-RNAseq/releases/tag/1.0.2) - 2017-04-11
A couple of tweaks to help the pipeline in production:

* Trimming FastQ files and intermediate BAM files now not saved by default
  * This is configurable in the config or with `--saveTrimmed` / `--saveAlignedIntermediates`
* featureCounts merge process uses `.collect()` for better consistency

## [1.0.1](https://github.com/SciLifeLab/NGI-RNAseq/releases/tag/1.0.1) - 2017-04-10
This release includes a bugfix for the last major release relating to the strandedness of `RSEQC`.

* Single end reverse is now correctly `+-,-+.`
* Single end forward is now correctly `++, --`
* PE forward is now correctly `-1++,1--,2+-,2-+`

## [1.0](https://github.com/SciLifeLab/NGI-RNAseq/releases/tag/1.0) - 2017-04-05
The pipeline has now been validated for use in our production work.
This version includes some new features:

* The output from featureCounts is now merged into a single table and supplied along side the individual reports.
* markDuplicates JVM memory is now automatically scaled based on the process memory
* an `html` file with results documentation is now generated and supplied amongst the results
* It's now possible to configure the pipeline for different stranded libraries with just a simple CL flag.
* Additional support and documentation for other platforms than Uppmax. Inluding C3SE.
* + Numerous minor tweaks and improvements.

## [0.3](https://github.com/SciLifeLab/NGI-RNAseq/releases/tag/0.3) - 2016-12-13
In order to properly validate this pipeline and take it into production we need to tag a stable release.
I've tagged specific software versions in the `uppmax.config` file.

## [0.2](https://github.com/SciLifeLab/NGI-RNAseq/releases/tag/0.2) - 2016-10-14
First (semi-) stable release of the new NGI-RNAseq pipeline, as we head towards deployment in production.
