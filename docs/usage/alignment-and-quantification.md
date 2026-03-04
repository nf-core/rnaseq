---
order: 3
---

# Alignment and quantification

The pipeline supports several alignment and quantification strategies. The default (`--aligner star_salmon`) uses STAR for alignment and Salmon for quantification, which is suitable for most use cases.

## Which aligner should I use?

| Strategy                           | Parameter                   | Best for                                 | Notes                                                              |
| ---------------------------------- | --------------------------- | ---------------------------------------- | ------------------------------------------------------------------ |
| **STAR + Salmon** (default)        | `--aligner star_salmon`     | Most users                               | Full BAM-level QC, splice-aware alignment, accurate quantification |
| **STAR + RSEM**                    | `--aligner star_rsem`       | When RSEM-specific features needed       | Uses ENCODE3 STAR settings                                         |
| **HISAT2**                         | `--aligner hisat2`          | Memory-constrained environments          | No quantification — alignment and QC only                          |
| **Salmon** (pseudoalignment)       | `--pseudo_aligner salmon`   | Fast quantification, no alignment needed | Can run alongside or instead of alignment                          |
| **Kallisto** (pseudoalignment)     | `--pseudo_aligner kallisto` | Alternative to Salmon pseudoalignment    | Can run alongside or instead of alignment                          |
| **Bowtie2 + Salmon** (prokaryotic) | `--aligner bowtie2_salmon`  | Bacterial/archaeal data                  | Automatically set by `-profile prokaryotic`                        |

:::tip
For GPU-accelerated STAR alignment, see [Advanced features — GPU acceleration](advanced-features.md#gpu-acceleration).
:::

## STAR + Salmon (default)

By default, the pipeline uses [STAR](https://github.com/alexdobin/STAR) to map the raw FastQ reads to the reference genome, project the alignments onto the transcriptome and to perform the downstream BAM-level quantification with [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html). STAR is fast but requires a lot of memory to run, typically around 38GB for the Human GRCh37 reference genome.

Two additional parameters `--extra_star_align_args` and `--extra_salmon_quant_args` allow you to append any custom parameters to the STAR align and Salmon quant commands, respectively. Note, the `--seqBias` and `--gcBias` are not provided to Salmon quant by default so you can provide these via `--extra_salmon_quant_args '--seqBias --gcBias'` if required.

## STAR + RSEM

Selecting `--aligner star_rsem` uses STAR for alignment with RSEM-specific settings based on the ENCODE3 STAR parameters.

:::note
Selecting `star_rsem` automatically applies the same STAR settings as [rsem-calculate-expression](https://deweylab.github.io/RSEM/rsem-calculate-expression.html) with the `--star` option. These are based on the ENCODE3 STAR settings, and are as follows:

<details>
<summary>View STAR parameters</summary>

```raw
--outSAMunmapped Within
--outFilterType BySJout
--outFilterMultimapNmax 20
--outFilterMismatchNmax 999
--outFilterMismatchNoverLmax 0.04
--alignIntronMin 20
--alignIntronMax 1000000
--alignMatesGapMax 1000000
--alignSJoverhangMin 8
--alignSJDBoverhangMin 1
--sjdbScore 1
```

</details>

**NOTE**: The pipeline parameter [`--extra_star_align_args`](#custom-star-parameters) cannot be used with aligner option `star_rsem`. It is possible to set, via [custom tool arguments](configuration.md#custom-tool-arguments), custom settings to the process `STAR_ALIGN` via `ext.args`. However, we discourage this and consider adjusting the STAR aligner while using RSEM to be an unsupported option in this pipeline. If you need to pass in custom star alignment options, we recommend using the aligner option `star_salmon`. After this warning, if you still wish to set custom STAR settings and use RSEM for quantification, then reference [this configuration file](https://github.com/nf-core/rnaseq/blob/master/subworkflows/local/align_star/nextflow.config).
:::

## HISAT2

Use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) (`--aligner hisat2`) if you have memory limitations. HISAT2 provides alignment and QC but does not perform quantification due to the lack of an appropriate option to calculate accurate expression estimates from HISAT2-derived genomic alignments — this may change in future releases (see [#822](https://github.com/nf-core/rnaseq/issues/822)). HISAT2 has been made available for those who have a preference for the alignment, QC and other types of downstream analysis compatible with its output.

## Pseudoalignment with Salmon or Kallisto

You have the option to pseudoalign and quantify your data directly with [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) or [Kallisto](https://pachterlab.github.io/kallisto/) by specifying `salmon` or `kallisto` to the `--pseudo_aligner` parameter. The selected pseudoaligner will then be run in addition to the standard alignment workflow defined by `--aligner`, mainly because it allows you to obtain QC metrics with respect to the genomic alignments. However, you can provide the `--skip_alignment` parameter if you would like to run Salmon or Kallisto in isolation.

By default, the pipeline will use the genome fasta and gtf file to generate the transcripts fasta file, and then to build the Salmon index. You can override these parameters using the `--transcript_fasta` and `--salmon_index` parameters, respectively. You can also supply additional arguments to Kallisto via `--extra_kallisto_quant_args`.

### Library type inference

The library preparation protocol (library type) used by Salmon quantification is inferred by the pipeline based on the information provided in the samplesheet, however, you can override it using the `--salmon_quant_libtype` parameter. You can find the available options in the [Salmon documentation](https://salmon.readthedocs.io/en/latest/library_type.html). Similarly, strandedness is taken from the sample sheet or calculated automatically, and passed to Kallisto on a per-library basis, but you can apply a global override by setting the Kallisto strandedness parameters in `--extra_kallisto_quant_args` like `--extra_kallisto_quant_args '--fr-stranded'` see the [Kallisto documentation](https://pachterlab.github.io/kallisto/manual).

### Decoy-aware transcriptome index

When running Salmon in mapping-based mode via `--pseudo_aligner salmon`, supplying a genome fasta via `--fasta` and not supplying a Salmon index, the entire genome of the organism is used by default for the decoy-aware transcriptome when creating the indices, as is recommended (see second bulleted option in [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode)). If you do not supply a FASTA file or an index, Salmon will index without those decoys, using only transcript sequences in the index. This second option is not usually recommended, but may be useful in limited circumstances. Note that Kallisto does not index with genomic sequences.

## Skip alignment and pseudoalignment

:::note
You can use `--skip_alignment --skip_pseudo_alignment` if you only want to run the pre-processing QC steps in the pipeline like FastQC and trimming. This will skip alignment, pseudoalignment and any post-alignment processing steps.
:::

Note that `--skip_alignment` and `--skip_pseudo_alignment` prevent both the execution of alignment/pseudoalignment steps and the building of their corresponding indices. For example, using `--skip_alignment` with `--aligner star_salmon` will skip both STAR alignment and index building.

## Quantification

The current options align with STAR and quantify using either Salmon (`--aligner star_salmon`) / RSEM (`--aligner star_rsem`). You also have the option to pseudoalign and quantify your data with Salmon or Kallisto by providing the `--pseudo_aligner salmon` or `--pseudo_aligner kallisto` parameter, respectively.

Since v3.0 of the pipeline, featureCounts is no longer used to perform gene/transcript quantification, however it is still used to generate QC metrics based on [biotype](http://www.ensembl.org/info/genome/genebuild/biotypes.html) information available within GFF/GTF genome annotation files. This decision was made primarily because of the limitations of featureCounts to appropriately quantify gene expression data. See [Zhao et al., 2015](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0141910#pone-0141910-t001) and [Soneson et al., 2015](https://f1000research.com/articles/4-1521/v1).

For similar reasons, quantification will not be performed if using `--aligner hisat2` due to the lack of an appropriate option to calculate accurate expression estimates from HISAT2 derived genomic alignments — this may change in future releases (see [#822](https://github.com/nf-core/rnaseq/issues/822)). HISAT2 has been made available for those who have a preference for the alignment, QC and other types of downstream analysis compatible with its output.
