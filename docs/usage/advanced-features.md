---
order: 5
---

# Advanced features

This page covers specialised use cases including prokaryotic RNA-seq, UMI handling, 3' digital gene expression assays, and GPU-accelerated alignment.

## Prokaryotic RNA-seq

The pipeline includes a dedicated profile for analysing bacterial and archaeal RNA-seq data:

:::note
This pipeline is primarily designed for eukaryotic RNA-seq. The prokaryotic profile provides basic support but does not include specialist features like operon detection or TSS mapping.
:::

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --fasta genome.fasta \
    --gff annotation.gff3 \
    --outdir results \
    -profile prokaryotic,docker
```

### Input requirements

| Input        | Parameter          | Requirements                                                                                                       |
| ------------ | ------------------ | ------------------------------------------------------------------------------------------------------------------ |
| Samplesheet  | `--input`          | Standard CSV format (see [Samplesheet input](samplesheet.md))                                                      |
| Genome FASTA | `--fasta`          | Genomic sequence file (`.fasta`, `.fa`, `.fna`, optionally gzipped)                                                |
| Annotation   | `--gff` or `--gtf` | Must contain **CDS features** with `gene_id` attributes. GFF3 format (`.gff3`, `.gff`) is typical for prokaryotes. |

**Key points:**

- **Use GFF3 format**: Prokaryotic annotations are typically distributed as GFF3 (not GTF). The pipeline accepts both via `--gff` or `--gtf`.
- **CDS features required**: The annotation must contain CDS (coding sequence) features. The pipeline extracts transcripts from these.
- **Matching contig names**: Chromosome/contig names in your FASTA must exactly match those in your GFF/GTF (for example, if your FASTA has `>NC_003197.2`, your GFF must use `NC_003197.2` in column 1).
- **No transcript FASTA needed**: The pipeline generates the transcript FASTA automatically using GFFREAD.

The `-profile prokaryotic` configures the pipeline with settings optimised for prokaryotic data:

- **Bowtie2 alignment**: Uses Bowtie2 instead of STAR as the default aligner. Prokaryotic genomes lack introns, so splice-aware aligners are unnecessary.
- **GFFREAD transcript extraction**: Uses GFFREAD instead of RSEM to generate transcript sequences. Prokaryotic annotations typically use CDS features rather than exon features, and RSEM requires exon features.
- **CDS-based counting**: Sets `featurecounts_feature_type` to `CDS` for QC metrics.
- **Skips eukaryote-specific QC**: Disables RSeQC and dupRadar, which assume eukaryotic gene structures.

You can override the aligner while keeping other prokaryotic settings:

```bash
# Use STAR instead of Bowtie2 with prokaryotic settings
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --fasta genome.fasta \
    --gff annotation.gff3 \
    --outdir results \
    --aligner star_salmon \
    -profile prokaryotic,docker
```

When using STAR with `-profile prokaryotic`, the pipeline automatically configures STAR to use CDS features instead of exons for both index building and alignment.

### Annotation sources

For prokaryotic genomes, we recommend using annotations from [Ensembl Bacteria](https://bacteria.ensembl.org/) when available, as these follow consistent formatting conventions. NCBI/RefSeq annotations can also work but require additional attention:

- **Contig name matching**: Ensure chromosome/contig names in your FASTA match those in your GTF/GFF. This is a common source of errors when mixing files from different sources.
- **Empty transcript IDs**: Some bacterial GTFs from RefSeq contain gene entries with empty `transcript_id` fields, which can cause Salmon errors. If you encounter issues, check your annotation file and consider filtering problematic entries.

### Use GFFREAD independently

The `--gffread_transcript_fasta` parameter can also be used independently of `-profile prokaryotic` for any situation where RSEM fails to extract transcripts correctly (for example, non-standard annotation formats):

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --fasta genome.fasta \
    --gtf annotation.gtf \
    --gffread_transcript_fasta \
    --outdir results \
    -profile docker
```

### Manual configuration (advanced)

If you prefer not to use the profile, you can manually configure the pipeline for prokaryotic data. The following parameters are set by `-profile prokaryotic`:

| Parameter                      | Value            | Purpose                     |
| ------------------------------ | ---------------- | --------------------------- |
| `--aligner`                    | `bowtie2_salmon` | Splice-unaware alignment    |
| `--gffread_transcript_fasta`   | `true`           | Handle CDS-only annotations |
| `--featurecounts_feature_type` | `CDS`            | QC counting on CDS features |
| `--skip_rseqc`                 | `true`           | Skip eukaryote-specific QC  |
| `--skip_dupradar`              | `true`           | Skip eukaryote-specific QC  |

Get in touch with us on the #rnaseq channel in the [nf-core Slack workspace](https://nf-co.re/join) if you are having problems or need any advice.

## UMI support

The pipeline supports Unique Molecular Identifiers to increase the accuracy of the quantification. UMIs are short sequences used to uniquely tag each molecule in a sample library and facilitate the accurate identification of read duplicates. They must be added during library preparation and prior to sequencing, therefore require appropriate arrangements with your sequencing provider.

To take UMIs into consideration during a workflow run, specify the `--with_umi` parameter. The pipeline currently supports UMIs which are embedded within a read's sequence and UMIs whose sequence is given inside the read's name. Consult your kit's manual and/or contact your sequencing provider regarding the exact specification.

The `--umitools_grouping_method` parameter affects [how similar, but non-identical UMIs](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html#method) are treated. `directional`, the default setting, is most accurate, but computationally very demanding. Consider `percentile` or `unique` if processing many samples.

### Examples

| UMI type     | Source                                                                                                                                                                                                                                                                                                                                                                                                                                                         | Pipeline parameters                                                                                                                                                         |
| ------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| In read name | [Illumina BCL convert >3.7.5](https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl_convert/bcl-convert-v3-7-5-software-guide-1000000163594-00.pdf)                                                                                                                                                                                                                                                | `--with_umi --skip_umi_extract --umitools_umi_separator ":"`                                                                                                                |
| In sequence  | [Lexogen QuantSeq® 3' mRNA-Seq V2 FWD](https://www.lexogen.com/quantseq-3mrna-sequencing) + [UMI Second Strand Synthesis Module](https://faqs.lexogen.com/faq/how-can-i-add-umis-to-my-quantseq-libraries)                                                                                                                                                                                                                                                    | `--with_umi --umitools_extract_method "regex" --umitools_bc_pattern "^(?P<umi_1>.{6})(?P<discard_1>.{4}).*"`                                                                |
| In sequence  | [Lexogen CORALL® Total RNA-Seq V1](https://www.lexogen.com/corall-total-rna-seq/)<br> > _mind [Appendix H](https://www.lexogen.com/wp-content/uploads/2020/04/095UG190V0130_CORALL-Total-RNA-Seq_2020-03-31.pdf) regarding optional trimming_                                                                                                                                                                                                                 | `--with_umi --umitools_extract_method "regex" --umitools_bc_pattern "^(?P<umi_1>.{12}).*"`<br>Optional: `--clip_r2 9 --three_prime_clip_r2 12`                              |
| In sequence  | Takara Bio [SMART-Seq Total RNA Pico Input with UMIs](https://www.takarabio.com/documents/User%20Manual/SMART/SMART-Seq%20Total%20RNA%20Pico%20Input%20with%20UMIs%20%28ZapR%20Mammalian%29%20User%20Manual.pdf) and [SMARTer® Stranded Total RNA-Seq Kit v3](https://www.takarabio.com/documents/User%20Manual/SMARTer%20Stranded%20Total%20RNA/SMARTer%20Stranded%20Total%20RNA-Seq%20Kit%20v3%20-%20Pico%20Input%20Mammalian%20User%20Manual-a_114950.pdf) | `--with_umi --umitools_extract_method "regex" --umitools_bc_pattern2 "^(?P<umi_1>.{8})(?P<discard_1>.{6}).*"`                                                               |
| In sequence  | [Watchmaker mRNA Library Prep Kit](https://watchmakergenomics.com/wp-content/uploads/2023/11/M223_mRNA-Library-Prep-Kit-_UG_WMUG214_v1-1-0823.pdf) with [Twist UMI Adapter System](https://www.twistbioscience.com/sites/default/files/resources/2023-03/DOC-001337_TechNote-ProcessingSequencingDataUtilizingUMI-REV1-singles.pdf)                                                                                                                            | `--with_umi --umitools_extract_method "regex" --umitools_bc_pattern "^(?P<umi_1>.{5})(?P<discard_1>.{2}).*" --umitools_bc_pattern2 "^(?P<umi_2>.{5})(?P<discard_2>.{2}).*"` |

> _No warranty for the accuracy or completeness of the parameters is implied_

## 3' digital gene expression assays

Some bulk RNA-seq library preparation protocols capture only a 3' tag from each transcript, for example [3'Pool-seq](https://pubmed.ncbi.nlm.nih.gov/31959126/), [DRUG-seq](https://pubs.acs.org/doi/10.1021/acschembio.1c00920), [BRB-seq](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1671-x) or Lexogen's commercial [QuantSeq 3' mRNA-seq FWD](https://www.lexogen.com/quantseq-3mrna-sequencing/) protocol. The following parameters have been validated for `QuantSeq 3' mRNA-seq FWD` data, and provide useful starting points for other 3' RNA-seq protocols:

### Custom STAR parameters

Lexogen provides an example analysis workflow [on their website](https://www.lexogen.com/quantseq-data-analysis/), which includes the _ENCODE standard options_ for the [STAR aligner](https://github.com/alexdobin/STAR). In addition, Lexogen also decreases the tolerance for mismatches and clips poly(A) tails. To apply these settings, add the following parameters when running the pipeline:

```bash
--extra_star_align_args "--alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNoverLmax 0.1 --clip3pAdapterSeq AAAAAAAA"
```

### Custom Salmon arguments

[Salmon's default quantitation algorithm](https://www.nature.com/articles/nmeth.4197) takes into account transcript length. Because 3' tag protocols do not capture full transcripts, this feature needs to be deactivated by specifying:

```bash
--extra_salmon_quant_args "--noLengthCorrection"
```

### QuantSeq analysis with UMIs

If unique molecular identifiers were used to prepare the library, add the following arguments as well, to extract the UMIs and deduplicated alignments:

```bash
--with_umi
--umitools_extract_method regex
--umitools_bc_pattern "^(?P<umi_1>.{6})(?P<discard_1>.{4}).*"
```

## GPU acceleration

### Sentieon acceleration for STAR

The STAR aligner can be accelerated through its Sentieon implementation using the parameter `--use_sentieon_star`.

Sentieon is a commercial solution to process genomics data, requiring a paid license. Sentieon's tooling contains an accelerated version of the [`STAR` aligner](https://support.sentieon.com/manual/usages/general/?highlight=star#star-binary), which nf-core/rnaseq supports. In order to use those functions, the user will need to supply a license for Sentieon.

Sentieon supply license in the form of a string-value (a url) or a file. It should be base64-encoded and stored in a nextflow secret named `SENTIEON_LICENSE_BASE64`. If a license string (url) is supplied, then the nextflow secret should be set like this:

```bash
nextflow secrets set SENTIEON_LICENSE_BASE64 $(echo -n <sentieon_license_string> | base64 -w 0)
```

:::note
`<sentieon_license_string>` is formatted as `IP:Port` for example: `12.12.12.12:8990`
:::

If a license file is supplied, then the nextflow secret should be set like this:

```bash
nextflow secrets set SENTIEON_LICENSE_BASE64 $(cat <sentieon_license_file.lic> | base64 -w 0)
```

:::note
If you're looking for documentation on how the nf-core Sentieon GitHub Actions and Sentieon License Server are set up: [Here be dragons.](https://github.com/nf-core/ops/blob/main/pulumi/sentieon_license_server/README.md)

For detailed instructions on how to test the modules and subworkflows separately, see [here](https://github.com/nf-core/modules/blob/master/modules/nf-core/sentieon/README.md).
:::

### Parabricks GPU acceleration for STAR

The STAR aligner can also be GPU-accelerated using NVIDIA Parabricks via the `--use_parabricks_star` parameter. Parabricks runs STAR alignment on NVIDIA GPUs, significantly reducing wall-clock time for large datasets.

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta genome.fa \
    --gtf annotation.gtf \
    --use_parabricks_star \
    -profile docker
```

Container GPU flags (`--gpus all` for Docker, `--nv` for Singularity/Apptainer) are automatically applied to GPU tasks based on the container engine. These flags are scoped to GPU tasks only, so non-GPU steps will run normally on CPU-only nodes in mixed clusters.

The container GPU flags can be overridden via `--gpu_container_options` (for example, `--gpu_container_options '--gpus 1'`).

#### Requirements

- One or more NVIDIA GPUs (with appropriate drivers installed)
- Docker or Singularity (Conda/Mamba is **not** supported for this module)
- The Parabricks container (`nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1`) will be pulled automatically

#### Behaviour differences

When using Parabricks, the pipeline automatically handles mark duplicates during the alignment step (via `pbrun rna_fq2bam`), so the separate Picard MarkDuplicates step is skipped.

#### Known differences from native STAR

Parabricks `rna_fq2bam` is based on STAR 2.7.2a. The following native STAR flags have no pbrun equivalent and are therefore not applied:

- `--outFilterType BySJout` — pbrun uses its own splice junction filtering defaults
- `--sjdbScore 1` — affects junction scoring priority
- `--quantTranscriptomeBan Singleend` — Salmon handles mixed single/paired-end transcriptome records gracefully
- `--runRNGseed 0` — pbrun uses deterministic primary alignment selection

These differences are unlikely to materially affect downstream quantification results, but be aware of them for reproducibility purposes. All other STAR parameters (multi-mapping limits, intron sizes, mate gap, splice junction overhangs, and others) have pbrun equivalents and are applied consistently.
