---
order: 2
---

# Reference genomes

Refer to the [nf-core website](https://nf-co.re/usage/reference_genomes) for general usage docs and guidelines regarding reference genomes.

## Provide reference files (recommended)

The minimum reference genome requirements for this pipeline are a FASTA file (genome and/or transcriptome) and GTF file. All other files required to run the pipeline can be generated from these files. For example, the latest reference files for human can be derived from Ensembl like:

```bash
latest_release=$(curl -s 'http://rest.ensembl.org/info/software?content-type=application/json' | grep -o '"release":[0-9]*' | cut -d: -f2)
wget -L ftp://ftp.ensembl.org/pub/release-${latest_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -L ftp://ftp.ensembl.org/pub/release-${latest_release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${latest_release}.gtf.gz
```

Specify these files to the workflow with the `--fasta` and `--gtf` parameters.

:::note
**Consistent reference resource usage**

When supplying reference files as discussed below, it is important to be consistent in the reference resource used (Ensembl, GENCODE, UCSC, and others), since differences in conventions between these resources can make their files incompatible. For example, UCSC prefixes chromosomes with `chr`, while Ensembl does not, so a GTF file from Ensembl cannot be supplied alongside a genome FASTA from UCSC. GENCODE also attaches version identifiers to gene and transcript names (for example, `ENSG00000254647.1`) while Ensembl does not.
:::

Notes:

- The pipeline supports compressed reference files, that is, standard files with the `.gz` extension and indices folders with the `tar.gz` extension.
- If `--gff` is provided as input then this will be converted to a GTF file, or the latter will be used if both are provided.
- If `--gene_bed` is not provided then it will be generated from the GTF file.
- If `--additional_fasta` is provided then the features in this file (such as ERCC spike-ins) will be automatically concatenated onto both the reference FASTA file as well as the GTF annotation before building the appropriate indices. Note: if you need the pipeline to build a pseudo-aligner index (Salmon/Kallisto), `--additional_fasta` cannot be used together with `--transcript_fasta` because the pipeline cannot append additional sequences to a user-provided transcriptome. Either omit `--transcript_fasta` and let the pipeline generate it, or provide a pre-built index that already contains the spike-ins.
- When using `--aligner star_rsem`, the pipeline will build separate STAR and RSEM indices. STAR performs alignment with RSEM-compatible parameters, then RSEM quantifies from the resulting BAM files using `--alignments` mode.
- If the `--skip_alignment` option is used along with `--transcript_fasta`, the pipeline can technically run without providing the genomic FASTA (`--fasta`). However, this approach is **not recommended** with `--pseudo_aligner salmon`, as any dynamically generated Salmon index will lack decoys. To ensure optimal indexing with decoys, include the genomic FASTA (`--fasta`) with Salmon, unless a pre-existing decoy-aware Salmon index is supplied. For more details on the benefits of decoy-aware indexing, refer to the [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode).

:::tip
For prokaryotic (bacterial/archaeal) data, see [Advanced features — Prokaryotic RNA-seq](advanced-features.md#prokaryotic-rna-seq).
:::

## Reference genome

Provide the most complete reference genome for your species, without additional loci (haplotypes) or patches. For model organisms such as mouse or human, this is the "primary assembly", which includes the reference chromosomes and some additional scaffolds. For the human assembly GRCh38 (hg38), use the `GRCh38.primary_assembly.genome.fa.gz` file from GENCODE or the `Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` file from Ensembl. These files cover the largest portion of the reference genome without including multiple copies of the same sequence, which would result in heavy mapping quality penalties.

For most other species (such as fly, cow, dog), no primary assembly is published. This reflects inadequately characterised genomic variation and a lower degree of curation, meaning that there are no established alternative loci (haplotypes), and that the toplevel file is equivalent to a primary assembly. Therefore, while the toplevel assembly can be used for these organisms, it is advisable to verify the absence of N-padded haplotype or patch regions first.

## Gene annotation

Gene annotations are updated more frequently than the reference genome sequence, so you must choose an appropriate annotation version (such as an Ensembl release). We recommend using sources with well-defined, versioned releases such as ENSEMBL or GENCODE. Generally, it is best to use the most recent release for the latest gene annotations. However, if you are combining your data with older datasets, use the annotation version previously used for consistency.

Once you have chosen a release, select the annotation file that matches your reference genome. For the human primary assembly, use the comprehensive annotation (such as `gencode.{release}.primary_assembly.annotation.gtf.gz` from GENCODE or `Homo_sapiens.GRCh38.{release}.gtf.gz` from Ensembl). For other species, like fly, use the annotation matching the toplevel assembly (such as `Drosophila_melanogaster.BDGP6.46.{release}.gtf.gz` from Ensembl).

Ensure that the annotation files use gene IDs as the primary identifier, not the gene name/symbol. For example, the Ensembl ID `ENSG00000254647` corresponds to the `INS` gene, which encodes the insulin protein. While gene names are more familiar, it is crucial to retain and use the primary identifiers as they are unique and easier to map between annotation versions or sources.

To take advantage of all the quality control modules implemented in the pipeline, the gene annotation must include a `gene_biotype` field which describes the function of each feature (protein coding, long non-coding, and others). This is usually the case for annotations from GENCODE or Ensembl but not always for annotations from other sources. If your annotation does not include this field, set the `--skip_biotype_qc` option to avoid running the steps that rely on it.

:::note
**GTF vs GFF**

GFF (General Feature Format) is a tab-separated text file format for representing genomic annotations, while GTF (General Transfer Format) is a specific implementation of this format corresponding to GFF version 2. The pipeline can accept both GFF and GTF but any GFF files will be converted to GTF so if a GTF is available for your annotation of choice it is better to provide that directly.

More information and links to further resources are [available from Ensembl](https://www.ensembl.org/info/website/upload/gff.html).
:::

## Reference transcriptome

In addition to the reference genome sequence and annotation, you can provide a reference transcriptome FASTA file. These files can be obtained from GENCODE or Ensembl. However, these sequences only cover the reference chromosomes and can cause inconsistencies if you are using a primary or toplevel genome assembly and annotation.

We recommend not providing a transcriptome FASTA file and instead allowing the pipeline to create it from the provided genome and annotation. Similar to aligner indexes, you can save the created transcriptome FASTA and BED files to a central location for future pipeline runs. This helps avoid redundant computation and having multiple copies on your system. Ensure that all genome, annotation, transcriptome, and index versions match to maintain consistency.

:::warning
If you are using `--additional_fasta` to add spike-in sequences (such as ERCC) and need the pipeline to build a pseudo-aligner index (Salmon/Kallisto), you **must not** provide `--transcript_fasta`. The pipeline needs to generate the transcriptome itself so that it includes the spike-in sequences. This combination will cause the pipeline to exit with an error unless you also provide a pre-built index (`--salmon_index` or `--kallisto_index`) that already contains the spike-in sequences.
:::

## Indices

By default, indices are generated dynamically by the workflow for tools such as STAR and Salmon. Since indexing is an expensive process in time and resources, ensure that it is only done once by retaining the indices generated from each batch of reference files by specifying `--save_reference`.

Once you have the indices from a workflow run, save them somewhere central and reuse them in subsequent runs using custom config files or command line parameters such as `--star_index '/path/to/STAR/index/'`.

Remember to note the genome and annotation versions as well as the versions of the software used for indexing, as an index created with one version is not always compatible with other versions.

## GENCODE

If you are using [GENCODE](https://www.gencodegenes.org/) reference genome files, specify the `--gencode` parameter because the format of these files is slightly different to ENSEMBL genome files:

- The `--featurecounts_group_type` parameter will automatically be set to `gene_type` as opposed to `gene_biotype`, respectively.
- If you are running Salmon, the `--gencode` flag will also be passed to the index building step to overcome parsing issues resulting from the transcript IDs in GENCODE fasta files being separated by vertical pipes (`|`) instead of spaces (see [this issue](https://github.com/COMBINE-lab/salmon/issues/15)).

As well as the standard annotations, GENCODE also provides "basic" annotations, which include only representative transcripts, but we do not recommend using these.

## Large chromosomes (plant genomes)

Genomes with very large chromosomes (>500 Mb), such as plant genomes, can encounter failures in the RSeQC `inner_distance` module due to a known limitation in the underlying bx-python library. The bx-python BitSet implementation has a maximum capacity of approximately 537 million bases, which can be exceeded by chromosomes in organisms like wheat, barley, and other plants.

If you encounter an error message similar to `IndexError: [coordinate] is larger than the size of this BitSet (536870912)`, work around this by excluding the `inner_distance` module from RSeQC analysis:

```bash
--rseqc_modules 'bam_stat,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication'
```

This removes `inner_distance` from the default list of RSeQC modules while retaining all other quality control metrics. Note that the inner_distance metric is only relevant for paired-end data and provides information about fragment size distribution.

For more information, see the upstream issues:

- [nf-core/rnaseq#608](https://github.com/nf-core/rnaseq/issues/608)
- [bxlab/bx-python#67](https://github.com/bxlab/bx-python/issues/67)

## iGenomes (not recommended)

If the `--genome` parameter is provided (such as `--genome GRCh37`) then the FASTA and GTF files (and existing indices) will be automatically obtained from AWS-iGenomes unless these have already been downloaded locally in the path specified by `--igenomes_base`.

However this is no longer recommended because:

- Gene annotations in iGenomes are extremely out of date. This can be particularly problematic for RNA-seq analysis, which relies on accurate gene annotation.
- Some iGenomes references (such as GRCh38) point to annotation files that use gene symbols as the primary identifier. This can cause issues for downstream analysis, such as the nf-core [differential abundance](https://nf-co.re/differentialabundance) workflow where a conventional gene identifier distinct from symbol is expected.

Notes:

- As of v3.7 of the pipeline, if you are using a genome downloaded from AWS iGenomes and using `--aligner star_salmon` (default) the version of STAR to use for the alignment will be auto-detected (see [#808](https://github.com/nf-core/rnaseq/issues/808)).

## Filter GTF

By default, the input GTF file will be filtered to ensure that sequence names correspond to those in the genome fasta file (where supplied), and to remove rows with empty transcript identifiers. Filtering can be bypassed completely where you are confident it is not necessary, using the `--skip_gtf_filter` parameter. If you want to skip only the `transcript_id` checking component of the GTF filtering script used in the pipeline, this can be disabled specifically using the `--skip_gtf_transcript_filter` parameter.
