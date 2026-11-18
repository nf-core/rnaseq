# hg19 Reference Build Scripts

These scripts build a consistent set of gene annotation references for three repositories — rnaseq (RNA quantification), pasta (epigenomic TSS scoring), and leming (expression prediction from epigenomics) — starting from the same NCBI source file.

All three outputs derive from the NCBI GRCh37.p13 annotation release 105.20220307 GTF, so gene names, transcript coordinates, and chromosome naming stay consistent across the set.

## Source

**NCBI GRCh37.p13 annotation release 105.20220307** (March 2022), the last NCBI annotation freeze for hg19.

```
https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/105.20220307/
  GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz
```

MD5 (gzipped): `7877d9e1383dc194c4b22d82ae07d759`

Each script downloads this file, verifies its checksum, and converts RefSeq accession chromosome names (NC_000001.10, etc.) to UCSC names (chr1, etc.) using the companion assembly report.

## Why this source

**Why NCBI rather than UCSC?** UCSC distributes RefSeq annotations in genePred format (ncbiRefSeqCurated.txt.gz, ncbiRefSeqSelect.txt.gz), not GTF. The rnaseq pipeline requires GTF for STAR and Salmon, and converting genePred back to GTF loses attributes like `gene_biotype`. The NCBI GTF has these attributes natively.

**Why not lift over from GRCh38?** NCBI provides a GRCh37-native annotation, so liftover is unnecessary. The 2022 GRCh37 freeze and the latest GRCh38 release share ~99.97% of BestRefSeq transcripts. The 500 or so GRCh38-only transcripts are mostly recently characterized lncRNAs and do not include well-known cancer genes.

**Why the 2022 freeze despite being older?** It is the last official NCBI annotation on this assembly, and it is comprehensive for protein-coding genes. The 98 gene symbol renames between 2022 NCBI and 2024 UCSC (DDX58 to RIGI, GBA to GBA1, various FAM/KIAA renames) are the main practical difference. By using a single source, we avoid these naming mismatches across pipelines.

**Why NCBI-derived RefSeq Select rather than the UCSC ncbiRefSeqSelect table?** The UCSC table uses 2024 gene symbols while the NCBI GTF uses 2022 symbols. Filtering the Select tag from the same GTF that produces the curated set guarantees name consistency. The two sources agree on transcript identity for 98.3% of genes; the discordant entries are documented in `assets/refseq_select_discordance_report.tsv` and `assets/refseq_select_gene_discordance.tsv`.

## Scripts

### `build_refseq_curated_gtf.sh <output_dir>` — for rnaseq

Filters the GTF to `source == "BestRefSeq"` on primary + mitochondrial chromosomes. This excludes Gnomon-predicted transcripts (XM_, XR_ accessions) and retains only manually curated RefSeq entries (NM_, NR_).

Output: a gzipped GTF suitable as the `--gtf` parameter for nf-core/rnaseq. The GTF includes `gene_biotype` attributes, which means the pipeline's biotype QC (previously skipped via `skip_biotype_qc = true`) can be re-enabled.

**28,401 genes, 85,432 transcripts.**

### `build_refseq_curated_indices.sh --gtf <path> --outdir <path>` — for rnaseq

Builds STAR and Salmon indices by running the nf-core/rnaseq pipeline end-to-end with `--save_reference` on a single sample. The pipeline only builds indices when alignment and pseudo-alignment are enabled, so the full pipeline must run (not just the genome preparation step). The sample results are incidental; the indices under `<outdir>/genome/index/` are the output.

### `build_refseq_select_genepred.sh <output_dir>` — for pasta

Filters the GTF to features tagged `RefSeq Select` (one representative transcript per gene), then converts to UCSC extended genePred format (16 columns). The output is a drop-in replacement for `pasta/data/ncbiRefSeqSelect.txt.gz`.

Genes are additionally filtered to those with at least one BestRefSeq curated transcript on a standard chromosome (chr1-22, X, Y, M). This removes two uncharacterized loci (LOC100288966, LOC389831) that have Select transcripts only on non-standard contigs.

One transcript override is applied: **CLDN18 uses NM_001002026 (CLDN18.2) instead of the RefSeq Select transcript NM_016369 (CLDN18.1)**. CLDN18.1 and CLDN18.2 are driven by separate promoters 11.3kb apart. The ADC therapeutic target (zolbetuximab) binds CLDN18.2 specifically, and the CLDN18.2 TSS is the relevant one for predicting target expression from epigenomic data. This override was validated by correlating PASTA combined z-scores with RNA-seq TPM across 1,248 cell line and tissue samples: CLDN18 had the lowest K4 promoter correlation among ADC target genes (r=0.22), consistent with PASTA scoring the wrong promoter. The `TX_OVERRIDES` dict in the script can hold additional per-gene overrides if others are needed later.

**19,159 transcripts (one per gene).**

### `build_gene_locations_rds.sh <output_dir>` — for leming

Extracts all BestRefSeq transcript coordinates per gene, merges overlapping transcripts into distinct loci (via `GenomicRanges::reduce()`), and picks the locus containing the TSS as the primary range. The output is an R named list matching the structure leming expects: `gene` (GRanges), `chr`, `strand`, `tss`, `all_ranges` (GRanges).

This replaces a file previously derived from `TxDb.Hsapiens.UCSC.hg19.knownGene`, which used Gencode-based transcript models and HGNC symbols from whatever version of `org.Hs.eg.db` was installed at build time.

Standard chromosomes only (chr1-22, X, Y). 30 genes have multiple non-overlapping loci, mostly pseudoautosomal region genes with copies on both X and Y.

**28,350 genes.**

### `adc_pasta_rna_correlation.R <output_dir>` — validation

Queries PASTA vd5_pebsnr z-scores and joins them with RNA-seq TPM for 48 ADC target genes across 1,248 matched samples. Outputs per-gene Pearson correlations (K4, K27, combined) and expression summary statistics.

The ADC genes have a median K4 correlation of 0.81 (genome-wide median is 0.45), confirming that RefSeq Select TSSes are well-placed for these targets. Results are in `assets/adc_pasta_rna_correlation.tsv`; genome-wide results are in `assets/all_genes_pasta_rna_correlation.tsv`.

## Quantification comparison

Sample PB23125 was run through the full nf-core/rnaseq pipeline with both the old Gencode/UCSC-based GTF (`Homo_sapiens_hg19.genes.gtf`) and the new NCBI BestRefSeq curated GTF. Salmon gene-level TPM was compared across 12,750 expressed genes (TPM > 0.1 in either reference). 98 genes changed symbol between old and new (GNB2L1→RACK1, PVRL4→NECTIN4, etc.) and were matched by genomic position.

**Log2(TPM+1) Pearson correlation: r = 0.966.**

333 genes showed >4-fold TPM change. Most are pseudogenes or sit in segmental duplication regions where small changes in transcript models shift multimapping reads. Among ADC-relevant genes, shifts were smaller but worth tracking:

| Gene | Old TPM | New TPM | Change | Old tx count | New tx count |
|------|---------|---------|--------|-------------|-------------|
| CDH3 | 112 | 134 | +20% | 1 | 3 |
| EGFR | 48 | 39 | −19% | 4 | 9 |
| LY75 | 6.3 | 5.5 | −13% | 1 | 1 |

These shifts trace to transcript set differences, not alignment changes. The new GTF averages 3.01 transcripts per gene (85,432 total) versus 2.05 (54,070 total) in the old GTF. With more isoforms available, Salmon's EM algorithm distributes ambiguous reads differently. CDH3 gained two transcripts and TPM rose as reads spread across shorter isoforms with higher per-read TPM contribution. EGFR gained five long isoforms (several >9 kb) that competed for shared reads, pulling the gene total down.

The richer transcript set is an improvement: reads that the old reference forced into wrong isoforms or left unassigned now have correct targets. Any gene-level TPM thresholds calibrated on old-reference output (e.g., expression cutoffs in ADC target selection) should be rechecked against the new baseline.

## Gene name consistency

All three outputs use NCBI 2022 gene symbols. Downstream code that joins across these files (e.g., correlating PASTA scores with RNA counts, or overlaying leming's gene locations with rnaseq expression) will match on gene name without reconciliation.

The rnaseq pipeline's existing Salmon TPM outputs use different gene names (from the old Gencode-based GTF). After switching to the new curated GTF, re-quantified samples will use the NCBI 2022 names. Any code that hardcodes renamed symbols (the Rmd notebooks manually rename PVRL4 to NECTIN4, for instance) should be reviewed — the new GTF already uses NECTIN4 natively.
