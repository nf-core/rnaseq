#!/usr/bin/env Rscript
#
# Correlate PASTA combined z-scores with RNA TPM for ADC target genes.
#
# Inputs:
#   - PASTA vd5_pebsnr table (queried via precedeR::queryPBMC)
#   - RNA TPM parquet from leming/RNA/
#   - Library mapping from leming/RNA/
#   - ADC gene list from leming/workflows/
#
# Outputs (written to output_dir):
#   - adc_pasta_rna_correlation.tsv   Per-gene correlation + RNA summary stats
#   - adc_pasta_rna_joined.parquet    Full joined table for downstream use
#
# Usage:
#   Rscript adc_pasta_rna_correlation.R <output_dir>

suppressPackageStartupMessages({
    library(precedeR)
    library(arrow)
    library(dplyr)
    library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript adc_pasta_rna_correlation.R <output_dir>")
}
output_dir <- args[1]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

leming_dir <- "/efs/home/knguyen/src/leming"
libs_file <- file.path(leming_dir, "RNA/pb_libraries.rna_id.20260514.tsv")
rna_file <- file.path(leming_dir, "RNA/combined_RNA.salmon.gene_tpm.averaged_long.20260514.parquet")
adc_file <- file.path(leming_dir, "workflows/leming_v2/assets/ADC_gene_list.txt")

# --- Load inputs ---
cat("Loading inputs...\n")
libs <- read.delim(libs_file, stringsAsFactors = FALSE)
adc <- read.csv(adc_file, stringsAsFactors = FALSE)
adc_genes <- adc$GENE_NAME

assay_set_ids <- unique(libs$assay_set_id)
cat("  Assay set IDs:", length(assay_set_ids), "\n")
cat("  ADC genes:", length(adc_genes), "\n")

# --- Query PASTA ---
cat("Querying PASTA vd5_pebsnr...\n")
quotes <- paste0("'", assay_set_ids, "'")
gene_quotes <- paste0("'", adc_genes, "'")

sql <- paste0(
    "SELECT assay_set_id, gene, ",
    "combined_zscore_k4, combined_zscore_k27, combined_zscore_k4_k27, ",
    "log2_mean_count_k4, log2_mean_count_k27 ",
    "FROM pasta.vd5_pebsnr WHERE assay_set_id IN (",
    paste(quotes, collapse = ", "), ") ",
    "AND gene IN (", paste(gene_quotes, collapse = ", "), ")"
)

pasta <- queryPBMC(sql)
cat("  PASTA rows:", nrow(pasta), "\n")

missing_pasta <- setdiff(adc_genes, unique(pasta$gene))
if (length(missing_pasta) > 0) {
    cat("  ADC genes missing from PASTA:", paste(missing_pasta, collapse = ", "), "\n")
}

# --- Load RNA ---
cat("Loading RNA TPM...\n")
rna <- read_parquet(rna_file)
rna_adc <- rna %>% filter(gene_name %in% adc_genes)
cat("  RNA rows for ADC genes:", nrow(rna_adc), "\n")

# --- Join ---
cat("Joining PASTA with RNA...\n")
k4_libs <- libs %>%
    filter(Assay == "H3K4me3") %>%
    select(assay_set_id, rna_id) %>%
    distinct()

joined <- pasta %>%
    inner_join(k4_libs, by = "assay_set_id") %>%
    inner_join(rna_adc, by = c("rna_id" = "rna_id", "gene" = "gene_name")) %>%
    mutate(log2_tpm1 = log2(tpm + 1))

cat("  Joined rows:", nrow(joined), "\n")
cat("  Samples:", length(unique(joined$assay_set_id)), "\n")

# --- Per-gene statistics ---
cat("Computing per-gene statistics...\n")
results <- joined %>%
    group_by(gene) %>%
    summarise(
        n_samples = n(),
        cor_k4 = cor(combined_zscore_k4, log2_tpm1, use = "complete.obs"),
        cor_k27 = cor(combined_zscore_k27, log2_tpm1, use = "complete.obs"),
        cor_k4_k27 = cor(combined_zscore_k4_k27, log2_tpm1, use = "complete.obs"),
        mean_tpm = mean(tpm, na.rm = TRUE),
        sd_tpm = sd(tpm, na.rm = TRUE),
        q00_tpm = quantile(tpm, 0.00, na.rm = TRUE),
        q05_tpm = quantile(tpm, 0.05, na.rm = TRUE),
        q25_tpm = quantile(tpm, 0.25, na.rm = TRUE),
        q50_tpm = quantile(tpm, 0.50, na.rm = TRUE),
        q75_tpm = quantile(tpm, 0.75, na.rm = TRUE),
        q95_tpm = quantile(tpm, 0.95, na.rm = TRUE),
        q100_tpm = quantile(tpm, 1.00, na.rm = TRUE),
        pct_below_1tpm = mean(tpm < 1) * 100,
        pct_below_5tpm = mean(tpm < 5) * 100,
        .groups = "drop"
    ) %>%
    left_join(adc %>% select(GENE_NAME, TARGET_NAME), by = c("gene" = "GENE_NAME")) %>%
    arrange(cor_k4)

# --- Write outputs ---
out_tsv <- file.path(output_dir, "adc_pasta_rna_correlation.tsv")
write.table(results, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Wrote", out_tsv, "\n")

out_parquet <- file.path(output_dir, "adc_pasta_rna_joined.parquet")
write_parquet(joined, out_parquet)
cat("  Wrote", out_parquet, "\n")

cat("\nDone.\n")
