#!/usr/bin/env bash
#
# Build gene_locations.hg19.single_range.rds from the NCBI GRCh37.p13
# BestRefSeq curated GTF.
#
# Each gene gets a single primary range (the locus containing its TSS) plus
# an all_ranges field with every distinct locus (merged overlapping transcripts).
# This is a drop-in replacement for the TxDb.Hsapiens.UCSC.hg19.knownGene-derived
# version used by leming.
#
# Requires: R (with GenomicRanges/IRanges from Bioconductor), curl, md5sum
#
# Usage:
#   bash build_gene_locations_rds.sh <output_dir>
#
# Outputs:
#   <output_dir>/gene_locations.hg19.single_range.rds  - R list (GRanges-based)
#   <output_dir>/build_manifest.txt                    - provenance metadata

set -euo pipefail

# ---------------------------------------------------------------------------
# Source
# ---------------------------------------------------------------------------
NCBI_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/105.20220307/GCF_000001405.25_GRCh37.p13"
GTF_URL="${NCBI_BASE}/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz"
REPORT_URL="${NCBI_BASE}/GCF_000001405.25_GRCh37.p13_assembly_report.txt"
EXPECTED_GZ_MD5="7877d9e1383dc194c4b22d82ae07d759"

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <output_dir>" >&2
    exit 1
fi

OUTDIR="$1"
mkdir -p "${OUTDIR}"
WORKDIR=$(mktemp -d "${OUTDIR}/build_geneloc_XXXXXX")
trap 'rm -rf "${WORKDIR}"' EXIT

echo "Working directory: ${WORKDIR}"
echo "Output directory:  ${OUTDIR}"
echo ""

# ---------------------------------------------------------------------------
# 1. Download
# ---------------------------------------------------------------------------
echo "[1/5] Downloading NCBI GRCh37.p13 GTF and assembly report..."
curl -fSL --retry 3 -o "${WORKDIR}/source.gtf.gz" "${GTF_URL}"
curl -fSL --retry 3 -o "${WORKDIR}/assembly_report.txt" "${REPORT_URL}"

# ---------------------------------------------------------------------------
# 2. Verify MD5
# ---------------------------------------------------------------------------
echo "[2/5] Verifying MD5 checksum..."
ACTUAL_MD5=$(md5sum "${WORKDIR}/source.gtf.gz" | awk '{print $1}')
if [[ "${ACTUAL_MD5}" != "${EXPECTED_GZ_MD5}" ]]; then
    echo "ERROR: MD5 mismatch (expected ${EXPECTED_GZ_MD5}, got ${ACTUAL_MD5})" >&2
    exit 1
fi
echo "  MD5 OK: ${ACTUAL_MD5}"

# ---------------------------------------------------------------------------
# 3. Build chromosome name map
# ---------------------------------------------------------------------------
echo "[3/5] Building chromosome name map..."
awk -F'\t' '
    !/^#/ && ($8 == "Primary Assembly" || $8 == "non-nuclear") {
        gsub(/\r/, "")
        print $7 "\t" $10
    }
' "${WORKDIR}/assembly_report.txt" > "${WORKDIR}/chr_map.txt"

# ---------------------------------------------------------------------------
# 4. Extract BestRefSeq transcripts
# ---------------------------------------------------------------------------
echo "[4/5] Extracting BestRefSeq transcript coordinates..."
gunzip -c "${WORKDIR}/source.gtf.gz" \
  | awk -F'\t' '$2 == "BestRefSeq" && $3 == "transcript"' \
  > "${WORKDIR}/curated_transcripts.gtf"

NUM_TX=$(wc -l < "${WORKDIR}/curated_transcripts.gtf")
echo "  ${NUM_TX} BestRefSeq transcripts found"

# ---------------------------------------------------------------------------
# 5. Build RDS in R
# ---------------------------------------------------------------------------
echo "[5/5] Building gene_locations RDS..."

cat > "${WORKDIR}/build_rds.R" <<'REOF'
suppressPackageStartupMessages(library(GenomicRanges))

args <- commandArgs(trailingOnly = TRUE)
tx_gtf_path  <- args[1]
chr_map_path <- args[2]
out_path     <- args[3]

# --- Load chromosome name map ---
chr_df <- read.delim(chr_map_path, header = FALSE, col.names = c("refseq", "ucsc"),
                     stringsAsFactors = FALSE)
chr_map <- setNames(chr_df$ucsc, chr_df$refseq)

# --- Parse transcript GTF (vectorized) ---
cat("Reading transcript GTF...\n")
lines <- readLines(tx_gtf_path)
fields <- strsplit(lines, "\t", fixed = TRUE)

gtf <- data.frame(
    chrom      = vapply(fields, `[`, character(1), 1),
    start      = as.integer(vapply(fields, `[`, character(1), 4)),
    end        = as.integer(vapply(fields, `[`, character(1), 5)),
    strand     = vapply(fields, `[`, character(1), 7),
    attributes = vapply(fields, `[`, character(1), 9),
    stringsAsFactors = FALSE
)

# Extract gene name — '; gene "' avoids matching gene_id
m <- regexec('; gene "([^"]*)"', gtf$attributes)
gtf$gene <- vapply(regmatches(gtf$attributes, m), function(x) {
    if (length(x) >= 2) x[2] else NA_character_
}, character(1))

# Map chromosomes and keep only standard chroms (chr1-22, X, Y, M)
gtf$ucsc_chrom <- chr_map[gtf$chrom]
standard_chroms <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
gtf <- gtf[!is.na(gtf$ucsc_chrom) & gtf$ucsc_chrom %in% standard_chroms, ]
cat("  ", nrow(gtf), "transcripts on standard chromosomes\n")

# --- Build per-gene GRanges, then split ---
tx_gr <- GRanges(
    seqnames = gtf$ucsc_chrom,
    ranges   = IRanges(start = gtf$start, end = gtf$end),
    strand   = "*",
    gene     = gtf$gene,
    tx_strand = gtf$strand
)

by_gene <- split(tx_gr, tx_gr$gene)
cat("  ", length(by_gene), "unique genes\n")

cat("Building per-gene loci...\n")
result <- vector("list", length(by_gene))
names(result) <- names(by_gene)

for (gname in names(by_gene)) {
    g <- by_gene[[gname]]
    chrom <- as.character(seqnames(g))[1]
    strand_val <- g$tx_strand[1]

    loci <- reduce(g)

    if (strand_val == "+") {
        tss <- min(start(g))
    } else {
        tss <- max(end(g))
    }

    tss_hits <- which(start(loci) <= tss & end(loci) >= tss)
    primary_idx <- if (length(tss_hits) == 0) 1L else tss_hits[1]
    primary <- loci[primary_idx]

    result[[gname]] <- list(
        gene       = primary,
        chr        = chrom,
        strand     = strand_val,
        tss        = tss,
        all_ranges = loci
    )
}

result <- result[sort(names(result))]
saveRDS(result, out_path)

multi <- sum(sapply(result, function(e) length(e$all_ranges) > 1))
cat("\nSummary:\n")
cat("  Total genes:       ", length(result), "\n")
cat("  Multi-locus genes: ", multi, "\n")
cat("  Chromosomes:       ", length(unique(sapply(result, function(e) e$chr))), "\n")
cat("  Plus strand:       ", sum(sapply(result, function(e) e$strand == "+")), "\n")
cat("  Minus strand:      ", sum(sapply(result, function(e) e$strand == "-")), "\n")
REOF

Rscript --vanilla "${WORKDIR}/build_rds.R" \
    "${WORKDIR}/curated_transcripts.gtf" \
    "${WORKDIR}/chr_map.txt" \
    "${OUTDIR}/gene_locations.hg19.single_range.rds"

# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
OUTPUT_MD5=$(md5sum "${OUTDIR}/gene_locations.hg19.single_range.rds" | awk '{print $1}')
NUM_GENES=$(Rscript --vanilla -e 'cat(length(readRDS(commandArgs(TRUE))))' \
    "${OUTDIR}/gene_locations.hg19.single_range.rds")

cat > "${OUTDIR}/build_manifest.txt" <<EOF
# hg19 Gene Locations RDS - Build Manifest
# Generated: $(date -u '+%Y-%m-%dT%H:%M:%SZ')

source_url:         ${GTF_URL}
source_md5_gz:      ${EXPECTED_GZ_MD5}
assembly_report:    ${REPORT_URL}
annotation_release: 105.20220307
genome_build:       GRCh37.p13 (hg19)

filter_criteria:    source "BestRefSeq" + feature "transcript" in NCBI GTF
chromosomes:        Primary Assembly + non-nuclear (chrM)
chromosome_names:   Converted from RefSeq accessions to UCSC (chr*) via assembly report
gene_range:         Per gene: reduce(transcript spans) -> distinct loci; primary locus = TSS-containing
output_format:      R named list, each element: gene (GRanges), chr, strand, tss, all_ranges (GRanges)

output_file:        gene_locations.hg19.single_range.rds
output_md5:         ${OUTPUT_MD5}
total_genes:        ${NUM_GENES}
EOF

echo ""
echo "Done."
echo "  RDS:      ${OUTDIR}/gene_locations.hg19.single_range.rds"
echo "  Manifest: ${OUTDIR}/build_manifest.txt"
