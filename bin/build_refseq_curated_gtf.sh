#!/usr/bin/env bash
#
# Build a curated RefSeq GTF for hg19 from NCBI annotation release 105.20220307.
#
# Downloads the official NCBI GRCh37.p13 GTF, verifies its checksum, filters to
# BestRefSeq (curated) entries on primary-assembly + mitochondrial chromosomes,
# and converts RefSeq accession names to UCSC-style chromosome names.
#
# The chromosome mapping is derived from the NCBI assembly report for
# GCF_000001405.25 (GRCh37.p13).
#
# Usage:
#   bash build_refseq_curated_gtf.sh <output_dir>
#
# Outputs:
#   <output_dir>/hg19.ncbiRefSeqCurated.gtf.gz   - final curated GTF
#   <output_dir>/build_manifest.txt               - provenance metadata

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
WORKDIR=$(mktemp -d "${OUTDIR}/build_refseq_XXXXXX")
trap 'rm -rf "${WORKDIR}"' EXIT

echo "Working directory: ${WORKDIR}"
echo "Output directory:  ${OUTDIR}"
echo ""

# ---------------------------------------------------------------------------
# 1. Download the GTF and assembly report
# ---------------------------------------------------------------------------
echo "[1/5] Downloading NCBI GRCh37.p13 GTF..."
curl -fSL --retry 3 -o "${WORKDIR}/source.gtf.gz" "${GTF_URL}"

echo "[1/5] Downloading assembly report..."
curl -fSL --retry 3 -o "${WORKDIR}/assembly_report.txt" "${REPORT_URL}"

# ---------------------------------------------------------------------------
# 2. Verify MD5
# ---------------------------------------------------------------------------
echo "[2/5] Verifying MD5 checksum of downloaded GTF..."
ACTUAL_MD5=$(md5sum "${WORKDIR}/source.gtf.gz" | awk '{print $1}')
if [[ "${ACTUAL_MD5}" != "${EXPECTED_GZ_MD5}" ]]; then
    echo "ERROR: MD5 mismatch" >&2
    echo "  Expected: ${EXPECTED_GZ_MD5}" >&2
    echo "  Got:      ${ACTUAL_MD5}" >&2
    exit 1
fi
echo "  MD5 OK: ${ACTUAL_MD5}"

# ---------------------------------------------------------------------------
# 3. Build the RefSeq-accession -> UCSC chromosome name map
#    from the assembly report. Keep primary assembly + non-nuclear (chrM).
# ---------------------------------------------------------------------------
echo "[3/5] Building chromosome name map from assembly report..."
awk -F'\t' '
    !/^#/ && ($8 == "Primary Assembly" || $8 == "non-nuclear") {
        gsub(/\r/, "")
        print $7 "\t" $10
    }
' "${WORKDIR}/assembly_report.txt" > "${WORKDIR}/chr_map.txt"

NUM_MAPPED=$(wc -l < "${WORKDIR}/chr_map.txt")
echo "  Mapped ${NUM_MAPPED} sequences (primary + mito)"

# ---------------------------------------------------------------------------
# 4. Filter to BestRefSeq on mapped chromosomes, rename chromosomes
# ---------------------------------------------------------------------------
echo "[4/5] Filtering to BestRefSeq (curated) and renaming chromosomes..."

# Build an awk-friendly lookup from the map file
gunzip -c "${WORKDIR}/source.gtf.gz" | awk -F'\t' -v OFS='\t' '
    # First pass: load the chromosome map
    NR == FNR {
        map[$1] = $2
        next
    }
    # Header lines
    /^#/ {
        print
        next
    }
    # Keep only BestRefSeq on mapped chromosomes
    $2 == "BestRefSeq" && ($1 in map) {
        $1 = map[$1]
        print
    }
' "${WORKDIR}/chr_map.txt" - > "${WORKDIR}/filtered_raw.gtf"

# Propagate gene_biotype from gene lines to transcript/exon/CDS lines.
# NCBI puts gene_biotype only on gene features; featureCounts biotype QC
# needs it on exon lines.
awk -F'\t' -v OFS='\t' '
    /^#/ { print; next }
    $3 == "gene" {
        match($9, /gene_biotype "([^"]*)"/, m)
        match($9, /gene_id "([^"]*)"/, gid)
        if (m[1] && gid[1]) biotype[gid[1]] = m[1]
        print
        next
    }
    {
        match($9, /gene_id "([^"]*)"/, gid)
        if (gid[1] && (gid[1] in biotype) && $9 !~ /gene_biotype/) {
            sub(/;[[:space:]]*$/, "; gene_biotype \"" biotype[gid[1]] "\"; ", $9)
        }
        print
    }
' "${WORKDIR}/filtered_raw.gtf" > "${WORKDIR}/filtered.gtf"

TOTAL_LINES=$(grep -cv '^#' "${WORKDIR}/filtered.gtf" || true)
NUM_GENES=$(awk -F'\t' '$3 == "gene"' "${WORKDIR}/filtered.gtf" | wc -l)
NUM_TRANSCRIPTS=$(awk -F'\t' '$3 == "transcript"' "${WORKDIR}/filtered.gtf" | wc -l)
echo "  Curated GTF: ${TOTAL_LINES} feature lines, ${NUM_GENES} genes, ${NUM_TRANSCRIPTS} transcripts"

# ---------------------------------------------------------------------------
# 5. Compress and write manifest
# ---------------------------------------------------------------------------
echo "[5/5] Compressing and writing manifest..."
gzip -c "${WORKDIR}/filtered.gtf" > "${OUTDIR}/hg19.ncbiRefSeqCurated.gtf.gz"

OUTPUT_MD5=$(md5sum "${OUTDIR}/hg19.ncbiRefSeqCurated.gtf.gz" | awk '{print $1}')

cat > "${OUTDIR}/build_manifest.txt" <<EOF
# hg19 NCBI RefSeq Curated GTF - Build Manifest
# Generated: $(date -u '+%Y-%m-%dT%H:%M:%SZ')

source_url:       ${GTF_URL}
source_md5_gz:    ${EXPECTED_GZ_MD5}
source_md5_actual:${ACTUAL_MD5}
assembly_report:  ${REPORT_URL}
annotation_release: 105.20220307
genome_build:     GRCh37.p13 (hg19)

filter_criteria:  source == "BestRefSeq" (NCBI curated transcripts)
chromosomes:      Primary Assembly + non-nuclear (chrM)
chromosome_names: Converted from RefSeq accessions to UCSC (chr*) via assembly report

output_file:      hg19.ncbiRefSeqCurated.gtf.gz
output_md5_gz:    ${OUTPUT_MD5}
total_features:   ${TOTAL_LINES}
total_genes:      ${NUM_GENES}
total_transcripts:${NUM_TRANSCRIPTS}
EOF

echo ""
echo "Done."
echo "  GTF:      ${OUTDIR}/hg19.ncbiRefSeqCurated.gtf.gz"
echo "  Manifest: ${OUTDIR}/build_manifest.txt"
