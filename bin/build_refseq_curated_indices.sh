#!/usr/bin/env bash
#
# Build STAR and Salmon indices for hg19 ncbiRefSeq Curated GTF.
#
# Runs the nf-core/rnaseq pipeline with a single sample and --save_reference.
# The pipeline must run alignment/pseudo-alignment (not skipped) for the index-
# building logic to trigger, so a small sample is processed end-to-end. The
# primary output of interest is the indices, not the sample results.
#
# Prerequisites:
#   - The curated GTF must already exist (run build_refseq_curated_gtf.sh first)
#   - nextflow, aws CLI configured
#
# Usage:
#   bash build_refseq_curated_indices.sh \
#       --gtf s3://precede-references/Homo_sapiens/hg19/v1/hg19.ncbiRefSeqCurated.gtf.gz \
#       --outdir s3://precede-references/rnaseq/refseq_curated \
#       [--input samplesheet.csv] \
#       [--work-dir s3://your-bucket/work]
#
# Outputs (under <outdir>/genome/index/):
#   star/    - STAR genome index
#   salmon/  - Salmon index (with genome decoys, k=31)

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
FASTA="s3://precede-references/Homo_sapiens/hg19/v0/Homo_sapiens_hg19.fa"
GTF=""
OUTDIR=""
INPUT=""
WORK_DIR=""
PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
usage() {
    cat >&2 <<USAGE
Usage: $0 --gtf <path> --outdir <path> [options]

Required:
  --gtf <path>       Path to the curated GTF (local or S3)
  --outdir <path>    Output directory for indices (local or S3)

Optional:
  --fasta <path>     Genome FASTA (default: precede hg19 reference)
  --input <path>     Samplesheet CSV (default: tests/precede_bio_test.csv)
  --work-dir <path>  Nextflow work directory (default: auto)
USAGE
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --gtf)      GTF="$2";      shift 2 ;;
        --fasta)    FASTA="$2";    shift 2 ;;
        --outdir)   OUTDIR="$2";   shift 2 ;;
        --input)    INPUT="$2";    shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        -h|--help)  usage ;;
        *)          echo "Unknown argument: $1" >&2; usage ;;
    esac
done

if [[ -z "$GTF" || -z "$OUTDIR" ]]; then
    echo "ERROR: --gtf and --outdir are required" >&2
    usage
fi

if [[ -z "$INPUT" ]]; then
    INPUT="${PIPELINE_DIR}/tests/precede_bio_test.csv"
fi

# ---------------------------------------------------------------------------
# Build the Nextflow command
# ---------------------------------------------------------------------------
# The pipeline only builds indices when alignment/pseudo-alignment is enabled,
# so we cannot use --skip_alignment or --skip_pseudo_alignment. Instead we run
# the full pipeline with a single sample. The indices are published under
# <outdir>/genome/index/ via --save_reference.
#
# We override --star_index and --salmon_index to false so the precede_bio
# profile's pre-built indices are ignored and new ones are generated.

NF_ARGS=(
    nextflow run "${PIPELINE_DIR}"
    -profile precede_bio
    --input "${INPUT}"
    --outdir "${OUTDIR}"
    --fasta "${FASTA}"
    --gtf "${GTF}"
    --star_index false
    --salmon_index false
    --save_reference
    --skip_deseq2_qc
)

if [[ -n "$WORK_DIR" ]]; then
    NF_ARGS+=(-w "${WORK_DIR}")
fi

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
echo "=============================================="
echo "  Build indices for hg19 RefSeq Curated GTF"
echo "=============================================="
echo ""
echo "Pipeline: ${PIPELINE_DIR}"
echo "FASTA:    ${FASTA}"
echo "GTF:      ${GTF}"
echo "Input:    ${INPUT}"
echo "Output:   ${OUTDIR}"
echo ""
echo "Note: The pipeline runs end-to-end with one sample because index"
echo "building is conditional on alignment/pseudo-alignment being enabled."
echo ""
echo "Command:"
echo "  ${NF_ARGS[*]}"
echo ""

"${NF_ARGS[@]}"

echo ""
echo "=============================================="
echo "  Done. Indices published to:"
echo "    STAR:   ${OUTDIR}/genome/index/star/"
echo "    Salmon: ${OUTDIR}/genome/index/salmon/"
echo ""
echo "  To use in future runs, set:"
echo "    --star_index   ${OUTDIR}/genome/index/star"
echo "    --salmon_index ${OUTDIR}/genome/index/salmon"
echo "=============================================="
