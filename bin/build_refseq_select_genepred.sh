#!/usr/bin/env bash
#
# Build an ncbiRefSeqSelect-compatible genePred file from the NCBI GRCh37.p13
# GTF, filtered by the "RefSeq Select" tag.
#
# This produces a file format-compatible with the UCSC ncbiRefSeqSelect.txt.gz
# track, but derived from the NCBI annotation release so that gene names are
# consistent with the BestRefSeq curated GTF used by the rnaseq pipeline.
#
# Requires: Python 3.6+, curl, md5sum, gzip
#
# Usage:
#   bash build_refseq_select_genepred.sh <output_dir>
#
# Outputs:
#   <output_dir>/ncbiRefSeqSelect.txt.gz   - genePred file (drop-in for pasta)
#   <output_dir>/build_manifest.txt         - provenance metadata

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
WORKDIR=$(mktemp -d "${OUTDIR}/build_select_XXXXXX")
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
# 4. Convert GTF RefSeq Select entries to genePred format
# ---------------------------------------------------------------------------
echo "[4/5] Converting RefSeq Select GTF entries to genePred format..."
gunzip -c "${WORKDIR}/source.gtf.gz" > "${WORKDIR}/source.gtf"

cat > "${WORKDIR}/gtf_to_genepred.py" <<'PYEOF'
import sys
import re
from collections import defaultdict

gtf_path = sys.argv[1]
chr_map_path = sys.argv[2]
out_path = sys.argv[3]

# Load chromosome map, restricted to standard chroms
standard = set(['chr%d' % i for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])
chr_map = {}
with open(chr_map_path) as f:
    for line in f:
        refseq, ucsc = line.strip().split('\t')
        if ucsc in standard:
            chr_map[refseq] = ucsc

# Parse attributes from a GTF line
attr_pat = re.compile(r'(\w+) "([^"]*)"')

def parse_attrs(attr_str):
    return dict(attr_pat.findall(attr_str))

# Transcript overrides: use a different BestRefSeq transcript instead of RefSeq Select.
# key = gene name, value = transcript accession (without version) to use instead.
TX_OVERRIDES = {
    'CLDN18': 'NM_001002026',  # CLDN18.2 (gastric/ADC target) instead of CLDN18.1
}

# First pass: collect gene names that have BestRefSeq curated transcripts
# on primary/standard chromosomes (matching what the rnaseq GTF and leming RDS include)
curated_genes = set()
with open(gtf_path) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[1] == 'BestRefSeq' and fields[2] == 'transcript':
            if fields[0] not in chr_map:
                continue
            attrs = parse_attrs(fields[8])
            gene = attrs.get('gene', '')
            if gene:
                curated_genes.add(gene)
print('BestRefSeq curated genes (primary chroms): %d' % len(curated_genes))

# Collect transcript and exon data for RefSeq Select transcripts
transcripts = {}  # tx_id -> {chrom, strand, txStart, txEnd, gene, cdsStart, cdsEnd, exons}
exons = defaultdict(list)  # tx_id -> [(start, end)]
cds_regions = defaultdict(list)  # tx_id -> [(start, end)]
skipped_genes = set()
override_targets = set(TX_OVERRIDES.values())

# Second pass: RefSeq Select transcripts (skipping genes that have overrides)
with open(gtf_path) as f:
    for line in f:
        if line.startswith('#'):
            continue
        if 'RefSeq Select' not in line:
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        if chrom not in chr_map:
            continue
        feature = fields[2]
        start = int(fields[3]) - 1  # GTF is 1-based, genePred is 0-based
        end = int(fields[4])
        strand = fields[6]
        attrs = parse_attrs(fields[8])
        tx_id = attrs.get('transcript_id', '')
        if not tx_id:
            continue

        if feature == 'transcript':
            gene = attrs.get('gene', attrs.get('gene_id', ''))
            if gene not in curated_genes:
                skipped_genes.add(gene)
                continue
            if gene in TX_OVERRIDES:
                continue
            transcripts[tx_id] = {
                'chrom': chr_map[chrom],
                'strand': strand,
                'txStart': start,
                'txEnd': end,
                'gene': gene,
            }
        elif feature == 'exon':
            exons[tx_id].append((start, end))
        elif feature == 'CDS':
            cds_regions[tx_id].append((start, end))

# Third pass: collect override transcripts from BestRefSeq
with open(gtf_path) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[1] != 'BestRefSeq':
            continue
        chrom = fields[0]
        if chrom not in chr_map:
            continue
        attrs = parse_attrs(fields[8])
        tx_id = attrs.get('transcript_id', '')
        tx_base = tx_id.rsplit('.', 1)[0]
        if tx_base not in override_targets:
            continue
        feature = fields[2]
        start = int(fields[3]) - 1
        end = int(fields[4])
        strand = fields[6]

        if feature == 'transcript':
            gene = attrs.get('gene', attrs.get('gene_id', ''))
            transcripts[tx_id] = {
                'chrom': chr_map[chrom],
                'strand': strand,
                'txStart': start,
                'txEnd': end,
                'gene': gene,
            }
            print('Override: %s -> %s (%s)' % (gene, tx_id, chr_map[chrom]))
        elif feature == 'exon':
            exons[tx_id].append((start, end))
        elif feature == 'CDS':
            cds_regions[tx_id].append((start, end))

# Write genePred output
# Columns: bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd,
#          exonCount, exonStarts, exonEnds, score, name2, cdsStartStat,
#          cdsEndStat, exonFrames
with open(out_path, 'w') as f:
    for tx_id in sorted(transcripts.keys()):
        tx = transcripts[tx_id]
        tx_exons = sorted(exons.get(tx_id, []))
        if not tx_exons:
            continue

        tx_cds = sorted(cds_regions.get(tx_id, []))
        if tx_cds:
            cds_start = tx_cds[0][0]
            cds_end = tx_cds[-1][1]
            cds_start_stat = 'cmpl'
            cds_end_stat = 'cmpl'
        else:
            # Non-coding: cdsStart == cdsEnd == txEnd (UCSC convention)
            cds_start = tx['txEnd']
            cds_end = tx['txEnd']
            cds_start_stat = 'none'
            cds_end_stat = 'none'

        exon_count = len(tx_exons)
        exon_starts = ','.join(str(e[0]) for e in tx_exons) + ','
        exon_ends = ','.join(str(e[1]) for e in tx_exons) + ','

        # Compute exon frames for CDS
        if tx_cds:
            frames = []
            if tx['strand'] == '+':
                running_frame = 0
                for es, ee in tx_exons:
                    # Does this exon overlap CDS?
                    cds_in_exon_start = max(es, cds_start)
                    cds_in_exon_end = min(ee, cds_end)
                    if cds_in_exon_start < cds_in_exon_end:
                        frames.append(str(running_frame))
                        cds_len = cds_in_exon_end - cds_in_exon_start
                        running_frame = (3 - (cds_len - running_frame) % 3) % 3
                    else:
                        frames.append('-1')
            else:
                running_frame = 0
                for es, ee in reversed(tx_exons):
                    cds_in_exon_start = max(es, cds_start)
                    cds_in_exon_end = min(ee, cds_end)
                    if cds_in_exon_start < cds_in_exon_end:
                        frames.insert(0, str(running_frame))
                        cds_len = cds_in_exon_end - cds_in_exon_start
                        running_frame = (3 - (cds_len - running_frame) % 3) % 3
                    else:
                        frames.insert(0, '-1')
            exon_frames = ','.join(frames) + ','
        else:
            exon_frames = ','.join(['-1'] * exon_count) + ','

        # bin field: use 0 as placeholder (pasta doesn't use it)
        row = [
            '0',
            tx_id,
            tx['chrom'],
            tx['strand'],
            str(tx['txStart']),
            str(tx['txEnd']),
            str(cds_start),
            str(cds_end),
            str(exon_count),
            exon_starts,
            exon_ends,
            '0',
            tx['gene'],
            cds_start_stat,
            cds_end_stat,
            exon_frames,
        ]
        f.write('\t'.join(row) + '\n')

print('Wrote %d transcripts' % len([t for t in transcripts if exons.get(t)]))
if skipped_genes:
    print('Skipped %d genes with no BestRefSeq curated transcript: %s' % (
        len(skipped_genes), ', '.join(sorted(skipped_genes))))
PYEOF

python3 -I "${WORKDIR}/gtf_to_genepred.py" "${WORKDIR}/source.gtf" "${WORKDIR}/chr_map.txt" "${WORKDIR}/output.txt"

NUM_TX=$(wc -l < "${WORKDIR}/output.txt")
echo "  ${NUM_TX} transcripts written"

# ---------------------------------------------------------------------------
# 5. Compress and write manifest
# ---------------------------------------------------------------------------
echo "[5/5] Compressing and writing manifest..."
gzip -c "${WORKDIR}/output.txt" > "${OUTDIR}/ncbiRefSeqSelect.txt.gz"

OUTPUT_MD5=$(md5sum "${OUTDIR}/ncbiRefSeqSelect.txt.gz" | awk '{print $1}')

cat > "${OUTDIR}/build_manifest.txt" <<EOF
# hg19 NCBI RefSeq Select genePred - Build Manifest
# Generated: $(date -u '+%Y-%m-%dT%H:%M:%SZ')

source_url:         ${GTF_URL}
source_md5_gz:      ${EXPECTED_GZ_MD5}
assembly_report:    ${REPORT_URL}
annotation_release: 105.20220307
genome_build:       GRCh37.p13 (hg19)

filter_criteria:    tag "RefSeq Select" in NCBI GTF, with transcript overrides (see below)
chromosomes:        Standard only (chr1-22, chrX, chrY, chrM)
chromosome_names:   Converted from RefSeq accessions to UCSC (chr*) via assembly report
output_format:      UCSC extended genePred (16 columns), compatible with ncbiRefSeqSelect track

output_file:        ncbiRefSeqSelect.txt.gz
output_md5_gz:      ${OUTPUT_MD5}
total_transcripts:  ${NUM_TX}

transcript_overrides:
  CLDN18:     NM_001002026 (CLDN18.2, gastric/ADC target) replaces NM_016369 (CLDN18.1, lung)
EOF

echo ""
echo "Done."
echo "  genePred: ${OUTDIR}/ncbiRefSeqSelect.txt.gz"
echo "  Manifest: ${OUTDIR}/build_manifest.txt"
