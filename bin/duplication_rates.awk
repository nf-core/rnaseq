#!/usr/bin/awk -f
# AWK script to calculate duplication rates from featureCounts output
# Usage: awk -f duplication_rates.awk -v output_prefix="sample" with_dups.txt no_dups.txt

BEGIN { 
    OFS="\t"
    print "ID", "Length", "Counts", "CountsNodup", "DupRate", "RPK", "RPKM"
}

# First pass: read with_dups file
FNR==NR && NR>2 { 
    with_dups[$1] = $7
    lengths[$1] = $6
    next 
}

# Second pass: read no_dups file and calculate metrics
FNR!=NR && NR>2 { 
    gene = $1
    length = lengths[gene]
    counts_with_dups = with_dups[gene]
    counts_no_dups = $7
    
    if (counts_with_dups > 0) {
        dup_rate = (counts_with_dups - counts_no_dups) / counts_with_dups
    } else {
        dup_rate = 0
    }
    
    if (length > 0) {
        rpk = counts_no_dups / (length / 1000)
    } else {
        rpk = 0
    }
    
    rpkm = rpk  # Simplified - could add total reads normalization if needed
    
    print gene, length, counts_with_dups, counts_no_dups, dup_rate, rpk, rpkm
}