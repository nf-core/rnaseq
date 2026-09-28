include { Rseqc } from '../bam_rseqc/types'

record BamQcPreseq {
    lc_extrap: Path
    log:       Path
}

record BamQcFeaturecounts {
    counts:  Path
    summary: Path
}

record BamQcBiotype {
    tsv:  Path
    rrna: Path
}

record BamQcDupradar {
    scatter2d:       Path
    boxplot:         Path
    hist:            Path
    dupmatrix:       Path
    intercept_slope: Path
    multiqc:         List<Path>
}

record BamQcRnaseq {
    id:            String
    preseq:        BamQcPreseq?
    featurecounts: BamQcFeaturecounts?
    biotype:       BamQcBiotype?
    qualimap:      Path?
    dupradar:      BamQcDupradar?
    rseqc:         Rseqc?
}
