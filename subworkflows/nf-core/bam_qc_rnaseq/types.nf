// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
include { Rseqc } from '../bam_rseqc/types'

record BamQcPreseq {
    id:        String
    meta:      Map
    lc_extrap: Path
    log:       Path
}

record BamQcFeaturecounts {
    id:      String
    meta:    Map
    counts:  Path
    summary: Path
}

record BamQcBiotype {
    tsv:  Path
    rrna: Path
}

record BamQcDupradar {
    id:              String
    meta:            Map
    scatter2d:       Path
    boxplot:         Path
    hist:            Path
    dupmatrix:       Path
    intercept_slope: Path
    multiqc:         List<Path>
    session_info:    Path
}

record BamQcRnaseq {
    id:            String
    meta:          Map
    preseq:        BamQcPreseq?
    featurecounts: BamQcFeaturecounts?
    biotype:       BamQcBiotype?
    qualimap:      Path?
    dupradar:      BamQcDupradar?
    rseqc:         Rseqc?
    mqc_files:     List<Path>
}
