
process makeSTARindex {
    tag "$fasta"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file fasta 
    file gtf

    output:
    file "star"

    script:
    def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir star
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --sjdbGTFfile $gtf \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        $avail_mem
    """
}

process makeHisatSplicesites {
    tag "$gtf"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file gtf

    output:
    file "${gtf.baseName}.hisat2_splice_sites.txt"

    script:
    """
    hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
    """
}

process makeHISATindex {
    tag "$fasta"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file fasta
    file indexing_splicesites
    file gtf

    output:
    file "${fasta.baseName}.*.ht2"

    script:
    if( !task.memory ){
        log.info "[HISAT2 index build] Available memory not known - defaulting to 0. Specify process memory requirements to change this."
        avail_mem = 0
    } else {
        log.info "[HISAT2 index build] Available memory: ${task.memory}"
        avail_mem = task.memory.toGiga()
    }
    if( avail_mem > params.hisatBuildMemory ){
        log.info "[HISAT2 index build] Over ${params.hisatBuildMemory} GB available, so using splice sites and exons in HISAT2 index"
        extract_exons = "hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt"
        ss = "--ss $indexing_splicesites"
        exon = "--exon ${gtf.baseName}.hisat2_exons.txt"
    } else {
        log.info "[HISAT2 index build] Less than ${params.hisatBuildMemory} GB available, so NOT using splice sites and exons in HISAT2 index."
        log.info "[HISAT2 index build] Use --hisatBuildMemory [small number] to skip this check."
        extract_exons = ''
        ss = ''
        exon = ''
    }
    """
    $extract_exons
    hisat2-build -p ${task.cpus} $ss $exon $fasta ${fasta.baseName}.hisat2_index
    """
}

process convertGFFtoGTF {
    tag "$gff"

    input:
    file gff

    output:
    file "${gff.baseName}.gtf"

    script:
    """
    gffread  $gff -T -o ${gff.baseName}.gtf
    """
}

process makeBED12 {
    tag "$gtf"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file gtf

    output:
    file "${gtf.baseName}.bed"

    script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
    """
    gtf2bed $gtf > ${gtf.baseName}.bed
    """
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    set val(name), file(reads)

    output:
    file "*_fastqc.{zip,html}"

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (!params.saveTrimmed && filename == "where_are_my_files.txt") filename
            else if (params.saveTrimmed && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    set val(name), file(reads)
    file wherearemyfiles

    output:
    file "*fq.gz"
    file "*trimming_report.txt"
    file "*_fastqc.{zip,html}" 
    file "where_are_my_files.txt"


    script:
    c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
    tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
    if (params.singleEnd) {
        """
        trim_galore --fastqc --gzip $c_r1 $tpc_r1 $reads
        """
    } else {
        """
        trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        """
    }
}

process star {
    tag "$prefix"
    publishDir "${params.outdir}/STAR", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") == -1) "logs/$filename"
            else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
            else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    file reads
    file index
    file gtf
    file wherearemyfiles

    output:
    set file("*Log.final.out"), file ('*.bam')
    file "*.out"
    file "*Log.out"
    file "${prefix}Aligned.sortedByCoord.out.bam.bai"
    file "*SJ.out.tab"
    file "where_are_my_files.txt"

    script:
    prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    def star_mem = task.memory ?: params.star_memory ?: false
    def avail_mem = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 100000000}" : ''
    seqCenter = params.seqCenter ? "--outSAMattrRGline ID:$prefix 'CN:$params.seqCenter'" : ''
    """
    STAR --genomeDir $index \\
        --sjdbGTFfile $gtf \\
        --readFilesIn $reads  \\
        --runThreadN ${task.cpus} \\
        --twopassMode Basic \\
        --outWigType bedGraph \\
        --outSAMtype BAM SortedByCoordinate $avail_mem \\
        --readFilesCommand zcat \\
        --runDirPerm All_RWX \\
            --outFileNamePrefix $prefix $seqCenter
        
    samtools index ${prefix}Aligned.sortedByCoord.out.bam
    """
}

process hisat2Align {
    tag "$prefix"
    publishDir "${params.outdir}/HISAT2", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
            else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
            else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    file reads
    file hs2_indices
    file alignment_splicesites 
    file wherearemyfiles 

    output:
    file "${prefix}.bam"
    file "${prefix}.hisat2_summary.txt"
    file "where_are_my_files.txt"

    script:
    index_base = hs2_indices[0].toString() - ~/.\d.ht2/
    prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    seqCenter = params.seqCenter ? "--rg-id ${prefix} --rg CN:${params.seqCenter.replaceAll('\\s','_')}" : ''
    def rnastrandness = ''
    if (forward_stranded && !unstranded){
        rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (reverse_stranded && !unstranded){
        rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
    }
    if (params.singleEnd) {
        """
        hisat2 -x $index_base \\
                -U $reads \\
                $rnastrandness \\
                --known-splicesite-infile $alignment_splicesites \\
                -p ${task.cpus} \\
                --met-stderr \\
                --new-summary \\
                --summary-file ${prefix}.hisat2_summary.txt $seqCenter \\
                | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
        """
    } else {
        """
        hisat2 -x $index_base \\
                -1 ${reads[0]} \\
                -2 ${reads[1]} \\
                $rnastrandness \\
                --known-splicesite-infile $alignment_splicesites \\
                --no-mixed \\
                --no-discordant \\
                -p ${task.cpus} \\
                --met-stderr \\
                --new-summary \\
                --summary-file ${prefix}.hisat2_summary.txt $seqCenter \\
                | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
        """
    }
}

process hisat2_sort {
    tag "${hisat2_bam.baseName}"
    publishDir "${params.outdir}/HISAT2", mode: 'copy',
        saveAs: { filename ->
            if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
            else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") "aligned_sorted/$filename"
            else null
        }

    input:
    file hisat2_bam
    file wherearemyfiles

    output:
    file "${hisat2_bam.baseName}.sorted.bam"
    file "${hisat2_bam.baseName}.sorted.bam.bai"
    file "where_are_my_files.txt"

    script:
    def avail_mem = task.memory ? "-m ${task.memory.toBytes() / task.cpus}" : ''
    """
    samtools sort \\
        $hisat2_bam \\
        -@ ${task.cpus} $avail_mem \\
        -o ${hisat2_bam.baseName}.sorted.bam
    samtools index ${hisat2_bam.baseName}.sorted.bam
    """
}


// === from here 


/*
 * STEP 4 - RSeQC analysis
 */
process rseqc {
    tag "${bam_rseqc.baseName - '.sorted'}"
    publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
            else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else filename
        }

    when:
    !params.skip_qc && !params.skip_rseqc

    input:
    file bam_rseqc
    file index
    file bed12

    output:
    file "*.{txt,pdf,r,xls}"

    script:
    """
    infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.infer_experiment.txt
    junction_annotation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    bam_stat.py -i $bam_rseqc 2> ${bam_rseqc.baseName}.bam_stat.txt
    junction_saturation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12 2> ${bam_rseqc.baseName}.junction_annotation_log.txt
    inner_distance.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    read_distribution.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc.baseName}.read_duplication
    """
}


/*
 * Step 4.1 Rseqc create BigWig coverage
 */

process createBigWig {
    tag "${bam.baseName - 'sortedByCoord.out'}"
    publishDir "${params.outdir}/bigwig", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_genebody_coverage

    input:
    file bam
    file index

    output:
    file "*.bigwig"

    script:
    """
    bamCoverage -b $bam -p ${task.cpus} -o ${bam.baseName}.bigwig
    """
}
/*
 * Step 4.2 Rseqc genebody_coverage
 */
process genebody_coverage {
    tag "${bigwig.baseName}"
       publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
            else if (filename.indexOf("geneBodyCoverage.r") > 0)           "geneBodyCoverage/rscripts/$filename"
            else if (filename.indexOf("geneBodyCoverage.txt") > 0)         "geneBodyCoverage/data/$filename"
            else if (filename.indexOf("log.txt") > -1) false
            else filename
        }

    when:
    !params.skip_qc && !params.skip_genebody_coverage

    input:
    file bigwig
    file bed12

    output:
    file "*.{txt,pdf,r}"

    script:
    """
    geneBody_coverage2.py \\
        -i $bigwig \\
        -o ${bigwig.baseName}.rseqc.txt \\
        -r $bed12
    """
}

/*
 * STEP 5 - preseq analysis
 */
process preseq {
    tag "${bam_preseq.baseName - '.sorted'}"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_preseq

    input:
    file bam_preseq

    output:
    file "${bam_preseq.baseName}.ccurve.txt"

    script:
    """
    preseq lc_extrap -v -B $bam_preseq -o ${bam_preseq.baseName}.ccurve.txt
    """
}


/*
 * STEP 6 Mark duplicates
 */
process markDuplicates {
    tag "${bam.baseName - '.sorted'}"
    publishDir "${params.outdir}/markDuplicates", mode: 'copy',
        saveAs: {filename -> filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"}

    when:
    !params.skip_qc && !params.skip_dupradar

    input:
    file bam

    output:
    file "${bam.baseName}.markDups.bam"
    file "${bam.baseName}.markDups_metrics.txt"
    file "${bam.baseName}.markDups.bam.bai"

    script:
    if( !task.memory ){
        log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
        avail_mem = 3
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    picard -Xmx${avail_mem}g MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${bam.baseName}.markDups.bam \\
        METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT
    samtools index ${bam.baseName}.markDups.bam
    """
}


/*
 * STEP 7 - dupRadar
 */
process dupradar {

    tag "${bam_md.baseName - '.sorted.markDups'}"
    publishDir "${params.outdir}/dupradar", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
            else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
            else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
            else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
            else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
            else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
            else "$filename"
        }

    when:
    !params.skip_qc && !params.skip_dupradar

    input:
    file bam_md
    file gtf

    output:
    file "*.{pdf,txt}" 

    script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
    def dupradar_direction = 0
    if (forward_stranded && !unstranded) {
        dupradar_direction = 1
    } else if (reverse_stranded && !unstranded){
        dupradar_direction = 2
    }
    def paired = params.singleEnd ? 'single' :  'paired'
    """
    dupRadar.r $bam_md $gtf $dupradar_direction $paired ${task.cpus}
    """
}


/*
 * STEP 8 Feature counts
 */
process featureCounts {
    tag "${bam_featurecounts.baseName - '.sorted'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("biotype_counts") > 0) "biotype_counts/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else "$filename"
        }

    input:
    file bam_featurecounts
    file gtf
    file biotypes_header

    output:
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt" 
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" 
    file "${bam_featurecounts.baseName}_biotype_counts*mqc.{txt,tsv}"

    script:
    def featureCounts_direction = 0
    def extraAttributes = params.fcExtraAttributes ? "--extraAttributes ${params.fcExtraAttributes}" : ''
    if (forward_stranded && !unstranded) {
        featureCounts_direction = 1
    } else if (reverse_stranded && !unstranded){
        featureCounts_direction = 2
    }
    // Try to get real sample name
    sample_name = bam_featurecounts.baseName - 'Aligned.sortedByCoord.out'
    """
    featureCounts -a $gtf -g gene_id -o ${bam_featurecounts.baseName}_gene.featureCounts.txt $extraAttributes -p -s $featureCounts_direction $bam_featurecounts
    featureCounts -a $gtf -g gene_biotype -o ${bam_featurecounts.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
    cut -f 1,7 ${bam_featurecounts.baseName}_biotype.featureCounts.txt | tail -n +3 | cat $biotypes_header - >> ${bam_featurecounts.baseName}_biotype_counts_mqc.txt
    mqc_features_stat.py ${bam_featurecounts.baseName}_biotype_counts_mqc.txt -s $sample_name -f rRNA -o ${bam_featurecounts.baseName}_biotype_counts_gs_mqc.tsv
    """
}

/*
 * STEP 9 - Merge featurecounts
 */
process merge_featureCounts {
    tag "${input_files[0].baseName - '.sorted'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy'

    input:
    file input_files

    output:
    file 'merged_gene_counts.txt'

    script:
    //if we only have 1 file, just use cat and pipe output to csvtk. Else join all files first, and then remove unwanted column names.
    def single = input_files instanceof Path ? 1 : input_files.size()
    def merge = (single == 1) ? 'cat' : 'csvtk join -t -f "Geneid,Start,Length,End,Chr,Strand,gene_name"'
    """
    $merge $input_files | csvtk cut -t -f "-Start,-Chr,-End,-Length,-Strand" | sed 's/Aligned.sortedByCoord.out.markDups.bam//g' > merged_gene_counts.txt
    """
}


/*
 * STEP 10 - stringtie FPKM
 */
process stringtieFPKM {
    tag "${bam_stringtieFPKM.baseName - '.sorted'}"
    publishDir "${params.outdir}/stringtieFPKM", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
            else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
            else if (filename.indexOf("ballgown") > 0) "ballgown/$filename"
            else "$filename"
        }

    input:
    file bam_stringtieFPKM
    file gtf

    output:
    file "${bam_stringtieFPKM.baseName}_transcripts.gtf"
    file "${bam_stringtieFPKM.baseName}.gene_abund.txt"
    file "${bam_stringtieFPKM}.cov_refs.gtf"
    file ".command.log" into stringtie_log
    file "${bam_stringtieFPKM.baseName}_ballgown"

    script:
    def st_direction = ''
    if (forward_stranded && !unstranded){
        st_direction = "--fr"
    } else if (reverse_stranded && !unstranded){
        st_direction = "--rf"
    }
    """
    stringtie $bam_stringtieFPKM \\
        $st_direction \\
        -o ${bam_stringtieFPKM.baseName}_transcripts.gtf \\
        -v \\
        -G $gtf \\
        -A ${bam_stringtieFPKM.baseName}.gene_abund.txt \\
        -C ${bam_stringtieFPKM}.cov_refs.gtf \\
        -e \\
        -b ${bam_stringtieFPKM.baseName}_ballgown
    """
}

/*
 * STEP 11 - edgeR MDS and heatmap
 */
process sample_correlation {
    tag "${input_files[0].toString() - '.sorted_gene.featureCounts.txt' - 'Aligned'}"
    publishDir "${params.outdir}/sample_correlation", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_edger

    input:
    file input_files
    val num_bams
    file mdsplot_header
    file heatmap_header 

    output:
    file "*.{txt,pdf,csv}"

    when:
    num_bams > 2 && (!params.sampleLevel)

    script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
    """
    edgeR_heatmap_MDS.r $input_files
    cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv >> tmp_file
    mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
    cat $heatmap_header log2CPM_sample_distances_mqc.csv >> tmp_file
    mv tmp_file log2CPM_sample_distances_mqc.csv
    """
}

/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' 

    script:
    """
    echo $workflow.manifest.version &> v_ngi_rnaseq.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    STAR --version &> v_star.txt
    hisat2 --version &> v_hisat2.txt
    stringtie --version &> v_stringtie.txt
    preseq &> v_preseq.txt
    read_duplication.py --version &> v_rseqc.txt
    echo \$(bamCoverage --version 2>&1) > v_deeptools.txt
    featureCounts -v &> v_featurecounts.txt
    picard MarkDuplicates --version &> v_markduplicates.txt  || true
    samtools --version &> v_samtools.txt
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * Pipeline parameters to go into MultiQC report
 */
process workflow_summary_mqc {

    when:
    !params.skip_multiqc

    output:
    file 'workflow_summary_mqc.yaml' 

    exec:
    def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nfcore-rnaseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nfcore/rnaseq Workflow Summary'
    section_href: 'https://github.com/nf-core/rnaseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()
}

/*
 * STEP 12 MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file multiqc_config 
    file (fastqc:'fastqc/*') 
    file ('trimgalore/*')
    file ('alignment/*')
    file ('rseqc/*')
    file ('rseqc/*')
    file ('preseq/*')
    file ('dupradar/*')
    file ('featureCounts/*')
    file ('featureCounts_biotype/*')
    file ('stringtie/stringtie_log*')
    file ('sample_correlation_results/*') 
    file ('software_versions/*')
    file ('workflow_summary/*')

    output:
    file "*multiqc_report.html" 
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m picard -m preseq -m rseqc -m featureCounts -m hisat2 -m star -m cutadapt -m fastqc
    """
}

/*
 * STEP 13 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}
