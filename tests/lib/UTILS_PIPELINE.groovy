class UTILS_PIPELINE {
    // These are the files to exclude when we want to snapshot
    static List<String> exclusionRegexesForUnstableFileContents = [
        // To exclude the pipeline_software_mqc_versions.yml file that contains the Nextflow version
        /nf_core_.*_software_mqc_versions\.yml/,

        // To exclude this folder that somehow is a file on stub tests
        /multiqc_plots/,

        // To exclude from qualimap css folders
        /.*\.(css|gif|js)/,

        // To exclude bbsplit files
        /.*\.stats\.txt/,

        // To exclude FASTQC reports
        /.*_raw\.(html|zip)/,
        /.*_fastqc\.(html|zip)/,

        // To exclude from the MultiQC reports
        /cutadapt_filtered_reads_plot-(cnt|pct)\.(pdf|svg)/,
        /cutadapt_trimmed_sequences_.*\.(pdf|svg)/,
        /dupradar-section-plot\.(pdf|svg)/,
        /fail_mapped_samples_table.*/,
        /fail_strand_check_table.*/,
        /fastqc-status-check-.*\.(pdf|svg)/,
        /fastqc_adapter_content_plot\.(pdf|png|svg)/,
        /fastqc_overrepresented_sequences_plot(.*)?\.(pdf|svg)/,
        /fastqc_per_base_.*_plot(.*)?\.(pdf|png|svg)/,
        /fastqc_per_sequence_.*\.(pdf|svg)/,
        /fastqc_sequence_(counts|duplication_levels)_plot(.*)?-(cnt|pct)\.(pdf|svg)/,
        /fastqc_sequence_length_distribution_plot\.(pdf|png|svg)/,
        /fastqc_top_overrepresented_sequences_table(.*)?\.(pdf|png|svg|txt)/,
        /featurecounts_biotype_plot-(cnt|pct)\.(pdf|png|svg)/,
        /general_stats_table\.(pdf|png|svg)/,
        /hisat2_(pe|se)_plot-(cnt|pct)\.(pdf|png|svg)/,
        /hisat2_pe_plot\.txt/,
        /hisat2_se_plot.*\.(png|svg)/,
        /junction_saturation_known\.txt/,
        /junction_saturation_novel\.txt/,
        /kallisto_alignment.*.(pdf|png|svg)/,
        /kallisto_alignment\.txt/,
        /multiqc_data\.json/,
        /multiqc_dupradar-section-plot\.txt/,
        /multiqc_fail_strand_check_table\.txt/,
        /multiqc_general_stats\.txt/,
        /multiqc_hisat2\.txt/,
        /multiqc_kallisto\.txt/,
        /multiqc_picard_dups\.txt/,
        /multiqc_report\.html/,
        /multiqc_rsem\.txt/,
        /multiqc_rseqc_bam_stat\.txt/,
        /multiqc_rseqc_infer_experiment\.txt/,
        /multiqc_rseqc_junction_annotation\.txt/,
        /multiqc_rseqc_read_distribution\.txt/,
        /multiqc_salmon\.txt/,
        /(multiqc_)?(salmon|star_salmon|star_rsem)_deseq2_clustering-plot(.*)?\.(pdf|png|txt)/,
        /multiqc_(salmon|star_salmon)_deseq2_pca-plot.*\.txt/,
        /multiqc_samtools_flagstat\.txt/,
        /multiqc_samtools_stats\.txt/,
        /multiqc_software_versions\.txt/,
        /multiqc_sortmerna\.txt/,
        /multiqc_sources\.txt/,
        /multiqc_star\.txt/,
        /picard_deduplication-(cnt|pct)\.(pdf|png|svg)/,
        /picard_deduplication\.txt/,
        /qualimap_gene_coverage_profile_Counts\.(pdf|svg)/,
        /qualimap_gene_coverage_profile_Normalised\.(pdf|svg)/,
        /qualimap_genomic_origin-(cnt|pct)\.(pdf|png|svg)/,
        /qualimap_genomic_origin\.txt/,
        /qualimap_rnaseq_genome_results\.txt/,
        /rsem_assignment_plot\.txt/,
        /rsem_assignment_plot-(cnt|pct)\.(pdf|png|svg)/,
        /rsem_multimapping_rates\.(pdf|svg|txt)/,
        /rseqc_bam_stat\.(pdf|png|svg|txt)/,
        /rseqc_infer_experiment_plot\.(pdf|svg)/,
        /rseqc_inner_.*/,
        /rseqc_junction_.*/,
        /rseqc_read_distribution_.*/,
        /rseqc_read_dups\.txt/,
        /rseqc_read_dups_plot\.(pdf|svg)/,
        /rseqc_read_dups_plot\.txt/,
        /(salmon|star_rsem)_deseq2_(clustering|pca)-plot\.(pdf|png|svg)/,
        /(salmon|star_rsem)_deseq2_(pca|pca)-plot\.(pdf|png|svg)/,
        /salmon_plot\.(pdf|png|svg|txt)/,
        /samtools-flagstat.*/,
        /samtools-idxstats-mapped-reads-plot_Normalised_Counts-(cnt|pct)\.(pdf|svg)/,
        /samtools-idxstats-mapped-reads-plot_Normalised_Counts-log\.(pdf|svg)/,
        /samtools-idxstats-mapped-reads-plot_Observed_over_Expected_Counts-(cnt|pct)\.(pdf|svg)/,
        /samtools-idxstats-mapped-reads-plot_Observed_over_Expected_Counts-log\.(pdf|svg)/,
        /samtools-idxstats-mapped-reads-plot_Raw_Counts-(cnt|pct)\.(pdf|svg)/,
        /samtools-idxstats-mapped-reads-plot_Raw_Counts-log\.(pdf|svg)/,
        /samtools-stats-dp\.(pdf|png|svg|txt)/,
        /(samtools|star)_alignment_plot.*/,
        /sortmerna-detailed-plot-(cnt|pct)\.(pdf|svg)/,
        /sortmerna-detailed-plot\.txt/,
        /(star_rsem|star_salmon)_deseq2_(clustering|pca)-plot(.*)?\.(pdf|png|svg)/,
        /star_summary_table\.(pdf|png|svg|txt)/,

        // To exclude from deseq2_qc
        /RAP1_IAA_30M_REP1\.txt/,
        /RAP1_UNINDUCED_REP1\.txt/,
        /RAP1_UNINDUCED_REP2\.txt/,
        /WT_REP1\.txt/,
        /WT_REP2\.txt/,
        /deseq2.*/,

        // To exclude from kallisto
        /abundance\.(h5|tsv)/,
        /kallisto_quant\.log/,
        /run_info\.json/,

        // To exclude from salmon quant
        /fld\.gz/,
        /meta_info\.json/,
        /flenDist\.txt/,
        /salmon_quant\.log/,
        /quant\.genes\.sf/,
        /quant\.sf/,

        // To exclude from kallisto|salmon aligners
        /(kallisto|salmon)\.merged.*/,

        // To exclude from star_rsem
        /.*.(genes|isoforms)\.results/,
        /.*.(cnt|model|theta)/,

        // To exclude bigwig
        /.*\.(forward|reverse)\.bigWig/,

        // To exclude dupradar
        /.*_(duprateExpBoxplot|duprateExpDens|expressionHist)\.(pdf|png|svg)/,

        // To exclude featurecounts
        /.*\.featureCounts\.txt\.summary/,

        // To exclude star salmon
        /.*\.Log\.?(final|progress)?\.out/,

        // To exclude Picard Markduplicates metrics
        /.*\.markdup\.sorted\.MarkDuplicates\.metrics\.txt/,

        // To exclude Qualimap files
        /Junction\sAnalysis\.png/,
        /Reads\sGenomic\sOrigin\.png/,
        /qualimapReport\.html/,
        /rnaseq_qc_results\.txt/,

        // To exclude rseqc
        /.*\.DupRate_plot\.(pdf|r)/,
        /.*\.inner_distance.*/,
        /.*\.junction.*/,
        /.*\.(pos|seq)\.DupRate\.xls/,
        /.*\.read_distribution\.txt/,
        /.*\.splice_(events|junction)\.(pdf|png|svg)/,
        /.*\.bam_stat\.txt/,

        // To exclude from samtools stats
        /.*\.sorted\.bam\.(flagstat|idxstats|stats)/,

        // To exclude from hisat2
        /.*hisat2\.summary/,

        // To exclude from stringtie
        /t_data\.ctab/,
        /.*\.coverage\.gtf/,
        /.*\.gene\.abundance\.txt/,
        /.*\.transcripts\.gtf/,

        // To exclude from sortmerna
        /.*\.sortmerna\.log/,

        // To exclude log from star rsem
        /RAP1_IAA_30M_REP1\.log/,
        /RAP1_UNINDUCED_REP1\.log/,
        /RAP1_UNINDUCED_REP2\.log/,
        /WT_REP1\.log/,
        /WT_REP2\.log/,

        // To exclude markdup
        /.*\.markdup\.sorted\.bam(\.bai)?/,

        // To exclude trimgalore
        /.*\.fastq\.gz_trimming_report\.txt/
    ]
}
