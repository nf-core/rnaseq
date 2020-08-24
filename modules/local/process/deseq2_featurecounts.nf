conda (params.conda ? "${baseDir}/environment.yml" : null)

// /*
//  * STEP 7.4: Differential analysis with DESeq2
//  */
// process CONSENSUS_PEAKS_DESEQ2 {
//     tag "${antibody}"
//     label 'process_medium'
//     publishDir "${params.outdir}/bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}/deseq2", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith('.igv.txt')) null
//                       else filename
//                 }
//
//     when:
//     params.macs_gsize && replicatesExist && multipleGroups && !params.skip_consensus_peaks && !params.skip_diff_analysis
//
//     input:
//     tuple val(antibody), val(replicatesExist), val(multipleGroups), path(counts) from ch_macs_consensus_counts
//     path deseq2_pca_header from ch_deseq2_pca_header
//     path deseq2_clustering_header from ch_deseq2_clustering_header
//
//     output:
//     path '*.tsv' into ch_macs_consensus_deseq_mqc
//     path '*igv.txt' into ch_macs_consensus_deseq_comp_igv
//     path '*.{RData,results.txt,pdf,log}'
//     path 'sizeFactors'
//     path '*vs*/*.{pdf,txt}'
//     path '*vs*/*.bed'
//
//     script:
//     prefix = "${antibody}.consensus_peaks"
//     bam_ext = params.single_end ? '.mLb.clN.sorted.bam' : '.mLb.clN.bam'
//     vst = params.deseq2_vst ? '--vst TRUE' : ''
//     """
//     featurecounts_deseq2.r \\
//         --featurecount_file $counts \\
//         --bam_suffix '$bam_ext' \\
//         --outdir ./ \\
//         --outprefix $prefix \\
//         --outsuffix '' \\
//         --cores $task.cpus \\
//         $vst
//
//     sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
//     sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
//     cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv
//
//     sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
//     sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
//     cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv
//
//     find * -type f -name "*.FDR0.05.results.bed" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}/deseq2/"{}"\\t255,0,0" \\; > ${prefix}.igv.txt
//     """
// }
