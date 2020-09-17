/*
 * This file holds several functions used to perform JSON parameter validation, help and summary rendering for the nf-core pipeline template.
 */

import groovy.json.JsonSlurper

class Schema {
    /*
     * This method tries to read a JSON params file
     */
    private static LinkedHashMap params_get(String path) {
        def params_map = new LinkedHashMap()
        try {
            params_map = params_try(path)
        } catch (Exception e) {
            println "Could not read parameters settings from JSON. $e"
            params_map = new LinkedHashMap()
        }
        return params_map
    }

    /*
    Method to actually read in JSON file using Groovy.
    Group (as Key), values are all parameters
        - Parameter1 as Key, Description as Value
        - Parameter2 as Key, Description as Value
        ....
    Group
        -
    */
    private static LinkedHashMap params_try(String path) throws Exception {
        def json = new File(path).text
        def Map json_params = (Map) new JsonSlurper().parseText(json).get('definitions')
        /* Tree looks like this in nf-core schema
         * definitions <- this is what the first get('definitions') gets us
             group 1
               title
               description
                 properties
                   parameter 1
                     type
                     description
                   parameter 2
                     type
                     description
             group 2
               title
               description
                 properties
                   parameter 1
                     type
                     description
        */
        def params_map = new LinkedHashMap()
        json_params.each { key, val ->
            def Map group = json_params."$key".properties // Gets the property object of the group
            def title = json_params."$key".title
            def sub_params = new LinkedHashMap()
            group.each { innerkey, value ->
                sub_params.put("$innerkey", [ "$value.type", "$value.description" ])
            }
            params_map.put("$title", sub_params)
        }
        return params_map
    }

    private static Integer params_max_chars(params_map) {
        Integer max_chars = 0
        for (group in params_map.keySet()) {
            def params = params_map.get(group)  // This gets the parameters of that particular group
            for (par in params.keySet()) {
                if (par.size() > max_chars) {
                    max_chars = par.size()
                }
            }
        }
        return max_chars
    }

    private static String params_beautify(params_map) {
        String output = ""
        def max_chars = params_max_chars(params_map) + 1
        for (group in params_map.keySet()) {
            output += group + "\n"
            def params = params_map.get(group)  // This gets the parameters of that particular group
            for (par in params.keySet()) {
                def type = "[" + params.get(par)[0] + "]"
                def description = params.get(par)[1]
                output+= "    \u001B[1m" +  par.padRight(max_chars) + "\u001B[1m" + type.padRight(10) + description + "\n"
            }
            output += "\n"
        }
        return output
    }

    private static String params_help(path, command) {
          String output = "Typical pipeline command:\n\n"
          output += "    ${command}\n\n"
          output += params_beautify(params_get(path))
    }

    private static LinkedHashMap params_summary(workflow, params, run_name) {
        def Map summary = [:]
        if (workflow.revision)             summary['Pipeline Release'] = workflow.revision
        summary['Run Name']                = run_name ?: workflow.runName
        summary['Design File']             = params.input
        summary['Genome']                  = params.genome ?: 'Not supplied'
        summary['Fasta File']              = params.fasta
        if (params.additional_fasta)       summary["Additional Fasta File"] = params.additional_fasta
        if (params.gtf)                    summary['GTF Annotation'] = params.gtf
        if (params.gff)                    summary['GFF3 Annotation'] = params.gff
        if (params.gene_bed)               summary['Gene BED Annotation'] = params.gene_bed
        if (params.gencode)                summary['GENCODE'] = params.gencode
        if (params.stringtie_ignore_gtf)   summary['StringTie Ignore GTF'] = params.stringtie_ignore_gtf
        if (params.fc_group_features_type) summary['Biotype GTF field'] = params.gencode ? "gene_type" : params.fc_group_features_type
        summary['Save Reference Files']    = params.save_reference ? 'Yes' : 'No'
        summary['Remove Ribosomal RNA']    = params.remove_ribo_rna
        if (params.skip_trimming) {
            summary['Trimming Step']       = 'Skipped'
        } else {
            summary['Trimming']            = "5'R1: $params.clip_r1 bp / 5'R2: $params.clip_r2 bp / 3'R1: $params.three_prime_clip_r1 bp / 3'R2: $params.three_prime_clip_r2 bp / NextSeq Trim: $params.trim_nextseq bp"
        }
        if (params.with_umi) {
            summary["With UMI"]                           = params.with_umi
            summary["umi_tools extract-method"]           = params.umitools_extract_method
            summary["umi_tools bc-pattern"]               = params.umitools_bc_pattern
            summary["umi_tools extract extra parameters"] = params.umitools_extract_extra
            summary["umi_tools dedup extra parameters"]   = params.umitools_dedup_extra
        }
        if (params.save_trimmed)           summary['Save Trimmed'] = 'Yes'
        if (params.aligner == 'star') {
            summary['Aligner'] = "STAR"
            if (params.star_index)         summary['STAR Index'] = params.star_index
        } else if (params.aligner == 'hisat2') {
            summary['Aligner'] = "HISAT2"
            if (params.hisat2_index)         summary['HISAT2 Index'] = params.hisat2_index
            if (params.splicesites)          summary['Splice Sites'] = params.splicesites
        }
        if (params.pseudo_aligner == 'salmon') {
            summary['Pseudo Aligner']    = "Salmon"
            if (params.transcript_fasta) summary['Transcript Fasta'] = params.transcript_fasta
        }
        if (params.seq_center)           summary['Sequencing Center'] = params.seq_center
        if (params.save_align_intermeds) summary['Save Align Intermeds'] =  'Yes'
        if (params.skip_fastqc)          summary['Skip FastQC'] = 'Yes'
        if (params.skip_preseq)          summary['Skip Preseq'] = 'Yes'
        if (params.skip_multiqc)         summary['Skip MultiQC'] = 'Yes'
        summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
        if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
        summary['Output dir']       = params.outdir
        summary['Launch dir']       = workflow.launchDir
        summary['Working dir']      = workflow.workDir
        summary['Script dir']       = workflow.projectDir
        summary['User']             = workflow.userName
        if (workflow.profile.contains('awsbatch')) {
            summary['AWS Region']   = params.awsregion
            summary['AWS Queue']    = params.awsqueue
            summary['AWS CLI']      = params.awscli
        }
        summary['Config Profile'] = workflow.profile
        if (params.config_profile_description) summary['Config Profile Descr']   = params.config_profile_description
        if (params.config_profile_contact)     summary['Config Profile Contact'] = params.config_profile_contact
        if (params.config_profile_url)         summary['Config Profile URL']     = params.config_profile_url
        summary['Config Files'] = workflow.configFiles.join(', ')
        if (params.email || params.email_on_fail) {
            summary['E-mail Address']    = params.email
            summary['E-mail on failure'] = params.email_on_fail
            summary['MultiQC maxsize']   = params.max_multiqc_email_size
        }
        return summary
    }

    static String params_summary_multiqc(summary) {
        String yaml_file_text  = """
        id: 'nf-core-rnaseq-summary'
        description: " - this information is collected when the pipeline is started."
        section_name: 'nf-core/rnaseq Workflow Summary'
        section_href: 'https://github.com/nf-core/rnaseq'
        plot_type: 'html'
        data: |
            <dl class=\"dl-horizontal\">
            ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
            </dl>
        """.stripIndent()
        return yaml_file_text
    }
}
