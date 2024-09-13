class UTILS {
    // Remove Nextflow version from pipeline_software_mqc_versions.yml
    public static Object removeNextflowVersion(pipeline_software_mqc_versions) {
        def softwareVersions = path(pipeline_software_mqc_versions).yaml
        if (softwareVersions.containsKey("Workflow")) softwareVersions.Workflow.remove("Nextflow")
        return softwareVersions
    }

    // Recursively list all files in a directory and its sub-directories, matching or not matching supplied suffixes
    public static getAllFilesFromDir(dir, List<String> includeRegexes = null, List<String> excludeRegexes = null) {
        def output = []
        new File(dir).eachFileRecurse() { file ->
            boolean matchesInclusion = (includeRegexes == null || includeRegexes.any  { regex -> file.name.toString() ==~ regex })
            boolean matchesExclusion = (excludeRegexes == null || !excludeRegexes.any { regex -> file.name.toString() ==~ regex })

            // Conditionally add either full path or just the file name
            if (matchesInclusion && matchesExclusion) {
                output.add(file)
            }
        }
        return output.sort { it.name }
    }

    // Static (global) exclusion regexes list
    static List<String> exclusionRegexesForUnstableFileNames = [/.*\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}.*/]
    static List<String> snapshottablePatterns = [/.*\.(txt|json|tsv)$/]
    static List<String> exclusionRegexesForUnstableFileContents = [
        // To exclude files with timestamps in the format YYYY-MM-DD_HH-MM-SS
        /\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}/,
        // To exclude "goleft_indexcov-bin-plot.txt"
        /goleft_indexcov-bin-plot\.txt/,
        // To exclude "multiqc_data.json"
        /multiqc_data\.json/,
        // To exclude "multiqc_sources.txt"
        /multiqc_sources\.txt/,
        // To exclude params files with timestamps
        /params_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}\.json/,
        // To exclude "*.deduped.cram.metrics.multiqc.tsv"
        /.*\.deduped\.cram\.metrics\.multiqc\.tsv/,
        // To exclude "goleft_bin.txt"
        /goleft_bin\.txt/,
        // To exclude "goleft_roc.txt"
        /goleft_roc\.txt/,
        // To exclude any "multiqc_sample*.txt" files
        /multiqc_sample.*\.txt/,
        // To exclude "multiqc_software_versions.txt"
        /multiqc_software_versions\.txt/
    ]

}
