class UTILS {
    // Remove Nextflow version from pipeline_software_mqc_versions.yml
    public static Object removeNextflowVersion(pipeline_software_mqc_versions) {
        def softwareVersions = path(pipeline_software_mqc_versions).yaml
        if (softwareVersions.containsKey("Workflow")) softwareVersions.Workflow.remove("Nextflow")
        return softwareVersions
    }

    // Recursively list all files in a directory and its sub-directories, matching or not matching supplied suffixes
    public static getAllFilesFromDir(outdir, boolean includeDir = true, List<String> excludeRegexes = null) {
        def output = []
        new File(outdir).eachFileRecurse() { file ->
            boolean matchesInclusion = includeDir     ? true : file.isFile()
            boolean matchesExclusion = excludeRegexes ? excludeRegexes.any { regex -> file.name.toString() ==~ regex } : false

            // Add files (or folders if includeDir is set to true) to the list that don't match excludeRegexes
            if (matchesInclusion && !matchesExclusion) output.add(file)
        }
        return output.sort { it.path }
    }
    // Static (global) exclusion regexes list
    static List<String> exclusionRegexesForUnstableFileNames = [/.*\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}.*/]
}
