// Helper functions for pipeline tests

class UTILS {

    // Function to remove Nextflow version from software_versions.yml
    public static String removeNextflowVersion(outputDir) {
        def softwareVersions = path("$outputDir/pipeline_info/software_versions.yml").yaml
        if (softwareVersions.containsKey("Workflow")) {
            softwareVersions.Workflow.remove("Nextflow")
        }
        return softwareVersions
    }

    // Function to filter lines from a file and return a new file
    public static File filterLines(String inFilePath, int linesToSkip) {
        if (linesToSkip >= 0) {
            File inputFile = new File(inFilePath)
            File outputFile = new File(inFilePath + ".filtered")
            def lineCount = 0
            inputFile.eachLine { line ->
                lineCount++
                if (lineCount > linesToSkip) {
                    outputFile.append(line + '\n')
                }
            }
            return outputFile
        } else {
            File inputFile = new File(inFilePath)
            File outputFile = new File(inFilePath + ".filtered")
            def lines = inputFile.readLines()
            def totalLines = lines.size()
            lines.take(totalLines + linesToSkip).each { line ->
                outputFile.append(line + '\n')
            }
            return outputFile
        }
    }
}
