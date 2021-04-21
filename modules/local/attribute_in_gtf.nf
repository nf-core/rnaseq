/*
 * Check if attribute is available in column 9 of GTF file
 */
process ATTRIBUTE_IN_GTF {

    input:
    val gtf
    val attribute

    output:
    val attribute_in_gtf
    
    exec:
    def hits = 0
    def gtf_file = file(gtf)
    gtf_file.eachLine { line ->
        def attributes = line.split('\t')[-1].split()
        if (attributes.contains(attribute)) {
            hits ++
        }
    }

    attribute_in_gtf = false
    if (hits) {
        attribute_in_gtf = true
    }
}
