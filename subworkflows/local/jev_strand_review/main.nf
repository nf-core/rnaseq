// Shadow judgments never feed strand settings, failure checks, or sample filters.
include { jev } from 'plugin/nf-jev'
include { reviewStrandEvidence; strandReviewMultiqc } from './review'

workflow JEV_STRAND_REVIEW {
    take:
    strand_data
    protocol
    model
    confidence_floor
    report_asset
    output_dir

    main:
    def template = new groovy.json.JsonSlurper().parseText(report_asset.text)
    reviews = strand_data.map { row ->
        reviewStrandEvidence(row, protocol, model, confidence_floor, { state, questions -> jev(state, questions) })
    }

    audit = reviews
        .map { record -> groovy.json.JsonOutput.toJson(record) }
        .collectFile(name: 'strand_review.jsonl', newLine: true, storeDir: output_dir)

    // Hash local IDs for filenames; the readable ID remains inside the report only.
    reports = reviews
        .collectFile { record ->
            def key = java.security.MessageDigest.getInstance('SHA-256').digest(record.sample.toString().bytes).encodeHex().toString()
            ["${key}_jev_strand_review_mqc.json", strandReviewMultiqc(template, [record])]
        }
        .map { report ->
            def document = new groovy.json.JsonSlurper().parseText(report.text)
            [[id: document.data.keySet().first()], report]
        }

    emit:
    multiqc = reports
    decisions = audit
}
