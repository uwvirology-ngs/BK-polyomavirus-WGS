/*
 * Read in the samplesheet and emit resulting metamaps to main workflow
 */
workflow READ_SAMPLESHEET {

    take:
    samplesheet

    main:
    channel.fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row -> build_metamap(row) }
        .set { reads }

    emit:
    reads
}

/*
 * Transforms each samplesheet row into a meta map
 */
def build_metamap(LinkedHashMap row) {
    
    // construct metadata map
    def meta = [:]
    meta.id = row.sample
    meta.single_end = false

    def reads_1 = file("${projectDir}/${row.fastq_1}")
    def reads_2 = file("${projectDir}/${row.fastq_2}")

    if (params.cloud_compute) {
        reads_1 = file(row.fastq_1)
        reads_2 = file(row.fastq_2)
    }

    // require twin fastq files for paired-end reads
    if (!reads_1.exists() || !reads_2.exists()) {
        exit 1, "ERROR: samplesheet requires two fastqs for paired-end reads."
    }

    return [ meta, [ reads_1, reads_2 ] ]
}