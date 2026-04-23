/*
 * Parameters
 */
params {
    samplesheet: Path
    _adapter_fasta: Path
}

/*
 * Read inparams. the samplesheet and emit resulting metamaps to main workflow
 */
workflow READ_SAMPLESHEET {

    main:
    channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> build_metamap(row) }
        .map { meta, fastqs -> [ meta, fastqs, params._adapter_fasta ] }
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

    // require twin fastq files for paired-end reads
    if (!file(row.fastq_1).exists() || !file(row.fastq_2).exists()) {
        exit 1, "ERROR: samplesheet requires two fastqs for paired-end reads."
    }

    return [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
}