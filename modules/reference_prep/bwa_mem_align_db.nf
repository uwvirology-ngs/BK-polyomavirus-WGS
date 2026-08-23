/*
 * Step 1 of the reference_prep subworkflow
 *
 * Align reads to the multifasta database for downstream evaluation 
 * of coverage and depth stats.
 */
process BWA_MEM_ALIGN_DB {

    label 'process_high'
    container 'quay.io/epil02/revica-strm:0.0.5'

    input:
    tuple val(meta), path(fastqs)
    path(db)

    output:
    tuple val(meta), path("*.bam"), emit: alignment_to_db

    script:
    """
    bwa index ${db}

    # read unmapped (0x4), supplementary alignment (0x800)
    FLAG=2052

    bwa mem ${db} ${fastqs} -t $task.cpus \\
    | samtools view -b -F \$FLAG -@ 2 > "${meta.id}_aligned_to_db.bam"
    """
}
