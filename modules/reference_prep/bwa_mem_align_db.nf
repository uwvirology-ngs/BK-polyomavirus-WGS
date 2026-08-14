process BWA_MEM_ALIGN_DB {

    label 'process_high'
    container 'quay.io/epil02/revica-strm:0.0.5'

    input:
    tuple val(meta), path(fastqs)
    path(db)

    output:
    tuple val(meta), path("*.bam"), path(db), emit: alignment_to_db

    script:
    """
    # this flag combines 0x4 and 0x800, which means:
    # 1. read is unmapped
    # 2. read is supplementary alignment (chimeric, not representative alignment)
    # we can't let these reads survive.
    FLAG=2052

    bwa index ${db}

    bwa mem \
        ${db} \
        ${fastqs} \
        -t $task.cpus \
        | samtools view -bS -F \$FLAG -@ 2 > ${meta.id}.bam
    """
}
