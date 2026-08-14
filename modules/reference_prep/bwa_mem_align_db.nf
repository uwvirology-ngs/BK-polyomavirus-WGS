process BWA_MEM_ALIGN_DB {

    label 'process_high'
    container 'quay.io/epil02/revica-strm:0.0.5'

    input:
    tuple val(meta), path(fastqs)
    path(db)

    output:
    tuple val(meta), path("*covstats.tsv"), emit: covstats

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

    ## run pandepth to get coverage/depth reporting
    pandepth -i ${meta.id}.bam -o ${meta.id} -t ${task.cpus}
    gunzip -f ${meta.id}.chr.stat.gz

    ## replace abbreviated ref names in pandepth with originals from db
    ## we do this because BWA_MEM only records the alignment ref before the first space, which is
    ## usually just the acc number. We need the rest of the fasta header for downstream analyses

    prep_pandepth_output.py ${meta.id}.chr.stat $db ${meta.id}_covstats.tsv
    """
}
