process PANDEPTH {

    label 'process_low'
    container 'quay.io/epil02/revica-strm:0.0.5'

    input: 
    tuple val(meta), path(bam), path(db)

    output:
    tuple val(meta), path("*.stat"),    emit: pandepth
    tuple val(meta), path("*covstats.tsv"), emit: covstats

    script:
    """
    ## run pandepth to get coverage/depth reporting
    pandepth -i ${bam} -o ${meta.id} -t ${task.cpus}
    gunzip -f ${meta.id}.chr.stat.gz

    ## replace abbreviated ref names in pandepth with originals from db
    ## we do this because BWA_MEM only records the alignment ref before the first space, which is
    ## usually just the acc number. We need the rest of the fasta header for downstream analyses

    prep_pandepth_output.py ${meta.id}.chr.stat $db ${meta.id}_covstats.tsv
    """
}