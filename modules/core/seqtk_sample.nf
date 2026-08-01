process SEQTK_SAMPLE {
    label 'twist'
    label 'process_high'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"),    emit: sampled_reads

    script:
    """
    fastq_1=\$(echo ${reads} | cut -d' ' -f1)
    fastq_2=\$(echo ${reads} | cut -d' ' -f2)
    
    seqtk sample -s 100 \$fastq_1 ${params.sampleTo} | gzip > ${meta.id}_R1.sample.fastq.gz
    seqtk sample -s 100 \$fastq_2 ${params.sampleTo} | gzip > ${meta.id}_R2.sample.fastq.gz
    """
}