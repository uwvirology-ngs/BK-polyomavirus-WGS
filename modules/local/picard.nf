process PICARD {

    container 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_unaligned.bam"),  emit: unaligned_bam

    script:
    """
    picard FastqToSam \\
        O="${meta.id}_unaligned.bam" \\
        F1=${meta.id}_R1.fastq.gz \\
        F2=${meta.id}_R2.fastq.gz \\
        SM=${meta.id} \\
        LB=Library1 \\
        PU=Unit1 \\
        PL=Illumina
    """
}