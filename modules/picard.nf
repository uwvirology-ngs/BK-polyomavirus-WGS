process PICARD {

    container 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'

    input:
    tuple val(sample_name), path(read1), path(read2)

    output:
    tuple val(sample_name), path("${sample_name}_unaligned.bam"),   emit: unaligned_bam

    script:
    """
    picard FastqToSam \\
        O="${sample_name}_unaligned.bam" \\
        F1=${read1} \\
        F2=${read2} \\
        SM=${sample_name} \\
        LB=Library1 \\
        PU=Unit1 \\
        PL=Illumina
    """
}