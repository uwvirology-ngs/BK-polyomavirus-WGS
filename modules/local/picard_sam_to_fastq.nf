process PICARD_SAM_TO_FASTQ {

    container 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'

    input:
    tuple val(meta), path(unaligned_bam_umi_extracted)

    output:
    tuple val(meta), path("${meta.id}_umi_extracted.fastq"),    emit: umi_extracted_fastq

    script:
    """
    picard SamToFastq \\
        I="${unaligned_bam_umi_extracted}" \\
        F="${meta.id}_umi_extracted.fastq" \\
        INTERLEAVE=true
    """
}