process FGBIO {

    container 'community.wave.seqera.io/library/fgbio:1.5.1--9658b275f6ef481d'

    input:
    tuple val(sample_name), path(unaligned_bam)

    output:
    tuple val(sample_name), path("${sample_name}_umi_extracted.bam"),   emit: umi_extracted_bam

    script:
    """
    fgbio ExtractUmisFromBam \\
        --input="${unaligned_bam}" \\
        --output="${sample_name}_umi_extracted.bam" \\
        --read-structure=5M2S+T 5M2S+T \\
        --molecular-index-tags=ZA ZB \\
        --single-tag=RX
    """
}