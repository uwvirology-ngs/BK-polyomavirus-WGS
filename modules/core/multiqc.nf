/*
 * Runs MULTIQC on collected FASTQC outputs
 */
process MULTIQC {

    label 'process_low'
    container 'community.wave.seqera.io/library/multiqc:1.35--c5313f6c74976ae9'

    input:
    path(fastqc)

    output:
    path("*"),  emit: multiqc

    script:
    """
    multiqc .
    """
}
