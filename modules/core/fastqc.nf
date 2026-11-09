/*
 * Runs FASTQC
 */
process FASTQC {

    label 'process_low'
    container 'community.wave.seqera.io/library/fastqc:0.12.1--9971ea336a9eddae'

    input:
    tuple val(meta), path(reads)

    output:
    path("*"),  emit: fastqc

    script:
    """
    fastqc ${reads}
    """
}
