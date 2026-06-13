/*
 * Step 8 of the Twist Bioscience UMI Protocol
 *
 * Calls duplex consensus sequences, accounting for information from both strands, 
 * from reads generated from the same double-stranded source molecule. 
 */
process FGBIO_CALL_DUPLEX_CONSENSUS_READS {

    label 'twist'

    input:
    tuple val(meta), path(grouped_bam), path(ref)

    output:
    tuple val(meta), path("*.bam"), path(ref),  emit: unaligned_consensus_bam

    script:
    """
    fgbio CallDuplexConsensusReads \\
        --input=${grouped_bam} \\
        --output="${ref.simpleName}_unaligned_consensus.bam" \\
        --error-rate-pre-umi=45 \\
        --error-rate-post-umi=30 \\
        --min-input-base-quality=30 \\
        --min-reads 2 1 1
    """
}