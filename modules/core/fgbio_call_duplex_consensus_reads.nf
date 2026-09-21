/*
 * Step 8 of the Twist Bioscience UMI Protocol
 *
 * Calls duplex consensus sequences, accounting for information from both strands, 
 * from reads generated from the same double-stranded source molecule. 
 */
process FGBIO_CALL_DUPLEX_CONSENSUS_READS {

    label 'twist'
    label 'process_high'

    input:
    tuple val(meta), val(ref_info), path(ref), path(grouped_bam)

    output:
    tuple val(meta), val(ref_info), path(ref), path("*.bam"),   emit: unaligned_consensus_bam

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2
    
    """
    fgbio CallDuplexConsensusReads \\
        -Xmx${avail_mem}g \\
        --input=${grouped_bam} \\
        --output="${ref.simpleName}_unaligned_consensus.bam" \\
        --error-rate-pre-umi=45 \\
        --error-rate-post-umi=30 \\
        --min-input-base-quality=30 \\
        --min-reads ${params.fgbio_consensus_min_reads}
    """
}
