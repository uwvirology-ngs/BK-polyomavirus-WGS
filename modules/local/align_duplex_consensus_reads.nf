/*
 * Step 9 of the Twist Bioscience UMI Protocol
 *
 * After consensus duplex collapse, convert reads back to FASTQ format and realign. 
 */
process ALIGN_DUPLEX_CONSENSUS_READS {

    label 'twist'
    label 'process_high'

    input:
    tuple val(meta), path(unaligned_consensus_bam), path(ref), val(ref_info)

    output:
    tuple val(meta), path("*.bam"), path(ref), val(ref_info),  emit: aligned_consensus_bam

    script:
    """
    picard SamToFastq \\
        I=${unaligned_consensus_bam} \\
        F="${ref.simpleName}_consensus_interleaved.fastq" \\
        INTERLEAVE=true

    bwa index ${ref}

    bwa mem -p -t 8 ${ref} "${ref.simpleName}_consensus_interleaved.fastq" |
    samtools sort -@ 8 -o "${ref.simpleName}_aligned_consensus.bam"
    """
}