/*
 * Step 10 of the Twist Bioscience UMI Protocol
 *
 * Merge the aligned and unaligned consensus BAM files to produce finalized reads.
 */ 
process PICARD_MERGE_CONSENSUS_BAMS {

    label 'twist'

    input:
    tuple val(meta), path(unaligned_consensus_bam), path(aligned_consensus_bam), path(ref)

    output:
    tuple val(meta), path("*.bam"), path(ref)

    script:
    """
    samtools dict ${ref} > ${ref.baseName}.dict

    picard MergeBamAlignment \\
        UNMAPPED=${unaligned_consensus_bam} \\
        ALIGNED=${aligned_consensus_bam} \\
        O="${ref.simpleName}_merged_no_read_group_consensus.bam" \\
        R=${ref} \\
        CLIP_ADAPTERS=false \\
        VALIDATION_STRINGENCY=SILENT \\
        CREATE_INDEX=true \\
        EXPECTED_ORIENTATIONS=FR \\
        MAX_GAPS=-1 \\
        SORT_ORDER=coordinate \\
        ALIGNER_PROPER_PAIR_FLAGS=false \\
        ATTRIBUTES_TO_RETAIN=X0 \\
        ATTRIBUTES_TO_RETAIN=ZS \\
        ATTRIBUTES_TO_RETAIN=ZI \\
        ATTRIBUTES_TO_RETAIN=ZM \\
        ATTRIBUTES_TO_RETAIN=ZC \\
        ATTRIBUTES_TO_RETAIN=ZN \\
        ATTRIBUTES_TO_RETAIN=ad \\
        ATTRIBUTES_TO_RETAIN=bd \\
        ATTRIBUTES_TO_RETAIN=cd \\
        ATTRIBUTES_TO_RETAIN=ae \\
        ATTRIBUTES_TO_RETAIN=be \\
        ATTRIBUTES_TO_RETAIN=ce
    """
}