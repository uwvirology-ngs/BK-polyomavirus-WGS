/*
 * Step 10 of the Twist Bioscience UMI Protocol
 *
 * Merge the aligned and unaligned consensus BAM files to produce finalized reads.
 */ 
process PICARD_MERGE_CONSENSUS_BAMS {

    label 'twist'
    label 'process_high'

    input:
    tuple val(meta), path(unaligned_consensus_bam), path(aligned_consensus_bam), path(ref), val(ref_info)

    output:
    tuple val(meta), path("*_consensus_aligned_merged.bam"), path(ref), val(ref_info), emit: final_consensus_bam

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2
    
    """
    picard CreateSequenceDictionary \\
        R=${ref}
        O="${ref.baseName}.dict"

    picard MergeBamAlignment \\
        -Xmx${avail_mem}g \\
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

    picard AddOrReplaceReadGroups \\
        I="${ref.simpleName}_merged_no_read_group_consensus.bam" \\
        O="${ref.simpleName}_consensus_aligned_merged.bam" \\
        RGID=${ref.simpleName} \\
        RGLB=sample_lib \\
        RGPL=Illumina RGSM=sample_name \\
        RGPU=NA
    """
}