/*
 * Step 5 of the Twist Bioscience UMI Protocol
 *
 * Merge each aligned BAM file to the unaliged BAM file from step two 
 * which contain the UMI bases in as metadata. 
 */ 
process PICARD_MERGE_BAM_ALIGNMENT {

    label 'twist'
    label 'process_high'

    input:
    tuple val(meta), path(aligned_bam), path(ref), val(ref_info), path(umi_extracted_bam)

    output:
    tuple val(meta), path("*.bam"), path(ref), val(ref_info),  emit: merged_bam

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2
    
    """
    samtools dict ${ref} > ${ref.baseName}.dict

    picard MergeBamAlignment \\
        -Xmx${avail_mem}g \\
        UNMAPPED=${umi_extracted_bam} \\
        ALIGNED=${aligned_bam} \\
        O="${ref.simpleName}_merged.bam" \\
        R=${ref} \\
        CLIP_ADAPTERS=false \\
        VALIDATION_STRINGENCY=SILENT \\
        CREATE_INDEX=true \\
        EXPECTED_ORIENTATIONS=FR \\
        MAX_GAPS=-1 \\
        SO=coordinate \\
        ALIGNER_PROPER_PAIR_FLAGS=false
    """
}