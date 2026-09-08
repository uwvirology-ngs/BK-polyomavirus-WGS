/*
 * Step 1 of the variant_calling subworkflow
 *
 * Ensures all sample are assigned to the same read group.
 */
process PICARD_ADDORREPLACEREADGROUPS {

    label 'process_low'
    container 'quay.io/biocontainers/picard:3.3.0--hdfd78af_0'

    input:
    tuple val(meta), path(input_bam), path(ref), val(ref_info)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), path(ref), val(ref_info),    emit: rg_bam

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2

    """
    picard AddOrReplaceReadGroups -Xmx${avail_mem}g \\
        --INPUT ${input_bam} \\
        --OUTPUT "${ref.baseName}_rg.bam" \\
        --RGID 4 \\
        --RGLB Library1 \\
        --RGPL Illumina \\
        --RGPU Unit1 \\
        --RGSM ${meta.id} \\
        --CREATE_INDEX true \\
        --SORT_ORDER coordinate
    """
}
