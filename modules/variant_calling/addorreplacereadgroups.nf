/*
 * Assigns all samples to the same read group for downstream processing
 */
process PICARD_ADDORREPLACEREADGROUPS {

    label 'process_low'
    container 'quay.io/biocontainers/picard:3.3.0--hdfd78af_0'

    input:
    tuple val(meta), path(merged_bam), path(ref), val(ref_info)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), path(ref), val(ref_info), emit: bam

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2

    """
    picard AddOrReplaceReadGroups \\
        -Xmx${avail_mem}g \\
        --SORT_ORDER coordinate \\
        --RGID foo \\
        --RGLB bar \\
        --RGPL illumina \\
        --RGPU unit1 \\
        --RGSM ${meta.id} \\
        --CREATE_INDEX true \\
        --INPUT ${merged_bam} \\
        --OUTPUT ${ref.baseName}_rg.bam || true
    """
}
