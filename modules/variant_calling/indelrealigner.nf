/*
 * Step 3 of the variant_calling subworkflow
 *
 * Realigns genomic regions to correct potential misalignments around indels.
 */
process GATK_INDELREALIGNER {

    label 'process_medium'
    container 'community.wave.seqera.io/library/gatk_samtools:5773d856edb307d7'

    input:
    tuple val(meta), val(ref_info), path(ref), path(bam), path(bai), path(intervals)

    output:
    tuple val(meta), val(ref_info), path(ref), path("*.bam"), path("*.bai"),    emit: realigned_bam

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2

    """
    samtools faidx "${ref}"
    samtools dict "${ref}" > "${ref.baseName}.dict"

    gatk3 -Xmx${avail_mem}g -T IndelRealigner \\
        -R ${ref} \\
        -targetIntervals ${intervals} \\
        -I ${bam} \\
        -o ${ref.baseName}.bam \\
        -maxReads 500000
    """
}
