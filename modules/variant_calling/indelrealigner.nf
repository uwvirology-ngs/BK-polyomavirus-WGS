process GATK_INDELREALIGNER {

    label 'process_single'
    container 'community.wave.seqera.io/library/gatk_samtools:5773d856edb307d7'

    input:
    tuple val(meta), path(bam), path(bai), path(ref), val(ref_info), path(intervals)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), path(ref), val(ref_info),    emit: bam

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2

    """
    samtools faidx "${ref}"
    samtools dict "${ref}" > "${ref.baseName}.dict"

    gatk3 \\
        -Xmx${avail_mem}g \\
        -T IndelRealigner \\
        -R ${ref} \\
        -I ${bam} \\
        --targetIntervals ${intervals} \\
        -o ${ref.baseName}.bam \\
        -maxReads 500000
    """
}
