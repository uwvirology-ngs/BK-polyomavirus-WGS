process GATK_REALIGNERTARGETCREATOR {

    label 'process_low'
    container 'community.wave.seqera.io/library/gatk_samtools:5773d856edb307d7'

    input:
    tuple val(meta), path(bam), path(bai), path(ref), val(ref_info)

    output:
    tuple val(meta), path(bam), path(bai), path(ref), val(ref_info), path("*.intervals"), emit: intervals

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
        -T RealignerTargetCreator \\
        -nt ${task.cpus} \\
        -I ${bam} \\
        -R ${ref} \\
        -o ${meta.id}.intervals
    """
}
