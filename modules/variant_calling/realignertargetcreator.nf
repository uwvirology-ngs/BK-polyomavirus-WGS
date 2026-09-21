/*
 * Step 2 of the variant_calling subworkflow
 *
 * Identifies genomic intervals which require local realignment around indels
 * by the IndelRealigner tool used in the next step.
 */
process GATK_REALIGNERTARGETCREATOR {

    label 'process_low'
    container 'community.wave.seqera.io/library/gatk_samtools:5773d856edb307d7'

    input:
    tuple val(meta), val(ref_info), path(ref), path(bam), path(bai)

    output:
    tuple val(meta), val(ref_info), path(ref), path(bam), path(bai), path("*.intervals"),   emit: intervals

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2

    """
    samtools faidx "${ref}"
    samtools dict "${ref}" > "${ref.baseName}.dict"

    gatk3 -Xmx${avail_mem}g -T RealignerTargetCreator \\
        -I ${bam} \\
        -R ${ref} \\
        -o "${ref.baseName}.intervals" \\
        -nt ${task.cpus}
    """
}
