process GATK_REALIGNERTARGETCREATOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/gatk_samtools:5773d856edb307d7'
    // gatk=3.8, samtools=1.23.1

    input:
    tuple val(meta), path(bam), path(bai), path(ref), val(ref_info)

    output:
    tuple val(meta), path(bam), path(bai), path(ref), val(ref_info), path("*.intervals"), emit: intervals
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def known = ""
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK RealignerTargetCreator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    samtools faidx "${ref}"
    samtools dict "${ref}" > "${ref.baseName}.dict"

    gatk3 \\
        -Xmx${avail_mem}M \\
        -T RealignerTargetCreator \\
        -nt ${task.cpus} \\
        -I ${bam} \\
        -R ${ref} \\
        -o ${prefix}.intervals \\
        ${known} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(echo \$(gatk3 --version))
    END_VERSIONS
    """
}
