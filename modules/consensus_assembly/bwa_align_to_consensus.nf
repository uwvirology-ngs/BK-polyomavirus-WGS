/*
 * Step 2 in the consensus assembly workflow. 
 * 
 * Aligns reads to the initial consensus sequence with bwa. 
 */
process BWA_ALIGN_TO_CONSENSUS {

    label 'process_high'
    container 'quay.io/epil02/revica-strm:0.0.5'

    input:
    tuple val(meta), val(ref_info), path(reads), path(consensus_fa)

    output:
    tuple val(meta), val(ref_info), path("*.bam"),  emit: realigned_bam

    script:
    """
    bwa index ${consensus_fa}

    # read unmapped (0x4), supplementary alignment (0x800)
    FLAG=2052

    bwa mem -p -t 8 ${consensus_fa} ${reads} \\
    | samtools sort -@ ${task.cpus} \\
    | samtools view -b -F \$FLAG -@ ${task.cpus} > "${ref_info.acc}_realigned.bam"
    """
}
