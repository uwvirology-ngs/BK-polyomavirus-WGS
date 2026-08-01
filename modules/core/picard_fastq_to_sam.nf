/*
 * Step 1 of the Twist Bioscience UMI Protocol
 * 
 * Convert the FASTQ files to unaligned BAMs in preparation for 
 * extraction of the UMIs. 
 */
process PICARD_FASTQ_TO_SAM {

    label 'twist'
    label 'process_high'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_unaligned.bam"),  emit: unaligned_bam

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2
    
    """
    picard FastqToSam \\
        -Xmx${avail_mem}g \\
        O="${meta.id}_unaligned.bam" \\
        F1=${meta.id}_R1.sample.fastq.gz \\
        F2=${meta.id}_R2.sample.fastq.gz \\
        SM=${meta.id} \\
        LB=Library1 \\
        PU=Unit1 \\
        PL=Illumina
    """
}