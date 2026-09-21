/*
 * Step 4 of the Twist Bioscience UMI Protocol
 * 
 * Align the UMI extracted FastQ files to the reference genome using BWA.
 */
process BWA_ALIGN_FASTQ {

    label 'twist'
    label 'process_high'

    input: 
    tuple val(meta), val(ref_info), path(ref), path(reads)

    output: 
    tuple val(meta), val(ref_info), path(ref), path("*.bam"),   emit: aligned_umi_extracted_bam

    script: 
    """
    bwa index ${ref}

    bwa mem -p -t 8 ${ref} ${reads} \\
    | samtools sort -@ 8 -o "${ref.simpleName}_aligned_umi_extracted.bam"
    """
}