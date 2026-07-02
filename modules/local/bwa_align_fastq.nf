/*
 * Step 4 of the Twist Bioscience UMI Protocol
 * 
 * Align the UMI extracted FastQ files to the reference genome using BWA.
 */
process BWA_ALIGN_FASTQ {

    label 'twist'
    label 'process_high'

    input: 
    tuple val(meta), path(reads), path(ref), val(ref_info)

    output: 
    tuple val(meta), path("*.bam"), path(ref),  emit: aligned_umi_extracted_bam

    script: 
    """
    bwa index ${ref}

    bwa mem -p -t 8 ${ref} ${reads} \\
    | samtools sort -@ 8 -o "${ref.simpleName}_aligned_umi_extracted.bam"
    """
}