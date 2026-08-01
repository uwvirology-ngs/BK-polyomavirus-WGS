/*
 * Step 3 of the Twist Bioscience UMI Protocol
 *
 * With the UMI sequences extracted, convert the unaligned BAM files back to 
 * FASTQ format in preparation for reference selection and alignment.
 * 
 * Generates output FASTQs in interleaved and paired format for continuation with the 
 * Twist protocol and input into REVICA-STRM for reference selection, respectively. 
 */
process PICARD_SAM_TO_FASTQ {

    label 'twist'
    label 'process_high'

    input:
    tuple val(meta), path(unaligned_bam_umi_extracted)

    output:
    tuple val(meta), path("${meta.id}_umi_extracted_i*"),   emit: umi_extracted_fastq_interleaved
    tuple val(meta), path("${meta.id}_umi_extracted_R*"),   emit: umi_extracted_fastq_paired

    script:
    def avail_mem = 8
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }
    avail_mem -= 2
    
    """
    picard SamToFastq \\
        -Xmx${avail_mem}g \\
        I="${unaligned_bam_umi_extracted}" \\
        F="${meta.id}_umi_extracted_interleaved.fastq" \\
        INTERLEAVE=true

    # generate paired fastqs for REVICA-STRM input
    seqtk seq -1 "${meta.id}_umi_extracted_interleaved.fastq" > "${meta.id}_umi_extracted_R1.fastq"
    seqtk seq -2 "${meta.id}_umi_extracted_interleaved.fastq" > "${meta.id}_umi_extracted_R2.fastq"
    """
}