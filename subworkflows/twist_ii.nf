/*
 * Modules
 */
include { BWA_ALIGN_FASTQ            } from '../modules/core/bwa_align_fastq.nf'
include { PICARD_MERGE_BAM_ALIGNMENT } from '../modules/core/picard_merge_bam_alignment'
include { FGBIO_GROUP_READS_BY_UMI   } from '../modules/core/fgbio_group_reads_by_umi.nf'

/*
 * Aligns the umi_extracted fastq files to the selected reference genome, 
 * then adds BAM headers with UMI sequences and filters to mapped reads
 * with acceptable quality scores grouped by UMI. 
 */
workflow TWIST_II_ALIGNMENT {

    take:
    input_ch    // channel: [ val(meta), val(ref_info), path(ref), path(reads), path(umi_extracted_bam) ]

    main:
    input_ch
        .multiMap { meta, ref_info, ref, reads, umi_extracted_bam -> 
            bwa_align_input:    tuple(meta, ref_info, ref, reads)
            umi_extracted_bam:  tuple(meta, umi_extracted_bam)
        } 
        .set { input }

    BWA_ALIGN_FASTQ (
        input.bwa_align_input
    ) 

    // channel: [ val(meta), val(ref_info), path(ref), path(aligned_bam), path(umi_extracted_bam) ]
    merged_bams_ch = BWA_ALIGN_FASTQ.out.aligned_umi_extracted_bam
        .join(input.umi_extracted_bam, by: 0)

    PICARD_MERGE_BAM_ALIGNMENT (
        merged_bams_ch
    )

    FGBIO_GROUP_READS_BY_UMI (
        PICARD_MERGE_BAM_ALIGNMENT.out.merged_bam
    )

    emit:
    aligned_umi_extracted_bam   = BWA_ALIGN_FASTQ.out.aligned_umi_extracted_bam
    merged_bam                  = PICARD_MERGE_BAM_ALIGNMENT.out.merged_bam
    grouped_bam                 = FGBIO_GROUP_READS_BY_UMI.out.grouped_bam
    grouped_hist                = FGBIO_GROUP_READS_BY_UMI.out.grouped_hist
}
