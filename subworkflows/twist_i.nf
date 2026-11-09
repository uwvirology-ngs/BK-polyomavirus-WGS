/*
 * Modules
 */
include { PICARD_FASTQ_TO_SAM         } from '../modules/core/picard_fastq_to_sam.nf'
include { FGBIO_EXTRACT_UMIS_FROM_BAM } from '../modules/core/fgbio_extract_umis_from_bam.nf'
include { PICARD_SAM_TO_FASTQ         } from '../modules/core/picard_sam_to_fastq.nf'

/*
 * Extracts UMI sequences from reads, storing in an unaligned BAM file.
 */
workflow TWIST_I_UMI_EXTRACTION {

    take:
    input_ch    // channel: [ val(meta), path(reads) ]

    main:
    PICARD_FASTQ_TO_SAM (
        input_ch
    )

    FGBIO_EXTRACT_UMIS_FROM_BAM (
        PICARD_FASTQ_TO_SAM.out.unaligned_bam
    )

    PICARD_SAM_TO_FASTQ (
        FGBIO_EXTRACT_UMIS_FROM_BAM.out.umi_extracted_bam
    ) 

    emit:
    unaligned_bam                   = PICARD_FASTQ_TO_SAM.out.unaligned_bam
    umi_extracted_bam               = FGBIO_EXTRACT_UMIS_FROM_BAM.out.umi_extracted_bam
    umi_extracted_fastq_interleaved = PICARD_SAM_TO_FASTQ.out.umi_extracted_fastq_interleaved
    umi_extracted_fastq_paired      = PICARD_SAM_TO_FASTQ.out.umi_extracted_fastq_paired
}
