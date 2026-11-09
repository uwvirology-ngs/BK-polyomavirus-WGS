/*
 * Modules
 */
include { FGBIO_CALL_DUPLEX_CONSENSUS_READS } from '../modules/core/fgbio_call_duplex_consensus_reads.nf'
include { ALIGN_DUPLEX_CONSENSUS_READS      } from '../modules/core/align_duplex_consensus_reads.nf'
include { PICARD_MERGE_CONSENSUS_BAMS       } from '../modules/core/picard_merge_consensus_bams.nf'

/*
 * Collapses the BAM file by UMI and, if set, filters for duplex consensus sequences. 
 */
workflow TWIST_III_COLLAPSE {

    take:
    input_ch    // channel: [ val(meta), val(ref_info), path(ref), path(grouped_bam) ]

    main:
    FGBIO_CALL_DUPLEX_CONSENSUS_READS (
        input_ch
    )

    ALIGN_DUPLEX_CONSENSUS_READS (
        FGBIO_CALL_DUPLEX_CONSENSUS_READS.out.unaligned_consensus_bam
    )

    consensus_ch = FGBIO_CALL_DUPLEX_CONSENSUS_READS.out.unaligned_consensus_bam
        .join(ALIGN_DUPLEX_CONSENSUS_READS.out.aligned_consensus_bam, by: [0, 1])
        .map { meta, ref_info, ref, bam1, _ref, bam2 -> tuple(meta, ref_info, ref, bam1, bam2) }

    PICARD_MERGE_CONSENSUS_BAMS (
        consensus_ch
    )

    emit:
    unaligned_consensus_bam = FGBIO_CALL_DUPLEX_CONSENSUS_READS.out.unaligned_consensus_bam
    aligned_consensus_bam   = ALIGN_DUPLEX_CONSENSUS_READS.out.aligned_consensus_bam
    final_consensus_bam     = PICARD_MERGE_CONSENSUS_BAMS.out.final_consensus_bam
}
