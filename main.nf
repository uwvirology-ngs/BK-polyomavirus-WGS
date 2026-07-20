#!/usr/bin/env nextflow

/*
 * Parameters
 */
params {
    input: Path
    ref: Path
    db: Path
}

/*
 * Subworkflows
 */
include { READ_SAMPLESHEET } from './subworkflows/local/read_samplesheet.nf'
include { REFERENCE_PREP   } from './subworkflows/revica/reference_prep'

/*
 * Modules
 */
include { PICARD_FASTQ_TO_SAM               } from './modules/local/picard_fastq_to_sam.nf'
include { FGBIO_EXTRACT_UMIS_FROM_BAM       } from './modules/local/fgbio_extract_umis_from_bam.nf'
include { PICARD_SAM_TO_FASTQ               } from './modules/local/picard_sam_to_fastq.nf'
include { BWA_ALIGN_FASTQ                   } from './modules/local/bwa_align_fastq.nf'
include { PICARD_MERGE_BAM_ALIGNMENT        } from './modules/local/picard_merge_bam_alignment.nf'
include { FGBIO_GROUP_READS_BY_UMI          } from './modules/local/fgbio_group_reads_by_umi.nf'
include { FGBIO_CALL_DUPLEX_CONSENSUS_READS } from './modules/local/fgbio_call_duplex_consensus_reads.nf'
include { ALIGN_DUPLEX_CONSENSUS_READS      } from './modules/local/align_duplex_consensus_reads.nf'
include { PICARD_MERGE_CONSENSUS_BAMS       } from './modules/local/picard_merge_consensus_bams.nf'
include { BUILD_ALIGNMENT_SUMMARY           } from './modules/local/build_alignment_summary.nf'

workflow {
    
    main:
    READ_SAMPLESHEET(
        params.input
    )

    PICARD_FASTQ_TO_SAM (
        READ_SAMPLESHEET.out.reads
    )

    FGBIO_EXTRACT_UMIS_FROM_BAM (
        PICARD_FASTQ_TO_SAM.out.unaligned_bam
    )

    PICARD_SAM_TO_FASTQ (
        FGBIO_EXTRACT_UMIS_FROM_BAM.out.umi_extracted_bam
    )

    // if a reference genome is not provided, select with revica-strm and combine
    // with umi-extracted FASTQs. otherwise, use provided reference.
    if (params.ref == null) {    
        REFERENCE_PREP (
            PICARD_SAM_TO_FASTQ.out.umi_extracted_fastq_paired,
            file(params.db),
            false
        )

        bwa_align_ch = PICARD_SAM_TO_FASTQ.out.umi_extracted_fastq_interleaved
            .combine(REFERENCE_PREP.out.ref, by: 0)
            .map {meta, reads, ref_info, ref -> tuple(meta, reads, ref, ref_info)}
    } else {
        bwa_align_ch = PICARD_SAM_TO_FASTQ.out.umi_extracted_fastq_interleaved
            .map {meta, reads -> tuple(meta, reads, file(params.ref), null)}
    }

    BWA_ALIGN_FASTQ(
        bwa_align_ch
    )

    // combine each aligned BAM with the correct UMI-extracted bam from Twist step 2 
    // and record the reference genome used in metadata
    merged_bams_ch = BWA_ALIGN_FASTQ.out.aligned_umi_extracted_bam
        .combine(FGBIO_EXTRACT_UMIS_FROM_BAM.out.umi_extracted_bam, by: 0)
        .map { meta, bam1, ref, ref_info, bam2 -> tuple(meta + [pair_name: ref.baseName], bam1, ref, ref_info, bam2) }

    PICARD_MERGE_BAM_ALIGNMENT(
        merged_bams_ch
    )

    FGBIO_GROUP_READS_BY_UMI(
        PICARD_MERGE_BAM_ALIGNMENT.out.merged_bam
    )

    FGBIO_CALL_DUPLEX_CONSENSUS_READS(
        FGBIO_GROUP_READS_BY_UMI.out.grouped_bam
    )

    ALIGN_DUPLEX_CONSENSUS_READS(
        FGBIO_CALL_DUPLEX_CONSENSUS_READS.out.unaligned_consensus_bam
    )

    consensus_ch = FGBIO_CALL_DUPLEX_CONSENSUS_READS.out.unaligned_consensus_bam
        .join(ALIGN_DUPLEX_CONSENSUS_READS.out.aligned_consensus_bam)
        .map { meta, bam1, ref1, ref_info1, bam2, _ref2, _ref_info2 -> tuple(meta, bam1, bam2, ref1, ref_info1) }

    PICARD_MERGE_CONSENSUS_BAMS (
        consensus_ch
    )

    summary_ch = PICARD_MERGE_BAM_ALIGNMENT.out.merged_bam
        .join(PICARD_MERGE_CONSENSUS_BAMS.out.final_consensus_bam)
        .map { meta, bam1, ref1, ref_info1, bam2, _ref2, _ref_info2 -> tuple(meta, bam1, bam2, ref1, ref_info1) }

    BUILD_ALIGNMENT_SUMMARY (
        summary_ch
    )

    publish:
    unaligned_bam                   = PICARD_FASTQ_TO_SAM.out.unaligned_bam
    umi_extracted_bam               = FGBIO_EXTRACT_UMIS_FROM_BAM.out.umi_extracted_bam
    umi_extracted_fastq_interleaved = PICARD_SAM_TO_FASTQ.out.umi_extracted_fastq_interleaved
    umi_extracted_fastq_paired      = PICARD_SAM_TO_FASTQ.out.umi_extracted_fastq_paired
    aligned_umi_extracted_bam       = BWA_ALIGN_FASTQ.out.aligned_umi_extracted_bam
    merged_bam                      = PICARD_MERGE_BAM_ALIGNMENT.out.merged_bam
    grouped_bam                     = FGBIO_GROUP_READS_BY_UMI.out.grouped_bam
    grouped_hist                    = FGBIO_GROUP_READS_BY_UMI.out.grouped_hist
    unaligned_consensus_bam         = FGBIO_CALL_DUPLEX_CONSENSUS_READS.out.unaligned_consensus_bam
    aligned_consensus_bam           = ALIGN_DUPLEX_CONSENSUS_READS.out.aligned_consensus_bam
    final_consensus_bam             = PICARD_MERGE_CONSENSUS_BAMS.out.final_consensus_bam
    alignment_summary               = BUILD_ALIGNMENT_SUMMARY.out.alignment_summary
}

output {
    unaligned_bam {
        path 'picard_fastq_to_sam'
    }
    umi_extracted_bam {
        path 'fgbio_extract_umis_from_bam'
    }
    umi_extracted_fastq_interleaved {
        path 'picard_sam_to_fastq/interleaved'
    }
    umi_extracted_fastq_paired {
        path 'picard_sam_to_fastq/paired'
    }
    aligned_umi_extracted_bam {
        path 'bwa_align_fastq'
    }
    merged_bam {
        path 'picard_merge_bam_alignment'
    }
    grouped_bam {
        path 'picard_group_reads_by_umi/bam'
    }
    grouped_hist {
        path 'picard_group_reads_by_umi/hist'
    }
    unaligned_consensus_bam {
        path 'fgbio_call_duplex_consensus_reads'
    }
    aligned_consensus_bam {
        path 'align_duplex_consensus_reads'
    }
    final_consensus_bam {
        path 'picard_merge_consensus_bams'
    }
    alignment_summary {
        path 'run_summary/alignment_summary'
    }
}