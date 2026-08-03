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
include { READ_SAMPLESHEET   } from './subworkflows/core/read_samplesheet.nf'
include { REFERENCE_PREP     } from './subworkflows/reference_prep/reference_prep.nf'
include { CONSENSUS_ASSEMBLY } from './subworkflows/consensus_assembly/consensus_assembly.nf'

/*
 * Modules
 */
include { SEQTK_SAMPLE                      } from './modules/core/seqtk_sample.nf'

include { PICARD_FASTQ_TO_SAM               } from './modules/core/picard_fastq_to_sam.nf'
include { FGBIO_EXTRACT_UMIS_FROM_BAM       } from './modules/core/fgbio_extract_umis_from_bam.nf'
include { PICARD_SAM_TO_FASTQ               } from './modules/core/picard_sam_to_fastq.nf'
include { BWA_ALIGN_FASTQ                   } from './modules/core/bwa_align_fastq.nf'
include { PICARD_MERGE_BAM_ALIGNMENT        } from './modules/core/picard_merge_bam_alignment.nf'
include { FGBIO_GROUP_READS_BY_UMI          } from './modules/core/fgbio_group_reads_by_umi.nf'
include { FGBIO_CALL_DUPLEX_CONSENSUS_READS } from './modules/core/fgbio_call_duplex_consensus_reads.nf'
include { ALIGN_DUPLEX_CONSENSUS_READS      } from './modules/core/align_duplex_consensus_reads.nf'
include { PICARD_MERGE_CONSENSUS_BAMS       } from './modules/core/picard_merge_consensus_bams.nf'
include { BUILD_ALIGNMENT_SUMMARY           } from './modules/core/build_alignment_summary.nf'
include { BUILD_RUN_SUMMARY                 } from './modules/core/build_run_summary.nf'

include { PICARD_ADDORREPLACEREADGROUPS     } from './modules/variant_calling/addorreplacereadgroups.nf'
include { GATK_REALIGNERTARGETCREATOR       } from './modules/variant_calling/realignertargetcreator.nf'
include { GATK_INDELREALIGNER               } from './modules/variant_calling/indelrealigner.nf'
include { CDS_VARIANTS                      } from './modules/variant_calling/cds_variants.nf'

workflow {
    
    main:
    READ_SAMPLESHEET(
        params.input
    )

    SEQTK_SAMPLE (
        READ_SAMPLESHEET.out.reads
    )

    PICARD_FASTQ_TO_SAM (
        SEQTK_SAMPLE.out.sampled_reads
    )

    FGBIO_EXTRACT_UMIS_FROM_BAM (
        PICARD_FASTQ_TO_SAM.out.unaligned_bam
    )

    PICARD_SAM_TO_FASTQ (
        FGBIO_EXTRACT_UMIS_FROM_BAM.out.umi_extracted_bam
    )

    // ----------------------------------- REFERENCE PREP ------------------------------------

    // select reference genome with revica-strm and combine with umi-extracted FASTQs
    REFERENCE_PREP (
        PICARD_SAM_TO_FASTQ.out.umi_extracted_fastq_paired,
        file(params.db),
        false
    )

    // ----------------------------------- Back to Twist ------------------------------------

    bwa_align_ch = PICARD_SAM_TO_FASTQ.out.umi_extracted_fastq_interleaved
        .combine(REFERENCE_PREP.out.ref, by: 0)
        .map {meta, reads, ref_info, ref -> tuple(meta, reads, ref, ref_info)}

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

    // ----------------------------------- VARIANT CALLING ----------------------------------- 

    PICARD_ADDORREPLACEREADGROUPS (
        PICARD_MERGE_BAM_ALIGNMENT.out.merged_bam,
    )

    GATK_REALIGNERTARGETCREATOR (
        PICARD_ADDORREPLACEREADGROUPS.out.bam
    )

    GATK_INDELREALIGNER (
        GATK_REALIGNERTARGETCREATOR.out.intervals
    )

    // tuple val(meta), path(bam), path(bai), path(ref), path(gff), region, val(save_mpileup)
    // needs ref and gff and save_mpileup
    variants_ch = GATK_INDELREALIGNER.out.bam
        .map { meta, bam, bai, ref, ref_info -> tuple(
            meta, bam, bai, ref, ref_info,
            "${projectDir}/assets/database/${ref_info.acc}.gff", 
            Utils.getGenomicRegion(ref_info.acc),
            false
        )}

    CDS_VARIANTS (
        variants_ch
    )

    // -------------------------------------- CONSENSUS --------------------------------------

    CONSENSUS_ASSEMBLY (
        REFERENCE_PREP.out.reads,
        REFERENCE_PREP.out.ref,
        false
    )

    // ---------------------------------------- TWIST ----------------------------------------

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

    alignment_summary_ch = PICARD_MERGE_BAM_ALIGNMENT.out.merged_bam
        .join(PICARD_MERGE_CONSENSUS_BAMS.out.final_consensus_bam)
        .map { meta, bam1, ref1, ref_info1, bam2, _ref2, _ref_info2 -> tuple(meta, bam1, bam2, ref1, ref_info1) }

    BUILD_ALIGNMENT_SUMMARY (
        alignment_summary_ch
    )

    run_summary_ch = BUILD_ALIGNMENT_SUMMARY.out.alignment_summary
        .map { _meta, tsv -> tsv }
        .collect()

    BUILD_RUN_SUMMARY (
        run_summary_ch
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
    run_summary                     = BUILD_RUN_SUMMARY.out.run_summary
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
        path 'summaries/alignment_summary'
    }
    run_summary {
        path 'summaries/run_summary'
    }
}