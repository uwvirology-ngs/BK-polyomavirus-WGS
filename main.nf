#!/usr/bin/env nextflow

/*
 * Parameters
 */
params {
    input: Path
    db: Path
}

/*
 * Subworkflows
 */
include { READ_SAMPLESHEET                  } from './subworkflows/read_samplesheet.nf'
include { TWIST_I_UMI_EXTRACTION            } from './subworkflows/twist_i.nf'
include { REFERENCE_PREP                    } from './subworkflows/reference_prep.nf'
include { TWIST_II_ALIGNMENT                } from './subworkflows/twist_ii.nf'
include { TWIST_III_COLLAPSE                } from './subworkflows/twist_iii.nf'
include { VARIANT_CALLING as VC_INITIAL     } from './subworkflows/variant_calling.nf'
include { VARIANT_CALLING as VC_COLLAPSED   } from './subworkflows/variant_calling.nf'
include { CONSENSUS_ASSEMBLY                } from './subworkflows/consensus_assembly.nf'

/*
 * Modules
 */
include { SEQTK_SAMPLE                      } from './modules/core/seqtk_sample.nf'
include { FASTQC as FASTQC_RAW              } from './modules/core/fastqc.nf'
include { FASTQC as FASTQC_UMI_EXTRACTED    } from './modules/core/fastqc.nf'
include { FGBIO_COLLECT_DUPLEX_SEQ_METRICS  } from './modules/core/fgbio_collect_duplex_seq_metrics'
include { UMI_QC                            } from './modules/core/umi_qc.nf'
include { BUILD_ALIGNMENT_SUMMARY           } from './modules/core/build_alignment_summary.nf'
include { BUILD_RUN_SUMMARY                 } from './modules/core/build_run_summary.nf'

workflow {
    
    main:
    READ_SAMPLESHEET (
        params.input
    )

    FASTQC_RAW (
        READ_SAMPLESHEET.out.reads
    )

    SEQTK_SAMPLE (
        READ_SAMPLESHEET.out.reads
    )

    // Extract UMI sequences from reads, storing in an unaligned BAM file
    TWIST_I_UMI_EXTRACTION (
        SEQTK_SAMPLE.out.sampled_reads
    )

    FASTQC_UMI_EXTRACTED (
        TWIST_I_UMI_EXTRACTION.out.umi_extracted_fastq_interleaved
    )

    // ----------------------------------- REFERENCE PREP ------------------------------------

    // select reference genome by mapping umi-extracted reads to a database
    REFERENCE_PREP (
        TWIST_I_UMI_EXTRACTION.out.umi_extracted_fastq_paired,
        file(params.db)
    )

    // ----------------------------------- Back to Twist ------------------------------------

    // channel: [ val(meta), val(ref_info), path(ref), path(reads), path(umi_extracted_bam) ]
    alignment_ch = REFERENCE_PREP.out.ref
        .combine(TWIST_I_UMI_EXTRACTION.out.umi_extracted_fastq_interleaved, by: 0)
        .combine(TWIST_I_UMI_EXTRACTION.out.umi_extracted_bam, by: 0)
    
    TWIST_II_ALIGNMENT (
        alignment_ch
    )

    // ---------------------------------- CONSENSUS ASSEMBLY ---------------------------------

    consensus_ch = TWIST_II_ALIGNMENT.out.aligned_umi_extracted_bam
        .combine(TWIST_I_UMI_EXTRACTION.out.umi_extracted_fastq_interleaved, by: 0)

    CONSENSUS_ASSEMBLY (
        consensus_ch
    )

    // ---------------------------------------- TWIST ----------------------------------------

    FGBIO_COLLECT_DUPLEX_SEQ_METRICS (
        TWIST_II_ALIGNMENT.out.grouped_bam
    )

    if (params.umiInput != null) {
        umi_qc_ch = TWIST_II_ALIGNMENT.out.grouped_bam
            .join(TWIST_II_ALIGNMENT.out.grouped_hist, by: [0, 1])

        UMI_QC (
            umi_qc_ch,
            file(params.umiInput)
        )
    }

    TWIST_III_COLLAPSE (
        TWIST_II_ALIGNMENT.out.grouped_bam
    )

    // ----------------------------------- VARIANT CALLING ---------------------------------- 

    VC_INITIAL (
        TWIST_II_ALIGNMENT.out.merged_bam,
        true
    )

    VC_COLLAPSED (
        TWIST_III_COLLAPSE.out.final_consensus_bam,
        false
    )

    // --------------------------------------- SUMMARY ---------------------------------------

    alignment_summary_ch = TWIST_II_ALIGNMENT.out.merged_bam
        .join(TWIST_III_COLLAPSE.out.final_consensus_bam, by: [0, 1])
        .map { meta, ref_info, ref, bam1, _ref, bam2 -> tuple(meta, ref_info, ref, bam1, bam2) }
        .join(CONSENSUS_ASSEMBLY.out.final_consensus, by: [0, 1])
        .combine(READ_SAMPLESHEET.out.reads, by: 0)
        .map { meta, ref_info, ref, bam1, bam2, consensus_fa, reads -> tuple(
            meta, ref_info, ref, bam1, bam2, consensus_fa, reads,
            Utils.getGenomicRegion(ref_info.acc), 
            Utils.getCDSLen(ref_info.acc)
        ) }

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
    unaligned_bam                   = TWIST_I_UMI_EXTRACTION.out.unaligned_bam
    umi_extracted_bam               = TWIST_I_UMI_EXTRACTION.out.umi_extracted_bam
    umi_extracted_fastq_interleaved = TWIST_I_UMI_EXTRACTION.out.umi_extracted_fastq_interleaved
    umi_extracted_fastq_paired      = TWIST_I_UMI_EXTRACTION.out.umi_extracted_fastq_paired
    aligned_umi_extracted_bam       = TWIST_II_ALIGNMENT.out.aligned_umi_extracted_bam
    merged_bam                      = TWIST_II_ALIGNMENT.out.merged_bam
    grouped_bam                     = TWIST_II_ALIGNMENT.out.grouped_bam
    grouped_hist                    = TWIST_II_ALIGNMENT.out.grouped_hist
    duplex_seq_metrics              = FGBIO_COLLECT_DUPLEX_SEQ_METRICS.out.duplex_seq_metrics
    umi_qc                          = params.umiInput ? UMI_QC.out.umi_qc_plots : channel.empty()
    unaligned_consensus_bam         = TWIST_III_COLLAPSE.out.unaligned_consensus_bam
    aligned_consensus_bam           = TWIST_III_COLLAPSE.out.aligned_consensus_bam
    final_consensus_bam             = TWIST_III_COLLAPSE.out.final_consensus_bam
    alignment_summary               = BUILD_ALIGNMENT_SUMMARY.out.alignment_summary
    run_summary                     = BUILD_RUN_SUMMARY.out.run_summary

    // reference_prep
    alignment_to_db     = REFERENCE_PREP.out.alignment_to_db
    samtools_coverage   = REFERENCE_PREP.out.samtools_coverage
    covstats            = REFERENCE_PREP.out.covstats
    covstats_pass       = REFERENCE_PREP.out.covstats_pass
    covstats_fail       = REFERENCE_PREP.out.covstats_fail
    refs_tsv            = REFERENCE_PREP.out.refs_tsv
    selected_refs       = REFERENCE_PREP.out.selected_refs

    // variant calling
    rg_bam              = VC_INITIAL.out.rg_bam
    intervals           = VC_INITIAL.out.intervals
    realigned_bam       = VC_INITIAL.out.realigned_bam
    mpileup             = VC_INITIAL.out.mpileup
    variants            = VC_INITIAL.out.variants
    vcf                 = VC_INITIAL.out.vcf

    rg_bam_II           = VC_COLLAPSED.out.rg_bam
    intervals_II        = VC_COLLAPSED.out.intervals
    realigned_bam_II    = VC_COLLAPSED.out.realigned_bam
    mpileup_II          = VC_COLLAPSED.out.mpileup
    variants_II         = VC_COLLAPSED.out.variants
    vcf_II              = VC_COLLAPSED.out.vcf

    // fastqc
    fastqc_I            = FASTQC_RAW.out.fastqc
    fastqc_II           = FASTQC_UMI_EXTRACTED.out.fastqc

    // consensus assembly
    final_consensus     = CONSENSUS_ASSEMBLY.out.final_consensus
}

output {

    fastqc_I {
        path 'fastqc/raw'
    }

    fastqc_II {
        path 'fastqc/umi_extracted'
    }

    unaligned_bam {
        path 'twist/picard_fastq_to_sam'
    }
    umi_extracted_bam {
        path 'twist/fgbio_extract_umis_from_bam'
    }
    umi_extracted_fastq_interleaved {
        path 'twist/picard_sam_to_fastq/interleaved'
    }
    umi_extracted_fastq_paired {
        path 'twist/picard_sam_to_fastq/paired'
    }
    aligned_umi_extracted_bam {
        path 'twist/bwa_align_fastq'
    }
    merged_bam {
        path 'twist/picard_merge_bam_alignment'
    }
    grouped_bam {
        path 'twist/picard_group_reads_by_umi/bam'
    }
    grouped_hist {
        path 'twist/picard_group_reads_by_umi/hist'
    }
    duplex_seq_metrics {
        path 'twist/duplex_seq_metrics'
    }
    umi_qc {
        path 'twist/umi_qc_plots'
    }
    unaligned_consensus_bam {
        path 'twist/fgbio_call_duplex_consensus_reads'
    }
    aligned_consensus_bam {
        path 'twist/align_duplex_consensus_reads'
    }
    final_consensus_bam {
        path 'twist/picard_merge_consensus_bams'
    }
    alignment_summary {
        path 'summaries/alignment_summary'
    }
    run_summary {
        path 'summaries/run_summary'
    }

    // reference_prep
    alignment_to_db {
        path 'reference_prep/align_to_db_bams'
    }
    samtools_coverage {
        path 'reference_prep/samtools_coverage'
    }
    covstats {
        path 'reference_prep/covstats/all'
    }
    covstats_pass {
        path 'reference_prep/covstats/pass'
    }
    covstats_fail {
        path 'reference_prep/covstats/fail'
    }
    refs_tsv {
        path 'reference_prep/refs_tsv'
    }
    selected_refs {
        path 'reference_prep/selected_refs'
    }

    // variant calling 
    rg_bam {
        path 'variant_calling/picard_addorreplacereadgroups'
    }
    intervals {
        path 'variant_calling/gatk_realignertargetcreator'
    }
    realigned_bam {
        path 'variant_calling/gatk_indelrealigner'
    }
    mpileup {
        path 'variant_calling/samtools_mpileup'
    }
    variants {
        path 'variant_calling/ivar_variants'
    }
    vcf {
        path 'variant_calling/bcftools'
    }

    // collapsed
    rg_bam_II {
        path 'variant_calling_collapsed/picard_addorreplacereadgroups'
    }
    intervals_II {
        path 'variant_calling_collapsed/gatk_realignertargetcreator'
    }
    realigned_bam_II {
        path 'variant_calling_collapsed/gatk_indelrealigner'
    }
    mpileup_II {
        path 'variant_calling_collapsed/samtools_mpileup'
    }
    variants_II {
        path 'variant_calling_collapsed/ivar_variants'
    }
    vcf_II {
        path 'variant_calling_collapsed/bcftools'
    }

    // consensus assembly
    final_consensus {
        path 'consensus_assembly/final_consensus'
    }
}
