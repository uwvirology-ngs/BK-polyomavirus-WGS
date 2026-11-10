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
include { TWIST_I_UMI_EXTRACTION            } from './subworkflows/twist_i.nf'
include { TWIST_II_ALIGNMENT                } from './subworkflows/twist_ii.nf'
include { TWIST_III_COLLAPSE                } from './subworkflows/twist_iii.nf'
include { READ_SAMPLESHEET                  } from './subworkflows/read_samplesheet.nf'
include { REFERENCE_PREP                    } from './subworkflows/reference_prep.nf'
include { VARIANT_CALLING as VC_INITIAL     } from './subworkflows/variant_calling.nf'
include { VARIANT_CALLING as VC_COLLAPSED   } from './subworkflows/variant_calling.nf'
include { CONSENSUS_ASSEMBLY                } from './subworkflows/consensus_assembly.nf'

/*
 * Modules
 */
include { SEQTK_SAMPLE                      } from './modules/core/seqtk_sample.nf'
include { FASTQC                            } from './modules/core/fastqc.nf'
include { MULTIQC                           } from './modules/core/multiqc.nf'
include { FGBIO_COLLECT_DUPLEX_SEQ_METRICS  } from './modules/core/fgbio_collect_duplex_seq_metrics'
include { UMI_QC                            } from './modules/core/umi_qc.nf'
include { SAMTOOLS_SORT                     } from './modules/core/samtools_sort.nf'
include { BUILD_ALIGNMENT_SUMMARY           } from './modules/core/build_alignment_summary.nf'
include { BUILD_RUN_SUMMARY                 } from './modules/core/build_run_summary.nf'

workflow {
    
    main:
    // --------------------------------------- FASTQ Input ----------------------------------------
    
    // Parse samplesheet and retrieve paired-end fastq files
    READ_SAMPLESHEET (
        params.input
    )

    // Downsample reads as needed
    SEQTK_SAMPLE (
        READ_SAMPLESHEET.out.reads
    )

    // -------------------------------------- UMI Extraction --------------------------------------

    // Extract UMI sequences from reads, storing in an unaligned BAM file
    TWIST_I_UMI_EXTRACTION (
        SEQTK_SAMPLE.out.sampled_reads
    )

    // ------------------------------------- Quality Control --------------------------------------

    // Run FASTQC on UMI-extracted reads
    FASTQC (
        TWIST_I_UMI_EXTRACTION.out.umi_extracted_fastq_interleaved
    )

    // Aggregate FASTQC results and report as HTML
    MULTIQC (
        FASTQC.out.fastqc.collect()
    )

    // -------------------------------------- Reference Prep --------------------------------------

    // Select reference genome(s) by mapping umi-extracted reads to a BKPyV sequence database
    REFERENCE_PREP (
        TWIST_I_UMI_EXTRACTION.out.umi_extracted_fastq_paired,
        file(params.db)
    )

    // ---------------------------------- Alignment to Reference ----------------------------------

    // channel: [ val(meta), val(ref_info), path(ref), path(reads), path(umi_extracted_bam) ]
    alignment_ch = REFERENCE_PREP.out.ref
        .combine(TWIST_I_UMI_EXTRACTION.out.umi_extracted_fastq_interleaved, by: 0)
        .combine(TWIST_I_UMI_EXTRACTION.out.umi_extracted_bam, by: 0)
    
    // Align umi-extracted reads to selected reference(s) and build BAM files grouped by UMI
    TWIST_II_ALIGNMENT (
        alignment_ch
    )

    // Generate coordinate-sorted BAM files for input into consensus assembly and variant calling
    SAMTOOLS_SORT (
        TWIST_II_ALIGNMENT.out.grouped_bam
    )

    // ---------------------------------- Quality Control (UMI) -----------------------------------

    // Collect QC metrics for duplex sequencing data from FGBIO output
    FGBIO_COLLECT_DUPLEX_SEQ_METRICS (
        TWIST_II_ALIGNMENT.out.grouped_bam
    )

    // If a UMI list is provided, build additional QC plots showing UMI distribution
    if (params.umiInput != null) {
        umi_qc_ch = TWIST_II_ALIGNMENT.out.grouped_bam
            .join(TWIST_II_ALIGNMENT.out.grouped_hist, by: [0, 1])

        UMI_QC (
            umi_qc_ch,
            file(params.umiInput)
        )
    }

    // --------------------------------------- UMI Collapse ---------------------------------------

    // Build final consensus BAM file
    TWIST_III_COLLAPSE (
        TWIST_II_ALIGNMENT.out.grouped_bam
    )

    // ------------------------------------ Consensus Assembly ------------------------------------

     // channel: [ val(meta), path(ref_info), path(ref), path(bam), path(reads) ]
    consensus_ch = SAMTOOLS_SORT.out.grouped_coord_sorted_bam
        .combine(TWIST_I_UMI_EXTRACTION.out.umi_extracted_fastq_interleaved, by: 0)

    // Assemble consensus genome(s) against each selected reference
    CONSENSUS_ASSEMBLY (
        consensus_ch
    )

    // ----------------------------------- Variant Calling ---------------------------------- 

    // Call variants against BAM file from initial alignment
    VC_INITIAL (
        SAMTOOLS_SORT.out.grouped_coord_sorted_bam,
        true
    )

    // Call variants against final consensus BAM file
    VC_COLLAPSED (
        TWIST_III_COLLAPSE.out.final_consensus_bam,
        false   // skip indel realignment
    )

    // --------------------------------------- SUMMARY ---------------------------------------

    // channel: [ val(meta), val(ref_info), path(ref), path(merged_bam), 
    //            path(consensus_bam), path(consensus_fa), path(reads), 
    //            val(genomic_region), val(genomic_region_len)           ]
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

    // Generate a summary tsv file for every alignment
    BUILD_ALIGNMENT_SUMMARY (
        alignment_summary_ch
    )

    // Compile the summary files into a run summary report
    BUILD_RUN_SUMMARY (
        BUILD_ALIGNMENT_SUMMARY.out.alignment_summary.collect()
    )

    publish:
    // core pipeline
    umi_extracted_bam   = TWIST_I_UMI_EXTRACTION.out.umi_extracted_bam
    grouped_bam         = TWIST_II_ALIGNMENT.out.grouped_bam
    grouped_hist        = TWIST_II_ALIGNMENT.out.grouped_hist
    final_consensus_bam = TWIST_III_COLLAPSE.out.final_consensus_bam

    // reference selection
    alignment_to_db     = REFERENCE_PREP.out.alignment_to_db
    samtools_coverage   = REFERENCE_PREP.out.samtools_coverage
    covstats_all        = REFERENCE_PREP.out.covstats
    covstats_pass       = REFERENCE_PREP.out.covstats_pass
    covstats_fail       = REFERENCE_PREP.out.covstats_fail

    // variant calling | umi-grouped bam (I) vs. final consensus bam (II)
    mpileup_I           = VC_INITIAL.out.mpileup
    variants_I          = VC_INITIAL.out.variants
    vcf_I               = VC_INITIAL.out.vcf

    mpileup_II          = VC_COLLAPSED.out.mpileup
    variants_II         = VC_COLLAPSED.out.variants
    vcf_II              = VC_COLLAPSED.out.vcf

    // consensus genome assembly
    final_consensus     = CONSENSUS_ASSEMBLY.out.final_consensus

    // qc and run summary
    fastqc              = FASTQC.out.fastqc
    multiqc             = MULTIQC.out.multiqc
    duplex_seq_metrics  = FGBIO_COLLECT_DUPLEX_SEQ_METRICS.out.duplex_seq_metrics
    umi_qc              = params.umiInput ? UMI_QC.out.umi_qc_plots : channel.empty()
    run_summary         = BUILD_RUN_SUMMARY.out.run_summary
}

output {

    // core pipeline
    umi_extracted_bam   { path 'core/fgbio_extract_umis_from_bam'             }
    grouped_bam         { path 'core/picard_group_reads_by_umi/bam'           }
    grouped_hist        { path 'core/picard_group_reads_by_umi/hist'          }
    final_consensus_bam { path 'core/picard_merge_consensus_bams'             }
    
    // reference selection
    alignment_to_db     { path 'reference_prep/align_to_db_bam'               }
    samtools_coverage   { path 'reference_prep/samtools_coverage'             }
    covstats_all        { path 'reference_prep/covstats/all'                  }
    covstats_pass       { path 'reference_prep/covstats/pass'                 }
    covstats_fail       { path 'reference_prep/covstats/fail'                 }

    // variant calling | umi-grouped bam (I) vs. final consensus bam (II)
    mpileup_I           { path 'variant_calling/initial/samtools_mpileup'     }
    variants_I          { path 'variant_calling/initial/ivar_variants'        }
    vcf_I               { path 'variant_calling/initial/bcftools'             }

    mpileup_II          { path 'variant_calling/collapsed/samtools_mpileup'   }
    variants_II         { path 'variant_calling/collapsed/ivar_variants'      }
    vcf_II              { path 'variant_calling/collapsed/bcftools'           }

    // consensus genome assembly
    final_consensus     { path 'consensus_assembly/final_consensus'           }
    
    // qc and run summary
    fastqc              { path 'qc/fastqc'                                    }
    multiqc             { path 'qc/multiqc'                                   }
    duplex_seq_metrics  { path 'qc/duplex_seq_metrics'                        }
    umi_qc              { path 'qc/umi_qc_plots'                              }
    run_summary         { path '.'                                            }
}
