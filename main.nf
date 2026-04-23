#!/usr/bin/env nextflow

/*
 * Subworkflows
 */
include { READ_SAMPLESHEET        } from './subworkflows/local/read_samplesheet.nf'
include { FASTQ_TRIM_FASTP_FASTQC } from './subworkflows/nf-core/fastq_trim_fastp_fastqc.nf'

/*
 * Modules
 */
include { PICARD } from './modules/local/picard.nf'
include { FGBIO  } from './modules/local/fgbio.nf'

workflow {
    
    main:
    READ_SAMPLESHEET()

    FASTQ_TRIM_FASTP_FASTQC (
        READ_SAMPLESHEET.out.reads,
        params.save_trimmed_fail,
        params.discard_trimmed_pass,
        params.save_merged,
        params.skip_fastp,
        params.skip_fastqc
    )

    PICARD (
        FASTQ_TRIM_FASTP_FASTQC.out.reads
    )

    FGBIO (
        PICARD.out.unaligned_bam
    )

    publish:
    unaligned_bam       = PICARD.out.unaligned_bam
    umi_extracted_bam   = FGBIO.out.umi_extracted_bam
}

output {
    unaligned_bam {
        path 'picard'
    }
    umi_extracted_bam {
        path 'fgbio'
    }
}