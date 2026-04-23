#!/usr/bin/env nextflow

/*
 * Subworkflows
 */
include { READ_SAMPLESHEET } from './subworkflows/read_samplesheet.nf'

/*
 * Modules
 */
include { PICARD } from './modules/picard.nf'
include { FGBIO  } from './modules/fgbio.nf'

workflow {
    
    main:
    READ_SAMPLESHEET()

    PICARD (
        READ_SAMPLESHEET.out.reads
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