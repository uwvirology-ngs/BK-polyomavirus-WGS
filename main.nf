#!/usr/bin/env nextflow

/*
 * Parameters
 */
params {
    samplesheet: Path
}

/*
 * Subworkflows
 */
include { READ_SAMPLESHEET } from './subworkflows/local/read_samplesheet.nf'

/*
 * Modules
 */
include { PICARD } from './modules/local/picard.nf'
include { FGBIO  } from './modules/local/fgbio.nf'

workflow {
    
    main:
    READ_SAMPLESHEET(
        params.samplesheet
    )

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