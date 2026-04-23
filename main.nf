#!/usr/bin/env nextflow

/*
 * Modules
 */
include { PICARD } from './modules/picard.nf'
include { FGBIO  } from './modules/fgbio.nf'

/*
 * Parameters
 */
params {
    samplesheet: Path
}

workflow {
    main:
    read_ch = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)] }
        .view()

    PICARD(read_ch)

    FGBIO(PICARD.out.unaligned_bam)

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