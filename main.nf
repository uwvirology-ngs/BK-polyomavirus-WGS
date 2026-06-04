#!/usr/bin/env nextflow

/*
 * Parameters
 */
params {
    samplesheet: Path
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
include { PICARD_FASTQ_TO_SAM         } from './modules/local/picard_fastq_to_sam.nf'
include { FGBIO_EXTRACT_UMIS_FROM_BAM } from './modules/local/fgbio_extract_umis_from_bam.nf'
include { PICARD_SAM_TO_FASTQ         } from './modules/local/picard_sam_to_fastq.nf'
include { BWA_ALIGN_FASTQ             } from './modules/local/bwa_align_fastq.nf'

workflow {
    
    main:
    READ_SAMPLESHEET(
        params.samplesheet
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

    // use revica-strm to select appropriate reference genomes
    REFERENCE_PREP (
        PICARD_SAM_TO_FASTQ.out.umi_extracted_fastqs,
        file(params.db),
        false
    )
    
    // join revica results to pair reads with each selected reference
    revica_ch = REFERENCE_PREP.out.reads.join(REFERENCE_PREP.out.ref)

    BWA_ALIGN_FASTQ(
        revica_ch
    )

    publish:
    unaligned_bam             = PICARD_FASTQ_TO_SAM.out.unaligned_bam
    umi_extracted_bam         = FGBIO_EXTRACT_UMIS_FROM_BAM.out.umi_extracted_bam
    umi_extracted_fastqs      = PICARD_SAM_TO_FASTQ.out.umi_extracted_fastqs
    aligned_umi_extracted_bam = BWA_ALIGN_FASTQ.out.aligned_umi_extracted_bam
}

output {
    unaligned_bam {
        path 'picard_fastq_to_sam'
    }
    umi_extracted_bam {
        path 'fgbio_extract_umis_from_bam'
    }
    umi_extracted_fastqs {
        path 'picard_sam_to_fastq'
    }
    aligned_umi_extracted_bam {
        path 'bwa_align_fastq'
    }
}