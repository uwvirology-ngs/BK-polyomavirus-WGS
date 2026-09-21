/*
 * Modules
 */
include { PICARD_ADDORREPLACEREADGROUPS } from '../../modules/variant_calling/addorreplacereadgroups.nf'
include { GATK_REALIGNERTARGETCREATOR   } from '../../modules/variant_calling/realignertargetcreator.nf'
include { GATK_INDELREALIGNER           } from '../../modules/variant_calling/indelrealigner.nf'
include { CDS_VARIANTS                  } from '../../modules/variant_calling/cds_variants.nf'

/*
 * Calls variants against the selected reference genome for a particular sample
 * using ivar variants after indel realignment with GATK. 
 */
workflow VARIANT_CALLING {

    take:
    input_bam_ch    // channel: [ val(meta), val(ref_info), path(ref), path(input_bam) ]
    realign_indels  // boolean: whether or not to realign indels before calling variants

    main:
    PICARD_ADDORREPLACEREADGROUPS (
        input_bam_ch
    )

    // indel realignment optional due to incompatibility with duplex consensus reads
    if (realign_indels) {
        GATK_REALIGNERTARGETCREATOR (
            PICARD_ADDORREPLACEREADGROUPS.out.rg_bam
        )

        GATK_INDELREALIGNER (
            GATK_REALIGNERTARGETCREATOR.out.intervals
        )

        // tuple val(meta), path(bam), path(bai), path(ref), path(gff), region, val(save_mpileup)
        // needs ref and gff and save_mpileup
        variants_ch = GATK_INDELREALIGNER.out.realigned_bam
            .map { meta, ref_info, ref, bam, bai -> tuple(
                meta, ref_info, ref, bam, bai,
                "${projectDir}/assets/database/${ref_info.acc}.gff", 
                Utils.getGenomicRegion(ref_info.acc)
            )}
    } else {
        variants_ch = PICARD_ADDORREPLACEREADGROUPS.out.rg_bam
            .map { meta, ref_info, ref, bam, bai -> tuple(
                meta, ref_info, ref, bam, bai,
                "${projectDir}/assets/database/${ref_info.acc}.gff", 
                Utils.getGenomicRegion(ref_info.acc)
            )}
    }

    CDS_VARIANTS (
        variants_ch
    )

    emit:
    rg_bam              = PICARD_ADDORREPLACEREADGROUPS.out.rg_bam
    intervals           = realign_indels ? GATK_REALIGNERTARGETCREATOR.out.intervals : channel.empty()
    realigned_bam       = realign_indels ? GATK_INDELREALIGNER.out.realigned_bam     : channel.empty()
    mpileup             = CDS_VARIANTS.out.mpileup
    variants            = CDS_VARIANTS.out.variants
}
