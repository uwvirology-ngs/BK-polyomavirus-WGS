/*
 * Step 1 of the bcftools variant_calling subworkflow
 *
 * Calls variants on the realigned bam file with bfctools
 */
process BCFTOOLS {

    label 'process_medium'
    container 'community.wave.seqera.io/library/bcftools:1.25--642ca8fc61ec2de8'

    input:
    tuple val(meta), val(ref_info), path(ref), path(input_bam)

    output:
    tuple val(meta), val(ref_info), path("*.vcf"),  emit: vcf 

    script:
    """
    bcftools mpileup -f ${ref} ${input_bam} | bcftools call -mv -Ov -o "${ref.baseName}_calls.vcf"
    """
}
