/*
 * Step 4 of the variant_calling subworkflow
 *
 * Calls variants on the realigned bam file with ivar variants.
 */
process CDS_VARIANTS {

    label 'process_medium'
    container 'quay.io/jefffurlong/ivar:1.4.4'

    input:
    tuple val(meta), val(ref_info), path(ref), path(bam), path(bai), path(gff), val(genomic_region)

    output:
    tuple val(meta), path("*.mpileup"), emit: mpileup
    tuple val(meta), path("*.tsv"),     emit: variants

    script:
    """
    samtools mpileup -B ${bam} -r ${genomic_region} -f ${ref} -Q 20 \\
        --max-depth 0 \\
        --count-orphans \\
        --disable-overlap-removal \\
        > "${ref.baseName}.mpileup"

    cat "${ref.baseName}.mpileup" | ivar variants \\
        -q ${params.ivar_variants_q} \\
        -t ${params.ivar_variants_t} \\
        -m ${params.ivar_variants_m} \\
        -G \\
        -r ${ref} \\
        -g ${gff} \\
        -p ${ref.baseName}
    """
}
