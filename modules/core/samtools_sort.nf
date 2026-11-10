/*
 * Resorts the UMI-grouped BAM file into coordinate-sorted order.
 */
process SAMTOOLS_SORT {

    label 'process_medium'
    label 'twist'

    input:
    tuple val(meta), val(ref_info), path(ref), path(grouped_bam)

    output:
    tuple val(meta), val(ref_info), path(ref), path("*.bam"),   emit: grouped_coord_sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} ${grouped_bam} -o "${ref.baseName}_grouped_coord_sorted.bam"
    """
}
