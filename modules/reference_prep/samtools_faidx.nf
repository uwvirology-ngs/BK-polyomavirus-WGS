/*
 * Step 4 of the reference_prep subworkflow
 *
 * Emits a copy of a particular reference genome for 
 * any set of reads assigned to it. 
 */
process SAMTOOLS_FAIDX {

    label 'process_single'
    label 'reference_prep'

    input:
    tuple val(meta), val(ref_info)
    path db

    output:
    tuple val(meta), val(ref_info), path("*.fa"),   emit: selected_refs

    script:
    """
    samtools faidx ${db} ${ref_info.acc} > "${meta.id}_${ref_info.acc}.fa"
    """
}
