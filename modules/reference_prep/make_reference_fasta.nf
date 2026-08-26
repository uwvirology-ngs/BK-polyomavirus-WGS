/*
 * TODO
 */
process MAKE_REFERENCE_FASTA {

    label 'process_single'
    label 'reference_prep'

    input:
    tuple val(meta), val(ref_info)
    path db

    output:
    tuple val(meta), val(ref_info), path("*.fa"),   emit: ref

    script:
    """
    samtools faidx ${db} ${ref_info.acc} > "${meta.id}_${ref_info.acc}.fa"
    """
}
