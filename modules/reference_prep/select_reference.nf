process SELECT_REFERENCE {
    tag "$meta.id"
    label 'process_single'
    label 'reference_prep'
    
    input:
    tuple val(meta), path(bbmap_db_covstats)

    output:
    tuple val(meta), path("*_refs.tsv"),            optional: true, emit: refs_tsv
    tuple val(meta), path("*_failed_assembly.tsv"), optional: true, emit: failed_assembly_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    select_reference.py \\
        -bbmap_covstats ${bbmap_db_covstats} \\
        -b ${prefix} \\
        -m ${params.ref_min_depth} \\
        -p ${params.ref_min_cov}
    """
}
