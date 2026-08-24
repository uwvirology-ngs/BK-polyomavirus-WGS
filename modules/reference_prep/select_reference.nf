/*
 * Step 4 of the reference_prep subworkflow
 * 
 * Selects reference genomes for which alignment against produces acceptable 
 * coverage and depth statistics by Pandepth.
 */
process SELECT_REFERENCE {

    label 'process_single'
    label 'reference_prep'
    
    input:
    tuple val(meta), path(covstats)

    output:
    tuple val(meta), path("*_refs.tsv"),            optional: true, emit: refs_tsv
    tuple val(meta), path("*_failed_assembly.tsv"), optional: true, emit: failed_assembly_summary

    script:
    """
    select_reference.py \\
        ${covstats} \\
        --sample_id ${meta.id} \\
        --min_coverage ${params.ref_min_cov} \\
        --min_depth ${params.ref_min_depth}
    """
}
