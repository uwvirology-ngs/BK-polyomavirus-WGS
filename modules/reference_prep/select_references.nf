/*
 * Step 3 of the reference_prep subworkflow
 * 
 * Selects reference genomes for which alignment against produces acceptable 
 * coverage and depth statistics by Pandepth.
 */
process SELECT_REFERENCES {

    label 'process_low'
    label 'reference_prep'
    
    input:
    tuple val(meta), path(pandepth)
    path(db)

    output:
    tuple val(meta), path("*covstats.tsv"),     emit: covstats
    tuple val(meta), path("*_refs.tsv"),        emit: refs_tsv, optional: true

    script:
    """
    select_references.py \\
        ${pandepth} ${db} ${meta.id} \\
        --min_coverage ${params.ref_min_cov} \\
        --min_depth ${params.ref_min_depth}
    """
}
