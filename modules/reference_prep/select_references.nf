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
    tuple val(meta), path(samtools_coverage)
    path(db)

    output:
    tuple val(meta), path("*_covstats.tsv"),        emit: covstats
    tuple val(meta), path("*_covstats_pass.tsv"),   emit: covstats_pass
    tuple val(meta), path("*_covstats_fail.tsv"),   emit: covstats_fail
    tuple val(meta), path("*_refs.tsv"),            emit: refs_tsv, optional: true

    script:
    """
    select_references.py ${samtools_coverage} ${db} ${meta.id} \\
        --min_coverage ${params.ref_min_cov} \\
        --min_depth ${params.ref_min_depth}
    """
}
