/* 
 * Step 3 of the reference_prep subworkflow
 *
 * Reformats pandepth output, adding the full reference genome headers.
 */
process REFORMAT_PANDEPTH {

    label 'process_low'
    container 'quay.io/epil02/revica-strm:0.0.5'

    input: 
    tuple val(meta), path(pandepth)
    path(db)

    output:
    tuple val(meta), path("*covstats.tsv"),     emit: covstats

    script:
    """
    reformat_pandepth.py ${pandepth} ${db} ${meta.id}_covstats.tsv
    """
}