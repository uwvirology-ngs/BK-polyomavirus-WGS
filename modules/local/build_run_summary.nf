/*
 * Concatenates all sample-reference pair summaries into a run-level .tsv
 */
process BUILD_RUN_SUMMARY {

    label 'twist'
    label 'process_low'

    input:
    path(summary)

    output:
    path("run_summary.tsv"),    emit: run_summary

    script:
    """
    awk 'FNR==1 && NR!=1 { next } { print }' *.tsv > run_summary.tsv
    """
}