/* 
 * Step 2 of the reference_prep subworkflow
 *
 * Generates coverage statistics with Samtools
 */
process SAMTOOLS_COVERAGE {

    label 'process_low'
    label 'reference_prep'

    input: 
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_coverage.tsv"),    emit: samtools_coverage

    script:
    """
    set -euo pipefail
    
    samtools sort -@ ${task.cpus} ${bam} \\
    | samtools coverage -d 0 - > ${meta.id}_coverage.tsv    
    """
}
