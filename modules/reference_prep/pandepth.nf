/* 
 * Step 2 of the reference_prep subworkflow
 *
 * Generates coverage and depth statistics with Pandepth.
 */
process PANDEPTH {

    label 'process_low'
    label 'reference_prep'

    input: 
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.chr.stat"),    emit: pandepth

    script:
    """
    pandepth \\
        -i ${bam} \\
        -o ${meta.id} \\
        -t ${task.cpus}
    
    gunzip -f ${meta.id}.chr.stat.gz
    """
}
