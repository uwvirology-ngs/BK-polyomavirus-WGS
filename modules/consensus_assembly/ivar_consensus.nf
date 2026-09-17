/*
 * Steps 1 and 3 of the consensus assembly workflow.
 *
 * Constructs a consensus genome with samtools mpileup and ivar consensus. 
 * docs: https://andersen-lab.github.io/ivar/html/manualpage.html
 */
process IVAR_CONSENSUS {

    label 'process_medium'
    container 'quay.io/biocontainers/ivar:1.4--h6b7c446_1'

    input:
    tuple val(meta), val(ref_info), path(bam), path(ref)
    tuple val(ivar_t), val(ivar_q), val(ivar_m)
    val(step)

    output:
    tuple val(meta), val(ref_info), path("*.fa"),                   emit: consensus_fa
    tuple val(meta), val(ref_info), path("*.mpileup"), path(ref),   emit: mpileup, optional: true

    script:
    """
    samtools mpileup \\
        --reference ${ref} \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 0 \\
        --min-BQ 0 \\
        -aa \\
        ${bam} \\
    | ivar consensus \\
        -t ${ivar_t} \\
        -q ${ivar_q} \\
        -m ${ivar_m} \\
        -n N \\
        -p "${meta.id}_${ref_info.acc}_${step}"

    # get rid of linebreaks except the header line
    awk '/^>/ {printf "%s\\n", \$0; next} {printf "%s", \$0} END {print ""}' ${meta.id}_${ref_info.acc}.fa > ${meta.id}_${ref_info.acc}_temp.fa
    
    # removing leading Ns
    sed '/^>/!s/^N\\+//' ${meta.id}_${ref_info.acc}_temp.fa > ${meta.id}_${ref_info.acc}_temp_frontNtrimmed.fa
    
    # remove trailing Ns
    sed '/^>/!s/N\\+\$//' ${meta.id}_${ref_info.acc}_temp_frontNtrimmed.fa > ${meta.id}_${ref_info.acc}.fa
    
    rm ${meta.id}_${ref_info.acc}_temp.fa
    rm ${meta.id}_${ref_info.acc}_temp_frontNtrimmed.fa
    """
}
