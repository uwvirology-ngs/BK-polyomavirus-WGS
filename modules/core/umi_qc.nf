/*
 * Generates visualizations to describe the sequencing of 
 * unique molecular identifiers. 
 */
process UMI_QC {

    container 'community.wave.seqera.io/library/pysam_matplotlib_openpyxl_pandas_seaborn:7f62bde10f25a4f0'
    label 'process_low'

    input:
    tuple val(meta), val(ref_info), path(ref), path(bam)
    path(umi_input)

    output:
    tuple val(meta), path("*.png"),     emit: umi_qc_plots

    script:
    """
    umi_qc.py ${bam} ${meta.id} ${umi_input}
    """
}
