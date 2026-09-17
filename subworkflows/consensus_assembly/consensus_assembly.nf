/*
 * Modules
 */
include { IVAR_CONSENSUS as IVAR_CONSENSUS_INI } from '../../modules/consensus_assembly/ivar_consensus.nf'
include { BWA_ALIGN_TO_CONSENSUS }               from '../../modules/consensus_assembly/bwa_align_to_consensus.nf'
include { IVAR_CONSENSUS as IVAR_CONSENSUS_FIN } from '../../modules/consensus_assembly/ivar_consensus.nf'

/*
 * Constructs a consensus genome sequence. Reads are realigned to 
 * the initial consensus assembly to counteract reference bias. 
 */
workflow CONSENSUS_ASSEMBLY {

    take:
    input_ch    // channel: [ val(meta), path(bam), path(ref), path(ref_info), path(reads) ]

    main:
    input_ch
        .multiMap { meta, bam, ref, ref_info, reads ->
            bam:    tuple(meta, ref_info, bam)
            ref:    tuple(meta, ref_info, ref)
            reads:  tuple(meta, ref_info, reads)
        }
        .set { input }

    initial_consensus_ch = input.bam
        .join(input.ref, by: [0, 1])

    IVAR_CONSENSUS_INI (
        initial_consensus_ch,
        tuple(
            params.ivar_consensus_ini_t, 
            params.ivar_consensus_ini_q, 
            params.ivar_consensus_ini_m
        ),
        "initial"
    )

    bwa_ch = input.reads
        .join(IVAR_CONSENSUS_INI.out.consensus_fa, by: [0, 1])

    BWA_ALIGN_TO_CONSENSUS (
        bwa_ch
    )

    final_consensus_ch = BWA_ALIGN_TO_CONSENSUS.out.realigned_bam
        .join(input.ref, by: [0, 1])

    
    IVAR_CONSENSUS_FIN (
        final_consensus_ch,
        tuple(
            params.ivar_consensus_fin_t, 
            params.ivar_consensus_fin_q, 
            params.ivar_consensus_fin_m
        ),
        "final"
    )

    emit:
    final_consensus = IVAR_CONSENSUS_FIN.out.consensus_fa
}
