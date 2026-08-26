/*
 * Modules
 */
include { BWA_MEM_ALIGN_DB     } from '../../modules/reference_prep/bwa_mem_align_db.nf'
include { PANDEPTH             } from '../../modules/reference_prep/pandepth.nf'
include { SELECT_REFERENCES    } from '../../modules/reference_prep/select_references.nf'                     
include { MAKE_REFERENCE_FASTA } from '../../modules/reference_prep/make_reference_fasta.nf'

/*
 * Selects the best reference among the provided database using coverage 
 * and depth statistics provided by Pandepth. Adds reference genome 
 * information (ref_info) as a secondary map adjacent to meta.
 */
workflow REFERENCE_PREP {        

    take:                                                                          
    reads_ch    // channel: [ val(meta), path(reads) ]
    db          // path:    file(params.db)
                                                                                   
    main:
    BWA_MEM_ALIGN_DB (
        reads_ch,
        db
    )

    PANDEPTH (
        BWA_MEM_ALIGN_DB.out.alignment_to_db
    )
                                                                            
    SELECT_REFERENCES (
        PANDEPTH.out.pandepth,
        db
    )

    // use meta and ref_info maps: [ [ meta.id, ... ], [ acc, tag, description ] ]              
    ref_info_ch = SELECT_REFERENCES.out.refs_tsv
        .flatMap { meta, refs_tsv -> Utils.add_ref_info_to_meta(meta, refs_tsv) }

    MAKE_REFERENCE_FASTA (                                                                     
        ref_info_ch,
        db
    )

    // duplicate the reads for each selected reference genome
    output_ch = reads_ch
        .cross(MAKE_REFERENCE_FASTA.out.ref)
        .multiMap { reads_tuple, ref_info_tuple -> 
            reads:  reads_tuple
            ref:    ref_info_tuple
        }

    emit:
    reads   = output_ch.reads   // channel: [ val(meta), path(reads) ]
    ref     = output_ch.ref     // channel: [ val(meta), val(ref_info), path(ref_fasta) ]
    
    alignment_to_db = BWA_MEM_ALIGN_DB.out.alignment_to_db
    covstats = SELECT_REFERENCES.out.covstats
    failed_assembly_summary = SELECT_REFERENCES.out.failed_assembly_summary
    reference_fasta = MAKE_REFERENCE_FASTA.out.ref
    pandepth = PANDEPTH.out.pandepth
    refs_tsv = SELECT_REFERENCES.out.refs_tsv
}
