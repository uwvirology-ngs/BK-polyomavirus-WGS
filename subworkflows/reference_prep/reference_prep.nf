/*
 * Modules
 */
include { BWA_ALIGN_TO_DB   } from '../../modules/reference_prep/bwa_align_to_db.nf'
include { PANDEPTH          } from '../../modules/reference_prep/pandepth.nf'
include { SELECT_REFERENCES } from '../../modules/reference_prep/select_references.nf'
include { SAMTOOLS_FAIDX    } from '../../modules/reference_prep/samtools_faidx.nf'

/*
 * Selects acceptable reference genomes among the provided database using 
 * coverage and depth statistics provided by Pandepth. Adds reference 
 * genome information (ref_info) as a secondary map adjacent to meta. 
 */
workflow REFERENCE_PREP {        

    take:                                                                          
    reads_ch    // channel: [ val(meta), path(reads) ]
    db          // path:    file(params.db)
                                                                                   
    main:
    BWA_ALIGN_TO_DB (
        reads_ch,
        db
    )

    PANDEPTH (
        BWA_ALIGN_TO_DB.out.alignment_to_db
    )
                                                                            
    SELECT_REFERENCES (
        PANDEPTH.out.pandepth,
        db
    )

    // use meta and ref_info maps: [ [ meta.id, ... ], [ acc, tag, description ] ]              
    ref_info_ch = SELECT_REFERENCES.out.refs_tsv
        .flatMap { meta, refs_tsv -> Utils.add_ref_info_to_meta(meta, refs_tsv) }

    SAMTOOLS_FAIDX (                                                                     
        ref_info_ch,
        db
    )

    // duplicate the reads for each selected reference genome
    output_ch = reads_ch
        .cross(SAMTOOLS_FAIDX.out.selected_refs)
        .multiMap { reads_tuple, ref_info_tuple -> 
            reads:  reads_tuple
            ref:    ref_info_tuple
        }

    emit:
    reads   = output_ch.reads   // channel: [ val(meta), path(reads) ]
    ref     = output_ch.ref     // channel: [ val(meta), val(ref_info), path(ref_fasta) ]
    
    // outfiles
    alignment_to_db = BWA_ALIGN_TO_DB.out.alignment_to_db
    pandepth        = PANDEPTH.out.pandepth
    covstats        = SELECT_REFERENCES.out.covstats
    covstats_pass   = SELECT_REFERENCES.out.covstats_pass
    covstats_fail   = SELECT_REFERENCES.out.covstats_fail
    refs_tsv        = SELECT_REFERENCES.out.refs_tsv
    selected_refs   = SAMTOOLS_FAIDX.out.selected_refs
}
