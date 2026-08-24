#!/usr/bin/env python3

import argparse

def select_reference(covstats_path, sample_id: str, min_depth: int, min_coverage: int) -> None:
    """
    DESCRIPTION
    """
    # read in covstats
    covstats = open(covstats_path, 'r').readlines()
    header = covstats[0]
    records = covstats[1:]
    
    # find covstats records which meet coverage and depth criteria
    passing_records = []
    records = [record.strip() for record in records if len(record) > 0]
    for record in records:
        reference_info = record.split('\t')
        coverage = float(reference_info[4])
        depth = float(reference_info[5])
        if coverage >= min_coverage and depth > min_depth:
            passing_records.append(reference_info)

    # sort by coverage and depth
    sorted_by_coverage_depth = sorted(passing_records, key = lambda x: (x[4], x[5]), reverse = True)
    init_ref_candidates = [entry[0] for entry in sorted_by_coverage_depth]

    # map reference tags to accession and header info
    init_ref_header: dict[str, list] = {}
    for candidate in init_ref_candidates: 
        meta = candidate.split(' ')
        acc = meta[0]   # accession
        tag = meta[1]   # header tag
        info = ' '.join(meta[2:])   # get header info
        if not tag in init_ref_header: 
            init_ref_header[tag] = [acc, info]

    # if reference(s) selected, output relevant info to file
    if init_ref_header:
        output_file_name = sample_id + '_refs.tsv'
        output_file = open(output_file_name, 'a+')
        for i in init_ref_header:
            # output format: reference accession <tab> reference header tag <tab> reference header info
            output_file.write(str(init_ref_header[i][0]) + "\t" + str(i) + "\t" + str(init_ref_header[i][1]) + "\n")
        output_file.close()        

    # no reference selected, output reference with highest covered percent and read distribution info
    else:
        output_file_name = sample_id + '_failed_assembly.tsv'
        output_file = open(output_file_name, 'w')

        output_text = str(sorted(records, key=lambda line: float(line.split('\t')[4]), reverse=True)[0])
        output_file.write(header)
        output_file.write('\n')
        output_file.write(output_text)
        output_file.write('\n')
        output_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("covstats")
    parser.add_argument("--sample_id", type = str)
    parser.add_argument("--min_depth", type = int)
    parser.add_argument("--min_coverage", type = int)
    args = parser.parse_args()

    select_reference(args.covstats, args.sample_id, args.min_depth, args.min_coverage)
