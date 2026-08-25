#!/usr/bin/env python3

import argparse
import pandas as pd

def select_reference(covstats_path, sample_id: str, min_depth: int, min_coverage: int) -> None:
    """
    DESCRIPTION
    """
    # read in covstats
    covstats: pd.DataFrame = pd.read_csv(covstats_path, sep = '\t')

    # filter for passing
    covstats = covstats[
        (covstats["coverage"] >= min_coverage) & (covstats["mean_depth"] >= min_depth)
    ]

    # sort covstats by coverage and depth
    covstats = covstats.sort_values(by = ["coverage", "mean_depth"], ascending = False)

    # build acc, tag, header_info cols
    meta = covstats["accession"].str.split(' ')
    covstats["acc"] = meta.str[0]
    covstats["tag"] = meta.str[1]
    covstats["header_info"] = meta.str[2:].str.join(' ')

    # write chosen files out 
    if not covstats.empty:
        outfile_name = sample_id + "_refs.tsv"
        refs_tsv = covstats[["acc", "tag", "header_info"]]
        refs_tsv.to_csv(outfile_name, sep = '\t', header = False, index = False)  

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("covstats")
    parser.add_argument("--sample_id", type = str)
    parser.add_argument("--min_depth", type = int)
    parser.add_argument("--min_coverage", type = int)
    args = parser.parse_args()

    select_reference(args.covstats, args.sample_id, args.min_depth, args.min_coverage)
