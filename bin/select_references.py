#!/usr/bin/env python3

import argparse
import pandas as pd

def build_covstats(samtools_cov: str, ref_db: str) -> pd.DataFrame:
    """
    Returns a pandas DataFrame reflecting read statistics from 
    Samtools Coverage output useful for reference selection.

    :param samtools_cov:    path to Samtools Coverage output
    :param ref_db:          path to reference genome database

    :returns:   coverage statistics as a pandas DataFrame
    """
    # read in Samtools Coverage output tsv
    covstats_df: pd.DataFrame = pd.read_csv(samtools_cov, sep='\t')

    # reformat and validate column types
    covstats_df = covstats_df.rename(columns = { "#rname": "accession" })
    covstats_df = covstats_df.astype({
        "accession": "string", 
        "startpos":  "int64",       "endpos":    "int64",
        "numreads":  "int64",       "covbases":  "int64", 
        "coverage":  "float64",     "meandepth": "float64", 
        "meanbaseq": "float64",     "meanmapq":  "float64"
    })

    # Samtools includes only the accession number in its output, so 
    # we reintroduce additional reference genome information here
    ref_info_map: dict[str, str] = dict()
    with open(ref_db, 'r') as db:
        for line in db:
            if line.startswith('>'):
                fasta_header = line[1:].strip()
                ref_acc = fasta_header.split(' ')[0]
                ref_info_map[ref_acc] = fasta_header

    covstats_df["fasta_header"] = covstats_df["accession"].map(ref_info_map)

    # include separate columns for the reference tag and description
    ref_info = covstats_df["fasta_header"].str.split(' ')
    covstats_df["tag"]          = ref_info.str[1]
    covstats_df["description"]  = ref_info.str[2:].str.join(' ')

    covstats_df = covstats_df.sort_values(
        by = ["coverage", "meandepth"], ascending = False
    )

    return covstats_df

def select_references(
        covstats: pd.DataFrame, min_coverage: int, min_depth: int
    ) -> tuple:
    """
    Selects acceptable reference genomes as determined by coverage and depth 
    statistics from Pandepth, then returns a pandas DataFrame with their
    information (accession, tag, description)

    :param  covstats:       formatted coverage and depth stats from Pandepth
    :param  min_coverage:   minimum acceptable coverage for selection
    :param  min_depth:      minimum acceptable depth for selection

    :returns:   A tuple of pandas DataFrames including covstats for passing vs
                failing references and a simple tsv with the accesion, tag, and 
                description of each acceptable reference genome.
    """
    # select reference genomes by minimum coverage and mean_depth per Pandepth
    filter = (covstats["coverage"] >= min_coverage) & (covstats["meandepth"] >= min_depth)
    covstats_pass = covstats.loc[ filter].copy()
    covstats_fail = covstats.loc[~filter].copy()

    # build simple tsv for shuttling accepted references to downstream processes
    refs_tsv = covstats_pass[["accession", "tag", "description"]].copy()

    return covstats_pass, covstats_fail, refs_tsv

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pandepth", type = str)
    parser.add_argument("database", type = str)
    parser.add_argument("sample_id", type = str)
    parser.add_argument("--min_coverage", type = int)
    parser.add_argument("--min_depth", type = int)
    
    args = parser.parse_args()

    # build covstats from pandepth output and write to disk
    covstats: pd.DataFrame = build_covstats(args.pandepth, args.database)
    covstats.to_csv(args.sample_id + "_covstats.tsv", sep = '\t', index = False)

    # select acceptable reference genomes and write their information to disk
    covstats_pass, covstats_fail, refs_tsv = select_references(
        covstats, args.min_coverage, args.min_depth
    )
    covstats_pass.to_csv(
        args.sample_id + "_covstats_pass.tsv", sep = '\t', header = False, index = False
    )
    covstats_fail.to_csv(
        args.sample_id + "_covstats_fail.tsv", sep = '\t', header = False, index = False
    )
    if not refs_tsv.empty:
        refs_tsv.to_csv(
            args.sample_id + "_refs.tsv", sep = '\t', header = False, index = False
        )
