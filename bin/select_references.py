#!/usr/bin/env python3

import argparse
import pandas as pd

def reformat_pandepth(pandepth: str, ref_db: str) -> pd.DataFrame:
    """
    Returns a pandas DataFrame reflecting Pandepth stats 
    as well as additonal reference genome information.

    :param pandepth:    filepath to Pandepth output
    :param ref_db:      filepath to reference genome database

    :returns:       pandas DataFrame reflecting Pandepth stats
    """
    # read in Pandepth tsv, removing metadata rows from the bottom
    covstats_df = pd.read_csv(pandepth, sep='\t')
    covstats_df = covstats_df[~covstats_df["#Chr"].str.startswith('##')]

    # reformat and validate Pandepth columns
    covstats_df = covstats_df.rename(columns = {
        "#Chr": "accession", "Length": "length",
        "CoveredSite": "covered_site", "TotalDepth": "total_depth",
        "Coverage(%)": "coverage", "MeanDepth": "mean_depth"
    })
    covstats_df = covstats_df.astype({
        "accession": "string", "length": "int64",
        "covered_site": "int64", "total_depth": "int64",
        "coverage": "float64", "mean_depth": "float64",
    })

    # Pandepth includes only the accession number in its output, so we 
    # add back the rest of the reference genome information here
    ref_info_map: dict[str, str] = dict()
    with open(ref_db, 'r') as db:
        for line in db:
            if line.startswith('>'):
                fasta_header = line[1:].strip()
                ref_acc = fasta_header.split(' ')[0]
                ref_info_map[ref_acc] = fasta_header

    covstats_df["fasta_header"] = covstats_df["accession"].map(ref_info_map)

    # build columns for the reference genome tag and description
    ref_info = covstats_df["fasta_header"].str.split(' ')
    covstats_df["tag"] = ref_info.str[1]
    covstats_df["description"] = ref_info.str[2:].str.join(' ')

    covstats_df = covstats_df.sort_values(
        by = ["coverage", "mean_depth"], ascending = False
    )

    return covstats_df

def select_references(
        covstats: pd.DataFrame, min_coverage: int, min_depth: int
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
    filter = (covstats["coverage"] >= min_coverage) & (covstats["mean_depth"] >= min_depth)
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
    covstats: pd.DataFrame = reformat_pandepth(args.pandepth, args.database)
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
