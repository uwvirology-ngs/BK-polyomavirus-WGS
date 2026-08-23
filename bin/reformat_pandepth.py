#!/usr/bin/env python3

import argparse
import pandas as pd

def clean_pandepth(infile: str, ref_db: str) -> pd.DataFrame:
    """
    Returns a pandas DataFrame reflecting Pandepth stats 
    as well as additonal reference genome information.

    :param infile:  filepath to Pandepth output
    :param ref_db:  filepath to reference genome database

    :returns:       pandas DataFrame reflecting Pandepth stats
    """
    # map reference genome accessions to full reference information
    ref_map: dict[str, str] = {}
    with open(ref_db, 'r') as db:
        for line in db:
            if line.startswith('>'):
                ref_name = line[1:].strip('\n')
                ref_acc = ref_name.split(' ')[0]
                ref_map[ref_acc] = ref_name

    # reformat pandepth output, add comprehensive reference genome information
    df_pandepth = pd.read_csv(infile, sep='\t')
    df_pandepth = df_pandepth.rename(columns = {
        "#Chr": "accession",
        "Length": "length",
        "CoveredSite": "covered_site",
        "TotalDepth": "total_depth",
        "Coverage(%)": "coverage",
        "MeanDepth": "mean_depth"
    })
    df_pandepth = df_pandepth[~df_pandepth["accession"].str.startswith('##')]
    df_pandepth = df_pandepth.astype({
        "accession": "string",
        "length": "int64",
        "covered_site": "int64",
        "total_depth": "int64",
        "coverage": "float64",
        "mean_depth": "float64",
    })
    df_pandepth["accession"] = df_pandepth["accession"].map(ref_map)

    df_pandepth = df_pandepth.sort_values(by = "coverage", ascending = False)

    return df_pandepth

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    arg_names: list[str] = ["infile", "ref", "outfile"]
    for a in arg_names:
        parser.add_argument(a)
    args = parser.parse_args()

    pandepth: pd.DataFrame = clean_pandepth(args.infile, args.ref)

    # write reformatted pandepth stats directly to the workdir
    pandepth.to_csv(args.outfile, sep = '\t', index = False)
