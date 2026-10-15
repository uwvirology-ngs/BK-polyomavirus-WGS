#!/usr/bin/env python3

from collections import Counter
import argparse
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", type = str)
    parser.add_argument("sample_id", type = str)
    parser.add_argument("umi_input", type = str)

    args = parser.parse_args()

    # build hash map counting occurrences of each UMI pair
    def count_umis_from_bam(bam_path):
        umi_counter = Counter()

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for record in bam:
                if record.is_unmapped:
                    continue
                try:
                    umi_pair = record.get_tag("RX")
                    umi_counter[umi_pair] += 1
                except KeyError:
                    continue

            return umi_counter

    counts = count_umis_from_bam(args.bam)

    # read in list of expected UMIs
    df = pd.read_excel(args.umi_input, usecols='A', skiprows=4, nrows=32, header=None)
    umi_list = df[0].tolist()

    # generate a matrix representing the frequency of each UMI pair
    def build_matrix(umi_counter):
        # initialize as a matrix of 0's
        matrix = pd.DataFrame(0, index=umi_list, columns=umi_list)

        # then fill with umi pair counts
        for umi_pair, count in umi_counter.items():
            if '-' in umi_pair:
                fwd, rev = umi_pair.split('-')
                if fwd in matrix.index and rev in matrix.columns:
                    matrix.loc[fwd, rev] = count

        return matrix

    matrix_df = build_matrix(counts)

    # visualize with maplotlib and seaborn
    plt.figure(figsize = (10, 8))
    sns.heatmap(matrix_df, cmap = "viridis", cbar = True)
    plt.title(f"{args.sample_id} UMI Frequencies")
    plt.xlabel("Rev UMI")
    plt.ylabel("Fwd UMI")

    plt.savefig(
        f"{args.sample_id}_matrix.png",
        dpi=300,
        bbox_inches="tight"
    )
    
    plt.close()
