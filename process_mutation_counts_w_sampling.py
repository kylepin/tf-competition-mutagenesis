"""
# Overview
To run: `python process_mutation_counts.py`

## Important note:
`process_mutation_counts.py` must be run before running this script. There are
files generated in that script that are required by this script.

Input files from external sources:
- hg19.fa
- new_refseq_exons_171007.bed
    file of exons to exclude
- facilitated_myc_mutations.tsv (already included)

Input files from other scripts from this repository:
- *Many*. `process_mutation_counts.py` should be run before running
    this file.

Output files:
- myc_control_sample_{i}_mutation_counts.tsv for i in 1,2..10
- myc_control_sample{i}_counts_table.tsv for i in 1,2..10
"""
import pandas as pd
from pybedtools import BedTool

from process_mutation_counts import (compute_mutation_counts, get_total_core_match_counts,
                                     reverse_complement)

INPUT_PATH = "./input"
OUTPUT_PATH = "./intermediate_and_output"

def compute_control_rates(myc_site_mutations_df, sample_num):
    """This function serves to compute control rates (output on the
        console) and create a table that can easily be transformed to be
        similar to that found in Table S3.

        This differs from the function of the same name in
        process_mutation_counts.py in that we do not need to adjust
        because we have ensured the number of control sites of each
        suffix match the corresponding number of MYC binding sites.
    """
    ctrl = pd.read_csv(
        f"{OUTPUT_PATH}/myc_control_sample_{sample_num}_mutation_counts.tsv",
        sep="\t")
    df = myc_site_mutations_df.merge(
        ctrl, on=["MSI", "suffix", "WT", "mutation", "mutation_pos"])

    df = df.rename(columns={"count_x": "myc_mutation_count", "count_y": "myc_match_count",
                            "count": "ctrl_mutation_count"})

    df = df.merge(pd.read_csv(
        f"{OUTPUT_PATH}/myc_control_sample_{sample_num}_match_counts.tsv",
        sep="\t"), on="suffix")
    df = df.rename(columns={"count": "ctrl_match_count"})
    df = df.sort_values(by=["MSI", "facilitated"])
    # prior to computing control rate, output the table so it can be transformed to Table S3
    df.to_csv(f"{OUTPUT_PATH}/myc_control_sample_{sample_num}_mutation_counts_table.tsv", sep="\t",
              index=False, header=True)

    # count cores regardless of nn
    df = df.groupby(["MSI", "facilitated", "suffix", "myc_match_count"]).agg(
        {"ctrl_mutation_count": "sum"}).reset_index()

    # count total control and MYC binding site
    # mutations in each quadrant ({MSS, MSI}x{facilitated, neutral})
    ctrl_rate_df = df.groupby(["MSI", "facilitated"]).agg(
        {"ctrl_mutation_count": "sum", "myc_match_count": "sum"})

    ctrl_rate_df = ctrl_rate_df["ctrl_mutation_count"] / ctrl_rate_df[
        "myc_match_count"]

    print(ctrl_rate_df)
    return ctrl_rate_df


def main():
    # intersect MYC binding site mutation counts with facilitated/neutral annotation
    # and the MYC binding site match counts
    mutation_counts = pd.read_csv(f"{OUTPUT_PATH}/myc_mutation_counts.tsv", sep="\t")
    facilitated_df = pd.read_csv(f"{INPUT_PATH}/facilitated_myc_mutations.tsv",
                                 sep="\t")
    myc_matches = pd.read_csv(f"{OUTPUT_PATH}/myc_match_counts.tsv", sep="\t")
    myc_mutation_df = mutation_counts.merge(
        facilitated_df, on=["suffix", "WT", "mutation", "mutation_pos"]).merge(
        myc_matches, on="suffix").groupby(
        ["MSI", "facilitated", "suffix", "WT", "mutation", "mutation_pos"]).sum().reset_index()

    # chunked to avoid putting it all in memory simultaneously..
    # create new files to hold all the control sites with matching suffixes (suffices?)
    i = 0
    for chunk in pd.read_csv(f"{OUTPUT_PATH}/myc_control_matches_no_exons.bed",
                             sep="\t", header=None, chunksize=1_000_000):
        chunk["upper"] = chunk[6].str.upper()  # column 6 is the pattern match
        chunk["suffix"] = chunk["upper"].str[2:]
        # column 5 is the stranded (+ or 0)
        chunk.loc[chunk[5] == "-", "suffix"] = chunk["upper"].apply(reverse_complement).str[2:]
        for suffix in ["CGCG", "CGAG", "CGTG", "CATG"]:
            chunk[chunk["suffix"] == suffix][[0, 1, 2, 3, 4, 5, 6]].to_csv(
                f"{OUTPUT_PATH}/myc_control_matches_nn{suffix}_no_exons.bed", sep="\t", index=False,
                header=None, mode="a")
        i += 1_000_000
        print(i)

    suffix_and_counts = list(pd.read_csv(f"{OUTPUT_PATH}/myc_match_counts.tsv",
                                         sep="\t").itertuples(index=False, name=None))

    # randomly generate 10 samples with the same numbers of each suffix/core as from the MYC binding
    # site set
    for sample_num in range(1, 11):
        print(f"Processing sample # {sample_num}")
        for suffix, myc_match_count in suffix_and_counts:
            print(suffix)
            print(myc_match_count)
            print(29+sample_num)
            f = f"{OUTPUT_PATH}/myc_control_sample_{sample_num}_matches_nn{suffix}_sampled.bed"
            BedTool(f"{OUTPUT_PATH}/myc_control_matches_nn{suffix}_no_exons.bed").sample(
                n=myc_match_count, seed=29+sample_num,
                output=f)

        # concatenate the sampled files for each of the 4 suffixes
        f = f"{OUTPUT_PATH}/myc_control_sample_{sample_num}_matches_no_exons.bed"
        with open(f, "w") as fout:
            for suffix in ["CATG", "CGCG", "CGAG", "CGTG"]:
                f = f"{OUTPUT_PATH}/myc_control_sample_{sample_num}_matches_nn{suffix}_sampled.bed"
                with open(f, "r") as fin:
                    fout.write(fin.read())

        # then do all the mutation count analysis using the sample control site set
        myc_control_matches = BedTool(
            f"{OUTPUT_PATH}/myc_control_sample_{sample_num}_matches_no_exons.bed")
        myc_control_matches.intersect(
            BedTool(f"{OUTPUT_PATH}/mutations.bed"), wa=True, wb=True,
            output=f"{OUTPUT_PATH}/myc_control_sample_{sample_num}_matches_mutations.bed")

        compute_mutation_counts(f"myc_control_sample_{sample_num}")
        get_total_core_match_counts(f"myc_control_sample_{sample_num}")

        # now compute the CTRL rates for MSS/MSI Myc-neutral/-facilitated
        compute_control_rates(myc_mutation_df, sample_num)


if __name__ == "__main__":
    main()
