"""
# Overview

Input files from external sources:
- hg19.fa
- new_refseq_exons_171007.bed
    file of exons to exclude
- TableS1.txt
    patient data for the mutation data (includes needed info on MSI/MSS status)
- somatic_mutations.txt
    mutation dataset
- facilitated_myc_mutations.tsv (already included)

Output files:
- myc_control_sample_{i}_nnCGTG_methylation_values.bed for i in 1,2..10
- myc_control_methylation_values.bed
- myc_methylation_values.bed

Also outputs a chart and MannWhitney U test results in the notebook
"""

import pandas as pd
from pybedtools import BedTool


def reverse_complement(sequence):
    complements = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
    }

    return "".join(complements[nuc] for nuc in reversed(sequence))


def get_mutations_file():
    """Load the mutation dataset and join with the patient data
    for the MSS/MSI annotations.
    """

    # We begin by setting up data on mutations annotated with MSI vs MSS
    patient_df = pd.read_csv("TableS1.txt", delimiter="\t")
    patient_df = patient_df[patient_df["Study cohort"] == "WGS"]
    mutations_df = pd.read_csv("somatic_mutations.txt", delimiter="\t")
    # include only single nucleotide mutations (df["mutation_type"] == "S"])
    mutations_df = mutations_df[mutations_df["mutation_type"] == "S"]

    # merge with patient_df table so that we can get the MSS/MSI annotation data
    # associated with each patient
    mutations_df = mutations_df.merge(patient_df, left_on="#sample", right_on="Sample code")
    mutations_df["sample"] = mutations_df["#sample"]
    mutations_df.drop(columns=["#sample"], inplace=True)

    # generate a file for use downstream
    mutations_df.to_csv("mutations.bed", sep="\t", index=False, header=False)


def compute_mutation_counts(prefix, positions=(4, 5), save_to_file=True, insert_missing=True):
    """Compute counts for mutation identities at certain positions,
    for particular MSS/MSI status.

    Result: "{prefix}_mutation_counts.tsv"
    """
    # load the dataframe of mutations within non-exon motif matches
    df = pd.read_csv(f"{prefix}_matches_mutations.bed", delimiter="\t", header=None)

    df.rename(
        columns={0: "chrom", 1: "motif_start", 2: "motif_end", 3: "name", 4: "length", 5: "strand",
                 6: "motif", 8: "mutation_coord",
                 10: "WT", 11: "mutation", 24: "MSI"}, inplace=True)

    # determine the length of the motif
    motif_length = df["length"].unique().item()

    df = df[["chrom", "motif_start", "strand", "motif", "mutation_coord", "WT", "mutation", "MSI"]]

    # identify the position of the mutation within the motif
    df["mutation_pos"] = df["mutation_coord"] - df["motif_start"]
    # account for reverse position counting if motif was on backward strand
    df.loc[df["strand"] == "-", "mutation_pos"] = motif_length + 1 - df["mutation_pos"]

    # we are only interested in mutations at certain positions within the motif
    df = df[df["mutation_pos"].isin(positions)]

    # WT is the original nucleotide at the position
    df.loc[df["strand"] == "-", "WT"] = df["WT"].apply(reverse_complement)
    # mutation is the after-mutation nucleotide at the position
    df.loc[df["strand"] == "-", "mutation"] = df["mutation"].apply(reverse_complement)
    # convert all motifs to upper-case
    df.loc[:, "motif"] = df["motif"].str.upper()
    # reverse complement of the motif if it is on the backward strand
    df.loc[df["strand"] == "-", "motif"] = df["motif"].apply(reverse_complement)
    # this serves to ensure that detected motifs match the WT nucleotide specified
    df["matches"] = df.apply(lambda row: row["motif"][row["mutation_pos"] - 1] == row["WT"], axis=1)
    df = df[df["matches"]]

    df = df.groupby(["MSI", "motif", "WT", "mutation", "mutation_pos"]).count()[
        "chrom"].reset_index()
    df = df.rename(columns={"chrom": "count"})

    # fill in zeroes for any columns for which there were no mutations
    if insert_missing:
        for msi in (0, 1):
            for motif in df["motif"].unique().tolist():
                for pos in positions:
                    wt = motif[pos - 1]
                    # print(wt)
                    for mutant_allele in list({"C", "A", "G", "T"} - {wt}):
                        if not (((df["MSI"] == msi) & (df["motif"] == motif) & (df["WT"] == wt)
                                 & (df["mutation"] == mutant_allele)
                                 & df["mutation_pos"] == pos).any()):
                            df.loc[len(df)] = [msi, motif, wt, mutant_allele, pos, 0]

        df = df.sort_values(by=["MSI", "motif", "WT", "mutation", "mutation_pos"])
        df = df.groupby(["MSI", "motif", "WT", "mutation", "mutation_pos"]).sum().reset_index()

    df["suffix"] = df["motif"].str[2:]

    # for control sites to this point, all the nn values will have separate counts
    # we must aggregate by the suffix (e.g. CGTG for nnCGTG).
    df = df.groupby(["MSI", "suffix", "WT", "mutation", "mutation_pos"])[
        "count"].sum().reset_index()

    if save_to_file:
        df.to_csv(f"{prefix}_mutation_counts.tsv", sep="\t", index=False, header=True)
    else:
        return df


def get_total_core_match_counts(prefix, save_to_file=True):
    """Retrieve total number of counts of the different MYC binding
    or control site patterns. For control sites this requires
    aggregating counts for the different values of nn for the same
    suffix.
    """
    df = pd.read_csv(f"{prefix}_matches_no_exons.bed", delimiter="\t", header=None)
    df[6] = df[6].str.upper()

    # we specifically do not remove non-standard chromosomes

    df.loc[df[5] == "-", 6] = df[6].apply(reverse_complement)

    df = df.groupby(6).count()[0].reset_index()
    df = df.rename(columns={6: "core", 0: "count"})

    df["suffix"] = df["core"].str[2:]
    df = df.groupby("suffix")["count"].sum().reset_index()

    if save_to_file:
        df.to_csv(f"{prefix}_match_counts.tsv", sep="\t", index=False, header=True)

    else:
        return df


def compute_control_rates(myc_site_mutations_df):
    """This function serves to compute control rates (output on the
    console) and create a table that can easily be transformed to be
    simliar to that found in Table S3.
    """
    ctrl = pd.read_csv(f"myc_control_mutation_counts.tsv", sep="\t")
    df = myc_site_mutations_df.merge(
        ctrl, on=["MSI", "suffix", "WT", "mutation", "mutation_pos"])

    df = df.rename(columns={"count_x": "myc_mutation_count", "count_y": "myc_match_count",
                            "count": "ctrl_mutation_count"})

    df = df.merge(pd.read_csv(f"myc_control_match_counts.tsv", sep="\t"),
                  on="suffix")
    df = df.rename(columns={"count": "ctrl_match_count"})
    df = df.sort_values(by=["MSI", "facilitated"])
    # prior to computing control rate, output the table so it can be transformed to Table S3
    df.to_csv(f"myc_control_mutation_counts_table.tsv", sep="\t", index=False, header=True)

    # multiply control site mutation counts by ratio of control sites to MYC binding sites
    # (per suffix)
    df["normalized_ctrl_mutation_count"] = (df["ctrl_mutation_count"] / df["ctrl_match_count"]
                                            * df["myc_match_count"])

    # count cores regardless of nn
    df = df.groupby(["MSI", "facilitated", "suffix", "myc_match_count"]).agg(
        {"normalized_ctrl_mutation_count": "sum"}).reset_index()

    # count total control and MYC binding site
    # mutations in each quadrant ({MSS, MSI}x{facilitated, neutral})
    ctrl_rate_df = df.groupby(["MSI", "facilitated"]).agg(
        {"normalized_ctrl_mutation_count": "sum", "myc_match_count": "sum"})

    ctrl_rate_df = ctrl_rate_df["normalized_ctrl_mutation_count"] / ctrl_rate_df["myc_match_count"]

    print(ctrl_rate_df)
    return ctrl_rate_df


def main():
    # put the mutations file in place
    get_mutations_file()

    # for both MYC binding sites and control sites, exclude exons
    for prefix in ("myc", "myc_control"):
        print(prefix)
        # intersect with exons
        BedTool(f"{prefix}_matches.bed").intersect(
            "new_refseq_exons_171007.bed", v=True, output=f"{prefix}_matches_no_exons.bed")

    # intersect MYC binding sites with mutations
    BedTool(f"{prefix}_matches_no_exons.bed").intersect(BedTool("mutations.bed"), wa=True, wb=True,
                                  output=f"myc_matches_mutations.bed")

    compute_mutation_counts("myc")
    get_total_core_match_counts("myc")

    # intersect MYC binding site mutation counts with facilitated/neutral annotation
    # and the MYC binding site match counts
    mutation_counts = pd.read_csv("myc_mutation_counts.tsv", sep="\t")
    facilitated_df = pd.read_csv("facilitated_myc_mutations.tsv", sep="\t")
    myc_matches = pd.read_csv("myc_match_counts.tsv", sep="\t")
    myc_mutation_df = mutation_counts.merge(
        facilitated_df, on=["suffix", "WT", "mutation", "mutation_pos"]).merge(
        myc_matches, on="suffix").groupby(
        ["MSI", "facilitated", "suffix", "WT", "mutation", "mutation_pos"]).sum().reset_index()

    # intersect control sites with mutations
    myc_control_matches = BedTool("myc_control_matches_no_exons.bed")
    myc_control_matches.intersect(BedTool("mutations.bed"), wa=True, wb=True,
                                  output=f"myc_control_matches_mutations.bed")

    compute_mutation_counts(f"myc_control")
    get_total_core_match_counts(f"myc_control")

    # compute the control rates for MSS/MSI Myc-neutral/-facilitated
    # and output a table
    compute_control_rates(myc_mutation_df)


if __name__ == "__main__":
    main()
