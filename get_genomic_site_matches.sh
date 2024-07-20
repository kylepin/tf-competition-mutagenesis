#!/bin/bash
# This should be run prior to process_mutation_counts
python fastaRegexFinder.py --fasta input/hg19.fa --regex "(?=(CAC(G[CAT]|AT)G))" > intermediate_and_output/myc_matches.bed
python fastaRegexFinder.py --fasta  input/hg19.fa --regex "(?=([CAGT]((?<!C)A|[CGT])C(G[CAT]|AT)G))" > intermediate_and_output/myc_control_matches.bed
