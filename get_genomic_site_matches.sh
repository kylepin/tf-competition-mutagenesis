#!/bin/bash
# This should be run prior to process_mutation_counts
/Users/kylepinheiro/opt/anaconda3/envs/wei_paper/bin/python fastaRegexFinder2.py --fasta  hg19.fa --regex "(?=(CAC(G[CAT]|AT)G))" > myc_matches.bed
/Users/kylepinheiro/opt/anaconda3/envs/wei_paper/bin/python fastaRegexFinder2.py --fasta  hg19.fa --regex "(?=([CAGT]((?<!C)A|[CGT])C(G[CAT]|AT)G))" > myc_control_matches.bed
