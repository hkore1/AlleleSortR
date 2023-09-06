# AlleleSortR
A program for sorting and removing non-monophyletic allele sequences from phased sequence alignments generated using the allele phasing pipeline described [here](https://github.com/hkore1/TargetAllelePhasing/tree/main).

It assumes that there are two allele sequences for each sample in every alignment. If homozygous sequences are treated as single sequences in your alignments, then you can 'phony-phase' them using the script 'phony_phase_homozygotes.py' from [here](https://github.com/hkore1/python_scripts/blob/main/phony_phase_homozygotes.py).

## Required packages:
AlleleSortR requires the following packages:
* optparse
* seqinr
* reshape2
* dplyr
* stringr
* tidyverse
* magrittr
* RColorBrewer

Batch install with:

`install.packages(c("optparse", "seqinr", "reshape2", "dplyr", "stringr", "tidyverse", "magrittr", "RColorBrewer"))`

## Overview
Steps:

1. Loop over each alignment:
   * Calculate pairwise distances between each allele sequence
   * For each sequence, take the best match(es) with the smallest distances and record these in a table
   * If the corresponding allele of a sequence is among the best matches, then remove this allele from the table
2. Combine the tables from each alignment
3. Calculate a co-occurrence matrix, with each cell containing a count of the number of times that each sequence is paired
   * This is plotted as a heatmap
4. Count the number of times that each allele sequence is split from its corresponding allele across all alignments

One iteration of this process is detailed in the below diagram:
![overview.pdf](overview.jpg)
