# AlleleSortR
A program for sorting and removing non-monophyletic allele sequences from phased sequence alignments generated using the allele phasing pipeline described [here](https://github.com/hkore1/TargetAllelePhasing/tree/main).

It assumes that there are two allele sequences for each sample in every alignment. If homozygous sequences are treated as single sequences in your alignments, then you can 'phony-phase' them with the script 'phony_phase_homozygotes.py' found [here](https://github.com/hkore1/python_scripts/blob/main/phony_phase_homozygotes.py).

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
The program works iteratively, where the user provides the number of iterations to run and the first iteration uses the provided input alignments. Subsequent iterations use the sorted/renamed alignments that were generated in the previous iteration.

Steps:

1. Loop over each alignment:
   * Calculate pairwise distances between each allele sequence
   * For each sequence, take the best match(es) with the smallest distances and record these in a table
   * If the corresponding allele of a sequence is among the best matches, then remove this allele from the table
2. Combine the tables from each alignment
3. Calculate a co-occurrence matrix, with each cell containing a count of the number of times that each sequence is paired
   * This is plotted as a heatmap
4. Count the number of times that each allele sequence is split from its corresponding allele across all alignments
5. Retain the names of allele sequences that are split in more than half of the total number of alignments
6. For each allele sequence that passes this filter, get the name of allele sequence that it most commonly pairs with across all alignments
   * These get associated in a one-to-one 'renaming table' that is output while the program is running
7. Sort and rename the alleles
   * For each alignment:
     - Re-run steps 1 and 2 as earlier outlined to determine the best matching sequence(s) for each allele
     - If the allele sequence in question does not match with its corresponding allele from the same biological sample & If the sequence it is associated with in the renaming table is among the best matches, then rename the sequence according to the renaming table. Renamed sequences follow the format 'ORIGINAL_SEQ_NAME_matched_to_NAME_2_MATCH'
     - Write out the alignment with newly renamed sequences.

There are options to tell the program whether to remove the renamed sequences (i.e., the "_matched_to_" sequences) between iterations and from the final alignments (). The user may want to retain these renamed sequences to infer potential relationships of a subgenome, however these may also simply be 'rogue' alleles that do not pair with their corresponding allele due to non-biological processes - in this case their inclusion will lead to the inclusion of noise.

The first iteration of this process is detailed in the below diagram:
![overview.pdf](overview.jpg)
