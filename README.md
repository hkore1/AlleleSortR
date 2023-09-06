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
