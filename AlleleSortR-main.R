#!/usr/bin/env Rscript

library(seqinr)
library(reshape2)
library(dplyr)
library(stringr)
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(optparse)


#################
### Arguments ###
#################

option_list = list(
  make_option(c("-a", "--input_alignment_directory"), type="character", default=NULL, 
              help="Path to input directory of fasta alignments", metavar="character"),
  make_option(c("-o", "--output_directory"), type="character", default="AlleleSortR_output", 
              help="Name for the output directory [default = %default]", metavar="character"),
  make_option(c("-n", "--number_of_samples"), type="integer", default=NULL, 
              help="The maximum number of sequences in any alignment ", metavar="integer"),
  make_option(c("-i", "--iterations"), type="integer", default=1, 
              help="The number of iterations to run the sorter for [default = %default]", metavar="integer"),
  make_option(c("-c", "--cut_between_iterations"), type="character", default="yes", 
              help="Strip previously sorted sequences between iterations? [default = %default]", metavar="yes or no"),
  make_option(c("-f", "--cut_on_final_iteration"), type="character", default="yes", 
              help="Strip sorted sequences in the final iteration? (i.e. so final alignments have no sorted seqs) [default= %default]", metavar="yes or no")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input_alignment_directory) | is.null(opt$number_of_samples)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input_alignment_directory) (number_of_samples)", call.=FALSE)
}

alignment_directory = opt$input_alignment_directory # Folder containing the input .fasta alignments 
output_dir = opt$output_directory # Where should the output be deposited
number_of_samples = opt$number_of_samples # The maximum number of samples in any alignment. (This is a bit arbitrary as it simply needs to be no less than the max, i.e. it can be any value greater than the max, but the greater the amount over will increase compute time of the first for loop.
iterations = opt$iterations # How many iterations to run the sorter for
cut_between_iterations =  opt$cut_between_iterations # Strip previously sorted sequences between iterations?
cut_on_final_iteration =  opt$cut_on_final_iteration # Strip sorted sequences in the final iteration? (i.e. so final alignments have no sorted seqs)



#################
### Functions ###
#################

# Define a function to find the most closely related samples for each sample from a distance matrix
find_closely_related_samples <- function(dist_matrix) {
  num_pairs <- nrow(dist_matrix)
  closely_related_samples <- matrix(NA, nrow = num_pairs, ncol = 2)
  
  for (i in 1:num_pairs) {
    distances <- dist_matrix[i, ]
    closest_sample_indices <- paste(names(distances[which(distances == min(distances, na.rm = T))]), collapse = "-")
    sample <- rownames(ali_dist)[i]
    closely_related_samples[i, ] <- c(sample, closest_sample_indices)
    
  }
  
  return(closely_related_samples)
}



######################################
### Print some feedback on startup ###
######################################

cat("\n\n\n*****************************\n")
cat("*** Beginning AlleleSortR ***\n")
cat("*****************************\n")
cat(paste0("\nUsing the following input directory: ", alignment_directory))
cat(paste0("\nOutputting to: ", output_dir))
cat(paste0("\nThe maximum number of samples in any alignment: ", number_of_samples))
cat(paste0("\nThe number of iterations to run: ", iterations))
cat(paste0("\nShould split alleles be cut between iterations? ", cut_between_iterations))
cat(paste0("\nShould split alleles be excluded from the output alignments of the last iteration? ", cut_on_final_iteration))



#########################
### Start the program ###
#########################

# Create output directory for renamed alignments
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

for (iter in 1:iterations){
  
  cat(paste0("\n\nRunning iteration ", iter, " out of ", max(iterations)))


  #############################################################################################
  ### Part 1 - Calculate distance matrices from the initial phased gene trees and summarise ###
  #############################################################################################
  
  cat("\n\nPerforming initial calculation of split alleles...")
  
  # Get the names of the alignment files
  ifelse(iter == 1,
         alignments <- list.files(path=alignment_directory, pattern="*.fasta", full.names=TRUE, recursive=FALSE),
         alignments <- list.files(path=paste0(output_dir, "_alignments_Iteration", iter-1), pattern="*.fasta", full.names=TRUE, recursive=FALSE))
  
  # Create the data frame to which all best matches for each sample across all alignments will be written.
  master_data <- data.frame()
  
  # Large for loop to calculate distance matrices and combine them
  for (i in 1:length(alignments)){
    
    cat(paste0("\nWorking on alignment: ", alignments[i]))
    
    # Read alignment
    alignment = alignments[i]
    ali <- read.alignment(alignment, format = "fasta")
    
    # Remove hyphens from sample names
    ali$nam <- gsub("-","",ali$nam)
    
    # If the sorter is being run for multiple iterations, the input alignments will have sequences that have 
    #   already been renamed. This excludes any already renamed samples from the calculation of master_data.
    if (cut_between_iterations == "yes"){
        ali <- as.alignment(nb = length(ali$nam[!grepl("_matched_to_", ali$nam)]),
                            nam = ali$nam[!grepl("_matched_to_", ali$nam)],
                            seq = ali$seq[!grepl("_matched_to_", ali$nam)],
                            com = ali$com[!grepl("_matched_to_", ali$nam)])
    }
    
    #Convert to pairwise dist matrix, with diagonal as NA
    ali_dist <- as.matrix(dist.alignment(ali))
    diag(ali_dist) <- NA
    
    # Find most closely related pairs for each pair in the distance matrix... 
    # This produces a data frame with the most similar sample(s) for each sample.
    matrix <- find_closely_related_samples(ali_dist)
    closely_related_samples <- as.data.frame(matrix)
    colnames(closely_related_samples) <- c("Sample", "BestMatches")
    
    # Remove samples where the best match column contains the corresponding allele of the same sample
    closely_related_samples_stripped = data.frame()
    for (i in 1:length(row.names(closely_related_samples))){
      if (grepl(gsub("_.*", "", closely_related_samples$Sample[i]), closely_related_samples$BestMatches[i]) != TRUE) {
        closely_related_samples_stripped <- rbind(closely_related_samples_stripped, closely_related_samples[i, ])
      }
    }
    
    # Split BestMatches to multiple columns
    closely_related_samples_split <- cbind(closely_related_samples_stripped$Sample,
                                           str_split_fixed(closely_related_samples_stripped$BestMatches, "-", number_of_samples)) 
    
    # Name the columns
    colnames(closely_related_samples_split) <- c("Sample",
                                                 c(1:number_of_samples))
    
    # Append the alignment info to the master_data
    master_data <- rbind(master_data, closely_related_samples_split)
  }
  
  
  ## Hidden option to save the initial master dataframe (can be useful for checkpointing 
  #   between the initial calculation of allele splits and renaming steps)
  #saveRDS(master_data, paste0(output_dir, "/master_data_Iteration", iter, ".rds"))
  #master_data <- readRDS(paste0(output_dir, "/master_data_Iteration", iter, ".rds"))
  
  # Convert empty cells to NA and filter out rows with no match
  master_data <- master_data %>% 
    mutate_all(na_if, "") %>%
    filter(!is.na(`1`))
  
  
  ## Visualise per sample relationships to gauge overall pairwise relationships
  
  # Melt to enable grouping of matching samples to target sample
  x <- melt(master_data, id.vars="Sample", na.rm = T)
  
  # If an allele sequence is always paired with its co-allele, it will not end up in the master_data (which is good),
  #   but this results in uneven sequences between x$Sample and x$value that break the frequency table and heatmap.
  #   This removes those sequences.
  all_names <- unique(c(x$Sample, x$value))
  dif_names <- c(setdiff(x$value, x$Sample), setdiff(x$Sample, x$value))
  while (length(dif_names) != 0){
    x <- x %>%
      filter(!Sample %in% dif_names) %>%
      filter(!value %in% dif_names)
    all_names <- unique(c(x$Sample, x$value))
    dif_names <- c(setdiff(x$value, x$Sample), setdiff(x$Sample, x$value))
  }
  
  # Convert to matrix
  y <- as.matrix(table(x$Sample, x$value))
  
  # Visualise as heatmap -- the convoluted nature of this code is to plot the heatmap without a dendrogram, but still use the
  #   order of the dendrogram generated in the 'z' object
  colPalette <- colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(256)
  pdf(NULL) # Null device so that stats::heatmap doesnt write a default plot in the run directory
  z <- stats::heatmap(y, cexRow = 0.2, cexCol = 0.2, symm = T)
  dev.off() # Turn off null 
  a <- y[,z$rowInd]
  b <- a[z$colInd,]

  pdf(file = paste0(output_dir, "/heatmap_Iteration", iter, ".pdf"), width = 20, height = 20)
  heatmap(scale(b), Rowv = NA, Colv = NA, cexRow = 0.2, cexCol = 0.2, col = colPalette)
  dev.off()
  
  
  ## Save the frequency table that was used to generate the heatmap
  
  # Convert frequency table to data frame
  pairwise_frequency_table <- as.data.frame.matrix(y, make.names = F)
  
  # Write frequency table to csv
  write.csv(pairwise_frequency_table, paste0(output_dir, "/cooccurence_matrix_Iteration", iter, ".csv"))
  
  
  
  ############################################################################
  ### Part 2 - Trim data to relevant samples and create table for renaming ###
  ############################################################################
  
  ## Count the number of times each allele is split from its pair ##
  highest_scoring <- master_data %>%
    group_by(Sample) %>%
    count()
  
  # Trim to most split samples, and convert to vector
  highest_scoring_names <- filter(highest_scoring, n > length(alignments)/2)$Sample
  
  
  ## Now we grab the best-pairing sample for each allele ##
  
  # First, trim the frequency table to include only alleles in highest_scoring_names
  pairwise_frequency_table_subset <- subset(pairwise_frequency_table, row.names(pairwise_frequency_table) %in% highest_scoring_names)
  
  # For each allele get the name of the highest-pairing sample across all genes
  renaming_table <- data.frame()
  for (i in 1:length(rownames(pairwise_frequency_table_subset))){
    sample_row = pairwise_frequency_table_subset[i, ]
    highest_paired <- gsub("_h.*", "", names(which.max(sample_row)))
    name_combination <- c(row.names(sample_row), highest_paired)
    renaming_table <- rbind(renaming_table, name_combination)
  }
  colnames(renaming_table) <- c("sample", "name2match")
  
  # In some cases, only one allele makes it through to the renaming table. This is because
  #  the count in highest_scoring might place each allele on either side of the outlier threshold.
  #  To deal with this, we just check that both alleles are in the renaming table, and if there is
  #  only one, then we remove it.
  renaming_table <- renaming_table[(duplicated(gsub("_.*", "", renaming_table$sample)) | duplicated(gsub("_.*", "", renaming_table$sample), fromLast = TRUE)), ]
  
  # Write renaming table to csv
  write.csv(renaming_table, paste0(output_dir, "/renaming_table_Iteration", iter, ".csv"), row.names = FALSE)
  
  # Write the count of split alleles to csv
  write.csv(highest_scoring, paste0(output_dir, "/allele_split_count_Iteration", iter, ".csv"), row.names = FALSE)
  
  
  
  ###############################################################################
  ### Part 3 - Re-iterate through the alignments and rename alleles sequences ###
  ###############################################################################
  
  cat("\n\nPerforming re-calculation of split alleles and re-naming...")
  
  output_dir2 = paste0(output_dir, "_alignments_Iteration", iter)
  
  if (exists("renaming_log") == FALSE){
    renaming_log <- data.frame() 
  }
  
  # Create output directory for renamed alignments
  if(!dir.exists(output_dir2)){
    dir.create(output_dir2)
  }
  
  # Loop through each alignment and, when corresponding alleles of a sample are not paired 
  #   with each other (based on pairwise distances between all samples), rename 
  for (i in 1:length(alignments)){
    
    # Read alignment
    alignment = alignments[i]
    ali <- read.alignment(alignment, format = "fasta")
    
    #Remove hyphens from sample names
    ali$nam <- gsub("-","",ali$nam)
    
    # If the sorter is being run for multiple iterations, the input alignments will have sequences that have 
    #   already been renamed. This trims any renamed samples from the alignments that are to be renamed if 
    #   cut_between_iterations =  "yes".
    if (cut_between_iterations == "yes"){
      ali <- as.alignment(nb = length(ali$nam[!grepl("_matched_to_", ali$nam)]),
                          nam = ali$nam[!grepl("_matched_to_", ali$nam)],
                          seq = ali$seq[!grepl("_matched_to_", ali$nam)],
                          com = ali$com[!grepl("_matched_to_", ali$nam)])
    }
    
    #Convert to pairwise dist matrix, with diagonal as NA
    ali_dist <- as.matrix(dist.alignment(ali))
    diag(ali_dist) <- NA
    
    # Find most closely related pairs for each pair in the distance matrix... 
    # This produces a data frame with the most similar sample(s) for each sample.
    matrix <- find_closely_related_samples(ali_dist)
    closely_related_samples <- as.data.frame(matrix)
    colnames(closely_related_samples) <- c("Sample", "BestMatches")
    
    # Remove samples where the best match column contains the corresponding allele of the same sample
    closely_related_samples_stripped = data.frame()
    for (i in 1:length(row.names(closely_related_samples))){
      if (grepl(gsub("_.*", "", closely_related_samples$Sample[i]), closely_related_samples$BestMatches[i]) != TRUE) {
        closely_related_samples_stripped <- rbind(closely_related_samples_stripped, closely_related_samples[i, ])
      }
    }
    
    # Remove samples where both alleles are in closely_related_samples_stripped
    closely_related_samples_stripped <- closely_related_samples_stripped[!(duplicated(gsub("_.*", "", closely_related_samples_stripped$Sample)) | duplicated(gsub("_.*", "", closely_related_samples_stripped$Sample), fromLast = TRUE)), ]
    
    ### The EXTRA IMPORTANT BIT!! - renaming ###
    # For each allele sequence in the alignment:
    # 1) Check if the sequence should be renamed (i.e. is in closely_related_samples_stripped$Sample)
    # 2) If yes, then rename the sequence as per the renaming table
    original_names <- c() # Empty vector for storing original names
    new_names <- c() # Empty vector for storing new names
    for (i in 1:length(renaming_table$sample)){
      if (renaming_table$sample[i] %in% closely_related_samples_stripped$Sample){
        seq_renamed <- ali$nam[match(renaming_table$sample[i], ali$nam)] # Keep track of sequence being renamed for logging
        renamed_seq <- paste0(gsub("_h.*","",renaming_table$sample[i]), "_matched_to_", renaming_table$name2match[i]) # Keep track of the name changed to for logging
        
        original_names <- c(original_names, seq_renamed)
        new_names <- c(new_names, renamed_seq)
      
        ali$nam[match(renaming_table$sample[i], ali$nam)] <- paste0(gsub("_h.*","",renaming_table$sample[i]), "_matched_to_", renaming_table$name2match[i]) # Does the renaming
        }
    }
    
    # Remove split alleles from final iteration if requested
    if (iter == iterations & cut_on_final_iteration == "yes"){
      ali <- as.alignment(nb = length(ali$nam[!grepl("_matched_to_", ali$nam)]),
                          nam = ali$nam[!grepl("_matched_to_", ali$nam)],
                          seq = ali$seq[!grepl("_matched_to_", ali$nam)],
                          com = ali$com[!grepl("_matched_to_", ali$nam)])
    }
    
    cat(paste0("\nWriting output alignment: ", paste0(output_dir2, "/", basename(alignment))))
    
    # Write alignment with renamed sequences (if any) to new file
    write.fasta(sequences = ali$seq,
                names = ali$nam,
                file.out = paste0(output_dir2, "/", basename(alignment)))
    
    # Only add to renaming log if there are names that have been changed
    if (length(original_names) != 0 & length(new_names) != 0){
      # Create the single alignment info data
      renamed_info <- data.frame("Original_names" = original_names,
                                 "New_names" = new_names,
                                 "Alignment" = basename(alignment),
                                 "Iteration" = iter)
      
      # Append the renaming information to the master_data
      renaming_log <- rbind(renaming_log, renamed_info)
      
    }
    
  }
  
  cat(paste0("\n\nCompleted iteration ", iter, "\n\n"))
  
}

write.csv(renaming_log, paste0(output_dir, "/renaming_log.csv"), row.names = FALSE)
