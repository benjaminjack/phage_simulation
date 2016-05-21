###########################################################################
## This script processes RNA abundances generated from kallisto.
##
## Author: Benjamin R. Jack
## Email: benjamin.r.jack@gmail.com
## May 2016
###########################################################################

rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(readxl)

# Each RNA seq rep is stored in a directory with all the sample ID information
# Actual abundunces are inside the directory saved as "abundances.tsv"
rna_dirs <- list.dirs("./data/rna_seq")

# Read in sample info
sample_info <- read_excel("./data/SampleList.xlsx") %>%
  select(`Sample #`, Strain, `Time point`) %>%
  rename(sample_id = `Sample #`, strain = Strain, time = `Time point`)

reps <- data.frame(dir_names = basename(rna_dirs[-(1:1)]), path = rna_dirs[-(1:1)]) %>%
  # Extract sample ID and biological replicate from file name
  separate(dir_names, c("sample_id"), sep="-", extra="drop", remove = FALSE) %>%
  separate(dir_names, c("drop", "b_rep"), sep=-2) %>%
  # Drop what's left of file name
  select(-drop) %>%
  # Join in sample_id info
  left_join(sample_info, by = "sample_id") %>%
  group_by(sample_id, time, strain, path, b_rep) %>%
  # Read in abundance data
  do(data = read_tsv(paste0(.$path, "/abundance.tsv"))) %>%
  unnest() %>%
  # Extract organism and gene info
  separate(target_id, c("org", "gene_name", "prot_id", "gene_type"), sep="\\|", extra="drop") %>%
  select(-path)

# Write out RNA data
write_csv(reps, "./data/all_rna.csv")
  
