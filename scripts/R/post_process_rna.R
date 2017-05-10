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
library(stringr)

# Each RNA seq rep is stored in a directory with all the sample ID information
# Actual abundunces are inside the directory saved as "abundances.tsv"
RNA_SRC <- "./data/rna_seq/kallisto_output/no_rrna_trna_merge/"

# Read in sample info
sample_info <- read_excel("./data/sample_list.xlsx") %>%
  filter(`Data Type` == "rna-seq") %>%
  select(`Sample`, Strain, `Time point`, `Biological replicate`, `File`) %>%
  rename(sample_id = `Sample`, strain = Strain, time = `Time point`, file = `File`, b_rep = `Biological replicate`)

reps <- sample_info %>%
  # Extract sample ID and biological replicate from file name
  mutate(path = paste0(RNA_SRC, str_replace(file, ".fastq", ""), "/abundance.tsv")) %>%
  # Read in abundance data
  mutate(data = purrr::map(path, read_tsv)) %>%
  unnest() %>%
  # Extract organism and gene info
  separate(target_id, c("org", "gene_name", "prot_id", "gene_type"), sep="\\|", extra="drop") %>%
  select(-path, -file) %>%
  # Replace org column wiht something more useful
  mutate(org = ifelse(org == "NC_001604.1", "phage", "ecoli"))

# Write out RNA data
write_csv(reps, "./data/rna_seq/post_processed/all_rna_no_ribo_trna_merge.csv")
  
