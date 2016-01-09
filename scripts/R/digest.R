###########################################################################
## Takes a fasta file with protein sequences and digests them in silico.
## Results are output in a TSV for use in the APEX protocol.
##
## Author: Benjamin R. Jack
## Email: benjamin.r.jack@gmail.com
## January 2016
###########################################################################

rm(list = ls())

library(cleaver)
library(Biostrings)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)

# Define function to read a fasta file into a data frame
cleave_fasta <- function(fastafile) {
  seqs <- readAAStringSet(fastafile)
  seq_name <- names(seqs)
  sequence <- paste(seqs)
  df <- data.frame(name = seq_name, sequence = sequence)  %>%
    # Extract protein ID
    mutate(protein_id = str_match(name, "^([A-Za-z0-9_\\.\\-:]+)\\s")[,2]) %>%
    select(-name) %>%
    # Cleave peptides using cleaver package
    mutate(sequence = as.character(sequence)) %>%
    mutate(peptide = cleave(sequence, "trypsin", missedCleavages=2)) %>% 
    unnest(peptide) %>%
    group_by(protein_id) %>%
    summarize(peptides = str_c(as.character(peptide), collapse=" "))
  return(df)
}

ecoli_phage <- cleave_fasta("E_coli_REL606_with_T7_bacteriophage_and_MQ_contam.fasta")

# Write out as TSV for APEX perl script
write.table(ecoli_phage, "peptides_test.tsv", quote=F, sep = "\t", row.names = F)
