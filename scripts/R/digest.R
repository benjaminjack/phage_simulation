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
    mutate(protein_id = str_match(name, "[(\\[protein_id=)|(ref|)]([A-Z]{2}_[0-9\\.]+)[(\\])|\\|]")[,2]) %>%
    select(-name) %>%
    # Cleave peptides using cleaver package
    mutate(sequence = as.character(sequence)) %>%
    mutate(peptide = cleave(sequence, "trypsin", missedCleavages=2)) %>% 
    unnest(peptide) %>%
    group_by(protein_id) %>%
    summarize(peptides = str_c(as.character(peptide), collapse=" "))
  return(df)
}

phage <- cleave_fasta("t7_proteome.txt")
ecoli <- cleave_fasta("ec_k12.fasta")
ecoli_phage <- bind_rows(phage, ecoli)

# Write out as TSV for perl script
write.table(ecoli_phage, "peptides_test.tsv", quote=F, sep = "\t", row.names = F)
