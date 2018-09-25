###
# This script takes raw TSV tables of RNA-seq counts, normalizes them, maps IDs, and combines them into a single
# CSV containing all gene count information. 
###

library(tidyverse)

# Read in ID map
labels <- read_csv("../../output/id_map.csv") %>% 
  mutate(accession = trimws(accession)) # Trim an whitespace that causes problems with joins

import_counts <- function(file_name, replicate, time_point) {
  read_tsv(file_name, col_names=FALSE) %>%
    filter(X1 == "V01146.1", X3 == "CDS") %>%
    select("start" = X4, "end" = X5, "gene" = X9, "count" = X10) %>%
    mutate(gene_number = str_match(gene, "Note=(possible )?gene ([0-9\\.\\-AB]*)")[,3],
           lengthkb = (end - start)/1000,
           rpk = count/lengthkb,
           scale_factor = sum(rpk)/1000000,
           tpm = rpk/scale_factor) %>%  # Compute TPM
    left_join(labels, by = c("gene_number" = "gene_number")) %>%  # Add labels
    mutate(gene_number = str_replace_all(gene_number, c("10A" = "10")))  # We cant distinguish between 10A and 10B, so relabeld 10A -> 10
} 

# Construct dataframe of input files, replicate numbers, and time points
# file_inputs <- bind_rows(
#   list(file = "../../data/rna_seq/37C/recounted/PA36_recounted_q10_9min_rep1.tsv", b_rep = 1, time = 9),
#   list(file = "../../data/rna_seq/37C/recounted/PA45_recounted_q10_9min_rep2.tsv", b_rep = 2, time = 9),
#   list(file = "../../data/rna_seq/37C/recounted/PA54_recounted_q10_9min_rep3.tsv", b_rep = 3, time = 9),
#   list(file = "../../data/rna_seq/37C/recounted/PA35_recounted_q10_5min_rep1.tsv", b_rep = 1, time = 5),
#   list(file = "../../data/rna_seq/37C/recounted/PA44_recounted_q10_5min_rep2.tsv", b_rep = 2, time = 5),
#   list(file = "../../data/rna_seq/37C/recounted/PA53_recounted_q10_5min_rep3.tsv", b_rep = 3, time = 5),
#   list(file = "../../data/rna_seq/37C/recounted/PA34_recounted_q10_1min_rep1.tsv", b_rep = 1, time = 1),
#   list(file = "../../data/rna_seq/37C/recounted/PA43_recounted_q10_1min_rep2.tsv", b_rep = 2, time = 1),
#   list(file = "../../data/rna_seq/37C/recounted/PA52_recounted_q10_1min_rep3.tsv", b_rep = 3, time = 1)
# )

file_inputs <- bind_rows(
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA36-wild10A-46-9min-rep1.counts.tsv", b_rep = 1, time = 9),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA45-wild10A-46-9min-rep2.counts.tsv", b_rep = 2, time = 9),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA54-wild10A-46-9min-rep3.counts.tsv", b_rep = 3, time = 9),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA35-wild10A-46-5min-rep1.counts.tsv", b_rep = 1, time = 5),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA44-wild10A-46-5min-rep2.counts.tsv", b_rep = 2, time = 5),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA53-wild10A-46-5min-rep3.counts.tsv", b_rep = 3, time = 5),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA34-wild10A-46-1min-rep1.counts.tsv", b_rep = 1, time = 1),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA43-wild10A-46-1min-rep2.counts.tsv", b_rep = 2, time = 1),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA52-wild10A-46-1min-rep3.counts.tsv", b_rep = 3, time = 1)
)

counts <- file_inputs %>% 
  mutate(read_counts = map(file, ~import_counts(.))) %>% # Run our import_counts function over all files
  unnest() %>%
  select(-file)

write_csv(counts, "../../output/transcripts_37C.csv")
