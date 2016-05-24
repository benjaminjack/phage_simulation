###########################################################################
## This script processes thermo MSF files that have been exported from
## from the Thermo software, after a database search has been conducted.
## The input files contain all peptides detected and the corresponding
## protein groups that they have been assigned.
##
## Author: Benjamin R. Jack
## Email: benjamin.r.jack@gmail.com
## May 2016
###########################################################################

# Remove any variables in memory
rm(list = ls())

# Load appropriate libraries
library(readr)
library(dplyr)
library(tidyr)
library(parsemsf)

normalize_areas <- function(prots) {
  # Replace an NaN (i.e. unobserved proteins, or proteins without enough peptides to quantitate) values with 0
  prots[is.na(prots)] <- 0

  # Normalize areas to total E. coli content
  prots <- prots %>%
    filter(!grepl('CONTAMINANT', Proteins)) %>% # Filter out contaminants
    ungroup() %>%
    mutate(org = ifelse(grepl('NP', Proteins), 'phage', 'ecoli')) %>% # Assign organism group
    group_by(strain, time, org) %>%
    # Compute sums for each organism type
    mutate(org_area = sum(area_mean)) %>%
    ungroup() %>%
    group_by(strain, time) %>%
    # Hack-y way of making a column that shows the sum of all ecoli protein areas
    mutate(ecoli_area = ifelse((org == "phage"), (sum(area_mean) - org_area), org_area)) %>%
    #Normalize to e. coli area
    mutate(area_norm = area_mean/ecoli_area) -> prots2

  return(prots2)
}

process_replicates <- function(files) {

  # Relabel protein groups
  relabel = c('NP_041997.1; NP_041998.1' = 'NP_041998.1',
              'NP_041975.1; NP_041977.1' = 'NP_041975.1',
              'NP_041997.1' = 'NP_041998.1',
              'NP_041977.1' = 'NP_041975.1'
              )

  out_df <- data.frame(file_names = files) %>%
    separate(file_names, into = c("strain", "time", "rep"), sep = "_", remove = F, extra = "drop") %>%
    mutate(strain = basename(strain)) %>%
    group_by(strain, time) %>%
    do(combine_tech_reps(as.character(.$file_names), normalize = F, relabel = relabel)) %>%
    group_by(strain, time) %>%
    do(normalize_areas(.))

  return(out_df)

}

rep1 <- process_replicates(list.files("./Rep1", full.names = T)) %>% mutate(b_rep = 1)
write_csv(rep1, "./post_processed/rep1.csv")

rep2 <- process_replicates(list.files("./Rep2", full.names = T)) %>% mutate(b_rep = 2)
write_csv(rep2, "./post_processed/rep2.csv")

rep3 <- process_replicates(list.files("./Rep3", full.names = T)) %>% mutate(b_rep = 3)
write_csv(rep3, "./post_processed/rep3.csv")

rep4 <- process_replicates(list.files("./Rep4", full.names = T)) %>% mutate(b_rep = 4)
write_csv(rep4, "./post_processed/rep4.csv")

rep5 <- process_replicates(list.files("./Rep5", full.names = T)) %>% mutate(b_rep = 5)
write_csv(rep5, "./post_processed/rep5.csv")
