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

# Set raw data directory
RAW_DATA_SRC <- "/Volumes/Seagate/data"

# Load appropriate libraries
library(readr)
library(dplyr)
library(tidyr)
library(parsemsf)
library(readxl)
library(stringr)

normalize_areas <- function(prots) {
  # Replace an NaN (i.e. unobserved proteins, or proteins without enough peptides to quantitate) values with 0
  prots[is.na(prots)] <- 0

  # Normalize areas to total E. coli content
  prots <- prots %>%
    filter(!grepl('CONTAMINANT', protein_desc)) %>% # Filter out contaminants
    ungroup() %>%
    mutate(org = ifelse(grepl('NP', protein_desc), 'phage', 'ecoli')) %>% # Assign organism group
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

process_replicates <- function(df) {

  # Relabel protein groups
  relabel = c('NP_041997.1; NP_041998.1' = 'NP_041998.1',
              'NP_041975.1; NP_041977.1' = 'NP_041975.1',
              'NP_041997.1' = 'NP_041998.1',
              'NP_041977.1' = 'NP_041975.1'
              )

  out_df <- df %>%
    group_by(strain, time) %>%
    do(quantitate(as.character(.$file), normalize = F, relabel = relabel)) %>%
    group_by(strain, time) %>%
    do(normalize_areas(.))

  return(out_df)

}

samples <- read_excel(paste0(RAW_DATA_SRC, "/sample_list.xlsx")) %>%
  filter(`Biological replicate` != 1, `Data Type` == "mass-spec") %>% # Skip replicate 1
  mutate(File = paste0(RAW_DATA_SRC, "/mass_spec/msf_1_4/", str_replace(File, ".raw", ".msf"))) %>%
  rename(b_rep = `Biological replicate`, time = `Time point`, strain = Strain, file = File) %>%
  select(b_rep, strain, time, file) %>%
  nest(-b_rep) %>%
  mutate(msf = purrr::map(data, process_replicates))
  
out <- samples %>% 
  select(b_rep, msf) %>% 
  unnest()

# Write out data
write_csv(out, "./data/proteomics/abundances/all_reps.csv")
