###########################################################################
## This script processes tab delimited files that have been exported from
## from the Thermo software, after a database search has been conducted.
## The input files contain all peptides detected and the corresponding
## protein groups that they have been assigned. Spectral counts are
## normalized using APEX.
##
## Author: Benjamin R. Jack
## Email: benjamin.r.jack@gmail.com
## January 2016
###########################################################################

# Remove any variables in memory
rm(list = ls())

# Load appropriate libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

parse_psms <- function(peps, apex) {
  
  # Correct protein accession groups for major capsid protein and helicase
  # There are actually two different versions of the capsid and helicase proteins due to read-through.
  # It is almost impossible to distinguish them via mass spec because they only differ by their n-termini,
  # so we combine them into one group each for the caspid and helicase.
  peps <- mutate(peps, `Protein Group Accessions` = str_replace(`Protein Group Accessions`, 'NP_041997.1 ;NP_041998.1', 'NP_041998.1')) %>%
    mutate(`Protein Group Accessions` = str_replace(`Protein Group Accessions`, 'NP_041975.1 ;NP_041977.1', 'NP_041975.1')) %>%
    mutate(`Protein Group Accessions` = str_replace(`Protein Group Accessions`, 'NP_041977.1', 'NP_041975.1'))
  
  # Split 'spectrum file' column into time, strain, and rep, and tech rep
  peps <- peps %>% 
    separate(`Spectrum File`, into=c('strain', 'time', 'rep', 'sample_date'), sep = '_') %>%
    separate(rep, into=c('b_rep', 't_rep'), sep=c(-2))
  
  # Replace modified cysteines and methionines with unmodified
  # In practice, this means replacing 'm' with 'M' and 'c' with 'C'
  peps <- peps %>% 
    mutate(Sequence = str_to_upper(Sequence)) # Make peptide sequences uppercase
  
  # Replace isoleucines (I) with leucines (L). They have the same mass and are thus,
  # indistinguishable via mass spec, so peptides are often represented twice, once for
  # each L/I variation.
  peps <- peps %>% 
    mutate(Sequence = str_replace_all(Sequence, 'I', 'L')) %>% # Replace I with L
    group_by(strain, time, b_rep, t_rep, `Protein Group Accessions`, `Search ID`, `First Scan`) %>%
    # Optional: Rank by PEP, keep peptide with highest score
    distinct(Sequence) # Remove duplicate peptides

  peps <- peps %>% ungroup() %>% # Ungroup to avoid errors in ifelse() call
    mutate(new_area = ifelse(!is.na(`Precursor Area`), `Precursor Area`, as.numeric(`Intensity`))) %>% # Define new area as Precursor Area or Intensity if Precursor Area is missing
    group_by(strain, time, b_rep, `Protein Group Accessions`) %>%
    mutate(PSM = n()) %>% # Count number of PSMs in each protein group
    group_by(strain, time, b_rep, t_rep, `Protein Group Accessions`) %>% 
    distinct(new_area) %>% # Remove peptides with identical sequences and areas
    group_by(strain, time, b_rep, t_rep, `Protein Group Accessions`, Sequence) %>% 
    summarize(PSM = unique(PSM), new_area = sum(new_area)) %>% # Sum areas identical peptides with different charges
    group_by(strain, time, b_rep, t_rep, `Protein Group Accessions`) %>%
    top_n(3, new_area) %>% # Take top 3 peptides by area
    group_by(strain, time, b_rep, `Protein Group Accessions`) %>% # Regroup by protein group
    summarize(PSM = unique(PSM), area = mean(new_area, na.rm = T)) %>% # Compute mean areas from top 3 peptides, while retaining PSMs
    rename(protein = `Protein Group Accessions`) # Rename protein group column to something more R friendly

  peps <- peps %>%
    filter(!grepl('CONTAMINANT', protein)) %>% # Filter out contaminants
    left_join(apex, by=c('protein')) %>%
    ungroup() %>%
    mutate(org = ifelse(grepl('NP', protein), 'phage', 'ecoli')) %>% # Assign organism group
    mutate(time = as.numeric(str_extract(time, '[0-9]+')) * 60) %>% # Convert time to seconds
    group_by(strain, time) %>%
    # APEX CALCULATIONS
    filter(Oi != 0) %>% # Remove proteins with Oi values of 0
    mutate(apex_pre = PSM/Oi) %>%
    mutate(apex_sum = sum(apex_pre, na.rm = T)) %>% # Remove NAs for the few Oi values that don't map
    mutate(apex = (PSM/(Oi*apex_sum))) %>%
    ungroup() %>%
    group_by(strain, time, org) %>%
    # Compute sums for each organism type
    mutate(org_area = sum(area), org_psm = sum(PSM), org_apex = sum(apex, na.rm = T)) %>%
    ungroup() %>%
    group_by(strain, time) %>%
    mutate(ecoli_area = ifelse((org == "phage"), (sum(area) - org_area), org_area),
           ecoli_psm = ifelse((org == "phage"), (sum(PSM) - org_psm), org_psm),
           ecoli_apex = ifelse((org == "phage"), (sum(apex, na.rm = T) - org_apex), org_apex)) %>%
    #Normalize to e. coli area
    mutate(area_norm = area/ecoli_area, psm_norm = PSM/ecoli_psm, apex_norm = apex/ecoli_apex)
  
  
#   areas <- peps %>% ungroup() %>% filter(org == e.coli)
#   %>%
#     group_by(strain, time) %>%
#     mutate(total_area = sum(area)) %>%
#     group_by(strain, time) %>%
#     mutate(area_norm = area/total_area)
  
  return(peps)
  
}

# Read in Oi values
apex <- read_tsv("../apex/weka_results/apex_values.oi") %>% 
  select(-3) %>%
  rename(protein = `# PROTEIN_ID`)

# Read in raw data

# Correct wonky Spectrum File formatting for rep 1
rep1 <- read_tsv('./T7_Phage_Rep1_psms.txt') %>% 
  mutate(`Spectrum File` = str_extract(`Spectrum File`, "(11-[0-9]{2})_([0-9]m)_([0-9])([abc])_([0-9]+)")) %>%
  parse_psms(apex)

# Replicate 2 has been mislabeled as replicate 1?
rep2 <- read_tsv('./T7_Phage_Rep2_psms.txt') %>% 
  mutate(`Spectrum File` = str_replace(`Spectrum File`, "(11-[0-9]{2})_([0-9]min)_([0-9])([abc])", "\\1_\\2_2\\4")) %>%
  parse_psms(apex)

# Read in replicate three
rep3 <- read_tsv('./T7 Phage Rep 3 PSMs.txt') %>%
  parse_psms(apex)

rep4 <- read_tsv('./T7_Phage_Rep4_psms.txt') %>% 
  parse_psms(apex)

all_reps <- bind_rows(rep1, rep2, rep3, rep4)

# Join ID Map for simpler tabasco IDs
tabasco <- read_excel('../id_map.xlsx') %>% mutate(Accession = trimws(Accession))
all_reps <- left_join(all_reps, tabasco, by=c('protein' = 'Accession'))

write_csv(all_reps, "../post_processed/rep1234_combined_apex.csv")

sims <- read_csv('../tabasco_runs/092115_A_avg.csv')
sims2 <- read_csv('../tabasco_runs/092115_B_avg.csv') %>% 
  rename(sim_count2 = sim_count) %>% 
  select(tabasco_id, time, sim_count2)

all_data <- left_join(all_reps, sims)
all_data <- left_join(all_data, sims2)
