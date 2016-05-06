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
    mutate(ecoli_area = ifelse((org == "phage"), (sum(area_mean) - org_area), org_area)) %>%
    #Normalize to e. coli area
    mutate(area_norm = area_mean/ecoli_area) -> prots2
  
    return(prots2)
}

files <- c("/Volumes/Seagate/proteomics/Rep1/Bull_11-44_9m_1a_08312014.msf", 
           "/Volumes/Seagate/proteomics/Rep1/Bull_11-44_9m_1b_08312014.msf")
rep1 <- combine_tech_reps(files, normalize = F, match_peps = F) %>% mutate(strain = '11-46', time = 9, rep = 1)
rep1 <- normalize_areas(rep1)
write_csv(rep1, "/Volumes/Seagate/proteomics/post_processed/rep1.csv")

files <- c("/Volumes/Seagate/proteomics/Rep4/11-44_9min_4a_20151113.msf", 
           "/Volumes/Seagate/proteomics/Rep4/11-44_9min_4b_20151113.msf")
rep2 <- combine_tech_reps(files, normalize = F, match_peps = T) %>% mutate(strain = '11-46', time = 9, rep = 4)
rep2 <- normalize_areas(rep2)
write_csv(rep2, "/Volumes/Seagate/proteomics/post_processed/rep4.csv")

files <- c("/Volumes/Seagate/proteomics/Rep5/11-44_9min_5a.msf", 
           "/Volumes/Seagate/proteomics/Rep5/11-44_9min_5b.msf",
           "/Volumes/Seagate/proteomics/Rep5/11-44_9min_5c.msf")
rep3 <- combine_tech_reps(files, normalize = F, match_peps = T) %>% mutate(strain = '11-46', time = 9, rep = 5)
rep3 <- normalize_areas(rep3)
write_csv(rep3, "/Volumes/Seagate/proteomics/post_processed/rep5.csv")

files <- c("/Volumes/Seagate/proteomics/Rep2/11-44_9min_1a_20150512.msf", 
           "/Volumes/Seagate/proteomics/Rep2/11-44_9min_1b_20150512.msf")
rep4 <- combine_tech_reps(files, normalize = F, match_peps = T) %>% mutate(strain = '11-46', time = 9, rep = 2)
rep4 <- normalize_areas(rep4)
write_csv(rep3, "/Volumes/Seagate/proteomics/post_processed/rep2.csv")

all_strains <- bind_rows(rep1, rep2, rep3, rep4)
