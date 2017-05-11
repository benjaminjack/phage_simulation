library(dplyr)
library(readr)
library(readxl)
library(stringr)
library(tidyr)
library(parsemsf)

# Set raw data directory
RAW_DATA_SRC <- "/Volumes/Seagate/data"

process_replicates <- function(df) {
  
  out_df <- df %>%
    group_by(strain, time) %>%
    mutate(tech_reps = n()) %>%
    group_by(strain, time, tech_reps, file) %>%
    do(map_peptides(as.character(.$file))) %>%
    group_by(strain, time, protein_desc, peptide_sequence, start, end) %>%
    # Some samples have more technical replicates, so I'm including a count that's been normalized
    # to the number of technical replicates
    summarize(spectral_count = n(), norm_count = n()/unique(tech_reps))

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

write_csv(out, "./data/proteomics/peptide_mappings/all_mappings.csv")
