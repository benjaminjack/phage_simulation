rm(list = ls())

library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(readr)

import_phage_data <- function(file_path, area_path) {
  # Read in excel file
  my_data <- read_excel(file_path)
  
  # Remove columns that are totally blank (NA in every row)
  my_data <- my_data[, colSums(is.na(my_data)) != nrow(my_data)]
  
  # Split 'Description' column
  my_data <- separate(my_data, 
                      Description, 
                      c('gene', 'protein_desc', 'protein_id', 'location'), 
                      sep="\\]\\s\\[[a-zA-Z_]+=")
  
  my_data <- gather(my_data, sample, area, 6:14, -rep) %>%
    separate(sample, c('strain', 'time'), sep=" ") %>%
    mutate(time = as.numeric(str_extract(time, '[0-9]+')) * 60, area = ifelse(is.na(area),0,area))
  
  areas <- read_excel(area_path)
  
  areas <- gather(areas, sample, area, 2:10) %>% 
    rename(org = `Total Area`) %>%
    separate(sample, c('strain', 'time'), sep=" ") %>%
    mutate(time = as.numeric(str_extract(time, '[0-9]+')) * 60) %>% 
    group_by(strain, time, org) %>% 
    summarize(total_area = sum(area)) %>% 
    filter(org == 'E. coli')

  my_data <- inner_join(my_data, areas)

  my_data <- mutate(my_data, n_area = area/total_area)
  
  return(my_data)
  
}

# Parse replicate 1
rep1 <- import_phage_data('../t7_exp_areas_rep1.xlsx', '../total_areas_rep1.xlsx')

# Parse replicate 2
rep2 <- import_phage_data('../t7_exp_areas_rep2.xlsx', '../total_areas_rep2.xlsx')

all_reps <- bind_rows(rep1, rep2)
# all_reps <- rep1

# Read in simulation data
sims <- read_csv('../092115_A_avg.csv')
sims2 <- read_csv('../092115_B_avg.csv') %>% rename(sim_count2 = sim_count) %>% select(tabasco_id, time, sim_count2)

# Read in ID map
id_map <- read_excel('../id_map.xlsx')

# Join tabasco gene IDs (e.g. gp3.5)
# sims <- full_join(sims, id_map)

all_reps <- inner_join(all_reps, id_map)

all_data <- left_join(all_reps, sims)
all_data <- left_join(all_data, sims2)