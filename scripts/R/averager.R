# This script averages together TABASCO simulation runs
rm(list=ls())

library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

read_mol_sim <- function(filename) {
  # Read in file lines
  myrun <- read_lines(filename)
  
  # Strip out wonky column headers that are not true columns
  myrun[1] <- str_replace(myrun[1], '\\tNum.+', '')
  
  # Collapse back into string and read in as TSV
  myrun <- str_c(myrun, collapse='\n')
  myrun <- read_tsv(myrun)
  
  return(myrun)
}

sim_files <- list.files(pattern = '^.+-Mol_sim1.txt$')

all_sims <- bind_rows(lapply(sim_files, read_mol_sim))

avg_sims <- mutate(all_sims, `>time` = round(`>time`)) %>% # Round time points so that we can properly group the data
  group_by(`>time`) %>%
  summarise_each(funs(mean), -`>time`) %>% # Compute mean counts for each gene and time point
  select(-contains('RNA')) # Select only proteins

# Convert to tidy data frame for easy plotting
avg_sims_tidy <- gather(avg_sims, gene, count, matches('*gp*')) %>%
  rename(tabasco_id = gene, sim_count = count, time = `>time`)

write_csv(avg_sims_tidy, "../092115_A_avg.csv")

# myplot <- ggplot(avg_sims_tidy, aes(`>time`, sim_count)) + 
#   geom_line() +
#   facet_wrap(~ gene) +
#   panel_border()