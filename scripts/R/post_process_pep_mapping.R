library(dplyr)
library(readr)
library(parsemsf)

process_replicates <- function(files) {

  out_df <- data.frame(file_names = files) %>%
    separate(file_names, into = c("strain", "time", "rep"), sep = "_", remove = F, extra = "drop") %>%
    mutate(strain = basename(strain)) %>%
    group_by(strain, time) %>%
    mutate(tech_reps = n()) %>%
    group_by(strain, time, rep, tech_reps, file_names) %>%
    do(map_peptides(as.character(.$file_names))) %>%
    group_by(strain, time, Proteins, Pep_seq, start, end) %>%
    # Some samples have more technical replicates, so I'm including a count that's been normalized
    # to the number of technical replicates
    summarize(spectral_count = n(), norm_count = n()/unique(tech_reps))

  return(out_df)

}

rep1 <- process_replicates(list.files("./Rep1", full.names = T)) %>% mutate(b_rep = 1)
write_csv(rep1, "./post_processed/peptide_mappings/rep1_pep_map.csv")

rep2 <- process_replicates(list.files("./Rep2", full.names = T)) %>% mutate(b_rep = 2)
write_csv(rep2, "./post_processed/peptide_mappings/rep2_pep_map.csv")

rep3 <- process_replicates(list.files("./Rep3", full.names = T)) %>% mutate(b_rep = 3)
write_csv(rep3, "./post_processed/peptide_mappings/rep3_pep_map.csv")

rep4 <- process_replicates(list.files("./Rep4", full.names = T)) %>% mutate(b_rep = 4)
write_csv(rep4, "./post_processed/peptide_mappings/rep4_pep_map.csv")

rep5 <- process_replicates(list.files("./Rep5", full.names = T)) %>% mutate(b_rep = 5)
write_csv(rep5, "./post_processed/peptide_mappings/rep5_pep_map.csv")
