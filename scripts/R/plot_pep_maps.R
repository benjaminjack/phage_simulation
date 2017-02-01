library(readr)
library(lmerTest)
library(tidyr)
library(dplyr)
library(cowplot)

rep2 <- read_csv("data/proteomics/peptide_mappings/rep2_pep_map.csv")
rep3 <- read_csv("data/proteomics/peptide_mappings/rep3_pep_map.csv")
rep4 <- read_csv("data/proteomics/peptide_mappings/rep4_pep_map.csv")
rep5 <- read_csv("data/proteomics/peptide_mappings/rep5_pep_map.csv")

all_reps <- bind_rows(rep2, rep3, rep4, rep5) %>%
  group_by(strain, time, Proteins, start, end) %>%
  summarize(mean_count = mean(spectral_count)) %>%
  filter(strain != "11-42") %>%
  group_by(time, Proteins, start, end) %>%
  spread(strain, mean_count) %>%
  na.omit() %>%
  mutate(norm_count = `11-44`/`11-46`)

pep_plot <- ggplot(all_reps %>% filter(Proteins == "NP_041998.1", time == "9min"),
       aes(x = start, xend = end, y = norm_count, yend = norm_count)) +
  geom_segment(size = 1) +
  geom_hline(aes(yintercept = 1)) +
  ylim(0, 2) + xlab("peptide position within major capsid protein") +
  ylab("recoded count relative to wildtype")

save_plot("./figures/peptides.pdf", pep_plot)

ggplot(all_reps %>% filter(Proteins == "NP_041998.1", time == "9min"),
       aes(x = start, width = end-start, y = norm_count)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_hline(aes(yintercept = 1))

capsid <- bind_rows(rep2, rep3, rep4, rep5) %>%
  filter(Proteins == "NP_041998.1", time == "9min")

# Create a model without and without the peptide as a random effect
mod1 <- lmer(spectral_count ~ strain + start + strain:start + (1 | Pep_seq), data = capsid, REML = F)
mod2 <- lmer(spectral_count ~ strain + start + (1 | Pep_seq), data = capsid, REML = F)

# Check to see which model is a better fit
anova(mod1,mod2)
