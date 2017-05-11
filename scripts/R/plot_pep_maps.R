library(readr)
library(lmerTest)
library(tidyr)
library(dplyr)
library(cowplot)

all_reps_csv <- read_csv("data/proteomics/peptide_mappings/all_mappings.csv")

all_reps <- all_reps_csv %>%
  group_by(strain, time, protein_desc, start, end) %>%
  summarize(mean_count = mean(spectral_count)) %>%
  filter(strain != "11-42") %>%
  group_by(time, protein_desc, start, end) %>%
  spread(strain, mean_count) %>%
  na.omit() %>%
  mutate(norm_count = `11-44`/`11-46`)

pep_plot <- ggplot(all_reps %>% filter(protein_desc %in% c("NP_041998.1", "NP_041997.1; NP_041998.1"), time == "9"),
       aes(x = start, xend = end, y = norm_count, yend = norm_count)) +
  geom_segment(size = 1) +
  geom_hline(aes(yintercept = 1)) +
  ylim(0, 2) + xlab("peptide position within major capsid protein") +
  ylab("recoded count relative to wild type")

save_plot("./figures/peptides.pdf", pep_plot)

ggplot(all_reps %>% filter(protein_desc %in% c("NP_041998.1", "NP_041997.1; NP_041998.1"), time == "9"),
       aes(x = start, width = end-start, y = norm_count)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_hline(aes(yintercept = 1))

capsid <- all_reps_csv %>%
  filter(protein_desc %in% c("NP_041998.1", "NP_041997.1; NP_041998.1"), time == "9")

# Create a model without and without the peptide as a random effect
mod1 <- lmer(spectral_count ~ strain + start + strain:start + (1 | peptide_sequence), data = capsid, REML = F)
mod2 <- lmer(spectral_count ~ strain + start + (1 | peptide_sequence), data = capsid, REML = F)

# Check to see which model is a better fit
anova(mod1,mod2)
