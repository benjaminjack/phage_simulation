library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(broom)

# Burst size and lysis time estimates

rep2 <- read_csv("data/proteomics/abundances/rep2.csv")
rep3 <- read_csv("data/proteomics/abundances/rep3.csv")
rep4 <- read_csv("data/proteomics/abundances/rep4.csv")
rep5 <- read_csv("data/proteomics/abundances/rep5.csv")

all_data <- bind_rows(rep2, rep3, rep4, rep5) %>%
  mutate(strain = str_replace_all(strain, c("11-44" = "atten", "11-42" = "evol", "11-46" = "wt"))) %>%
  mutate(Proteins = trimws(Proteins))


labels <- read_csv("data/id_map.csv") %>% mutate(Accession = trimws(Accession))
prots <- filter(all_data, org == "phage") %>%
  group_by(strain, time, b_rep) %>%
  # Join by each group seperately so that missing phage proteins are filled in with NAs
  do(full_join(., labels, by = c("Proteins" = "Accession")) %>% fill(strain, time, b_rep)) %>%
  mutate(area_norm = ifelse(is.na(area_norm), 0, area_norm)) %>% # Fill in missing proteins (NAs) with 0
  mutate(tabasco_id = str_replace_all(tabasco_id, c("gp10A" = "gp10"))) # Technically we can't distinguish 10A from 10B so we just compile them together as gp10

prots <- prots %>% filter(org=="phage", tabasco_id == "gp10")

# model <- lm(log(area_norm) ~ time, data = prots)

wildtype <- function(x)2^(43.2*(x/3600))
attenuated <- function(x)2^(35.7*(x/3600))
evolved <- function(x)2^(38.7*(x/3600))
wt_atten_ratio <- function(x)(2^(43.2*(x/3600)))/2^(35.7*(x/3600))
wt_evol_ratio <- function(x)(2^(43.2*(x/3600)))/2^(38.7*(x/3600))

time_points = seq(0, 600, 1)
df1 <- data.frame(x = time_points/60,
                particles = sapply(time_points, wildtype),
                strain = "wildtype")

df2 <- data.frame(x = time_points/60,
                particles = sapply(time_points, attenuated),
                strain = "recoded")

df <- bind_rows(df1, df2)

growth_curve <- ggplot(df, aes(x, particles, group = strain, color = strain)) +
  geom_vline(xintercept = 9, linetype = 2) +
  geom_line() +
  ylab("estimated phage particles") +
  xlab("time (minutes)") +
  scale_x_continuous(breaks = c(0, 1, 5, 9), limits = c(0, 10))

prots2 <- group_by(prots, strain, time) %>% summarize(area_mean = mean(area_norm))

ggplot(prots2, aes(x=time, y=area_mean, group=strain, color=strain)) + geom_point() + geom_line()

prots3 <- group_by(prots, time, b_rep) %>%
  summarize(wt_atten = area_norm[strain=="wt"]/area_norm[strain=="atten"], wt_evol = area_norm[strain=="wt"]/area_norm[strain=="evol"]) %>%
  gather(strain, ratio, -time, -b_rep)

prots3 <- group_by(prots2, time) %>%
  summarize(wt_atten = area_mean[strain=="wt"]/area_mean[strain=="atten"], wt_evol = area_mean[strain=="wt"]/area_mean[strain=="evol"]) %>%
  gather(strain, ratio, -time) %>%
  filter(strain == "wt_atten", time == "9min")

burst <- read_csv("data/burst.csv")

burst %>%
  select(-X1) %>%
  spread(strain, burst_size) %>%
  mutate(wt_atten_ratio = `11_46  ( /ul)`/`11_44  ( /ul)`) %>%
  summarize(mean_wt_atten_ratio = mean(wt_atten_ratio)) -> burst_ratio

ratio_data <- data_frame(wildtype_to_phage = c("capsid (gene 10)", "burst size", "estimated particles (9 min)", "estimated particles (12 min)"),
                         ratio = c(prots3$ratio[1], burst_ratio$mean_wt_atten_ratio[1], wildtype(540)/attenuated(540), wildtype(720)/attenuated(720)))

xlab_order <- rev(c("capsid (gene 10)", "burst size", "estimated particles (9 min)", "estimated particles (12 min)"))

ratio_plot <- ggplot(ratio_data, aes(x=factor(wildtype_to_phage, levels = xlab_order), y=ratio, fill = factor(wildtype_to_phage))) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(x = factor(wildtype_to_phage), y = ratio + 0.15, label = round(ratio, 1))) +
  ylab("wildtype-to-recoded ratio") +
  xlab("") +
  ylim(0, 3) +
  coord_flip() +
  scale_fill_manual(values = c("capsid (gene 10)" = "grey70", "burst size" = "tan", "estimated particles (9 min)" = "gold", "estimated particles (12 min)" = "gold")) + theme(legend.position="none")

# growth_plot <- plot_grid(growth_curve, ratio_plot, labels = c("A","B"), rel_widths = c(1.4, 1))

save_plot("./figures/ratios.pdf", ratio_plot, base_aspect_ratio = 1.7)
