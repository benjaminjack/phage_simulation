library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# Burst size plot and analysis

burst <- read_csv("data/burst_lysis/burst.csv")

strain_order <- c("11_44  ( /ul)", "11_42 ( /ul)", "11_46  ( /ul)")
strain_labels <- c("recoded", "evolved", "wild-type")
strain_colors <- c("11_44  ( /ul)" = "orange", "11_42 ( /ul)" = "lightblue", "11_46  ( /ul)" = "lightgreen")

burst_plot <- ggplot(burst, aes(x = strain, y = burst_size)) +
  stat_summary(fun.y = "mean", geom = "bar", aes(fill=strain)) +
  geom_line(aes(group=rep)) +
  geom_point() +
  scale_x_discrete(name = "strain", limits = strain_order, labels = strain_labels) +
  ylab("burst size (PFUs)") +
  scale_fill_manual(values = strain_colors) + theme(legend.position = "none")

# Do some t-tests
burst_wide <- dplyr::select(burst, rep, strain, burst_size) %>%
  spread(strain, burst_size)
# Evolved v recoded
t.test(burst_wide$`11_42 ( /ul)`, burst_wide$`11_44  ( /ul)`, paired = T)
# Evolved v wildtype
t.test(burst_wide$`11_42 ( /ul)`, burst_wide$`11_46  ( /ul)`, paired = T)
# Recoded v wildtype
t.test(burst_wide$`11_44  ( /ul)`, burst_wide$`11_46  ( /ul)`, paired = T)

# Lysis time plot and analysis
lysis_data <- read_csv("data/burst_lysis/lysis.csv") %>% filter(!is.na(time))

strain_order2 <- c("11-44", "11-42", "11-46")
strain_labels2 <- c("recoded", "evolved", "wild-type")
strain_colors2 <- c("11-44" = "orange", "11-42" = "lightblue", "11-46" = "lightgreen")

lysis_plot <- ggplot(lysis_data, aes(x = strain, y = time, group = strain)) +
  stat_summary(geom="bar", fun.y = "mean", aes(fill = strain)) +
  geom_jitter(width = 0.25) +
  scale_x_discrete(name = "strain", limits = strain_order2, labels = strain_labels2) +
  ylab("lysis time (minutes)") +
  scale_fill_manual(values = strain_colors2) + theme(legend.position = "none")


burst_lysis <- plot_grid(burst_plot, lysis_plot, labels = c("A", "B"))

save_plot("./figures/burst_lysis.pdf", burst_lysis, base_aspect_ratio = 2)
