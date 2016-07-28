library(readr)
library(ggplot2)

burst <- read_csv("data/burst.csv")

strain_order <- c("11_46  ( /ul)", "11_44  ( /ul)", "11_42 ( /ul)")
strain_labels <- c("wildtype", "attenuated", "evolved")

burst_plot <- ggplot(burst, aes(x = strain, y = burst_size)) + 
  geom_point() + 
  geom_line(aes(group=rep)) + 
  stat_summary(fun.y = "mean", geom = "crossbar", fun.ymin="mean", fun.ymax="mean") +
  scale_x_discrete(name = "Strain", limits = strain_order, labels = strain_labels) +
  ylab("Burst Size")

save_plot("./figures/burst_size.pdf", burst_plot)