library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

burst <- read_csv("data/burst.csv")

strain_order <- c("11_44  ( /ul)", "11_42 ( /ul)", "11_46  ( /ul)")
strain_labels <- c("recoded", "evolved", "wildtype")

burst_plot <- ggplot(burst, aes(x = strain, y = burst_size)) + 
  geom_point(color="grey85") + 
  geom_line(aes(group=rep), color="grey85") + 
  stat_summary(fun.y = "mean", geom = "errorbar", fun.ymin="mean", fun.ymax="mean", size=1, width=0.75) +
  scale_x_discrete(name = "strain", limits = strain_order, labels = strain_labels) +
  ylab("burst size")

save_plot("./figures/burst_size.pdf", burst_plot)

# Do some t-tests
burst_wide <- select(burst, rep, strain, burst_size) %>% 
  spread(strain, burst_size)
# Evolved v recoded
t.test(burst_wide$`11_42 ( /ul)`, burst_wide$`11_44  ( /ul)`, paired = T)
# Evolved v wildtype
t.test(burst_wide$`11_42 ( /ul)`, burst_wide$`11_46  ( /ul)`, paired = T)
# Recoded v wildtype
t.test(burst_wide$`11_44  ( /ul)`, burst_wide$`11_46  ( /ul)`, paired = T)

