rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)


k <- 2
s_1 <- 40
s_2 <- 40
s_3 <- 40
tau_1 <- 10
tau_3 <- 15

a_in <- function() min(k, s_1 / tau_1)
# These functions produce a non-linear curve because we vary tau_2, which is in the denominator
b_in <- function(tau_2) min(a_in(), s_2 / tau_2)
c_in <- function(tau_2) min(b_in(tau_2), s_3 / tau_3)

taus <- seq(10, 30, 0.01)

df <- data.frame(delay_B = taus,
                 A = a_in(),
                 B = sapply(taus, b_in),
                 C = sapply(taus, c_in))
df <- gather(df, protein, rate, A:C)

plot <- ggplot(df, aes(x = delay_B, y = rate)) + 
  geom_line() + 
  facet_wrap(~protein, labeller = as_labeller(function(protein) paste0("protein ", protein))) +
  panel_border() +
  xlab("time to translate protein B (s)") +
  ylab("rate of protein production") +
  theme(legend.position="none")

save_plot("figures/model.pdf", plot, base_aspect_ratio = 2.5)



