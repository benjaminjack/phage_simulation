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
b_in <- function(tau_2, c) min(a_in(), s_2 / tau_2) * c + (1 - c) * k
c_in <- function(tau_2, c) min(b_in(tau_2, c), s_3 / tau_3) * c + (1 - c) * k

taus <- seq(15, 30, 0.01)

c <- 0.2
df1 <- data.frame(delay_B = taus,
                 A = a_in(),
                 B = sapply(taus, b_in, c),
                 C = sapply(taus, c_in, c),
                 coupling = c)

c <- 0.4
df2 <- data.frame(delay_B = taus,
                  A = a_in(),
                  B = sapply(taus, b_in, c),
                  C = sapply(taus, c_in, c),
                  coupling = c)
c <- 0.6
df3 <- data.frame(delay_B = taus,
                  A = a_in(),
                  B = sapply(taus, b_in, c),
                  C = sapply(taus, c_in, c),
                  coupling = c)
c <- 0.8
df4 <- data.frame(delay_B = taus,
                  A = a_in(),
                  B = sapply(taus, b_in, c),
                  C = sapply(taus, c_in, c),
                  coupling = c)

df <- bind_rows(df1, df2, df3, df4)
df <- gather(df, protein, rate, A:C)

plot <- ggplot(df, aes(x = delay_B, y = rate, color = protein)) +
  geom_line() +
  facet_wrap(~coupling, nrow = 1) +
  panel_border() +
  xlab("time to translate protein B (s)") +
  ylab("rate of protein production")

save_plot("figures/model.pdf", plot, base_aspect_ratio = 3, base_height = 3)



