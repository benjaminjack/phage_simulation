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
b_in <- function(tau_2, c) min(a_in(), s_2 / tau_2) * c  + (1 - c) * min(k, s_2 /tau_2)
c_in <- function(tau_2, c) min(b_in(tau_2, c), s_3 / tau_3) * c  + (1 - c) * min(k, s_3 /tau_3)

delay_func <- function(x) {
  if (x < 0) {
    return(0)
  } else {
    return(1)
  }
}

a_dot <- function(t) a_in()*delay_func(t - tau_1)
b_dot <- function(tau_2, c, z, t) z*b_in(tau_2, c)*delay_func(t - tau_2 - tau_1) + (1 - z)*min(k, s_2/tau_2)*delay_func(t - tau_2)
c_dot <- function(tau_2, c, z, t) z*z*c_in(tau_2, c)*delay_func(t - tau_2 - tau_1 - tau_3) + (1 - z)*min(k, s_3/tau_3)*delay_func(t - tau_3) + z*(1 - z)*min(k, s_2/tau_2)*delay_func(t - tau_2 - tau_3)
# c_dot <- function(tau_2, c, z, t) z*b_dot(tau_2, c, z, t)*delay_func(t - tau_2 - tau_1 - tau_3) + (1 - z)*min(k, s_3/tau_3)*delay_func(t - tau_3)

taus <- seq(15, 30, 0.01)

c <- 0

make_df <- function(z, t) {
  data.frame(delay_B = taus,
             A = a_dot(t),
             B = sapply(taus, b_dot, c, z, t),
             C = sapply(taus, c_dot, c, z, t),
             coupling = z,
             time = t)
}

coupling <- list(0, 0.2, 0.4, 0.6, 0.8, 1)
times <- list(20, 30, 40, 50, 60, 70)



df <- purrr::cross_d(list(coupling = coupling, times = times)) %>% mutate(data = purrr::map2(coupling, times, make_df)) %>% unnest()
df <- gather(df, protein, rate, A:C)

plot <- ggplot(df, aes(x = delay_B, y = rate, color = protein)) +
  geom_line() +
  facet_grid(time~coupling) +
  panel_border() +
  xlab("time to translate protein B (s)") +
  ylab("rate of protein production")

# save_plot("figures/model.pdf", plot, base_aspect_ratio = 3, base_height = 3)



