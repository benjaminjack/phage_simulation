rm(list = ls())

library(deSolve)
library(dplyr)
library(tidyr)
library(cowplot)

my_function <- function(t, y, p, delayAC, delayB, M_init, R_init) {

  lag_M <- ifelse((t - delayAC) <= 0, 0, lagvalue(t - delayAC)[7])
  lag_M_b <- ifelse((t - delayAC - delayB) <= 0, 0, lagvalue(t - delayAC - delayB)[7])
  lag_M_c <- ifelse((t - 2 * delayAC - delayB) <= 0, 0, lagvalue(t - 2 * delayAC - delayB)[7])
  
  # k_beta <- min(p["k"], (1 - (y["beta"]/(y["M"]*p["s"])))*(p["s"]/delayB))
  k_beta <- p["k"]
  
  d_alpha <- p["k"] * (y["M"] *  p["R"] - lag_M *  p["R"])
  d_A <- p["k"] * lag_M * p["R"]
  
  d_beta <- k_beta * lag_M *  p["R"] - k_beta * lag_M_b *  p["R"]
  d_B <- k_beta * lag_M_b *  p["R"]
  
  d_chi <- p["k"] * (lag_M_b * p["R"] - lag_M_c * p["R"])
  d_C <- p["k"] * lag_M_c * p["R"]
  
  d_M <- p["k"] * (-(y["M"] * p["R"]) + lag_M * p["R"]) + 50
  
  return(list(c(d_alpha,
                d_A,
                d_beta,
                d_B,
                d_chi,
                d_C,
                d_M)))
  
}

parms <- c(k=1.9e-5, s=40, R = 10000)

times <- seq(0, 300, 0.5)

init <- c(alpha = 0, A = 0,
          beta = 0, B = 0,
          chi = 0, C = 0,
          M = 1000)

delayAC = 17
delayB = 50

out <- dede(y = init, times = times, func = my_function,
            p = parms, delayAC = delayAC, delayB = delayB,
            M_init = init["M"], R_init = parms["R"])

out_df <- as.data.frame(out)

out_df <- gather(out_df, variable, value, -time)

out_df1 <- filter(out_df, variable %in% c("A", "B", "C"))

out_df2 <- filter(out_df, variable %in% c("alpha", "beta", "chi"))

out_df3 <- filter(out_df, variable %in% c("R", "M"))

out_df4 <- filter(out_df, variable %in% c("q.k"))

protein <- ggplot(out_df1, aes(x = time, y = value, group = variable, color = variable)) + 
  geom_line(size=1.5) +
  ylim(0,100000) +
  xlim(0,250) +
  ylab("protein abundance") +
  scale_color_discrete(guide = guide_legend(title = "gene"))

save_plot("figures/dde_model2.pdf", protein, base_aspect_ratio = 1.4)

ggplot(out_df2, aes(x = time, y = value, group = variable, color = variable)) + geom_line()

ggplot(out_df3, aes(x = time, y = value, group = variable, color = variable)) + geom_line()

ggplot(out_df4, aes(x = time, y = value, group = variable, color = variable)) + geom_line()
