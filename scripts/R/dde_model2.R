rm(list = ls())

library(deSolve)
library(dplyr)
library(tidyr)
library(cowplot)

my_function <- function(t, y, p, delayAC, delayB, M_init, R_init) {

  lag_R <- ifelse((t - delayAC) <= 0, 0, lagvalue(t - delayAC)[7])
  lag_R_b <- ifelse((t - delayAC - delayB) <= 0, 0, lagvalue(t - delayAC - delayB)[7])
  lag_R_c <- ifelse((t - 2 * delayAC - delayB) <= 0, 0, lagvalue(t - 2 * delayAC - delayB)[7])
  
  # k_beta <- min(p["k"], (1 - (y["beta"]/(M_init*38)))*(38/delayB))
  k_beta <- p["k"]
  
  d_alpha <- p["k"] * (p["M"] * y["R"] - p["M"] * lag_R)
  d_A <- p["k"] * p["M"] * lag_R
  
  d_beta <- k_beta * p["M"] * lag_R - k_beta * p["M"] * lag_R_b
  d_B <- k_beta * p["M"] * lag_R_b
  
  d_chi <- p["k"] * (p["M"] * lag_R_b - p["M"] * lag_R_c)
  d_C <- p["k"] * p["M"] * lag_R_c
  
  d_R <- p["k"] * (-p["M"] * y["R"] + p["M"] * lag_R_c)
  
  return(list(c(d_alpha,
                d_A,
                d_beta,
                d_B,
                d_chi,
                d_C,
                d_R)))
  
}

parms <- c(k=1.9e-5, s=40, M = 1000)

times <- seq(0, 1000, 0.5)

init <- c(alpha = 0, A = 0,
          beta = 0, B = 0,
          chi = 0, C = 0,
          R = 10000)

delayAC = 17
delayB = 34

out <- dede(y = init, times = times, func = my_function,
            p = parms, delayAC = delayAC, delayB = delayB,
            M_init = init["M"], R_init = init["R"])

out_df <- as.data.frame(out)

out_df <- gather(out_df, variable, value, -time)

out_df1 <- filter(out_df, variable %in% c("A", "B", "C"))

out_df2 <- filter(out_df, variable %in% c("alpha", "beta", "chi"))

out_df3 <- filter(out_df, variable %in% c("R", "M"))

out_df4 <- filter(out_df, variable %in% c("q.k"))

ggplot(out_df1, aes(x = time, y = value, group = variable, color = variable)) + geom_line()

ggplot(out_df2, aes(x = time, y = value, group = variable, color = variable)) + geom_line()

ggplot(out_df3, aes(x = time, y = value, group = variable, color = variable)) + geom_line()

ggplot(out_df4, aes(x = time, y = value, group = variable, color = variable)) + geom_line()
