rm(list = ls())

library(deSolve)
library(dplyr)
library(tidyr)
library(cowplot)

my_function <- function(t, y, p, delayAC, delayB, M_init, R_init) {

  lag_M <- ifelse((t - delayAC) <= 0, 0, lagvalue(t - delayAC)[8])
  lag_M_b <- ifelse((t - delayAC - delayB) <= 0, 0, lagvalue(t - delayAC - delayB)[8])
  lag_M_c <- ifelse((t - 2 * delayAC - delayB) <= 0, 0, lagvalue(t - 2 * delayAC - delayB)[8])
  lag_R <- ifelse((t - delayAC) <= 0, 0, lagvalue(t - delayAC)[7])
  lag_R_b <- ifelse((t - delayAC - delayB) <= 0, 0, lagvalue(t - delayAC - delayB)[7])
  lag_R_c <- ifelse((t - 2 * delayAC - delayB) <= 0, 0, lagvalue(t - 2 * delayAC - delayB)[7])
  
  d_alpha <- p["k"] * (y["M"] * y["R"] - lag_M * lag_R)
  d_A <- p["k"] * lag_M * lag_R
  
  d_beta <- p["k"] * (lag_M * lag_R - lag_M_b * lag_R_b)
  d_B <- p["k"] * lag_M_b * lag_R_b
  
  d_chi <- p["k"] * (lag_M_b * lag_R_b - lag_M_c * lag_R_c)
  d_C <- p["k"] * lag_M_c * lag_R_c
  
  d_R <- p["k"] * (-(y["M"] * y["R"]) + lag_M_c * lag_R_c)
  d_M <- p["k"] * (-(y["M"] * y["R"]) + lag_M_c * lag_R_c)
  
  # d_R <- p["k"] * (-(y["M"] * y["R"]))
  # d_M <- p["k"] * (-(y["M"] * y["R"]))
  
  return(list(c(d_alpha,
                d_A,
                d_beta,
                d_B,
                d_chi,
                d_C,
                d_R,
                d_M)))
  
}

parms <- c(k=1.9e-8, q=10)

times <- seq(0, 300, 0.1)

init <- c(alpha = 0, A = 0,
          beta = 0, B = 0,
          chi = 0, C = 0,
          R = 100,
          M = 1000)

delayAC = 9
delayB = 30

out <- dede(y = init, times = times, func = my_function,
            p = parms, delayAC = delayAC, delayB = delayB,
            M_init = init["M"], R_init = init["R"])

out_df <- as.data.frame(out)

out_df <- gather(out_df, variable, value, -time)

out_df1 <- filter(out_df, variable %in% c("A", "B", "C"))

out_df2 <- filter(out_df, variable %in% c("alpha", "beta", "chi"))

out_df3 <- filter(out_df, variable %in% c("R", "M"))

ggplot(out_df1, aes(x = time, y = value, group = variable, color = variable)) + geom_line()

ggplot(out_df2, aes(x = time, y = value, group = variable, color = variable)) + geom_line()

ggplot(out_df3, aes(x = time, y = value, group = variable, color = variable)) + geom_line()
