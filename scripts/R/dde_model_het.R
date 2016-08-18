rm(list = ls())

library(deSolve)
library(dplyr)
library(tidyr)
library(cowplot)

my_function <- function(t, y, p) {
  
  Rf <- p["Rpool"] - (y["Re"] + y["Rh"])
  beta <- min(p["a_h"], (1 - y["Rh"] / (p["Nh"] * p["sh"])) * (p["sh"] / p["delay_h"]))
  lag_beta <- ifelse((t - p["delay_h"]) <= 0, 0, min(p["a_h"], (1 - lagvalue(t - p["delay_h"])[2] / (p["Nh"] * p["sh"])) * (p["sh"] / p["delay_h"])))
  lag_Rfe <- ifelse((t - p["delay_e"]) <= 0, 0, p["Rpool"] - lagvalue(t - p["delay_e"])[1] - lagvalue(t - p["delay_e"])[2])
  lag_Rfh <- ifelse((t - p["delay_h"]) <= 0, 0, p["Rpool"] - lagvalue(t - p["delay_h"])[1] - lagvalue(t - p["delay_h"])[2])
  
  # lag_Re <- ifelse((t - p["delay_e"]) <= 0, 0, lagvalue(t - p["delay_e"])[1])
  # lag_Rh <- ifelse((t - p["delay_h"]) <= 0, 0, lagvalue(t - p["delay_h"])[2])
  
  d_Re <- p["a_e"] * p["Ne"] * (Rf - lag_Rfe)
  d_Rh <- p["Nh"] * (beta * Rf - lag_beta * lag_Rfh)
  
  d_e <- p["a_e"] * p["Ne"] * lag_Rfe
  d_h <- p["a_h"] * p["Nh"] * lag_Rfh
  
  return(list(c(d_Re,
                d_Rh,
                d_e,
                d_h)))
  
}

parms <- c("Rpool" = 26300,
           "Ne" = 500,
           "Nh" = 500,
           "a_e" = 1.9e-5,
           "a_h" = 1.8e-5,
           "delay_e" = 17,
           "delay_h" = 34,
           "sh" = 38)

times <- seq(0, 300, 0.1)

init <- c("Re" = 0,
          "Rh" = 0,
          "e" = 0,
          "h" = 0)

out <- dede(y = init, times = times, func = my_function,
            p = parms)

out_df <- as.data.frame(out)

out_df <- gather(out_df, variable, value, -time)

out_df1 <- filter(out_df, variable %in% c("e", "h"))

ggplot(out_df1, aes(x = time, y = value, group = variable, color = variable)) + geom_line()
