library(readr)
library(dplyr)
library(ggplot2)
library(broom)

# Burst size and lysis time estimates

prots <- read_csv("rep123_combined.csv", col_types = list(strain = col_character()))

prots <- prots %>% filter(org=="phage", tabasco_id == "gp10A", b_rep != 1) %>%
  mutate(strain = str_replace_all(strain, c("11-44" = "atten", "11-42" = "evol", "11-46" = "wt")))

# model <- lm(log(area_norm) ~ time, data = prots)

wildtype <- function(x)2^(43.2*(x/3600))
attenuated <- function(x)2^(35.7*(x/3600))
evolved <- function(x)2^(38.7*(x/3600))
wt_atten_ratio <- function(x)(2^(43.2*(x/3600)))/2^(35.7*(x/3600))
wt_evol_ratio <- function(x)(2^(43.2*(x/3600)))/2^(38.7*(x/3600))

ggplot(data.frame(x=c(0,540)), aes(x)) + 
  stat_function(fun=wildtype, geom="line", color="red") + 
  stat_function(fun=attenuated, geom="line", color="green") + 
  geom_hline(y=wildtype(540)) + geom_hline(y=attenuated(540)) + 
  geom_vline(x=540) +
  geom_text(aes(x=250, y=60, label=wildtype(540)/attenuated(540))) +
  ylab("phage particles") +
  xlab("time (min)")

prots <- prots %>% filter(org=="phage", tabasco_id == "gp10A", b_rep != 1)

prots2 <- group_by(prots, strain, time) %>% summarize(area_mean = mean(area_norm))

ggplot(prots2, aes(x=time, y=area_mean, group=strain, color=strain)) + geom_point() + geom_line()

prots3 <- group_by(prots, time, b_rep) %>% 
  summarize(wt_atten = area_norm[strain=="wt"]/area_norm[strain=="atten"], wt_evol = area_norm[strain=="wt"]/area_norm[strain=="evol"]) %>%
  gather(strain, ratio, -time, -b_rep)

prots3 <- group_by(prots2, time) %>% 
  summarize(wt_atten = area_mean[strain=="wt"]/area_mean[strain=="atten"], wt_evol = area_mean[strain=="wt"]/area_mean[strain=="evol"]) %>%
  gather(strain, ratio, -time)

ggplot(prots3, aes(x=time, y=ratio, group=strain, color=strain)) + 
  stat_function(fun=wt_atten_ratio, geom="line", color="red") +
  stat_function(fun=wt_evol_ratio, geom="line", color="blue") +
  geom_point() + xlim(0,600)
