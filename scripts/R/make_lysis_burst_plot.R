library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Burst size plot and analysis

burst <- read_csv("data/burst.csv")

strain_order <- c("11_44  ( /ul)", "11_42 ( /ul)", "11_46  ( /ul)")
strain_labels <- c("recoded", "evolved", "wildtype")

burst_plot <- ggplot(burst, aes(x = strain, y = burst_size)) +
  stat_summary(fun.y = "mean", geom = "bar", fill="grey70") +
  geom_line(aes(group=rep)) +
  geom_point() +
  scale_x_discrete(name = "strain", limits = strain_order, labels = strain_labels) +
  ylab("burst size")

# Do some t-tests
burst_wide <- dplyr::select(burst, rep, strain, burst_size) %>%
  spread(strain, burst_size)
# Evolved v recoded
t.test(burst_wide$`11_42 ( /ul)`, burst_wide$`11_44  ( /ul)`, paired = T)
# Evolved v wildtype
t.test(burst_wide$`11_42 ( /ul)`, burst_wide$`11_46  ( /ul)`, paired = T)
# Recoded v wildtype
t.test(burst_wide$`11_44  ( /ul)`, burst_wide$`11_46  ( /ul)`, paired = T)

# Lysis time plot and analysis

lysis <- read_csv('data/lysis_clean.csv',
                  col_types = list(Date = col_date("%m/%d/%y")))

lysis <- mutate(lysis, `11-42` = ifelse(!is.na(Volume), `11-42`/Volume, `11-42`))
lysis <- mutate(lysis, `11-44` = ifelse(!is.na(Volume), `11-44`/Volume, `11-44`))
lysis <- mutate(lysis, `11-46` = ifelse(!is.na(Volume), `11-46`/Volume, `11-46`))

lysis <- gather(lysis, strain, concentration, `11-42`:`11-46`)

lysis_plot <- ggplot(lysis, aes(x = Time, y = concentration, group = Date, color = factor(Date))) +
  geom_line() +
  facet_wrap( ~ strain) +
  geom_point() +
  panel_border()

extract_midpoint <- function(model) {
  model$coefficients["e:(Intercept)"]
}

lysis_mid <- filter(lysis, !is.na(concentration)) %>%
  group_by(strain, Date) %>%
  nest() %>%
  mutate(mod = purrr::map(data, ~ drm(concentration ~ Time, data = ., fct = LL.3()))) %>%
  mutate(inflection = purrr::map(mod, extract_midpoint)) %>%
  dplyr::select(strain, Date, inflection) %>%
  unnest()

mod <- lm(data = lysis_mid, formula = inflection ~ strain)
anova(mod)

strain_order <- c("11-44", "11-42", "11-46")
strain_labels <- c("recoded", "evolved", "wildtype")

lysis_plot_mid <- ggplot(lysis_mid, aes(x = strain, y = inflection, group = strain)) +
  stat_summary(geom="bar", fun.y = "mean", fill = "grey70") +
  geom_point() +
  scale_x_discrete(name = "strain", limits = strain_order, labels = strain_labels) +
  ylab("lysis time (minutes)")

burst_lysis <- plot_grid(burst_plot, lysis_plot_mid, labels = c("A", "B"))

save_plot("./figures/burst_lysis.pdf", burst_lysis, base_aspect_ratio = 2)
