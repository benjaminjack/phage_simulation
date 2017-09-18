require(tidyverse)

cbPalette <- c("#555555", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sim_data <- read.csv("../runs/post_processed/T7_082017_all_counts.csv")
sim_data$class <- factor(sim_data$class)
 
make_plot <- function(t) {
  p <- ggplot(filter(sim_data, time == t, gene_number != 1, weight == 3, condition=="wt"),
         aes(y = log10(proteins), x = log10(transcripts), color=class)) +
    geom_point() +
    scale_color_manual("class", values = cbPalette[c(3, 4, 7)], drop = FALSE) +
    xlim(0, 4) + ylim(0, 5) +
    labs(title = paste0("time t = ", t, "s")) +
    theme_minimal() +
    theme(legend.position = c(.1, .8),
          legend.background = element_rect(fill = "white", color = "white"))
  
  #prot <- function(c) simulate_cvd(c, protan_transform(1))
  #edit_colors(p, prot)
  p
}

# function for saving images
plot_save<-function(t = 1200){
  # add 50000 to index so images are in the right order (10 comes after 9)
  file_path = paste0("plot-", 50000+t, ".png")
  ggsave(filename=file_path, width=5, height=3.75, dpi=300, make_plot(t))
}

system("rm plot-5*.png plot.gif")

map(5*(1:280), plot_save)

system("magick convert -delay 1 -loop 0 *.png plot.gif")
