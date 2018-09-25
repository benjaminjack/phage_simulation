make_plot <- function(file) {
  data <- read_tsv(file)
  data %>% filter(species %in% c("proteinX", "proteinY", "proteinZ")) -> data
  plot1 <- ggplot(data, aes(x = time, y = transcript, color=species)) + 
    geom_line() + 
    theme(legend.position="none")
  plot2 <- ggplot(data %>% filter(time > 499), aes(fill=species, x=species, y=transcript)) + 
    geom_bar(stat="identity")
  plot_grid(plot1, plot2, ncol=2, rel_widths = c(1, 1.5), align = "h", labels="AUTO")
}

make_plot2 <- function(file) {
  data <- read_tsv(file)
  data %>% filter(species %in% c("proteinX", "proteinY", "proteinZ")) -> data
  plot1 <- ggplot(data, aes(x = time, y = transcript, color=species)) + 
    geom_line() + 
    theme(legend.position="none")
  plot3 <- ggplot(data %>% filter(time > 499), aes(fill=species, x=species, y=transcript)) + 
    geom_bar(stat="identity")
  plot2 <- ggplot(data %>% filter(time > 100, time < 101), aes(fill=species, x=species, y=transcript)) + 
    geom_bar(stat="identity") + theme(legend.position="none")
  plot_grid(plot1, plot2, plot3, ncol=3, rel_widths = c(1, 1, 1.3), align = "h", labels="AUTO")
}

three_genes_no_deg <- make_plot("../data/simulation/three_genes/three_genes_no_deg_counts.tsv")
save_plot("three_genes_no_deg.pdf", three_genes_no_deg, ncol = 2, base_width = 4)

three_genes_w_deg <- make_plot("../data/simulation/three_genes/three_genes_w_deg_counts.tsv")
save_plot("three_genes_w_deg.pdf", three_genes_w_deg, ncol = 2, base_width = 4)

three_genes_prom <- make_plot("../data/simulation/three_genes/three_genes_prom_counts.tsv")
save_plot("three_genes_prom.pdf", three_genes_prom, ncol = 2, base_width = 4)

three_genes_prom_deg <- make_plot("../data/simulation/three_genes/three_genes_prom_deg_counts.tsv")
save_plot("three_genes_prom_deg.pdf", three_genes_prom_deg, ncol = 2, base_width = 4)

three_genes_rnase_site <- make_plot("../data/simulation/three_genes/three_genes_rnase_site_counts.tsv")
save_plot("three_genes_rnase_site.pdf", three_genes_rnase_site, ncol = 2, base_width = 4)

three_genes_rnase_prom <- make_plot2("../data/simulation/three_genes/three_genes_rnase_prom_counts.tsv")
save_plot("three_genes_rnase_prom.pdf", three_genes_rnase_prom, ncol = 2, base_width = 5.5)
