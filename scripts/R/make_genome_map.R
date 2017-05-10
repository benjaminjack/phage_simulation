library(tidyverse)
library(cowplot)
library(ggrepel)
library(stringr)

genes <- read_csv("./data/T7_gene_coords.csv")
id_map <- read_csv("./data/proteomics/id_map.csv")

make_transcript <- function(order, start_gene, stop_gene, df) {
  out_df <- filter(df, gene_number == start_gene | gene_number == stop_gene)
  data_frame("pos" = order, "start_transcript" = min(out_df$start), "stop_transcript" = max(out_df$stop))
}

genes <- left_join(id_map, genes, by = c("accession" = "protein_id")) %>% filter(!is.na(class))

transcript_coords <- bind_rows(
  make_transcript(1, 1, 11, genes),
  make_transcript(2, 2,11, genes),
  make_transcript(3, 5,11, genes)
)

gene_map <- ggplot(data = genes, aes(xmin = start, xmax = stop, ymin=0, ymax=1, fill=factor(class))) +
  # geom_segment(data = transcript_coords, aes(x = start_transcript, xend = stop_transcript, y = -pos, yend = -pos)) + 
  geom_text_repel(data = filter(genes, !str_detect(gene_number, "\\.")), 
                  aes(x = ((stop-start)/2 + start), y = 1, label = gene_number), nudge_y = 0.6, size = 4, segment.color="grey80") +
  geom_rect(color = "black", size = 0.3) +
  ylim(0,2) + 
  theme_nothing() + 
  theme(legend.position = "none")

save_plot("figures/gene_map.pdf", gene_map, base_width = 10, base_height = 1)
