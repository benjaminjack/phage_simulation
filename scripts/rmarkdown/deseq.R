library(tidyverse)
library(ggplot2)
library(gtools)
library(stringr)
library(cowplot)
library(DESeq2)

rep1 <- read_tsv("gene_counts_rep1.gff", col_names=FALSE) %>%
  filter(X3 == "CDS") %>%
  select("start" = X4, "end" = X5, "gene" = X9, "count" = X10) %>%
  mutate(gene_name = str_match(gene, "Note=([^;|%]*)")[,2],
         length = end - start)
rep1_filter <- filter(rep1, str_detect(gene, "cds*")) %>% 
  mutate(gene_name = str_match(gene, "Note=([^;|%]*)")[,2])

ggplot(rep1_filter, aes(x = factor(gene_name, levels=unique(mixedsort(gene_name))), y = count)) + 
  geom_bar(stat="identity")

rep1_df <- as.data.frame(rep1_filter %>% select(count, gene_name) %>% distinct(gene_name, .keep_all=TRUE))
row.names(rep1_df) <- rep1_df$gene_name
rep1_df <- select(rep1_df, -gene_name)

cts <- as.matrix(rep1_df %>% select(-gene_name))
coldata <- data.frame(condition = c("untreated"), type = c("single-read"), row.names = c("count"))

dds <- DESeqDataSetFromMatrix(cts, colData = coldata, design= ~1)
dds <- estimateSizeFactors(dds)

rep1_norm <- rownames_to_column(as.data.frame(counts(dds, normalized=TRUE)), var="gene_name")

ggplot(rep1_norm, aes(x = factor(gene_name, levels=unique(mixedsort(gene_name))), y = count)) + 
  geom_bar(stat="identity")


