library(argparse)
library(ggplot2)
library(Matrix)
library(dplyr)
library(scico)


parser <- ArgumentParser()
parser$add_argument("-i", "--input_dir", required = TRUE)
parser$add_argument("-n", "--sample_name", required = TRUE)
parser$add_argument("-o", "--out_dir", required = TRUE)
args <- parser$parse_args()


# Read in the count matrix that was output by `kb`.
kb_mtx <- readMM(args$input_dir)


umi_x_gene <- tibble(nCount = rowSums(kb_mtx),
                     nGene = rowSums(kb_mtx > 0))



plot <- ggplot(umi_x_gene, aes(nCount, nGene)) +
  geom_bin2d(bins = 50) +
  scale_fill_scico(palette = "devon", direction = -1, end = 0.95) +
  scale_x_log10() + scale_y_log10() + annotation_logticks() +
  labs(x = "Total UMI counts",
       y = "Number of genes detected",
       title = args$sample_name) +
  theme(plot.title = element_text(hjust = 0.5))


ggsave(args$out_dir)