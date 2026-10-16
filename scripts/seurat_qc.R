library(Seurat)
library("ggplot2")
library(argparse)


parser <- ArgumentParser()
parser$add_argument("-i", "--input_dir", required = TRUE)
parser$add_argument("-o", "--out_dir", required = TRUE)
args <- parser$parse_args()


seurat_sample <- readRDS(args$input_dir)


png()


plot <- VlnPlot(seurat_sample,
                features = c("nFeature_RNA",
                             "nCount_RNA",
                             "percent.mt"),
                ncol = 3)


ggsave(args$out_dir, width = 15, height = 15)


dev.off()