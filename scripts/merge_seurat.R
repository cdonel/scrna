library(Seurat)
library(stringr)
library(argparse)


parser <- ArgumentParser()
parser$add_argument("-i", "--input_dir", required = TRUE)
parser$add_argument("-o", "--out_dir", required = TRUE)
args <- parser$parse_args()


sra_ids <- readLines(args$input_dir)


seurat_samples <- c()


for (sample in sra_ids) {
     seurat_samples <- append(seurat_samples,
                              readRDS(str_glue("../data/seurat/seurat_objs/{sample}/{sample}.norm.rds")))
}




merged_seurat <- SeuratObject:::merge.Seurat(x = seurat_samples[[1]],
                                             y = seurat_samples[-1],
                                             add.cell.ids = sra_ids,
                                             project = "test")


# Normalise merged object
merged_seurat <- NormalizeData(merged_seurat,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)


# Merged feature selection
merged_seurat <- FindVariableFeatures(merged_seurat,
                                      election.method = "vst",
                                      nfeatures = 2000)


# Scale merged data
all.genes <- rownames(merged_seurat)
merged_seurat <- ScaleData(merged_seurat, features = all.genes)


saveRDS(merged_seurat,
        file = args$out_dir)