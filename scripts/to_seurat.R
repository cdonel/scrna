library(Seurat)
library(SeuratObject)
library(stringr)
library(argparse)


parser <- ArgumentParser()
parser$add_argument("-i", "--input_adata", required = TRUE)
parser$add_argument("-g", "--gene_ids", required = TRUE)
parser$add_argument("-n", "--sample_name", required = TRUE)
parser$add_argument("-o", "--out_dir", required = TRUE)
args <- parser$parse_args()


#Read gene names.
gene_ids <- read.delim(args$gene_ids, col.names="gene_ids") # nolint


#Read h5ad in using anndata.
#Var_names = gene ids
#Obs_names = cell barcode
#X = gene counts
h5ad <- anndata::read_h5ad(args$input_adata) # nolint


#Convert Anndata to dataframe.
sample_df <- h5ad$to_df()


#Convert EMBL ids to gene names.
colnames(sample_df) <- gene_ids$gene_ids


#Transpose
sample_df_t <- as.data.frame(t(sample_df))


#Create Seurat objects
seurat_sample <- CreateSeuratObject(counts = as.matrix(sample_df_t),
                                    min.cells = 3, min.features = 200,
                                    project = args$sample_name)


# Add metadata
seurat_sample[["percent.mt"]] <- PercentageFeatureSet(seurat_sample, pattern = "^MT-") # nolint


# Cell selection criteria
seurat_sample <- subset(seurat_sample,
                        subset = nFeature_RNA > 200 &
                          nFeature_RNA < 2500 & percent.mt < 5)


# Normalize
seurat_sample <- NormalizeData(seurat_sample,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)


# Feature selection
seurat_sample <- FindVariableFeatures(seurat_sample,
                                      selection.method = "vst",
                                      nfeatures = 2000)


# Scale data
all.genes <- rownames(seurat_sample)
seurat_sample <- ScaleData(seurat_sample, features = all.genes)


saveRDS(seurat_sample,
        file = args$out_dir)