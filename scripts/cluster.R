library(Seurat)
library(ggplot2)
library(stringr)
library(argparse)


parser <- ArgumentParser()
parser$add_argument("-i", "--input_dir", required = TRUE)
args <- parser$parse_args()


# Read in merged and integrated seurat object.
merged_seurat <- readRDS(args$input_dir)


# PCA merged object
merged_seurat <- RunPCA(merged_seurat,
                        features = VariableFeatures(object = merged_seurat))


# Integrate
intg_seurat <- IntegrateLayers(
  object = merged_seurat,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  assay = "RNA",
  verbose = TRUE
)


p1 <- ElbowPlot(intg_seurat)


ggsave("../results/elbow_plot.png",
       plot = p1,
       bg = "white")


# Using elbow plot determine number of dimensions
# used in finding nearest neighbors (NNs).
intg_seurat <- FindNeighbors(intg_seurat,
                             reduction = "integrated.cca",
                             dims = 1:10)


# Find clusters within NNs. Set resolution, usually between 0.5 and 1.2.
# Higher sample size = higher resolution.
intg_seurat <- FindClusters(intg_seurat,
                            resolution = 1,
                            cluster.name = "cca_clusters")


# Generate a UMAP. Dimensions used based on elbow plot.
intg_seurat <- RunUMAP(intg_seurat,
                       reduction = "integrated.cca",
                       reduction.name = "umap.cca",
                       dims = 1:10)


# Plot and save UMAP.
p2 <- DimPlot(intg_seurat,
              reduction = "umap.cca",
              group.by = c("cca_clusters"),
              combine = FALSE,
              label.size = 2)


ggsave("../results/umap_cluster.png",
       plot = p2,
       bg = "white")


merged_seurat <- JoinLayers(intg_seurat)


# Find markers in clusters.
merged_seurat.markers <- FindAllMarkers(merged_seurat, only.pos = TRUE)


write.csv(merged_seurat.markers, "../results/clusterfeatures.csv")