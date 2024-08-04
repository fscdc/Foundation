library(GSVA)
library(tidyr)
library(zellkonverter)
# library(Seurat)
library(GSVA)
library(GSEABase)
library(SingleCellExperiment)
library(pheatmap)
library(ggplot2)
library(limma)
library(stringr)

folder_path <- "/home/foundation/program/Foundation/record/temp-h5ad"

# get all .txt files in the folder
h5ad_files <- list.files(path = folder_path, pattern = "\\.h5ad$", full.names = TRUE)

s.sets <- getGmt("./c5.go.bp.v7.5.1.symbols.gmt")
options(timeout = 300)

for (file in h5ad_files) {

  file_name <- basename(file)

  pattern <- "^result-(.*)-(.*)-(.*)-(.*)-(.*).txt$"
  matches <- str_match(file_name, pattern)

  dataset <- matches[1,2]
  genelist_name <- matches[1,3]
  cluster_method <- matches[1,4]
  resolution <- matches[1,5]
  embedding_method <- matches[1,6]

  adata <- readH5AD(file)
  X_matrix <- assay(adata, "X")
  celltype <- colData(adata)$celltype
  cluster <- colData(adata)$louvain

  gsva_data <- gsva(X_matrix, s.sets, method = "gsva",kcdf="Gaussian",parallel.sz = 50)

  colormap <- colorRampPalette(c("#4575B4","white","#D73027"))(100)
  breaks = seq(-1, 1, length.out=100)
  gsva_data_used <- gsva_data[c(1:10),]

  average_scores_celltype <- apply(gsva_data_used, 1, function(row) {
    tapply(row, celltype, mean)
  })
  average_scores_cluster <- apply(gsva_data_used, 1, function(row) {
    tapply(row, cluster, mean)
  })

  gsva_figure <- pheatmap(average_scores_celltype, 
          cellheight = 8, cellwidth = 12,
          border = "white",
          fontsize = 10,
          main = "Celltype GSVA",
          show_rownames = T,
          scale = "row", # column
          cluster_cols = F,
          cluster_rows = T,
          color = colormap, breaks = breaks)
  gsva_path <- paste0("/home/foundation/program/Foundation/record/figures/enrich/",dataset, "-", genelist, "-", cluster_method, "-", resolution, "-", embedding_method, "-GSVA-Celltype.pdf")

  ggsave(gsva_path, plot = gsva_figure, width = 30, height = 20)

  gsva_figure <- pheatmap(average_scores_cluster, 
          cellheight = 8, cellwidth = 12,
          border = "white",
          fontsize = 10,
          main = "Cluster GSVA",
          show_rownames = T,
          scale = "row", # column
          cluster_cols = F,
          cluster_rows = T,
          color = colormap, breaks = breaks)
  gsva_path <- paste0("/home/foundation/program/Foundation/record/figures/enrich/",dataset, "-", genelist, "-", cluster_method, "-", resolution, "-", embedding_method, "-GSVA-cluster.pdf")

  ggsave(gsva_path, plot = gsva_figure, width = 30, height = 20)

}