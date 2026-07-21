#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

option_list <- list(
  make_option("--input", type = "character", help = "Input annotated Seurat RDS"),
  make_option("--outdir", type = "character", default = "outputs/02_Figure2_multiome_annotation/Figure2A_RNA_UMAP")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(opt$input)

myeloid_clusters <- c(
  "Homeostatic Microglia",
  "Transiently Activated Microglia",
  "IFN-Responsive Microglia",
  "Primed Microglia",
  "Mitochondrial-Activated Microglia",
  "IEG-High Microglia",
  "Infiltrating Macrophages",
  "Cycling Myeloid Progenitors"
)

obj$broad_celltype <- ifelse(
  obj$annotated_clusters %in% myeloid_clusters,
  "Myeloid",
  "Non-myeloid"
)

broad_colors <- c(
  "Myeloid" = "#d62728",
  "Non-myeloid" = "#9ecae1"
)

granular_colors <- c(
  "Homeostatic Microglia" = "darkorange",
  "Transiently Activated Microglia" = "#fcae91",
  "IFN-Responsive Microglia" = "#4B0082",
  "Primed Microglia" = "#cb181d",
  "Mitochondrial-Activated Microglia" = "#99000d",
  "IEG-High Microglia" = "#fb6a4a",
  "Infiltrating Macrophages" = "#67000d",
  "Cycling Myeloid Progenitors" = "#008080",
  "CD8+ T Cells" = "lightblue",
  "CNS Endothelial Cells" = "gainsboro",
  "Vascular Smooth Muscle Cells" = "darkgray"
)

p_broad <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "broad_celltype",
  cols = broad_colors,
  pt.size = 1.1
) +
  theme_void() +
  labs(title = "")

ggsave(
  file.path(opt$outdir, "Figure2A_RNA_UMAP_broad.pdf"),
  p_broad,
  width = 5,
  height = 4
)

p_granular <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "annotated_clusters",
  cols = granular_colors,
  pt.size = 1.1
) +
  theme_void() +
  labs(title = "")

ggsave(
  file.path(opt$outdir, "Figure2A_RNA_UMAP_granular.pdf"),
  p_granular,
  width = 9,
  height = 5
)

write.csv(
  obj@meta.data[, c("seurat_clusters", "annotated_clusters", "broad_celltype", "condition")],
  file.path(opt$outdir, "Figure2A_RNA_UMAP_metadata_used.csv")
)

sink(file.path(opt$outdir, "sessionInfo_Figure2A_RNA_UMAP.txt"))
sessionInfo()
sink()
