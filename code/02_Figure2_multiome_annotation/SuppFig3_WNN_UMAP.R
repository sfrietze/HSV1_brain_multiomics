#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

option_list <- list(
  make_option("--input_rds", type = "character"),
  make_option(
    "--outdir",
    type = "character",
    default = "outputs/02_Figure2_multiome_annotation"
  ),
  make_option(
    "--group_by",
    type = "character",
    default = "annotated_clusters"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

granular_palette <- c(
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


obj <- readRDS(opt$input_rds)

message("Assays: ", paste(Assays(obj), collapse = ", "))
message("Reductions before WNN: ", paste(Reductions(obj), collapse = ", "))

if (!all(c("pca", "lsi") %in% Reductions(obj))) {
  stop("The object must contain both pca and lsi reductions.")
}

if (!opt$group_by %in% colnames(obj@meta.data)) {
  stop(
    "Metadata column '", opt$group_by, "' not found. Available columns: ",
    paste(colnames(obj@meta.data), collapse = ", ")
  )
}

obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:10, 2:30),
  modality.weight.name = c("CombinedRNA.weight", "ATAC.weight")
)

obj <- RunUMAP(
  object = obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  seed.use = 123
)

p <- DimPlot(
  object = obj,
  reduction = "wnn.umap",
  group.by = opt$group_by,
  label = TRUE,
  repel = TRUE
) +
  scale_color_manual(values = granular_palette, drop = FALSE) +
  labs(
    title = NULL,
    x = "WNN UMAP 1",
    y = "WNN UMAP 2"
  ) +
  theme_classic()

ggsave(
  filename = file.path(opt$outdir, "SuppFig3_WNN_UMAP.pdf"),
  plot = p,
  width = 7,
  height = 6
)

ggsave(
  filename = file.path(opt$outdir, "SuppFig3_WNN_UMAP.png"),
  plot = p,
  width = 7,
  height = 6,
  dpi = 300
)


sink(file.path(opt$outdir, "sessionInfo_SuppFig3_WNN_UMAP.txt"))
sessionInfo()
sink()
