#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(ggrastr)
  library(ggrepel)
  library(patchwork)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--outdir", type = "character", default = "outputs/02_Figure2_multiome_annotation/Figure2C_ATAC_LSI_UMAP")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(opt$input)
DefaultAssay(obj) <- "ATAC"

cluster_order <- c(
  "Homeostatic Microglia",
  "Transiently Activated Microglia",
  "IFN-Responsive Microglia",
  "Primed Microglia",
  "Mitochondrial-Activated Microglia",
  "IEG-High Microglia",
  "Infiltrating Macrophages",
  "CNS Endothelial Cells",
  "CD8+ T Cells",
  "Cycling Myeloid Progenitors",
  "Vascular Smooth Muscle Cells"
)

obj$annotated_clusters <- factor(obj$annotated_clusters, levels = cluster_order)
obj <- obj[, !is.na(obj$annotated_clusters)]

obj@reductions$umap.atac <- NULL

set.seed(123)
obj <- RunUMAP(
  object = obj,
  reduction = "lsi",
  dims = 2:30,
  reduction.name = "umap.atac",
  umap.method = "uwot",
  metric = "cosine"
)

umap_df <- Embeddings(obj, reduction = "umap.atac") %>%
  as.data.frame()

colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cluster <- factor(obj$annotated_clusters, levels = cluster_order)

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

umap_df$myeloid_binary <- ifelse(
  umap_df$cluster %in% myeloid_clusters,
  "Myeloid",
  "Non-Myeloid"
)

umap_df$myeloid_binary <- factor(
  umap_df$myeloid_binary,
  levels = c("Non-Myeloid", "Myeloid")
)

binary_palette <- c(
  "Non-Myeloid" = "#9ecae1",
  "Myeloid" = "#d62728"
)

label_df <- umap_df %>%
  group_by(cluster) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  ) %>%
  mutate(
    label = recode(
      as.character(cluster),
      "Transiently Activated Microglia" = "Transiently\nActivated",
      "IFN-Responsive Microglia" = "IFN-Responsive",
      "Mitochondrial-Activated Microglia" = "Mitochondrial\nActivated",
      "IEG-High Microglia" = "IEG-High",
      "Homeostatic Microglia" = "Homeostatic",
      "Infiltrating Macrophages" = "Infiltrating\nMacrophages",
      "Cycling Myeloid Progenitors" = "Myeloid\nprogenitor",
      "CD8+ T Cells" = "CD8+ T cells",
      "CNS Endothelial Cells" = "CNS Endothelial Cells",
      "Vascular Smooth Muscle Cells" = "Vascular Smooth\nMuscle Cells"
    )
  )

xlims <- range(umap_df$UMAP_1)
ylims <- range(umap_df$UMAP_2)

p_binary <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = myeloid_binary)) +
  geom_point_rast(size = 1.4, alpha = 0.9) +
  scale_color_manual(values = binary_palette) +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_void(base_size = 12) +
  theme(legend.position = "none")

p_granular <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point_rast(size = 1.4, alpha = 0.9) +
  scale_color_manual(values = granular_palette) +
  geom_text_repel(
    data = label_df,
    aes(label = label),
    color = "black",
    size = 3.2,
    fontface = "bold",
    segment.color = NA,
    show.legend = FALSE
  ) +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_void(base_size = 12) +
  theme(legend.position = "none")

combined_plot <- p_binary + p_granular + plot_layout(widths = c(1, 1))

ggsave(
  file.path(opt$outdir, "Figure2C_ATAC_LSI_UMAP_binary.pdf"),
  p_binary,
  width = 5.5,
  height = 5,
  device = cairo_pdf
)

ggsave(
  file.path(opt$outdir, "Figure2C_ATAC_LSI_UMAP_granular.pdf"),
  p_granular,
  width = 5.5,
  height = 5,
  device = cairo_pdf
)

ggsave(
  file.path(opt$outdir, "Figure2C_ATAC_LSI_UMAP_combined.pdf"),
  combined_plot,
  width = 10,
  height = 5,
  device = cairo_pdf
)

write.csv(
  umap_df,
  file.path(opt$outdir, "Figure2C_ATAC_LSI_UMAP_coordinates.csv"),
  row.names = FALSE
)

sink(file.path(opt$outdir, "sessionInfo_Figure2C_ATAC_LSI_UMAP.txt"))
sessionInfo()
sink()
