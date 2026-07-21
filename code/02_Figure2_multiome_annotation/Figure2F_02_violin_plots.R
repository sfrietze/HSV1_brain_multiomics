#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--outdir", type = "character",
              default = "outputs/02_Figure2_multiome_annotation/Figure2F_violin_plots")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(opt$input)
DefaultAssay(obj) <- "RNA"

cluster_order <- c(
  "Homeostatic Microglia",
  "Transiently Activated Microglia",
  "IFN-Responsive Microglia",
  "Primed Microglia",
  "Mitochondrial-Activated Microglia",
  "IEG-High Microglia",
  "Infiltrating Macrophages"
)

cluster_colors <- c(
  "Homeostatic Microglia" = "darkorange",
  "Transiently Activated Microglia" = "#fcae91",
  "IFN-Responsive Microglia" = "#4B0082",
  "Primed Microglia" = "#cb181d",
  "Mitochondrial-Activated Microglia" = "#99000d",
  "IEG-High Microglia" = "#fb6a4a",
  "Infiltrating Macrophages" = "#67000d"
)

obj$annotated_clusters <- factor(obj$annotated_clusters, levels = cluster_order)
obj <- obj[, !is.na(obj$annotated_clusters)]

plot_violin <- function(gene, outfile) {
  df <- FetchData(obj, vars = c(gene, "annotated_clusters"))
  colnames(df) <- c("expression", "cluster")
  df$cluster <- factor(df$cluster, levels = cluster_order)

  p <- ggplot(df, aes(x = cluster, y = expression, fill = cluster)) +
    geom_violin(scale = "width", trim = TRUE, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    theme_classic(base_size = 11) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = gene, y = "Expression")

  ggsave(
    file.path(opt$outdir, outfile),
    p,
    width = 5.5,
    height = 3,
    device = cairo_pdf
  )
}

plot_violin("Gvin1", "Figure2F_Gvin1_violin.pdf")
plot_violin("Jund", "Figure2F_Jund_violin.pdf")

sink(file.path(opt$outdir, "sessionInfo_Figure2F_02_violin_plots.txt"))
sessionInfo()
sink()
