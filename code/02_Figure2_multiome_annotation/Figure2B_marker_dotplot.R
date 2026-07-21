#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(forcats)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--outdir", type = "character",
              default = "outputs/02_Figure2_multiome_annotation/Figure2B_marker_dotplot")
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
  "Infiltrating Macrophages",
  "CNS Endothelial Cells",
  "CD8+ T Cells",
  "Cycling Myeloid Progenitors",
  "Vascular Smooth Muscle Cells"
)

new_annotations <- c(
  "0"  = "Homeostatic Microglia",
  "1"  = "Transiently Activated Microglia",
  "2"  = "IFN-Responsive Microglia",
  "3"  = "Primed Microglia",
  "4"  = "Mitochondrial-Activated Microglia",
  "5"  = "IEG-High Microglia",
  "6"  = "Infiltrating Macrophages",
  "7"  = "CNS Endothelial Cells",
  "8"  = "CD8+ T Cells",
  "9"  = "Cycling Myeloid Progenitors",
  "10" = "Vascular Smooth Muscle Cells"
)

if ("seurat_clusters" %in% colnames(obj@meta.data)) {
  obj$annotated_clusters <- plyr::mapvalues(
    as.character(obj$seurat_clusters),
    from = names(new_annotations),
    to = unname(new_annotations)
  )
}

obj$annotated_clusters <- factor(obj$annotated_clusters, levels = cluster_order)
Idents(obj) <- "annotated_clusters"

gene_order <- c(
  "P2ry12", "Fcrls", "Gpr34",
  "Rgs1", "Ier5", "Sox4",
  "Ccl12", "Gm4951", "Ifi204",
  "Cst3", "Hexb", "Tanc2",
  "mt-Cytb", "mt-Nd1", "mt-Nd2",
  "Egr1", "Fosb", "Gm34455",
  "F13a1", "Cd163", "Ms4a4a",
  "Myrip", "Shank3", "Dipk2b",
  "Ncr1", "Gimap3", "Fasl",
  "Knl1", "Neil3", "Esco2",
  "Slc6a20a", "Ror1", "Carmn"
)

missing_genes <- setdiff(gene_order, rownames(obj))
write_csv(
  tibble(missing_gene = missing_genes),
  file.path(opt$outdir, "Figure2B_missing_genes.csv")
)

gene_order <- intersect(gene_order, rownames(obj))

cell_counts <- obj@meta.data %>%
  count(annotated_clusters, name = "n_cells") %>%
  mutate(annotated_clusters = factor(annotated_clusters, levels = cluster_order))

write_csv(
  cell_counts,
  file.path(opt$outdir, "Figure2B_cluster_cell_counts.csv")
)

dp <- DotPlot(
  obj,
  features = gene_order,
  assay = "RNA",
  dot.scale = 6
)$data

dot_df <- dp %>%
  group_by(features.plot) %>%
  mutate(
    avg.exp.mm = scales::rescale(avg.exp, to = c(0, 1)),
    pct.exp.plot = ifelse(avg.exp < 0.10 | pct.exp < 5, 0, pct.exp)
  ) %>%
  ungroup() %>%
  filter(pct.exp.plot > 0) %>%
  mutate(
    features.plot = factor(features.plot, levels = gene_order),
    id = factor(id, levels = cluster_order)
  )

write_csv(
  dot_df,
  file.path(opt$outdir, "Figure2B_dotplot_data_used.csv")
)

dot_plot <- ggplot(dot_df, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp.plot, color = avg.exp.mm)) +
  scale_color_gradientn(
    colors = c("blue", "grey90", "red"),
    limits = c(0, 1),
    name = "Scaled expression"
  ) +
  scale_size(
    range = c(0, 8),
    limits = c(0, 100),
    breaks = c(25, 50, 75, 100),
    name = "Percent expressed"
  ) +
  scale_y_discrete(limits = rev(cluster_order)) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right"
  )

cluster_colors <- c(
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

bar_plot <- ggplot(
  cell_counts,
  aes(x = n_cells, y = fct_rev(annotated_clusters), fill = annotated_clusters)
) +
  geom_col(width = 0.65, color = "black", linewidth = 0.25) +
  geom_text(aes(label = n_cells), hjust = -0.1, size = 2.8) +
  scale_fill_manual(values = cluster_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

ggsave(
  file.path(opt$outdir, "Figure2B_marker_dotplot.pdf"),
  dot_plot,
  width = 12,
  height = 5,
  useDingbats = FALSE
)

final_plot <- dot_plot + bar_plot + plot_layout(widths = c(5.5, 0.9))

ggsave(
  file.path(opt$outdir, "Figure2B_marker_dotplot_with_cell_counts.pdf"),
  final_plot,
  width = 12,
  height = 6,
  useDingbats = FALSE
)

sink(file.path(opt$outdir, "sessionInfo_Figure2B_marker_dotplot.txt"))
sessionInfo()
sink()
