#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Signac)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--outdir", type = "character",
              default = "outputs/02_Figure2_multiome_annotation/Figure2D_GEX_ATAC_heatmap"),
  make_option("--link_p", type = "double", default = 1e-3),
  make_option("--link_z", type = "double", default = 2),
  make_option("--top_n", type = "integer", default = 400)
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

obj <- readRDS(opt$input)

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
Idents(obj) <- "annotated_clusters"

DefaultAssay(obj) <- "ATAC"

links_df <- as.data.table(as.data.frame(Links(obj[["ATAC"]])))
links_df[, peak := gsub(":", "-", peak)]

links_df <- links_df[
  pvalue < opt$link_p &
    zscore > opt$link_z &
    !is.na(gene) &
    !is.na(peak)
]

links_df <- unique(links_df[, .(
  peak,
  gene,
  zscore,
  pvalue
)])

links_df <- links_df[
  peak %in% rownames(obj[["ATAC"]]) &
    gene %in% rownames(obj[["RNA"]])
]

stopifnot(nrow(links_df) > 0)

average_by_cluster <- function(mat, clusters, cluster_levels) {
  clusters <- factor(clusters, levels = cluster_levels)
  design <- sparse.model.matrix(~ 0 + clusters)
  colnames(design) <- gsub("^clusters", "", colnames(design))

  avg <- mat %*% design
  n_cells <- Matrix::colSums(design)
  avg <- t(t(avg) / n_cells[colnames(avg)])
  avg <- as.matrix(avg)
  avg[, cluster_levels, drop = FALSE]
}

peak_ids <- unique(links_df$peak)
gene_ids <- unique(links_df$gene)

atac_mat <- GetAssayData(obj, assay = "ATAC", slot = "data")[peak_ids, , drop = FALSE]
rna_mat  <- GetAssayData(obj, assay = "RNA",  slot = "data")[gene_ids, , drop = FALSE]

atac_avg <- average_by_cluster(
  mat = atac_mat,
  clusters = obj$annotated_clusters,
  cluster_levels = cluster_order
)

rna_avg <- average_by_cluster(
  mat = rna_mat,
  clusters = obj$annotated_clusters,
  cluster_levels = cluster_order
)

access_avg <- as.data.table(atac_avg, keep.rownames = "peak")
access_avg <- melt(
  access_avg,
  id.vars = "peak",
  variable.name = "cluster",
  value.name = "access"
)
access_avg[, cluster := as.character(cluster)]

access_avg <- merge(
  access_avg,
  links_df[, .(peak, gene, zscore, pvalue)],
  by = "peak",
  allow.cartesian = TRUE
)

rna_avg_dt <- as.data.table(rna_avg, keep.rownames = "gene")
rna_avg_dt <- melt(
  rna_avg_dt,
  id.vars = "gene",
  variable.name = "cluster",
  value.name = "expr"
)
rna_avg_dt[, cluster := as.character(cluster)]

plot_dt <- merge(
  access_avg,
  rna_avg_dt,
  by = c("gene", "cluster")
)

plot_dt[, cluster := factor(cluster, levels = cluster_order)]
plot_dt[, row_id := paste(gene, peak, sep = " | ")]

plot_dt[, rna_scaled := {
  m <- max(expr, na.rm = TRUE)
  if (!is.finite(m) || m == 0) 0 else expr / m
}, by = row_id]

plot_dt[, atac_scaled := {
  m <- max(access, na.rm = TRUE)
  if (!is.finite(m) || m == 0) 0 else access / m
}, by = row_id]

cor_per_row <- plot_dt[
  ,
  .(
    row_cor = if (
      length(unique(rna_scaled)) > 1 &&
        length(unique(atac_scaled)) > 1
    ) cor(atac_scaled, rna_scaled, method = "pearson") else NA_real_,
    gene = gene[1],
    peak = peak[1],
    zscore = zscore[1],
    pvalue = pvalue[1]
  ),
  by = row_id
]

cor_per_row <- cor_per_row[!is.na(row_cor)]
cor_per_row <- cor_per_row[order(-row_cor, -zscore)]

best_row_per_gene <- cor_per_row[
  order(-row_cor, -zscore),
  .SD[1],
  by = gene
]

keep_rows <- head(best_row_per_gene[order(-row_cor, -zscore)]$row_id, opt$top_n)

plot_dt_final <- plot_dt[row_id %in% keep_rows]

plot_dt_final <- merge(
  plot_dt_final,
  cor_per_row[, .(row_id, row_cor)],
  by = "row_id",
  all.x = TRUE
)

# deterministic row ordering by strongest combined RNA + ATAC signal
row_order_dt <- plot_dt_final[
  ,
  .SD[which.max(rna_scaled + atac_scaled)],
  by = row_id
][
  ,
  .(
    row_id,
    row_cor = row_cor[1],
    order_cluster = cluster,
    max_combined = rna_scaled + atac_scaled
  )
]

row_order_dt[, order_cluster := factor(order_cluster, levels = cluster_order)]

row_order_dt <- row_order_dt[
  order(order_cluster, -row_cor, -max_combined)
]

plot_dt_final[, row_id := factor(row_id, levels = row_order_dt$row_id)]

theme_heat <- theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

p_rna <- ggplot(plot_dt_final, aes(x = cluster, y = row_id, fill = rna_scaled)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = c("gray95", "darkorange"),
    limits = c(0, 1),
    name = "scRNA-seq\navg. expression"
  ) +
  labs(title = "scRNA-seq") +
  theme_heat

p_atac <- ggplot(plot_dt_final, aes(x = cluster, y = row_id, fill = atac_scaled)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = c("gray95", "steelblue"),
    limits = c(0, 1),
    name = "scATAC-seq\navg. accessibility"
  ) +
  labs(title = "scATAC-seq") +
  theme_heat

plot_combined <- cowplot::plot_grid(
  p_rna,
  p_atac,
  nrow = 1,
  rel_widths = c(1, 1)
)

ggsave(
  file.path(opt$outdir, "Figure2D_GEX_ATAC_heatmap.pdf"),
  plot_combined,
  width = 10,
  height = 7,
  device = cairo_pdf
)

ggsave(
  file.path(opt$outdir, "Figure2D_GEX_ATAC_heatmap.png"),
  plot_combined,
  width = 10,
  height = 7,
  dpi = 300
)

fwrite(
  links_df,
  file.path(opt$outdir, "Figure2D_signac_links_used.tsv"),
  sep = "\t"
)

fwrite(
  cor_per_row,
  file.path(opt$outdir, "Figure2D_peak_gene_correlations.tsv"),
  sep = "\t"
)

fwrite(
  plot_dt_final,
  file.path(opt$outdir, "Figure2D_heatmap_matrix_long.tsv"),
  sep = "\t"
)

sink(file.path(opt$outdir, "sessionInfo_Figure2D_GEX_ATAC_heatmap.txt"))
sessionInfo()
sink()

cat("Links retained:", nrow(links_df), "\n")
cat("Peak-gene rows plotted:", uniqueN(plot_dt_final$row_id), "\n")
cat("Genes plotted:", uniqueN(plot_dt_final$gene), "\n")
