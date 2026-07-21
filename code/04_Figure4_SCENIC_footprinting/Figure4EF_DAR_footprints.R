#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Signac)
  library(Seurat)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(zoo)
  library(patchwork)
  library(grid)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--outdir", type = "character",
              default = "outputs/04_Figure4_SCENIC_footprinting/Figure4EF_DAR_footprints"),
  make_option("--panel", type = "character",
              help = "ifn or macrophage")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

stopifnot(opt$panel %in% c("ifn", "macrophage"))

obj <- readRDS(opt$input)
DefaultAssay(obj) <- "ATAC"

cluster_map <- c(
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

if (!"annotated_clusters" %in% colnames(obj@meta.data)) {
  obj$annotated_clusters <- unname(cluster_map[as.character(obj$seurat_clusters)])
}

obj$annotated_clusters <- as.character(obj$annotated_clusters)
obj$condition <- as.character(obj$condition)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)

make_dar_object <- function(sub_obj, dar_peaks, pfm, use_existing_bias = FALSE, full_obj = NULL) {
  counts_mat <- GetAssayData(sub_obj, assay = "ATAC", slot = "counts")
  counts_mat_dar <- counts_mat[dar_peaks, , drop = FALSE]

  dar_assay <- CreateChromatinAssay(
    counts = counts_mat_dar,
    sep = c("-", "-"),
    genome = "mm10",
    fragments = Fragments(sub_obj)
  )

  dar_obj <- CreateSeuratObject(
    counts = dar_assay,
    assay = "ATAC",
    meta.data = sub_obj@meta.data
  )

  DefaultAssay(dar_obj) <- "ATAC"

  if (use_existing_bias && !is.null(full_obj)) {
    dar_obj[["ATAC"]]@bias <- full_obj[["ATAC"]]@bias
  } else {
    dar_obj <- InsertionBias(dar_obj, genome = BSgenome.Mmusculus.UCSC.mm10)
  }

  dar_obj <- AddMotifs(
    dar_obj,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    pfm = pfm
  )

  if (use_existing_bias && !is.null(full_obj)) {
    dar_obj[["ATAC"]]@bias <- full_obj[["ATAC"]]@bias
  }

  dar_obj
}

get_mean_fp <- function(obj, motif, group_col, keep_groups) {
  mat <- as.matrix(obj[["ATAC"]]@positionEnrichment[[motif]])

  pos <- suppressWarnings(as.numeric(colnames(mat)))
  if (all(is.na(pos))) {
    n <- ncol(mat)
    half <- floor(n / 2)
    pos <- seq(-half, half, length.out = n)
  }

  df <- data.frame(
    cell = rep(rownames(mat), each = ncol(mat)),
    position = rep(pos, times = nrow(mat)),
    value = as.vector(t(mat))
  )

  groups <- obj@meta.data[[group_col]]
  names(groups) <- rownames(obj@meta.data)

  df$group <- groups[df$cell]
  df <- df[!is.na(df$group), ]
  df$group <- factor(df$group, levels = keep_groups)

  df %>%
    filter(group %in% keep_groups) %>%
    group_by(group, position) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    mutate(feature = motif)
}

plot_fp <- function(df, colors, label_map = NULL, k = 7) {
  df <- df %>%
    arrange(group, position) %>%
    group_by(group) %>%
    mutate(smooth = zoo::rollmean(value, k = k, fill = NA, align = "center")) %>%
    ungroup()

  ttl <- unique(df$feature)
  if (!is.null(label_map) && ttl %in% names(label_map)) ttl <- label_map[[ttl]]

  ggplot(df, aes(x = position, y = smooth, color = group)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45) +
    coord_cartesian(xlim = c(-80, 80)) +
    scale_x_continuous(breaks = c(-80, -40, 0, 40, 80)) +
    scale_color_manual(values = colors) +
    labs(
      x = "Distance from motif (bp)",
      y = "Mean Tn5 insertion enrichment",
      title = ttl
    ) +
    theme_classic(base_size = 9) +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      axis.text = element_text(color = "black", size = 8),
      axis.title = element_text(color = "black", size = 9),
      plot.title = element_text(face = "bold", hjust = 0, size = 10),
      plot.margin = margin(4, 4, 4, 4)
    )
}


if (opt$panel == "ifn") {
  sub <- subset(
    obj,
    subset = annotated_clusters %in% c("Homeostatic Microglia", "IFN-Responsive Microglia")
  )
  DefaultAssay(sub) <- "ATAC"

  sub <- InsertionBias(sub, genome = BSgenome.Mmusculus.UCSC.mm10)
  sub <- AddMotifs(sub, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)

  Idents(sub) <- sub$annotated_clusters

  dars <- FindMarkers(
    sub,
    ident.1 = "IFN-Responsive Microglia",
    ident.2 = "Homeostatic Microglia",
    test.use = "LR",
    latent.vars = "nCount_ATAC",
    min.pct = 0.05,
    logfc.threshold = 0.25
  )

  dars <- dars[dars$p_val < 0.05 & dars$avg_log2FC > 0.25, , drop = FALSE]
  stopifnot(nrow(dars) > 0)

  fwrite(
    as.data.frame(dars) %>% tibble::rownames_to_column("peak"),
    file.path(opt$outdir, "Figure4E_IFN_vs_Homeostatic_DARs.tsv"),
    sep = "\t"
  )

  dar_obj <- make_dar_object(sub, rownames(dars), pfm, use_existing_bias = FALSE)

  motifs <- c("STAT1::STAT2", "Stat2", "IRF1", "IRF7")
  colors <- c(
    "Homeostatic Microglia" = "darkorange",
    "IFN-Responsive Microglia" = "#4B0082"
  )
  group_col <- "annotated_clusters"
  keep_groups <- names(colors)
  prefix <- "Figure4E_IFN_microglia"
}

if (opt$panel == "macrophage") {
  if (is.null(obj[["ATAC"]]@bias)) {
    obj <- InsertionBias(obj, genome = BSgenome.Mmusculus.UCSC.mm10)
  }

  sub <- subset(obj, subset = annotated_clusters == "Infiltrating Macrophages")
  DefaultAssay(sub) <- "ATAC"

  sub <- AddMotifs(sub, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)

  Idents(sub) <- sub$condition

  dars <- FindMarkers(
    sub,
    ident.1 = "HSV1",
    ident.2 = "Mock",
    test.use = "LR",
    latent.vars = "nCount_ATAC",
    min.pct = 0.05,
    logfc.threshold = 0.25
  )

  dars <- dars[dars$p_val < 0.05 & dars$avg_log2FC > 0.25, , drop = FALSE]
  stopifnot(nrow(dars) > 0)

  fwrite(
    as.data.frame(dars) %>% tibble::rownames_to_column("peak"),
    file.path(opt$outdir, "Figure4F_HSV1_vs_Mock_macrophage_DARs.tsv"),
    sep = "\t"
  )

  dar_obj <- make_dar_object(sub, rownames(dars), pfm, use_existing_bias = TRUE, full_obj = obj)

  motifs <- c("STAT1::STAT2", "IRF1", "IRF3", "IRF7")
  colors <- c("Mock" = "gray40", "HSV1" = "#b2182b")
  group_col <- "condition"
  keep_groups <- c("Mock", "HSV1")
  prefix <- "Figure4F_macrophage"
}

available <- names(dar_obj[["ATAC"]]@positionEnrichment)

motif_names <- unique(unlist(dar_obj[["ATAC"]]@motifs@motif.names))
motifs <- intersect(motifs, motif_names)

stopifnot(length(motifs) > 0)

dar_obj <- Footprint(
  object = dar_obj,
  motif.name = motifs,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

label_map <- c(
  "Stat2" = "STAT2"
)

plots <- list()

for (m in motifs) {
  df <- get_mean_fp(
    obj = dar_obj,
    motif = m,
    group_col = group_col,
    keep_groups = keep_groups
  )

  write.csv(
    df,
    file.path(opt$outdir, paste0(prefix, "_", gsub("[:\\.\\(\\)]", "", m), "_footprint_data.csv")),
    row.names = FALSE
  )

  p <- plot_fp(df, colors = colors, label_map = label_map)
  plots[[m]] <- p

  ggsave(
    file.path(opt$outdir, paste0(prefix, "_", gsub("[:\\.\\(\\)]", "", m), ".pdf")),
    p,
    width = 3.4,
    height = 2.7,
    useDingbats = FALSE
  )
}

grid <- wrap_plots(plots, ncol = 2)

# Add one shared legend to the right, matching manuscript panel layout
legend_df <- data.frame(
  position = rep(c(-80, 80), times = length(colors)),
  smooth = rep(seq_along(colors), each = 2),
  group = factor(rep(names(colors), each = 2), levels = names(colors))
)

legend_plot <- ggplot(legend_df, aes(position, smooth, color = group)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = colors) +
  theme_void(base_size = 9) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8)
  )

final_grid <- grid + patchwork::wrap_elements(
  ggplotGrob(legend_plot)
) + patchwork::plot_layout(widths = c(1, 0.28))

ggsave(
  file.path(opt$outdir, paste0(prefix, "_footprint_grid.pdf")),
  final_grid,
  width = 7.0,
  height = 5.0,
  useDingbats = FALSE
)

sink(file.path(opt$outdir, paste0("sessionInfo_", prefix, ".txt")))
sessionInfo()
sink()

cat("Done:", prefix, "\n")
