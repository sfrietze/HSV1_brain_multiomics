#!/usr/bin/env Rscript

# ============================================================
# Figure 1 / Supplementary Figure 1
# GeoMx DSP protein analysis
#
# This script performs DSP protein differential abundance,
# creates marker plots, heatmap, and normalization-control plots.
#
# Usage:
#   Rscript Figure1_DSP_protein_analysis.R \
#     --protein data/01_Figure1_DSP_protein/DSP_protein_counts.xlsx \
#     --metadata data/01_Figure1_DSP_protein/GeoMX_ROI_metadata.csv \
#     --out outputs/01_Figure1_DSP_protein
#
# expected inputs:
#   1) DSP protein nCounter Excel export
#   2) ROI metadata CSV containing Scan Name, Slide Name,
#      ROI (Label), antigen, and infected_region columns
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(pheatmap)
  library(ggrepel)
  library(cowplot)
  library(ggplotify)
  library(patchwork)
})

# -----------------------------
# Parse arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[idx + 1]
}

protein_file <- get_arg("--protein", "data/01_Figure1_DSP_protein/DSP_protein_counts.xlsx")
meta_file    <- get_arg("--metadata", "data/01_Figure1_DSP_protein/GeoMX_ROI_metadata.csv")
out_dir      <- get_arg("--out", "outputs/01_Figure1_DSP_protein")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Helper functions
clean_antigen <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x[x %in% c("", "NA", "NaN")] <- NA
  x <- ifelse(is.na(x), "Unknown", x)

  # normalize common variants
  x <- dplyr::case_when(
    grepl("positive", x, ignore.case = TRUE) ~ "HSV1 positive",
    grepl("negative", x, ignore.case = TRUE) ~ "HSV1 negative",
    TRUE ~ x
  )

  factor(x, levels = c("HSV1 negative", "HSV1 positive", "Unknown"))
}

safe_save <- function(plot, filename, width = 7, height = 5) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )
}

plot_violin_marker <- function(expr_mat, meta_df, protein_name) {
  if (!protein_name %in% rownames(expr_mat)) {
    warning("Protein not found: ", protein_name)
    return(
      ggplot() +
        theme_void() +
        labs(title = paste0(protein_name, " not found"))
    )
  }

  df <- data.frame(
    Expression = as.numeric(expr_mat[protein_name, ]),
    Antigen = factor(meta_df$antigen, levels = c("HSV1 negative", "HSV1 positive"))
  )

  ggplot(df, aes(x = Antigen, y = Expression, fill = Antigen)) +
    geom_violin(trim = FALSE, color = "black", linewidth = 0.3) +
    geom_jitter(width = 0.08, size = 2, color = "black", alpha = 0.8) +
    scale_fill_manual(values = c(
      "HSV1 negative" = "#74ADD1",
      "HSV1 positive" = "#B2182B"
    )) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(x = NULL, y = "log2 protein abundance", title = protein_name)
}

# Load DSP protein matrix
raw <- readxl::read_excel(protein_file, sheet = 1, col_names = FALSE) |>
  as.data.frame()

start_row <- which(raw[, 1] == "#Target Group")
stopifnot(length(start_row) == 1)

expr_raw <- raw[(start_row + 1):nrow(raw), ]

protein_names <- expr_raw[, 4]
expr_matrix <- expr_raw[, 5:ncol(raw)]
expr_matrix <- apply(expr_matrix, 2, as.numeric) |> as.data.frame()

rownames(expr_matrix) <- protein_names
expr_matrix <- expr_matrix[!is.na(rownames(expr_matrix)), ]

roi_row <- which(raw[, 1] == "ROI (Label)")
scan_row <- which(raw[, 1] == "Scan Name")

stopifnot(length(roi_row) == 1)
stopifnot(length(scan_row) == 1)

roi_labels <- as.character(unlist(raw[roi_row, 5:ncol(raw)]))
scan_names <- as.character(unlist(raw[scan_row, 5:ncol(raw)]))

colnames(expr_matrix) <- paste(scan_names, roi_labels, sep = "_")

expr_log <- log2(expr_matrix + 1)

protein_meta <- data.frame(
  ProteinROI = colnames(expr_log),
  stringsAsFactors = FALSE
)

protein_meta$ROI_number <- stringr::str_extract(protein_meta$ProteinROI, "[0-9]{3}$")

protein_meta$Slide_number <- ifelse(
  grepl("Slide 11", protein_meta$ProteinROI), "11",
  ifelse(grepl("Slide5|Slide 5", protein_meta$ProteinROI), "5", NA)
)

message("Protein matrix dimensions:")
print(dim(expr_log))

message("Protein names:")
print(head(rownames(expr_log), 20))

message("ROI metadata:")
print(head(protein_meta))
print(table(protein_meta$Slide_number, useNA = "ifany"))

writeLines(
  rownames(expr_log),
  file.path(out_dir, "available_proteins.txt")
)

# Parse ROI IDs and merge metadata
protein_meta <- data.frame(
  ProteinROI = colnames(expr_log),
  stringsAsFactors = FALSE
)

protein_meta$ROI_number <- stringr::str_extract(protein_meta$ProteinROI, "[0-9]{3}$")

protein_meta$Slide_number <- dplyr::case_when(
  grepl("Slide 11", protein_meta$ProteinROI) ~ "11",
  grepl("Slide5|Slide 5", protein_meta$ProteinROI) ~ "5",
  TRUE ~ NA_character_
)

metadata_raw <- readr::read_csv(meta_file, show_col_types = FALSE)

required_cols <- c("Scan Name", "Slide Name", "ROI (Label)", "antigen", "infected_region")
missing_cols <- setdiff(required_cols, colnames(metadata_raw))

metadata_clean <- metadata_raw %>%
  dplyr::filter(`Scan Name` %in% c("5_WTA_final", "11_WTA_final")) %>%
  dplyr::filter(!is.na(`ROI (Label)`)) %>%
  dplyr::mutate(
    ROI_number = sprintf("%03d", as.numeric(`ROI (Label)`)),
    Slide_number = ifelse(`Slide Name` == "5_WTA", "5", "11")
  )

metadata <- dplyr::left_join(
  protein_meta,
  metadata_clean,
  by = c("ROI_number", "Slide_number")
) %>%
  dplyr::filter(!is.na(infected_region))

expr_log <- expr_log[, metadata$ProteinROI, drop = FALSE]

metadata$antigen <- clean_antigen(metadata$antigen)

write.csv(
  metadata,
  file.path(out_dir, "Figure1_DSP_ROI_metadata_used.csv"),
  row.names = FALSE
)

# QC plots

pca <- prcomp(t(expr_log), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Antigen = metadata$antigen,
  Infection = metadata$infected_region,
  Slide = metadata$Slide_number
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Antigen, shape = Slide)) +
  geom_point(size = 3) +
  theme_classic(base_size = 14) +
  labs(title = "DSP protein PCA")

safe_save(p_pca, "QC_DSP_protein_PCA.pdf", width = 6, height = 5)

roi_mean_df <- data.frame(
  ROI = colnames(expr_log),
  MeanExpression = colMeans(expr_log),
  MedianExpression = apply(expr_log, 2, median),
  Antigen = metadata$antigen,
  Slide = metadata$Slide_number
)

p_roi_median <- ggplot(
  roi_mean_df %>% dplyr::filter(Slide == "5", Antigen %in% c("HSV1 negative", "HSV1 positive")),
  aes(x = Antigen, y = MedianExpression, fill = Antigen)
) +
  geom_boxplot(alpha = 0.85, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 2.5, color = "black") +
  scale_fill_manual(values = c(
    "HSV1 negative" = "#74ADD1",
    "HSV1 positive" = "#B2182B"
  )) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Per-ROI median log2 abundance")

# Standalone Supplementary Figure 1A output omitted; included in combined controls figure.
# Slide 5 antigen-positive vs antigen-negative analysis

hsv_idx <- which(metadata$Slide_number == "5")
hsv_expr <- expr_log[, hsv_idx, drop = FALSE]
hsv_meta <- metadata[hsv_idx, , drop = FALSE]

hsv_meta$antigen <- clean_antigen(hsv_meta$antigen)

keep_antigen <- hsv_meta$antigen %in% c("HSV1 negative", "HSV1 positive")
hsv_expr <- hsv_expr[, keep_antigen, drop = FALSE]
hsv_meta <- hsv_meta[keep_antigen, , drop = FALSE]

hsv_meta$antigen <- factor(hsv_meta$antigen, levels = c("HSV1 negative", "HSV1 positive"))

if (ncol(hsv_expr) < 3) {
  stop("Too few HSV slide ROIs after filtering.")
}

design <- model.matrix(~ antigen, data = hsv_meta)

fit <- limma::lmFit(hsv_expr, design)
fit <- limma::eBayes(fit)

results_raw <- limma::topTable(
  fit,
  coef = "antigenHSV1 positive",
  number = Inf
)

results_raw <- results_raw[order(results_raw$adj.P.Val), ]
results_raw$Protein <- rownames(results_raw)

write.csv(
  results_raw,
  file.path(out_dir, "Figure1E_DSP_antigen_positive_vs_negative_results_raw.csv"),
  row.names = FALSE
)

# Median-centered sensitivity analysis

hsv_expr_centered <- sweep(
  hsv_expr,
  2,
  apply(hsv_expr, 2, median)
)

fit_centered <- limma::lmFit(hsv_expr_centered, design)
fit_centered <- limma::eBayes(fit_centered)

results_centered <- limma::topTable(
  fit_centered,
  coef = "antigenHSV1 positive",
  number = Inf
)

results_centered <- results_centered[order(results_centered$adj.P.Val), ]
results_centered$Protein <- rownames(results_centered)

write.csv(
  results_centered,
  file.path(out_dir, "Supplementary_Figure1_centered_DSP_results.csv"),
  row.names = FALSE
)

common_proteins <- intersect(rownames(results_raw), rownames(results_centered))

delta_specific <- results_raw$logFC - results_centered$logFC

delta_df <- data.frame(
  Protein = rownames(results_raw),
  delta_logFC = delta_specific,
  logFC_raw = results_raw$logFC,
  logFC_centered = results_centered$logFC,
  FDR_raw = results_raw$adj.P.Val,
  FDR_centered = results_centered$adj.P.Val
)

delta_df <- delta_df[order(delta_df$delta_logFC), ]


write.csv(
  delta_df,
  file.path(out_dir, "Supplementary_Figure1_raw_vs_centered_comparison.csv"),
  row.names = FALSE
)

# Figure 1E volcano
volc_df <- results_raw %>%
  dplyr::mutate(
    negLogFDR = -log10(adj.P.Val),
    category = dplyr::case_when(
      adj.P.Val < 0.05 & logFC > 0.5 ~ "Up",
      adj.P.Val < 0.05 & logFC < -0.5 ~ "Down",
      TRUE ~ "NS"
    )
  )

top_labels <- volc_df %>% dplyr::arrange(adj.P.Val) %>% head(12)

p_volcano <- ggplot(volc_df, aes(x = logFC, y = negLogFDR)) +
  geom_point(aes(color = category), size = 3.2, alpha = 0.9) +
  scale_color_manual(values = c(
    "Up" = "firebrick",
    "Down" = "royalblue",
    "NS" = "grey70"
  )) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = Protein),
    size = 4,
    max.overlaps = 50
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "log2 FC (HSV1 antigen-positive vs antigen-negative)",
    y = "-log10 FDR",
    title = "DSP protein abundance"
  )

safe_save(p_volcano, "Figure1E_DSP_protein_volcano.pdf", width = 6.5, height = 5.5)


# Figure 1F heatmap
top_proteins <- rownames(results_raw)[1:min(25, nrow(results_raw))]
heat_data <- hsv_expr[top_proteins, , drop = FALSE]

ord <- order(hsv_meta$antigen)
heat_data <- heat_data[, ord, drop = FALSE]
heat_meta <- hsv_meta[ord, , drop = FALSE]

annotation_col <- data.frame(
  Antigen = factor(heat_meta$antigen, levels = c("HSV1 negative", "HSV1 positive"))
)
rownames(annotation_col) <- colnames(heat_data)

annotation_colors <- list(
  Antigen = c("HSV1 negative" = "#74ADD1", "HSV1 positive" = "#B2182B")
)

pheat_silent <- pheatmap::pheatmap(
  heat_data,
  scale = "row",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  border_color = NA,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  silent = TRUE
)

p_heatmap <- ggplotify::as.ggplot(pheat_silent)
safe_save(p_heatmap, "Figure1F_DSP_top_protein_heatmap.pdf", width = 6.5, height = 6.5)

# Figure 1G selected marker plots
markers_to_plot <- c("CD11b", "IBA1", "Ctsd", "MHC II")

marker_plots <- lapply(markers_to_plot, function(x) {
  plot_violin_marker(hsv_expr, hsv_meta, x)
})

p_markers <- cowplot::plot_grid(plotlist = marker_plots, ncol = 4)

safe_save(p_markers, "Figure1G_selected_DSP_markers.pdf", width = 10, height = 4)

# Individual marker PDFs omitted; markers are saved in the combined panel above.

# Supplementary Figure 1 normalization controls

centered_volc_df <- results_centered %>%
  dplyr::mutate(
    negLogFDR = -log10(adj.P.Val)
  )

p_centered_volcano <- ggplot(centered_volc_df, aes(x = logFC, y = negLogFDR)) +
  geom_point(size = 2.5, color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, color = "grey80") +
  theme_classic(base_size = 14) +
  labs(
    x = "Median-centered log2 FC",
    y = "-log10 FDR"
  )

proteins_to_label <- c("CD68", "CSF1R", "MSR1", "MHC II", "CD11b", "IBA1", "Ctsd")

p_raw_vs_centered <- ggplot(delta_df, aes(x = logFC_raw, y = logFC_centered)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = subset(delta_df, Protein %in% proteins_to_label),
    aes(label = Protein),
    size = 4,
    max.overlaps = 50
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Raw log2 FC",
    y = "Median-centered log2 FC"
  )

p_sfig1 <- (p_roi_median | p_centered_volcano) / p_raw_vs_centered

safe_save(p_sfig1, "Supplementary_Figure1_DSP_normalization_controls.pdf", width = 9, height = 8)

# 11. Summary and session info
summary_lines <- c(
  paste0("Number of proteins tested: ", nrow(results_raw)),
  paste0("Number of Slide 5 ROIs used: ", ncol(hsv_expr)),
  paste0("HSV1 negative ROIs: ", sum(hsv_meta$antigen == "HSV1 negative")),
  paste0("HSV1 positive ROIs: ", sum(hsv_meta$antigen == "HSV1 positive")),
  paste0("Significant proteins raw limma FDR < 0.05: ", sum(results_raw$adj.P.Val < 0.05, na.rm = TRUE)),
  paste0("Significant proteins centered limma FDR < 0.05: ", sum(results_centered$adj.P.Val < 0.05, na.rm = TRUE))
)

writeLines(summary_lines, con = file.path(out_dir, "Figure1_DSP_summary.txt"))

sink(file.path(out_dir, "sessionInfo_Figure1_DSP.txt"))
sessionInfo()
sink()

message("Done. Outputs written to: ", out_dir)
