suppressPackageStartupMessages({
  library(SpatialDecon)
  library(GeomxTools)
  library(Biobase)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

input_file <- "data/05_Figure5_spatial_transcriptomics/processed/GeoMx_WTA_processed_manuscript.rds"
output_dir <- "outputs/05_Figure5_spatial_transcriptomics/Figure5F_deconvolution"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==== Load processed GeoMx object ====
geoMx <- readRDS(input_file)

# ==== Expression matrix ====
expr_mat <- assayDataElement(geoMx, "q_norm")
rownames(expr_mat) <- fData(geoMx)$TargetName
colnames(expr_mat) <- sampleNames(geoMx)

# ==== Mouse Cell Atlas brain reference ====
brain_ref <- download_profile_matrix(
  species = "Mouse",
  age_group = "Adult",
  matrixname = "Brain_MCA"
)

# ==== Background estimate ====
if (!"NegProbe-WTX" %in% rownames(expr_mat)) {
  stop("NegProbe-WTX was not found in the q_norm expression matrix.")
}

per_observation_neg <- expr_mat["NegProbe-WTX", ]
bg <- sweep(expr_mat * 0, 2, per_observation_neg, "+")

# ==== Run SpatialDecon ====
annot <- pData(geoMx)
annot$HSV1_positive <- annot$antigen == "HSV1 positive"

res_full <- spatialdecon(
  norm = expr_mat,
  bg = bg,
  X = brain_ref,
  is_pure_tumor = annot$HSV1_positive,
  cell_counts = annot$Nuclei,
  n_tumor_clusters = 3,
  align_genes = TRUE
)

# ==== Collapse reference cell types ====
celltypes <- rownames(res_full$beta)

matching <- list(
  Astrocytes = grep("^Astro", celltypes, value = TRUE),
  Macrophages = grep("Macrophage", celltypes, value = TRUE),
  Microglia = grep("Microglia", celltypes, value = TRUE),
  Neurons = grep("Neuron", celltypes, value = TRUE),
  Oligodendrocytes = grep("Oligo", celltypes, value = TRUE)
)

matching <- matching[lengths(matching) > 0]

res_collapsed <- collapseCellTypes(
  res_full,
  matching
)

# ==== Prepare plotting data ====
beta_df <- as.data.frame(t(res_collapsed$beta))
beta_df$ROI <- rownames(beta_df)

beta_long <- beta_df %>%
  left_join(
    annot %>%
      as.data.frame() %>%
      mutate(ROI = rownames(.)) %>%
      dplyr::select(ROI, antigen),
    by = "ROI"
  ) %>%
  filter(
    antigen %in% c(
      "HSV1 negative",
      "HSV1 positive"
    )
  ) %>%
  pivot_longer(
    cols = all_of(names(matching)),
    names_to = "CellType",
    values_to = "Abundance"
  ) %>%
  mutate(
    antigen = factor(
      antigen,
      levels = c(
        "HSV1 negative",
        "HSV1 positive"
      ),
      labels = c(
        "negative",
        "positive"
      )
    ),
    CellType = factor(
      CellType,
      levels = c(
        "Astrocytes",
        "Macrophages",
        "Microglia",
        "Neurons",
        "Oligodendrocytes"
      )
    )
  )

# ==== Wilcoxon tests ====
wilcox_stats <- beta_long %>%
  group_by(CellType) %>%
  summarise(
    p_value = wilcox.test(
      Abundance ~ antigen,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(
      p_value,
      method = "fdr"
    ),
    label = paste0(
      "italic(padj) == ",
      format(
        signif(p_adj, 2),
        trim = TRUE,
        scientific = FALSE
      )
    )
  )

write.csv(
  wilcox_stats,
  file.path(
    output_dir,
    "Figure5F_deconvolution_wilcox_stats.csv"
  ),
  row.names = FALSE
)

beta_long <- beta_long %>%
  left_join(
    wilcox_stats,
    by = "CellType"
  )

# ==== Plot ====
figure5f <- ggplot(
  beta_long,
  aes(
    x = antigen,
    y = Abundance,
    fill = antigen
  )
) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    alpha = 0.65,
    color = "black",
    linewidth = 0.7
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.7
  ) +
  facet_wrap(
    vars(CellType, label),
    scales = "free_y",
    nrow = 1,
    labeller = labeller(
      label = label_parsed
    )
  ) +
  scale_fill_manual(
    values = c(
      negative = "gray55",
      positive = "#D85C41"
    )
  ) +
  labs(
    x = "HSV-1 infection",
    y = "Estimated Abundance"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    ),
    strip.text = element_text(
      size = 11
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    axis.title.x = element_text(
      size = 13
    ),
    axis.title.y = element_text(
      size = 13
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.8
    ),
    panel.spacing.x = grid::unit(
      1.1,
      "lines"
    )
  )

ggsave(
  file.path(
    output_dir,
    "Figure5F_deconvolution.pdf"
  ),
  figure5f,
  width = 10,
  height = 6.5,
  device = cairo_pdf
)

ggsave(
  file.path(
    output_dir,
    "Figure5F_deconvolution.png"
  ),
  figure5f,
  width = 10,
  height = 6.5,
  dpi = 300
)

message(
  "Saved Figure 5F to: ",
  normalizePath(output_dir)
)
