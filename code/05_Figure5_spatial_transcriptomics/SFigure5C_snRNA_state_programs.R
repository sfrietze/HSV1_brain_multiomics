#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(GeomxTools)
  library(Biobase)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

sn_file <- "data/03_add_final_annotations/combined_seurat_final_annotated.rds"
geomx_file <- "data/05_Figure5_spatial_transcriptomics/processed/GeoMx_WTA_processed_manuscript.rds"
outdir <- "outputs/05_Figure5_spatial_transcriptomics/SFigure5C_snRNA_state_programs"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

states <- c(
  "Homeostatic Microglia",
  "Transiently Activated Microglia",
  "IFN-Responsive Microglia",
  "Primed Microglia",
  "Mitochondrial-Activated Microglia",
  "IEG-High Microglia",
  "Infiltrating Macrophages"
)

sn <- readRDS(sn_file)
geo <- readRDS(geomx_file)

sn <- subset(sn, subset = annotated_clusters %in% states)
DefaultAssay(sn) <- "RNA"
Idents(sn) <- "annotated_clusters"

geomx_genes <- unique(fData(geo)$TargetName)

markers <- FindAllMarkers(
  sn,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25
)

sig <- markers %>%
  filter(
    gene %in% geomx_genes,
    p_val_adj < 0.05,
    avg_log2FC >= 0.25
  ) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 25) %>%
  ungroup()

write.csv(
  sig,
  file.path(outdir, "SupplementaryFigure5C_snRNA_state_marker_genes.csv"),
  row.names = FALSE
)

q <- assayDataElement(geo, "q_norm")
rownames(q) <- fData(geo)$TargetName
colnames(q) <- sampleNames(geo)

annot <- as.data.frame(pData(geo))
annot$ROI <- rownames(annot)

z <- t(scale(t(log2(q + 1))))
z <- z[complete.cases(z), , drop = FALSE]

scores <- bind_rows(lapply(states, function(state) {
  genes <- intersect(
    sig$gene[sig$cluster == state],
    rownames(z)
  )

  data.frame(
    ROI = colnames(z),
    CellState = state,
    Score = colMeans(z[genes, , drop = FALSE])
  )
})) %>%
  left_join(
    annot %>% select(ROI, antigen),
    by = "ROI"
  ) %>%
  filter(
    antigen %in% c("HSV1 negative", "HSV1 positive")
  ) %>%
  mutate(
    antigen = factor(
      antigen,
      levels = c("HSV1 negative", "HSV1 positive"),
      labels = c("negative", "positive")
    ),
    CellState = factor(CellState, levels = states)
  )

stats <- scores %>%
  group_by(CellState) %>%
  summarise(
    n_negative = sum(antigen == "negative"),
    n_positive = sum(antigen == "positive"),
    median_negative = median(Score[antigen == "negative"]),
    median_positive = median(Score[antigen == "positive"]),
    delta_median = median_positive - median_negative,
    p_value = wilcox.test(Score ~ antigen, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    FDR_label = case_when(
      p_adj < 0.001 ~ paste0("FDR = ", formatC(p_adj, format = "e", digits = 1)),
      TRUE ~ paste0("FDR = ", formatC(p_adj, format = "f", digits = 3))
    )
  )

write.csv(
  scores,
  file.path(outdir, "SupplementaryFigure5C_state_program_scores_by_ROI.csv"),
  row.names = FALSE
)

write.csv(
  stats,
  file.path(outdir, "SupplementaryFigure5C_statistics.csv"),
  row.names = FALSE
)

plot_data <- scores %>%
  left_join(
    stats %>% select(CellState, FDR_label),
    by = "CellState"
  ) %>%
  mutate(
    facet_label = paste0(CellState, "\n", FDR_label),
    facet_label = factor(
      facet_label,
      levels = paste0(
        states,
        "\n",
        stats$FDR_label[match(states, stats$CellState)]
      )
    )
  )

p <- ggplot(
  plot_data,
  aes(x = antigen, y = Score, fill = antigen)
) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    linewidth = 0.65,
    alpha = 0.75
  ) +
  geom_jitter(
    width = 0.10,
    size = 1.5,
    alpha = 0.75
  ) +
  facet_wrap(
    ~ facet_label,
    nrow = 2
  ) +
  scale_fill_manual(
    values = c(
      negative = "gray55",
      positive = "#D85C41"
    )
  ) +
  labs(
    x = "HSV-1 antigen",
    y = "snRNA-derived state score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.7
    ),
    strip.text = element_text(size = 9),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    panel.spacing = grid::unit(0.8, "lines")
  )

ggsave(
  file.path(
    outdir,
    "SupplementaryFigure5C_snRNA_state_programs.pdf"
  ),
  p,
  width = 10,
  height = 7,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "SupplementaryFigure5C_snRNA_state_programs.png"
  ),
  p,
  width = 10,
  height = 7,
  dpi = 300
)

print(as.data.frame(stats), row.names = FALSE)

message("Saved Supplementary Figure 5C to: ", normalizePath(outdir))
