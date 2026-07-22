#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

input_file <- "data/00_preprocessing/combined_seurat_prefilter_QC_metadata.csv"
output_dir <- "outputs/00_Multiome_QC"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

qc <- read.csv(input_file, check.names = TRUE)

qc$condition <- factor(qc$condition, levels = c("Mock", "HSV1"))

cols <- c(
  Mock = "#3B5BA9",
  HSV1 = "#FF2A1A"
)

## percent.mt vs nCount_RNA
p1 <- ggplot(
  qc,
  aes(
    x = nCount_RNA,
    y = percent.mt,
    colour = condition
  )
) +
  geom_point(
    size = 0.5,
    alpha = 0.75
  ) +
  scale_colour_manual(values = cols) +
  labs(
    x = "nCount_RNA",
    y = "percent.mt"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  )

## nFeature_RNA vs nCount_RNA
p2 <- ggplot(
  qc,
  aes(
    x = nCount_RNA,
    y = nFeature_RNA,
    colour = condition
  )
) +
  geom_point(
    size = 0.5,
    alpha = 0.75
  ) +
  scale_colour_manual(values = cols) +
  labs(
    x = "nCount_RNA",
    y = "nFeature_RNA"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    axis.line = element_blank(),
    legend.position = "right",
    legend.title = element_blank()
  )

panel_c <- p1 + p2

ggsave(
  file.path(output_dir, "SuppFig2C_RNA_QC_scatterplots.pdf"),
  panel_c,
  width = 9,
  height = 4
)

ggsave(
  file.path(output_dir, "SuppFig2C_RNA_QC_scatterplots.png"),
  panel_c,
  width = 9,
  height = 4,
  dpi = 300
)

pC <- panel_c
