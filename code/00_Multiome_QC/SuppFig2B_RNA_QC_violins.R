library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

input_file <- "data/00_preprocessing/combined_seurat_prefilter_QC_metadata.csv"
output_dir <- "outputs/00_Multiome_QC"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

qc_df <- read.csv(input_file, check.names = TRUE) %>%
  mutate(
    condition = factor(condition, levels = c("Mock", "HSV1"))
  )

qc_colors <- c(
  "Mock" = "#4863A0",
  "HSV1" = "#E41A1C"
)

make_violin <- function(df, variable, title, y_label) {
  ggplot(df, aes(x = condition, y = .data[[variable]], fill = condition)) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      color = "black",
      linewidth = 0.4
    ) +
    geom_jitter(
      width = 0.18,
      size = 0.18,
      alpha = 0.55,
      color = "black"
    ) +
    scale_fill_manual(values = qc_colors) +
    labs(
      title = title,
      x = NULL,
      y = y_label
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)
    )
}

pB1 <- make_violin(
  qc_df,
  "nFeature_RNA",
  "nFeature_RNA",
  "nFeature_RNA"
)

pB2 <- make_violin(
  qc_df,
  "nCount_RNA",
  "nCount_RNA",
  "nCount_RNA"
)

pB3 <- make_violin(
  qc_df,
  "percent.mt",
  "percent.mt",
  "percent.mt"
)

pB <- pB1 + pB2 + pB3

ggsave(
  file.path(output_dir, "SuppFig2B_RNA_QC_violins.pdf"),
  pB,
  width = 10,
  height = 3.5
)

ggsave(
  file.path(output_dir, "SuppFig2B_RNA_QC_violins.png"),
  pB,
  width = 10,
  height = 3.5,
  dpi = 300
)
