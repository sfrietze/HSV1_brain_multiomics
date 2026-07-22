#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

object_path <- "data/03_add_final_annotations/combined_seurat_final_annotated.rds"

output_dir <- "outputs/00_Multiome_QC"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mult_obj <- readRDS(object_path)

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

if (!"annotated_clusters" %in% colnames(mult_obj@meta.data) ||
    all(grepl("^[0-9]+$", unique(mult_obj$annotated_clusters)))) {

  mult_obj$annotated_clusters <- unname(
    new_annotations[as.character(mult_obj$seurat_clusters)]
  )
}

plot_order <- rev(unname(new_annotations))

cell_counts <- mult_obj@meta.data %>%
  dplyr::count(annotated_clusters, condition, name = "nuclei") %>%
  mutate(
    annotated_clusters = factor(
      annotated_clusters,
      levels = plot_order
    ),
    condition = factor(
      condition,
      levels = c("HSV1", "Mock")
    )
  )

qc_colors <- c(
  "HSV1" = "#F02B24",
  "Mock" = "#3B5BA9"
)

pF <- ggplot(
  cell_counts,
  aes(
    x = annotated_clusters,
    y = nuclei,
    fill = condition
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.75
  ) +
  scale_fill_manual(values = qc_colors) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top"
  ) +
  labs(
    x = "Cell state",
    y = "Number of nuclei",
    fill = "Condition"
  )

ggsave(
  file.path(output_dir, "SuppFig2F_cell_state_composition.pdf"),
  pF,
  width = 9,
  height = 5
)

ggsave(
  file.path(output_dir, "SuppFig2F_cell_state_composition.png"),
  pF,
  width = 9,
  height = 5,
  dpi = 300
)
