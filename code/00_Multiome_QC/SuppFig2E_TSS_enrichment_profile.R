#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
})

object_path <- "data/03_add_final_annotations/combined_seurat_final_annotated.rds"
output_dir <- "outputs/00_Multiome_QC"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mult_obj <- readRDS(object_path)
DefaultAssay(mult_obj) <- "ATAC"

mult_obj$condition <- factor(
  mult_obj$condition,
  levels = c("HSV1", "Mock")
)

if (!"TSS" %in% names(mult_obj[["ATAC"]]@positionEnrichment)) {
  mult_obj <- TSSEnrichment(
    mult_obj,
    assay = "ATAC",
    fast = FALSE
  )
}

pE <- TSSPlot(
  mult_obj,
  assay = "ATAC",
  group.by = "condition"
) +
  scale_color_manual(
    values = c(
      "HSV1" = "#F02B24",
      "Mock" = "#3B5BA9"
    )
  ) +
  labs(
    title = "TSS Enrichment",
    x = "Distance from TSS (bp)",
    y = "Mean TSS enrichment score",
    color = "group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    strip.background = element_rect(
      fill = "white",
      color = "black"
    ),
    strip.text = element_text(face = "bold"),
    panel.spacing = grid::unit(0.15, "lines")
  )

ggsave(
  file.path(output_dir, "SuppFig2E_TSS_enrichment_profile.pdf"),
  pE,
  width = 8,
  height = 4.5
)

ggsave(
  file.path(output_dir, "SuppFig2E_TSS_enrichment_profile.png"),
  pE,
  width = 8,
  height = 4.5,
  dpi = 300
)
