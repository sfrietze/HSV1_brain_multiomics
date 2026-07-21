suppressPackageStartupMessages({
  library(GeomxTools)
  library(Biobase)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

input <- "data/05_Figure5_spatial_transcriptomics/processed/GeoMx_WTA_processed_manuscript.rds"
dea_file <- "outputs/05_Figure5_spatial_transcriptomics/Figure5B_volcano/Figure5B_DEA_HSV1_Antigen_Pos_vs_Neg.csv"
output_dir <- "outputs/05_Figure5_spatial_transcriptomics/Figure5C_heatmap"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

geoMx <- readRDS(input)
dea_result <- read.csv(dea_file, check.names = FALSE)

pData(geoMx)$region_annotation <- ifelse(
  as.numeric(pData(geoMx)$`ROI (Label)`) <= 12,
  "TG",
  "TN"
)

assayDataElement(geoMx, "log_q") <-
  assayDataApply(geoMx, 2, log, base = 2, elt = "q_norm")

expr_mat <- as.matrix(assayDataElement(geoMx, "log_q"))

relaxed_deg_genes <- dea_result %>%
  filter(
    p.value < 0.05,
    abs(log2foldchange) >= 0.8
  ) %>%
  pull(gene)

relaxed_deg_genes <- intersect(
  relaxed_deg_genes,
  rownames(expr_mat)
)

expr_mat_relaxed <- expr_mat[
  relaxed_deg_genes,
  ,
  drop = FALSE
]

z_scaled_relaxed <- t(scale(t(expr_mat_relaxed)))

z_scaled_capped_relaxed <- pmax(
  pmin(z_scaled_relaxed, 2),
  -2
)

anno_df <- pData(geoMx) %>%
  as.data.frame()

heatmap_object <- Heatmap(
  z_scaled_capped_relaxed,
  name = "Z-score",
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = anno_df$condition,
  top_annotation = HeatmapAnnotation(
    Condition = anno_df$condition,
    Antigen = anno_df$antigen,
    Region = anno_df$region_annotation,
    col = list(
      Condition = c(
        "HSV" = "purple4",
        "Mock" = "gray60"
      ),
      Antigen = c(
        "HSV1 positive" = "darkred",
        "HSV1 negative" = "lightblue"
      ),
      Region = c(
        "TG" = "grey80",
        "TN" = "grey60"
      )
    )
  ),
  col = colorRamp2(
    c(-2, 0, 2),
    c("navy", "white", "firebrick3")
  )
)

pdf(
  file.path(output_dir, "Figure5C_DEG_heatmap.pdf"),
  width = 10,
  height = 8,
  useDingbats = FALSE
)

draw(heatmap_object)

dev.off()

png(
  file.path(output_dir, "Figure5C_DEG_heatmap.png"),
  width = 3000,
  height = 2400,
  res = 300
)

draw(heatmap_object)

dev.off()

write.csv(
  z_scaled_capped_relaxed,
  file.path(output_dir, "Figure5C_heatmap_source_data.csv")
)

message(
  "Saved Figure 5C to: ",
  normalizePath(output_dir)
)
