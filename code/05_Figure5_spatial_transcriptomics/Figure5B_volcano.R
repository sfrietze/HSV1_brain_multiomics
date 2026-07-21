suppressPackageStartupMessages({
  library(GeomxTools)
  library(Biobase)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

input <- "data/05_Figure5_spatial_transcriptomics/processed/GeoMx_WTA_processed_manuscript.rds"
output_dir <- "outputs/05_Figure5_spatial_transcriptomics/Figure5B_volcano"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

geoMx <- readRDS(input)

pData(geoMx)$infection_region_simple <- case_when(
  pData(geoMx)$antigen == "HSV1 positive" ~ "Infected",
  pData(geoMx)$antigen == "HSV1 negative" ~ "Uninfected",
  TRUE ~ NA_character_
)

geoMx_subset <- geoMx[, !is.na(pData(geoMx)$infection_region_simple)]

run_geoMx_DEA <- function(
  gset,
  group_col,
  group1,
  group2,
  assay = "q_norm"
) {
  expr <- assayDataElement(gset, elt = assay)
  meta <- pData(gset)

  group <- meta[[group_col]]
  group1_index <- which(group == group1)
  group2_index <- which(group == group2)

  results <- apply(expr, 1, function(gene_expression) {
    group1_values <- gene_expression[group1_index]
    group2_values <- gene_expression[group2_index]

    log2_fold_change <-
      log2(mean(group2_values, na.rm = TRUE) + 1) -
      log2(mean(group1_values, na.rm = TRUE) + 1)

    p_value <- suppressWarnings(
      wilcox.test(group1_values, group2_values)$p.value
    )

    c(
      log2foldchange = log2_fold_change,
      p.value = p_value
    )
  })

  results <- as.data.frame(t(results))
  results$FDR <- p.adjust(results$p.value, method = "fdr")
  results$gene <- rownames(results)

  results %>%
    mutate(
      direction = case_when(
        p.value < 0.05 & log2foldchange > 1 ~ "up",
        p.value < 0.05 & log2foldchange < -1 ~ "down",
        TRUE ~ "ns"
      )
    )
}

dea_result <- run_geoMx_DEA(
  gset = geoMx_subset,
  group_col = "infection_region_simple",
  group1 = "Uninfected",
  group2 = "Infected",
  assay = "q_norm"
)

write.csv(
  dea_result,
  file.path(output_dir, "Figure5B_DEA_HSV1_Antigen_Pos_vs_Neg.csv"),
  row.names = FALSE
)

label_genes <- c(
  "Trim30a",
  "Oasl2",
  "Ifit1",
  "Ifit3",
  "Irgm1",
  "Irgm2",
  "Samhd1",
  "Gbp3",
  "Rev3l",
  "Smim22",
  "Ctf1",
  "Hoxd9",
  "Map3k21",
  "Tmprss2",
  "Dcps",
  "Hikeshi",
  "Moxd1"
)

label_data <- dea_result %>%
  filter(gene %in% label_genes)

volcano_plot <- ggplot(
  dea_result,
  aes(
    x = log2foldchange,
    y = -log10(p.value)
  )
) +
  geom_point(
    aes(color = direction),
    shape = 16,
    alpha = 0.6,
    size = 2
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    size = 4,
    box.padding = 0.35,
    point.padding = 0.2,
    segment.color = "black",
    segment.size = 0.5,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_color_manual(
    values = c(
      "down" = "#67A9CF",
      "ns" = "#BDBDBD",
      "up" = "#C95F5F"
    ),
    breaks = c("down", "ns", "up"),
    labels = c("down", "ns", "up")
  ) +
  labs(
    x = expression(Log[2] ~ Fold ~ Change),
    y = expression(-Log[10] ~ p - value),
    color = "Direction"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(
  filename = file.path(
    output_dir,
    "Volcano_HSV1_Antigen_Pos_vs_Neg.pdf"
  ),
  plot = volcano_plot,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "Volcano_HSV1_Antigen_Pos_vs_Neg.png"
  ),
  plot = volcano_plot,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

message(
  "Saved Figure 5B to: ",
  normalizePath(output_dir)
)
