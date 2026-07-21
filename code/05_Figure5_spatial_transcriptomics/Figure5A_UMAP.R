suppressPackageStartupMessages({
  library(GeomxTools)
  library(Biobase)
  library(umap)
  library(ggplot2)
  library(dplyr)
})

input <- "data/05_Figure5_spatial_transcriptomics/processed/GeoMx_WTA_processed_manuscript.rds"
output_dir <- "outputs/05_Figure5_spatial_transcriptomics/Figure5A_UMAP"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

geoMx <- readRDS(input)

# Reproduce the UMAP calculation used for the manuscript
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

umap_input <- t(log2(assayDataElement(geoMx, "q_norm")))

umap_output <- umap::umap(
  umap_input,
  config = custom_umap
)

pData(geoMx)$UMAP1 <- umap_output$layout[, 1]
pData(geoMx)$UMAP2 <- umap_output$layout[, 2]

# Build plotting data
umap_df <- pData(geoMx) %>%
  as.data.frame() %>%
  mutate(
    condition_plot = recode(
      as.character(condition),
      "HSV" = "HSV-1",
      "Mock" = "Mock"
    ),
    antigen_plot = case_when(
      antigen == "HSV1 positive" ~ "Infected",
      antigen == "HSV1 negative" ~ "Uninfected",
      TRUE ~ NA_character_
    ),
    condition_plot = factor(
      condition_plot,
      levels = c("HSV-1", "Mock")
    ),
    antigen_plot = factor(
      antigen_plot,
      levels = c("Infected", "Uninfected")
    )
  )

# Plot
umap_plot <- ggplot(
  umap_df,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = condition_plot,
    shape = antigen_plot
  )
) +
  geom_point(size = 3.8, alpha = 0.75) +
  scale_color_manual(
    values = c(
      "HSV-1" = "#B33259",
      "Mock" = "#BDBDBD"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "Infected" = 16,
      "Uninfected" = 17
    ),
    drop = FALSE
  ) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "Condition",
    shape = "HSV-1 antigen"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key.height = grid::unit(0.65, "cm"),
    legend.spacing.y = grid::unit(0.35, "cm"),
    plot.margin = margin(8, 8, 8, 8)
  ) +
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(
        shape = 16,
        size = 4,
        alpha = 1
      )
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        color = "black",
        size = 4,
        alpha = 1
      )
    )
  ) +
  coord_cartesian(clip = "off")

# Save
ggsave(
  filename = file.path(output_dir, "UMAP_condition_antigen.pdf"),
  plot = umap_plot,
  width = 7,
  height = 5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(output_dir, "UMAP_condition_antigen.png"),
  plot = umap_plot,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)

message(
  "Saved Figure 5A to: ",
  normalizePath(output_dir)
)
