#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(data.table)
  library(dplyr)
  library(seqsetvis)
  library(chiptsne2)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(grid)
})

option_list <- list(
  make_option("--dar_regions", type = "character"),
  make_option("--bw_dir", type = "character"),
  make_option("--outdir", type = "character",
              default = "outputs/02_Figure2_multiome_annotation/Figure2E_DAR_ATAC_signal"),
  make_option("--view_size", type = "integer", default = 2500),
  make_option("--window_size", type = "integer", default = 50),
  make_option("--signal_scale", type = "double", default = 10000)
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cluster_order <- c(
  "Homeostatic Microglia",
  "IEG-High Microglia",
  "IFN-Responsive Microglia",
  "Mitochondrial-Activated Microglia",
  "Primed Microglia",
  "Transiently Activated Microglia",
  "Infiltrating Macrophages"
)

cluster_order_simple <- make.names(cluster_order)

cluster_colors_simple <- c(
  "Homeostatic.Microglia" = "#ff8c00",
  "IEG.High.Microglia" = "#fb6a4a",
  "IFN.Responsive.Microglia" = "#54278f",
  "Mitochondrial.Activated.Microglia" = "#8b0000",
  "Primed.Microglia" = "#000000",
  "Transiently.Activated.Microglia" = "#fcae91",
  "Infiltrating.Macrophages" = "#1f78b4"
)

dar_regions <- fread(opt$dar_regions)
dar_regions[, cluster_simple := make.names(cluster_assignment)]
dar_regions <- dar_regions[cluster_simple %in% cluster_order_simple]
dar_regions[, cluster_simple := factor(cluster_simple, levels = cluster_order_simple)]

olaps_union <- GRanges(
  seqnames = dar_regions$seqnames,
  ranges = IRanges(start = dar_regions$start, end = dar_regions$end)
)
names(olaps_union) <- dar_regions$region_id

region_metadata_df <- data.frame(
  cluster_simple = as.character(dar_regions$cluster_simple),
  row.names = dar_regions$region_id
)

bw_files <- list.files(opt$bw_dir, pattern = "\\.bw$", full.names = TRUE)
names(bw_files) <- tools::file_path_sans_ext(basename(bw_files))
bw_files <- bw_files[!grepl("_Mock$|_HSV1$|^NA$", names(bw_files))]
names(bw_files) <- gsub("_", " ", names(bw_files))
bw_files <- bw_files[cluster_order]
stopifnot(all(cluster_order %in% names(bw_files)))

fetch_cfg <- FetchConfig.from_files(
  bw_files,
  read_mode = "bigwig",
  view_size = opt$view_size,
  window_size = opt$window_size
)

ct2_obj <- ChIPtsne2.from_FetchConfig(
  fetch_config = fetch_cfg,
  query_gr = olaps_union,
  region_metadata = region_metadata_df
)

colData(ct2_obj)$sample_label <- gsub("\\.bw$", "", basename(colData(ct2_obj)$file))
colnames(assay(ct2_obj)) <- rownames(colData(ct2_obj))

rowData(ct2_obj)$cluster_simple <- factor(
  rowData(ct2_obj)$cluster_simple,
  levels = cluster_order_simple
)

ct2_obj <- sortRegions(
  ct2_obj,
  sort_strategy = "sort",
  group_VAR = "cluster_simple"
)

ht_data <- plotSignalHeatmap(
  ct2_obj,
  group_VARS = "cluster_simple",
  heatmap_fill_limits = c(0, 50000),
  max_rows = length(olaps_union),
  heatmap_colors = c("#577AB2", "#FFFFE0", "#BC412B"),
  return_data = TRUE
)

ht_plot_df <- ht_data %>%
  filter(position >= -1000, position <= 1000) %>%
  mutate(
    name = factor(as.character(name), levels = cluster_order),
    cluster_simple = factor(cluster_simple, levels = cluster_order_simple),
    value_plot = pmin(value / opt$signal_scale, 5),
    region = factor(region, levels = levels(ht_data$region))
  ) %>%
  filter(!is.na(name), !is.na(cluster_simple))

bar_df <- ht_plot_df %>%
  distinct(region, cluster_simple) %>%
  mutate(x = 1)

p_cluster_bar <- ggplot(bar_df, aes(x = x, y = region, fill = cluster_simple)) +
  geom_tile() +
  facet_grid(rows = vars(cluster_simple), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = cluster_colors_simple, guide = "none") +
  theme_void() +
  theme(
    strip.text.y = element_blank(),
    panel.spacing.y = unit(0.02, "lines"),
    plot.margin = margin(5, 0, 25, 5)
  )

p_heatmap <- ggplot(ht_plot_df, aes(x = position, y = region, fill = value_plot)) +
  geom_raster() +
  facet_grid(rows = vars(cluster_simple), cols = vars(name), scales = "free_y", space = "free_y") +
  scale_fill_gradientn(
    colors = viridisLite::viridis(100),
    limits = c(0, 5),
    oob = scales::squish,
    name = "normalized\nATAC signal"
  ) +
  scale_x_continuous(breaks = c(-1000, 0, 1000), labels = c("-1kb", "0", "+1kb")) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.4, color = "black"),
    axis.line.x = element_line(linewidth = 0.5, color = "black"),
    axis.line.y = element_blank(),
    strip.text.x = element_text(angle = 45, hjust = 0, size = 9),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.spacing.x = unit(0.18, "lines"),
    panel.spacing.y = unit(0.02, "lines"),
    legend.position = "bottom",
    panel.border = element_blank()
  ) +
  labs(x = "position from peak center (±1kb)")

p_heatmap_with_bar <- p_cluster_bar + p_heatmap + plot_layout(widths = c(0.08, 1))

line_data <- plotSignalLinePlot(
  ct2_obj,
  group_VAR = "sample_label",
  color_VAR = "cluster_simple",
  facet_VAR = NULL,
  extra_VARS = c("cluster_simple", "sample_label"),
  n_splines = 5,
  moving_average_window = 3,
  return_data = TRUE
)

line_plot_df <- line_data %>%
  mutate(
    sample = factor(as.character(name), levels = cluster_order),
    cluster_simple = factor(as.character(cluster_simple), levels = cluster_order_simple),
    value_plot = pmin(value / opt$signal_scale, 5)
  ) %>%
  filter(!is.na(sample), !is.na(cluster_simple))

line_plot <- ggplot(line_plot_df, aes(x = position, y = value_plot, color = cluster_simple)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey45") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey85") +
  geom_line(linewidth = 0.8) +
  facet_grid(rows = vars(cluster_simple), cols = vars(sample), scales = "free_y") +
  scale_color_manual(values = cluster_colors_simple, guide = "none") +
  scale_y_continuous(limits = c(0, 5), breaks = c(0, 3, 5)) +
  scale_x_continuous(breaks = c(-1000, 0, 1000), labels = c("-1kb", "0", "+1kb")) +
  theme_classic(base_size = 8) +
  theme(
    strip.text.x = element_text(angle = 45, hjust = 0, size = 8),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.4, color = "black"),
    axis.line.x = element_line(linewidth = 0.5, color = "black"),
    axis.line.y = element_line(linewidth = 0.5, color = "black"),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.18, "lines"),
    panel.spacing.y = unit(0.08, "lines"),
    legend.position = "none"
  ) +
  labs(
    x = "position from peak center (±1kb)",
    y = "normalized ATAC signal (RPM)"
  )

ggsave(
  file.path(opt$outdir, "Figure2E_DAR_ATAC_signal_heatmap.pdf"),
  p_heatmap_with_bar,
  width = 5.8,
  height = 4.2,
  units = "in"
)

ggsave(
  file.path(opt$outdir, "Figure2E_DAR_ATAC_signal_lineplot.pdf"),
  line_plot,
  width = 5.8,
  height = 4.2,
  units = "in"
)

write.csv(ht_plot_df, file.path(opt$outdir, "Figure2E_heatmap_signal_data.csv"), row.names = FALSE)
write.csv(line_plot_df, file.path(opt$outdir, "Figure2E_lineplot_signal_data.csv"), row.names = FALSE)

sink(file.path(opt$outdir, "sessionInfo_Figure2E_DAR_ATAC_signal.txt"))
sessionInfo()
sink()
