#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Signac)
  library(Seurat)
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(ggforce)
  library(ssvTracks)
  library(seqsetvis)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--bw_dir", type = "character"),
  make_option("--gtf", type = "character"),
  make_option("--outdir", type = "character",
              default = "outputs/02_Figure2_multiome_annotation/Figure2F_locus_tracks")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(opt$input)
DefaultAssay(obj) <- "ATAC"

track_order <- c(
  "Homeostatic_Microglia",
  "Transiently_Activated_Microglia",
  "IFN-Responsive_Microglia",
  "Primed_Microglia",
  "Mitochondrial-Activated_Microglia",
  "IEG-High_Microglia",
  "Infiltrating_Macrophages"
)

track_colors <- c(
  "Homeostatic_Microglia" = "darkorange",
  "Transiently_Activated_Microglia" = "#fcae91",
  "IFN-Responsive_Microglia" = "#4B0082",
  "Primed_Microglia" = "#cb181d",
  "Mitochondrial-Activated_Microglia" = "#99000d",
  "IEG-High_Microglia" = "#fb6a4a",
  "Infiltrating_Macrophages" = "#67000d"
)

bw_files <- file.path(opt$bw_dir, paste0(track_order, ".bw"))
names(bw_files) <- track_order
stopifnot(all(file.exists(bw_files)))

cfg_bw <- data.table(
  file = bw_files,
  id = names(bw_files),
  group = names(bw_files),
  color = track_colors[names(bw_files)],
  nspline = 50,
  nMovingAverage = 3,
  ceiling_value = 175000
)

make_region <- function(chr, start, end, genome_build = "mm10") {
  gr <- GRanges(chr, IRanges(start, end))
  genome(gr) <- genome_build
  gr
}

make_arc_track <- function(region_gr, link_gr, zscore_cutoff = 0, gene_filter = NULL) {
  link_df <- as.data.table(link_gr)

  link_df[, c("peak_chr", "peak_start", "peak_end") := tstrsplit(peak, "-", fixed = TRUE)]
  link_df[, `:=`(
    peak_start = as.integer(peak_start),
    peak_end = as.integer(peak_end)
  )]

  link_subset <- link_df[
    seqnames == as.character(seqnames(region_gr)) &
      (
        (start <= end(region_gr) & end >= start(region_gr)) |
          (peak_start <= end(region_gr) & peak_end >= start(region_gr))
      )
  ]

  if (!is.null(gene_filter)) {
    link_subset <- link_subset[gene %in% gene_filter]
  }

  link_subset <- link_subset[zscore >= zscore_cutoff]

  if (nrow(link_subset) == 0) {
    return(ggplot() + theme_void())
  }

  link_subset[, gene_mid := (start + end) / 2]
  link_subset[, peak_mid := (peak_start + peak_end) / 2]

  arc_data <- rbindlist(lapply(seq_len(nrow(link_subset)), function(i) {
    d <- link_subset[i]
    data.table(
      x = c(d$peak_mid, mean(c(d$peak_mid, d$gene_mid)), d$gene_mid),
      y = c(0, 1200 + 300 * i, 0),
      arc_id = i
    )
  }))

  ggplot() +
    geom_bezier(
      data = arc_data,
      aes(x = x, y = y, group = arc_id),
      color = "firebrick",
      linewidth = 0.7
    ) +
    scale_y_continuous(limits = c(-100, max(arc_data$y) + 500)) +
    theme_void()
}

plot_tracks <- function(region_gr, output_path, y_max, zscore_cutoff, gene_filter) {
  color_mapping <- setNames(cfg_bw$color, cfg_bw$id)

  p_signal <- track_chip(
    signal_files = cfg_bw,
    query_gr = region_gr,
    fetch_fun = ssvFetchBigwig,
    fill_VAR = "id",
    color_VAR = "id",
    fill_mapping = color_mapping,
    color_mapping = color_mapping,
    nspline = 50,
    nMovingAverage = 3,
    ceiling_value = y_max,
    y_label = "ATAC signal"
  )

  p_links <- make_arc_track(
    region_gr = region_gr,
    link_gr = Links(obj[["ATAC"]]),
    zscore_cutoff = zscore_cutoff,
    gene_filter = gene_filter
  )

  p_gene <- track_gene_reference(
    ref = opt$gtf,
    query_gr = region_gr,
    show_tss = TRUE,
    tss_arrow_size = 0.2
  )

  p <- assemble_tracks(
    list(p_signal, p_links, p_gene),
    query_gr = region_gr,
    rel_heights = c(4, 1, 1)
  )

  ggsave(output_path, p, width = 14, height = 4.5, device = cairo_pdf)
}

plot_tracks(
  region_gr = make_region("chr7", 106100000, 106245000),
  output_path = file.path(opt$outdir, "Figure2F_Gvin1_locus_tracks.pdf"),
  y_max = 50000,
  zscore_cutoff = 0,
  gene_filter = c("Gvin1", "Gm4070")
)

plot_tracks(
  region_gr = make_region("chr8", 70675000, 70723000),
  output_path = file.path(opt$outdir, "Figure2F_Jund_locus_tracks.pdf"),
  y_max = 175000,
  zscore_cutoff = 5,
  gene_filter = "Jund"
)

sink(file.path(opt$outdir, "sessionInfo_Figure2F_locus_tracks.txt"))
sessionInfo()
sink()
