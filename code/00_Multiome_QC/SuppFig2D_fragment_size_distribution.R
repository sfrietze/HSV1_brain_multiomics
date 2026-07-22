#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    paste(
      "Usage:",
      "Rscript code/00_Multiome_QC/SuppFig2D_fragment_size_distribution.R",
      "<Mock fragments.tsv.gz> <HSV1 fragments.tsv.gz>"
    )
  )
}

mock_fragment_file <- args[1]
hsv1_fragment_file <- args[2]

if (!file.exists(mock_fragment_file)) {
  stop("Mock fragment file not found: ", mock_fragment_file)
}

if (!file.exists(hsv1_fragment_file)) {
  stop("HSV1 fragment file not found: ", hsv1_fragment_file)
}

output_dir <- "outputs/00_Multiome_QC"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

frag_dt_mock <- fread(mock_fragment_file, skip = 52)
frag_dt_mock[, condition := "Mock"]

frag_dt_hsv1 <- fread(hsv1_fragment_file, skip = 52)
frag_dt_hsv1[, condition := "HSV1"]

frag_dt <- rbind(frag_dt_mock, frag_dt_hsv1)

frag_dt[, size := V3 - V2]

frag_dt <- frag_dt[
  size > 0 & size < 500
]

frag_dt[, size_bin := cut(
  size,
  breaks = seq(0, 500, by = 5),
  include.lowest = TRUE
)]

frag_dt_sum <- frag_dt[
  ,
  .N,
  by = .(condition, size_bin)
]

frag_dt_sum[, size := as.numeric(size_bin) * 5 - 2.5]

frag_dt_sum[, condition := factor(
  condition,
  levels = c("Mock", "HSV1")
)]

qc_colors <- c(
  "Mock" = "gray60",
  "HSV1" = "gray60"
)

pD <- ggplot(
  frag_dt_sum,
  aes(x = size, y = N, fill = condition)
) +
  geom_col(width = 5) +
  geom_vline(
    xintercept = c(60, 200),
    color = "gray80"
  ) +
  facet_wrap(
    ~condition,
    nrow = 1
  ) +
  scale_fill_manual(values = qc_colors) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  labs(
    x = "ATAC-seq fragment size (bp)",
    y = "Fragment count"
  )

ggsave(
  file.path(output_dir, "SuppFig2D_fragment_size_distribution.pdf"),
  pD,
  width = 10,
  height = 3.2
)

ggsave(
  file.path(output_dir, "SuppFig2D_fragment_size_distribution.png"),
  pD,
  width = 10,
  height = 3.2,
  dpi = 300
)
