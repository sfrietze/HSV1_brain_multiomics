#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

input_file <- "data/00_preprocessing/combined_seurat_regenerated_metadata.csv"
output_dir <- "outputs/00_Multiome_QC"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

qc <- read.csv(
  input_file,
  check.names = FALSE,
  row.names = 1
)

if (!"condition" %in% colnames(qc)) {
  stop(
    "No 'condition' column found. Available columns:\n",
    paste(colnames(qc), collapse = ", ")
  )
}

if (!"nCount_RNA" %in% colnames(qc)) {
  stop("No 'nCount_RNA' column found in Seurat metadata.")
}

qc$condition <- factor(
  qc$condition,
  levels = c("Mock", "HSV1")
)

qc <- qc[
  !is.na(qc$condition) &
    !is.na(qc$nCount_RNA) &
    qc$nCount_RNA > 0,
  ,
  drop = FALSE
]

make_rank_df <- function(df, condition_name) {

  umi <- df$nCount_RNA[df$condition == condition_name]
  umi <- sort(umi, decreasing = TRUE)

  data.frame(
    rank = seq_along(umi),
    total_umi = umi,
    condition = condition_name
  )
}

rank_df <- rbind(
  make_rank_df(qc, "Mock"),
  make_rank_df(qc, "HSV1")
)

rank_df$condition <- factor(
  rank_df$condition,
  levels = c("Mock", "HSV1")
)

cell_counts <- table(qc$condition)

message("Mock cells: ", cell_counts["Mock"])
message("HSV1 cells: ", cell_counts["HSV1"])


cols <- c(
  Mock = "#3B5BA9",
  HSV1 = "#FF2A1A"
)

p <- ggplot(
  rank_df,
  aes(
    x = rank,
    y = total_umi,
    colour = condition
  )
) +
  geom_line(linewidth = 0.9) +
  geom_hline(
    yintercept = 500,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "black"
  ) +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_log10(
    breaks = c(300, 1000, 3000),
    limits = c(200, NA),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_colour_manual(values = cols) +
  labs(
    title = "Barcode Rank (Knee) Plot",
    x = "Barcode Rank (log10)",
    y = "Total UMI Counts (log10)",
    colour = "Condition"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 10
    ),
    axis.text = element_text(
      colour = "black",
      size = 8
    ),
    axis.title = element_text(
      colour = "black",
      size = 9
    ),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = "right"
  )

ggsave(
  file.path(output_dir, "SuppFig2A_barcode_rank.pdf"),
  p,
  width = 5.2,
  height = 4
)

ggsave(
  file.path(output_dir, "SuppFig2A_barcode_rank.png"),
  p,
  width = 5.2,
  height = 4,
  dpi = 300
)

pA <- p
