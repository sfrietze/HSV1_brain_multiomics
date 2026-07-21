#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(data.table)
})

option_list <- list(
  make_option("--mock_fragments", type = "character"),
  make_option("--hsv1_fragments", type = "character"),
  make_option("--mock_peaks", type = "character"),
  make_option("--hsv1_peaks", type = "character"),
  make_option("--mock_metrics", type = "character"),
  make_option("--hsv1_metrics", type = "character"),
  make_option("--outdir", type = "character", default = "data/01_scATAC_create_assay")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "mock_fragments", "hsv1_fragments",
  "mock_peaks", "hsv1_peaks",
  "mock_metrics", "hsv1_metrics"
)

missing_args <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

set.seed(123)

mock_peaks <- read.table(opt$mock_peaks, col.names = c("seqnames", "start", "end"))
hsv1_peaks <- read.table(opt$hsv1_peaks, col.names = c("seqnames", "start", "end"))

combined_peaks <- reduce(c(
  makeGRangesFromDataFrame(mock_peaks),
  makeGRangesFromDataFrame(hsv1_peaks)
))

combined_peaks <- combined_peaks[width(combined_peaks) < 10000 & width(combined_peaks) > 20]

write.table(
  data.frame(
    seqnames = as.character(seqnames(combined_peaks)),
    start = start(combined_peaks),
    end = end(combined_peaks)
  ),
  file = file.path(opt$outdir, "merged_peak_set.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

mock_metrics <- fread(opt$mock_metrics)
hsv1_metrics <- fread(opt$hsv1_metrics)

mock_bcs <- mock_metrics[is_cell == 1]$barcode
hsv1_bcs <- hsv1_metrics[is_cell == 1]$barcode

frags_mock <- CreateFragmentObject(
  path = opt$mock_fragments,
  cells = mock_bcs
)

frags_hsv1 <- CreateFragmentObject(
  path = opt$hsv1_fragments,
  cells = hsv1_bcs
)

counts_mock <- FeatureMatrix(
  fragments = frags_mock,
  features = combined_peaks,
  cells = mock_bcs
)

counts_hsv1 <- FeatureMatrix(
  fragments = frags_hsv1,
  features = combined_peaks,
  cells = hsv1_bcs
)

mock_atac <- CreateChromatinAssay(
  counts = counts_mock,
  fragments = frags_mock
)

hsv1_atac <- CreateChromatinAssay(
  counts = counts_hsv1,
  fragments = frags_hsv1
)

merged_atac <- merge(
  hsv1_atac,
  y = mock_atac,
  add.cell.ids = c("HSV1", "Mock")
)

saveRDS(
  merged_atac,
  file.path(opt$outdir, "merged_atac_assay.rds")
)

sink(file.path(opt$outdir, "sessionInfo_01_scATAC_create_assay.txt"))
sessionInfo()
sink()
