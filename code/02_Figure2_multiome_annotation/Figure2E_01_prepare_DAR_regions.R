#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(data.table)
  library(dplyr)
  library(seqsetvis)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--outdir", type = "character",
              default = "data/02_Figure2_multiome_annotation/Figure2E_DAR_ATAC_signal"),
  make_option("--top_n", type = "integer", default = 1000)
)

opt <- parse_args(OptionParser(option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(opt$input)

DefaultAssay(obj) <- "ATAC"
Idents(obj) <- "annotated_clusters"

cluster_order <- c(
  "Homeostatic Microglia",
  "IEG-High Microglia",
  "IFN-Responsive Microglia",
  "Mitochondrial-Activated Microglia",
  "Primed Microglia",
  "Transiently Activated Microglia",
  "Infiltrating Macrophages"
)

cluster_dars <- FindAllMarkers(
  object = obj,
  assay = "ATAC",
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.25,
  test.use = "LR"
)

top_dars <- cluster_dars %>%
  filter(cluster %in% cluster_order) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = opt$top_n, with_ties = FALSE) %>%
  ungroup()

dar_dt <- as.data.table(top_dars)

dar_dt[, c("chr","start","end") := tstrsplit(gene, "-", fixed = TRUE)]
dar_dt[, `:=`(
  start = as.integer(start),
  end = as.integer(end)
)]

dar_gr_list <- split(dar_dt, by = "cluster", keep.by = FALSE)

dar_gr_list <- lapply(dar_gr_list, function(x) {
  GRanges(
    seqnames = x$chr,
    ranges = IRanges(x$start, x$end)
  )
})

dar_gr_list <- dar_gr_list[cluster_order]
dar_gr_list <- dar_gr_list[!vapply(dar_gr_list, is.null, logical(1))]

olaps_union <- trim(ssvOverlapIntervalSets(dar_gr_list))
names(olaps_union) <- paste0("region_", seq_along(olaps_union))

membership <- as.data.frame(mcols(olaps_union))

cluster_assignment <- apply(membership, 1, function(x) {
  cl <- names(x)[as.logical(x)]
  if (length(cl) == 1) {
    cl
  } else if (length(cl) > 1) {
    "Shared"
  } else {
    "Other"
  }
})

dar_regions <- data.frame(
  seqnames = as.character(seqnames(olaps_union)),
  start = start(olaps_union),
  end = end(olaps_union),
  region_id = names(olaps_union),
  cluster_assignment = cluster_assignment
)

write.csv(
  dar_regions,
  file.path(opt$outdir, "Figure2E_DAR_regions.csv"),
  row.names = FALSE,
  quote = TRUE
)

sink(file.path(opt$outdir,
  "sessionInfo_Figure2E_01_prepare_DAR_regions.txt"))
sessionInfo()
sink()
