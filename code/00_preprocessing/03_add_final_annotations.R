#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--outdir", type = "character", default = "data")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(opt$input)

new_annotations <- c(
  "0"  = "Homeostatic Microglia",
  "1"  = "Transiently Activated Microglia",
  "2"  = "IFN-Responsive Microglia",
  "3"  = "Primed Microglia",
  "4"  = "Mitochondrial-Activated Microglia",
  "5"  = "IEG-High Microglia",
  "6"  = "Infiltrating Macrophages",
  "7"  = "CNS Endothelial Cells",
  "8"  = "CD8+ T Cells",
  "9"  = "Cycling Myeloid Progenitors",
  "10" = "Vascular Smooth Muscle Cells"
)

obj$annotated_clusters <- unname(new_annotations[as.character(obj$seurat_clusters)])

if (any(is.na(obj$annotated_clusters))) {
  stop("Cluster annotation failed.")
}

Idents(obj) <- "annotated_clusters"

saveRDS(obj, file.path(opt$outdir, "combined_seurat_final_annotated.rds"))
write.csv(obj@meta.data, file.path(opt$outdir, "combined_seurat_final_annotated_metadata.csv"))

sink(file.path(opt$outdir, "sessionInfo_03_add_final_annotations.txt"))
sessionInfo()
sink()
