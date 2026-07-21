#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Signac)
  library(Seurat)
})

option_list <- list(
  make_option("--rna_rds", type = "character"),
  make_option("--atac_assay_rds", type = "character"),
  make_option("--outdir", type = "character", default = "data/02_scATAC_add_to_scRNA")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

rna_obj <- readRDS(opt$rna_rds)
merged_atac <- readRDS(opt$atac_assay_rds)

rna_obj[["ATAC"]] <- merged_atac
DefaultAssay(rna_obj) <- "ATAC"

rna_obj <- RunTFIDF(rna_obj)
rna_obj <- FindTopFeatures(rna_obj, min.cutoff = 20)
rna_obj <- RunSVD(rna_obj)

rna_obj <- LinkPeaks(
  object = rna_obj,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)

saveRDS(
  rna_obj,
  file.path(opt$outdir, "combined_multiome_with_ATAC_linked_peaks.rds")
)

sink(file.path(opt$outdir, "sessionInfo_02_scATAC_add_to_scRNA.txt"))
sessionInfo()
sink()
