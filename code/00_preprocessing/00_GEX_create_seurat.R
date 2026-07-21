#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

option_list <- list(
  make_option("--mock_path", type = "character"),
  make_option("--hsv1_path", type = "character"),
  make_option("--outdir", type = "character", default = "data"),
  make_option("--resolution", type = "double", default = 0.45)
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(123)

get_gex_counts <- function(path) {
  x <- Read10X(data.dir = path)

  if (inherits(x, "dgCMatrix")) {
    return(x)
  }

  if (is.list(x) && "Gene Expression" %in% names(x)) {
    return(x[["Gene Expression"]])
  }

  if (is.list(x) && "RNA" %in% names(x)) {
    return(x[["RNA"]])
  }

  stop("Could not identify Gene Expression matrix in: ", path)
}

mock_counts <- get_gex_counts(opt$mock_path)
hsv1_counts <- get_gex_counts(opt$hsv1_path)

mock <- CreateSeuratObject(
  counts = mock_counts,
  project = "Mock",
  min.cells = 0,
  min.features = 200
)

hsv1 <- CreateSeuratObject(
  counts = hsv1_counts,
  project = "HSV1",
  min.cells = 0,
  min.features = 200
)

mock$condition <- "Mock"
hsv1$condition <- "HSV1"

all_genes <- union(rownames(mock), rownames(hsv1))
mock <- mock[all_genes, ]
hsv1 <- hsv1[all_genes, ]

obj <- merge(
  mock,
  y = hsv1,
  add.cell.ids = c("Mock", "HSV1"),
  project = "Mouse_HSV1"
)

obj <- JoinLayers(obj, assay = "RNA")

obj[["CombinedRNA"]] <- CreateAssayObject(
  counts = LayerData(obj, assay = "RNA", layer = "counts")
)

DefaultAssay(obj) <- "CombinedRNA"

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

obj <- subset(
  obj,
  subset = nFeature_CombinedRNA > 200 &
    nFeature_CombinedRNA < 2500 &
    percent.mt < 10
)

obj <- NormalizeData(obj, assay = "CombinedRNA")
obj <- FindVariableFeatures(obj, assay = "CombinedRNA")
obj <- ScaleData(obj, vars.to.regress = "percent.mt", assay = "CombinedRNA")
obj <- RunPCA(obj, assay = "CombinedRNA")
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = opt$resolution)
obj <- RunUMAP(obj, dims = 1:10, seed.use = 123)

saveRDS(
  obj,
  file.path(opt$outdir, "combined_seurat_regenerated_clustered.rds")
)

write.csv(
  obj@meta.data,
  file.path(opt$outdir, "combined_seurat_regenerated_metadata.csv")
)

sink(file.path(opt$outdir, "sessionInfo_00_GEX_create_seurat.txt"))
sessionInfo()
sink()
