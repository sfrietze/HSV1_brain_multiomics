#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

option_list <- list(
  make_option(
    "--seurat",
    type = "character",
    default = "data/03_add_final_annotations/combined_seurat_final_annotated.rds"
  ),
  make_option(
    "--ding",
    type = "character",
    default = NA
  ),
  make_option(
    "--outdir",
    type = "character",
    default = "outputs/02_Figure2_multiome_annotation/SuppFig5_Ding_comparison"
  ),
  make_option(
    "--assay",
    type = "character",
    default = "CombinedRNA"
  ),
  make_option(
    "--condition_col",
    type = "character",
    default = "condition"
  ),
  make_option(
    "--cluster_col",
    type = "character",
    default = "annotated_clusters"
  ),
  make_option(
    "--hsv_label",
    type = "character",
    default = "HSV1"
  ),
  make_option(
    "--mock_label",
    type = "character",
    default = "Mock"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$seurat)) {
  stop("Seurat object not found: ", opt$seurat)
}

if (is.na(opt$ding) || !file.exists(opt$ding)) {
  stop("Provide the Ding et al. Table S3 file using --ding.")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

as_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

check_cols <- function(x, required, object_name) {
  missing <- setdiff(required, colnames(x))

  if (length(missing) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

message("Loading Seurat object: ", opt$seurat)
obj <- readRDS(opt$seurat)

if (!opt$assay %in% Assays(obj)) {
  stop(
    "Assay '",
    opt$assay,
    "' not found. Available assays: ",
    paste(Assays(obj), collapse = ", ")
  )
}

check_cols(
  obj@meta.data,
  c(opt$condition_col, opt$cluster_col),
  "Seurat metadata"
)

DefaultAssay(obj) <- opt$assay

micro_states <- c(
  "Homeostatic Microglia",
  "Transiently Activated Microglia",
  "IFN-Responsive Microglia",
  "Primed Microglia",
  "Mitochondrial-Activated Microglia",
  "IEG-High Microglia"
)

micro_cells <- rownames(obj@meta.data)[
  obj@meta.data[[opt$cluster_col]] %in% micro_states
]

if (length(micro_cells) == 0) {
  stop("No cells matched the expected microglial states.")
}

micro_obj <- subset(obj, cells = micro_cells)
Idents(micro_obj) <- opt$condition_col

message("Microglial cells by condition:")
print(table(Idents(micro_obj)))

message("Loading Ding Table S3: ", opt$ding)

ding_celltype <- read_xlsx(
  opt$ding,
  sheet = "celltypeDEGs.fig2"
)

check_cols(
  ding_celltype,
  c(
    "gene",
    "day",
    "cellType",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj"
  ),
  "Ding cell-type table"
)

ding_celltype <- ding_celltype %>%
  mutate(
    avg_log2FC = as_num(avg_log2FC),
    pct.1 = as_num(pct.1),
    pct.2 = as_num(pct.2),
    p_val_adj = as_num(p_val_adj)
  )

ding_d6 <- ding_celltype %>%
  filter(
    cellType == "Microglia",
    day == "iD6"
  )

ding_up <- ding_d6 %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0,
    pct.1 > pct.2
  ) %>%
  pull(gene) %>%
  unique()

message("Ding D6 microglia upregulated genes: ", length(ding_up))

our_de <- FindMarkers(
  micro_obj,
  ident.1 = opt$hsv_label,
  ident.2 = opt$mock_label,
  assay = opt$assay,
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0
) %>%
  rownames_to_column("gene")

our_up <- our_de %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0
  ) %>%
  pull(gene) %>%
  unique()

shared_genes <- intersect(ding_up, our_up)

message("Current study upregulated genes: ", length(our_up))
message("Shared upregulated genes: ", length(shared_genes))

if (length(shared_genes) == 0) {
  stop("No shared upregulated genes were identified.")
}

ding_shared <- ding_d6 %>%
  filter(gene %in% shared_genes) %>%
  transmute(
    Gene = gene,
    Ding_log2FC = avg_log2FC,
    Ding_pct_HSV = pct.1,
    Ding_pct_Mock = pct.2,
    Ding_padj = p_val_adj
  )

our_shared <- our_de %>%
  filter(gene %in% shared_genes) %>%
  transmute(
    Gene = gene,
    CurrentStudy_log2FC = avg_log2FC,
    CurrentStudy_pct_HSV = pct.1,
    CurrentStudy_pct_Mock = pct.2,
    CurrentStudy_padj = p_val_adj
  )

avg_shared <- AverageExpression(
  micro_obj,
  assays = opt$assay,
  features = shared_genes,
  group.by = opt$cluster_col,
  slot = "data"
)[[opt$assay]]

avg_shared <- avg_shared[
  shared_genes,
  micro_states,
  drop = FALSE
]

avg_shared_df <- as.data.frame(avg_shared) %>%
  rownames_to_column("Gene")

shared_table <- ding_shared %>%
  left_join(our_shared, by = "Gene") %>%
  left_join(avg_shared_df, by = "Gene") %>%
  arrange(desc(CurrentStudy_log2FC))

write_csv(
  shared_table,
  file.path(
    opt$outdir,
    "Supplementary_Table_S5A_conserved_HSV_induced_genes.csv"
  )
)

shared_scaled <- t(scale(t(avg_shared)))
shared_scaled[is.na(shared_scaled)] <- 0

shared_scaled <- shared_scaled[
  order(
    shared_scaled[, "IFN-Responsive Microglia"],
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

expression_colors <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

heatmap_a <- Heatmap(
  shared_scaled,
  name = "Scaled\nexpression",
  col = expression_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 90,
  column_title = "Microglial states",
  row_title = "Conserved HSV-induced genes",
  border = TRUE,
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2)
  )
)

pdf(
  file.path(
    opt$outdir,
    "SuppFig5A_Ding_conserved_HSV_induced_microglial_genes.pdf"
  ),
  width = 6.5,
  height = 10
)

draw(
  heatmap_a,
  heatmap_legend_side = "right"
)

dev.off()

ding_subpop <- read_xlsx(
  opt$ding,
  sheet = "subpopulation.fig3.4.6.7"
)

check_cols(
  ding_subpop,
  c(
    "gene",
    "celltype",
    "cluster",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj"
  ),
  "Ding subpopulation table"
)

ding_subpop <- ding_subpop %>%
  mutate(
    cluster = as.character(cluster),
    avg_log2FC = as_num(avg_log2FC),
    pct.1 = as_num(pct.1),
    pct.2 = as_num(pct.2),
    p_val_adj = as_num(p_val_adj)
  ) %>%
  filter(celltype == "Microglia")

get_signature <- function(cluster_id) {
  ding_subpop %>%
    filter(
      cluster == as.character(cluster_id),
      p_val_adj < 0.05,
      avg_log2FC > 0.5,
      pct.1 > pct.2
    ) %>%
    arrange(desc(avg_log2FC)) %>%
    pull(gene) %>%
    unique() %>%
    intersect(rownames(obj[[opt$assay]]))
}

ding_signatures <- list(
  "Homeostatic\ncluster 0" = get_signature(0),
  "Activated innate\ncluster 6" = get_signature(6),
  "Antiviral/IFN\ncluster 14" = get_signature(14)
)

message("Ding signature sizes:")
print(lengths(ding_signatures))

signature_table <- bind_rows(
  lapply(
    names(ding_signatures),
    function(signature_name) {
      tibble(
        Ding_signature = signature_name,
        Gene = ding_signatures[[signature_name]]
      )
    }
  )
)

write_csv(
  signature_table,
  file.path(
    opt$outdir,
    "Supplementary_Table_S5B_Ding_signature_genes.csv"
  )
)

avg_all <- AverageExpression(
  micro_obj,
  assays = opt$assay,
  group.by = opt$cluster_col,
  slot = "data"
)[[opt$assay]]

avg_all <- avg_all[
  ,
  micro_states,
  drop = FALSE
]

signature_matrix <- sapply(
  ding_signatures,
  function(genes) {
    colMeans(
      avg_all[
        genes,
        micro_states,
        drop = FALSE
      ]
    )
  }
)

signature_matrix <- t(signature_matrix)

colnames(signature_matrix) <- c(
  "Homeo",
  "Transient",
  "IFN",
  "Primed",
  "Mito",
  "IEG"
)

signature_scaled <- t(scale(t(signature_matrix)))
signature_scaled[is.na(signature_scaled)] <- 0

scaled_colors <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

raw_colors <- colorRamp2(
  c(
    min(signature_matrix),
    mean(range(signature_matrix)),
    max(signature_matrix)
  ),
  c("#2166AC", "white", "#B2182B")
)

heatmap_b1 <- Heatmap(
  signature_scaled,
  name = "Row-scaled\nprogram\nexpression",
  col = scaled_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45,
  column_title = "Program mapping",
  border = TRUE,
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2)
  )
)

heatmap_b2 <- Heatmap(
  signature_matrix,
  name = "Mean\nprogram\nexpression",
  col = raw_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45,
  column_title = "Mean expression",
  border = TRUE
)

pdf(
  file.path(
    opt$outdir,
    "SuppFig5B_Ding_signature_mapping.pdf"
  ),
  width = 10,
  height = 3.5
)

draw(
  heatmap_b1 + heatmap_b2,
  heatmap_legend_side = "right"
)

dev.off()

summary_table <- tibble(
  metric = c(
    "Ding D6 microglia upregulated genes",
    "Current study upregulated genes",
    "Shared upregulated genes",
    "Ding homeostatic cluster 0 signature genes",
    "Ding activated innate cluster 6 signature genes",
    "Ding antiviral/IFN cluster 14 signature genes"
  ),
  value = c(
    length(ding_up),
    length(our_up),
    length(shared_genes),
    lengths(ding_signatures)
  )
)

write_csv(
  summary_table,
  file.path(
    opt$outdir,
    "SuppFig5_Ding_comparison_summary.csv"
  )
)

session_lines <- sub(
  "[[:space:]]+$",
  "",
  capture.output(sessionInfo())
)

writeLines(
  session_lines,
  file.path(
    opt$outdir,
    "sessionInfo_SuppFig5_Ding_comparison.txt"
  )
)

message("Done: ", normalizePath(opt$outdir))
