#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ReactomePA)
  library(ComplexHeatmap)
  library(circlize)
  library(dendextend)
  library(RColorBrewer)
  library(grid)
})

option_list <- list(
  make_option("--indir", type = "character", default = "data/03_Figure3_DEG_pathways/pathway_inputs"),
  make_option("--outdir", type = "character", default = "outputs/Figure3"),
  make_option("--cap", type = "double", default = 8),
  make_option("--k", type = "integer", default = 8)
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

sym2entrez <- function(symbols) {
  symbols <- unique(na.omit(symbols))
  if (length(symbols) == 0) return(data.frame())
  suppressMessages(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  ) %>% distinct(ENTREZID, .keep_all = TRUE)
}

run_all_ora <- function(gene_symbols, universe_symbols, cluster, direction) {
  genes <- sym2entrez(gene_symbols)$ENTREZID
  universe <- sym2entrez(universe_symbols)$ENTREZID
  if (length(genes) < 10) return(NULL)

  out <- bind_rows(
    tryCatch({
      enrichKEGG(gene = genes, universe = universe, organism = "mmu") %>%
        as.data.frame() %>% mutate(gene_set_db = "KEGG")
    }, error = function(e) NULL),
    tryCatch({
      enrichPathway(gene = genes, universe = universe, organism = "mouse") %>%
        as.data.frame() %>% mutate(gene_set_db = "Reactome")
    }, error = function(e) NULL),
    tryCatch({
      enrichGO(gene = genes, universe = universe, OrgDb = org.Mm.eg.db, ont = "BP") %>%
        as.data.frame() %>% mutate(gene_set_db = "GO_BP")
    }, error = function(e) NULL)
  )

  if (is.null(out) || nrow(out) == 0) return(NULL)

  out %>%
    mutate(
      cluster = cluster,
      direction = direction,
      cluster_dir = paste(cluster, direction, sep = "_"),
      Pathway = paste0(Description, " [", gene_set_db, "]"),
      logp = -log10(p.adjust)
    ) %>%
    filter(!is.na(p.adjust), is.finite(logp), logp > 0)
}

clusters <- c(
  "Infiltrating_Macrophages",
  "Transiently_Activated_Microglia",
  "Homeostatic_Microglia"
)

ora_all <- list()

for (cl in clusters) {
  up <- read.csv(file.path(opt$indir, paste0(cl, "_DE_sig_UP.csv")))
  down <- read.csv(file.path(opt$indir, paste0(cl, "_DE_sig_DOWN.csv")))
  full <- read.csv(file.path(opt$indir, paste0(cl, "_DE_full.csv")))

  ora_all[[paste0(cl, "_UP")]] <- run_all_ora(up$gene, full$gene, cl, "UP")
  ora_all[[paste0(cl, "_DOWN")]] <- run_all_ora(down$gene, full$gene, cl, "DOWN")
}

ifn_file <- file.path(opt$indir, "IFN-Responsive_Microglia_MARKERS_genes.txt")
background_file <- file.path(opt$indir, "Infiltrating_Macrophages_DE_full.csv")

if (file.exists(ifn_file)) {
  ifn_genes <- readLines(ifn_file)
  background <- read.csv(background_file)

  ora_all[["IFN_Responsive_Microglia_UP"]] <- run_all_ora(
    gene_symbols = ifn_genes,
    universe_symbols = background$gene,
    cluster = "IFN_Responsive_Microglia",
    direction = "UP"
  )
}

combined <- bind_rows(ora_all) %>%
  mutate(logp = as.numeric(logp)) %>%
  filter(!is.na(logp), is.finite(logp))

write.csv(combined, file.path(opt$outdir, "Figure3E_ORA_all_results.csv"), row.names = FALSE)

heat_df <- combined %>%
  dplyr::select(Pathway, cluster_dir, logp) %>%
  group_by(Pathway, cluster_dir) %>%
  summarise(logp = max(logp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = cluster_dir, values_from = logp, values_fill = list(logp = 0))

mat <- as.matrix(heat_df[, -1])
rownames(mat) <- heat_df$Pathway
storage.mode(mat) <- "numeric"

desired_col_order <- c(
  "Infiltrating_Macrophages_UP",
  "Homeostatic_Microglia_UP",
  "Transiently_Activated_Microglia_UP",
  "Infiltrating_Macrophages_DOWN",
  "Homeostatic_Microglia_DOWN",
  "Transiently_Activated_Microglia_DOWN",
  "IFN_Responsive_Microglia_UP"
)

desired_col_order <- desired_col_order[desired_col_order %in% colnames(mat)]
mat_ordered <- mat[, desired_col_order, drop = FALSE]
mat_ordered <- pmin(mat_ordered, opt$cap)

row_dend <- hclust(dist(mat_ordered), method = "ward.D2")
row_dend_cut <- cutree(row_dend, k = opt$k)

dend <- as.dendrogram(row_dend)
dend <- dend %>%
  dendextend::set("branches_k_color", k = opt$k) %>%
  dendextend::set("branches_lwd", 1.2) %>%
  dendextend::ladderize()

cluster_ids <- as.factor(row_dend_cut)
cluster_colors <- setNames(brewer.pal(opt$k, "Set2"), levels(cluster_ids))

row_ha <- rowAnnotation(
  Cluster = cluster_ids,
  col = list(Cluster = cluster_colors),
  show_annotation_name = TRUE,
  annotation_legend_param = list(title = "Cluster")
)

pdf(file.path(opt$outdir, "Figure3E_pathway_heatmap.pdf"), width = 10, height = 10)

Heatmap(
  mat_ordered,
  name = "-log10(p.adj)",
  cluster_rows = dend,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  row_dend_side = "left",
  row_split = opt$k,
  left_annotation = row_ha,
  col = colorRamp2(c(0, 2, opt$cap), c("white", "#41b6c4", "#045a8d")),
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45,
  heatmap_legend_param = list(
    title_position = "topcenter",
    legend_direction = "horizontal"
  )
)

dev.off()

write.csv(
  data.frame(Pathway = names(row_dend_cut), Cluster = row_dend_cut),
  file.path(opt$outdir, "Figure3E_pathway_clusters.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(mat_ordered) %>% rownames_to_column("Pathway"),
  file.path(opt$outdir, "Figure3E_heatmap_matrix.csv"),
  row.names = FALSE
)

sink(file.path(opt$outdir, "sessionInfo_Figure3E_pathway_heatmap.txt"))
sessionInfo()
sink()
