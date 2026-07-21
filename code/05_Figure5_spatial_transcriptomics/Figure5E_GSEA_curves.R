suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
})

dea_file <- "outputs/05_Figure5_spatial_transcriptomics/Figure5B_volcano/Figure5B_DEA_HSV1_Antigen_Pos_vs_Neg.csv"
output_dir <- "outputs/05_Figure5_spatial_transcriptomics/Figure5E_GSEA_curves"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==== Load DEA ranking ====
dea_result <- read_csv(dea_file, show_col_types = FALSE)

ranking_df <- dea_result %>%
  filter(
    !is.na(gene),
    !is.na(log2foldchange)
  ) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(desc(log2foldchange))

gene_list_symbol <- setNames(
  ranking_df$log2foldchange,
  ranking_df$gene
)

# ==== Convert SYMBOL to ENTREZ ====
mapped_genes <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = names(gene_list_symbol),
  columns = "ENTREZID",
  keytype = "SYMBOL"
) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

mapped_genes$log2foldchange <- gene_list_symbol[mapped_genes$SYMBOL]

mapped_genes <- mapped_genes %>%
  filter(!is.na(log2foldchange)) %>%
  arrange(desc(log2foldchange))

gene_list_entrez <- setNames(
  mapped_genes$log2foldchange,
  mapped_genes$ENTREZID
)

gene_list_entrez <- sort(
  gene_list_entrez,
  decreasing = TRUE
)

# ==== Retrieve pathway gene sets ====
hallmark_sets <- msigdbr(
  species = "Mus musculus",
  collection = "H"
) %>%
  transmute(
    pathway = gs_name,
    gene = as.character(ncbi_gene)
  )

reactome_sets <- msigdbr(
  species = "Mus musculus",
  collection = "C2",
  subcollection = "CP:REACTOME"
) %>%
  transmute(
    pathway = gs_name,
    gene = as.character(ncbi_gene)
  )

# Check available pathway names if database naming changes
pathway_names <- unique(c(
  hallmark_sets$pathway,
  reactome_sets$pathway
))

find_pathway <- function(pattern) {
  hits <- grep(
    pattern,
    pathway_names,
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(hits) == 0) {
    stop("Could not find pathway matching: ", pattern)
  }

  hits[[1]]
}

cholesterol_id <- find_pathway(
  "^HALLMARK_CHOLESTEROL_HOMEOSTASIS$"
)

il6_id <- find_pathway(
  "^HALLMARK_IL6_JAK_STAT3_SIGNALING$"
)

mitochondrial_id <- find_pathway(
  "MITOCHONDRIAL_TRANSLATION$"
)

antimicrobial_id <- find_pathway(
  "ANTIMICROBIAL_PEPTIDES$"
)

selected_term2gene <- bind_rows(
  hallmark_sets %>%
    filter(pathway == cholesterol_id) %>%
    mutate(pathway = "Cholesterol homeostasis"),

  hallmark_sets %>%
    filter(pathway == il6_id) %>%
    mutate(pathway = "IL6 JAK STAT3 signaling"),

  reactome_sets %>%
    filter(pathway == mitochondrial_id) %>%
    mutate(pathway = "Mitochondrial translation"),

  reactome_sets %>%
    filter(pathway == antimicrobial_id) %>%
    mutate(pathway = "Antimicrobial peptides")
) %>%
  filter(!is.na(gene)) %>%
  distinct(pathway, gene)

cat("\nPathways used for Figure 5E:\n")
print(
  selected_term2gene %>%
    count(pathway, name = "gene_set_size"),
  row.names = FALSE
)

# ==== Run GSEA on the four displayed pathways ====
gsea_figure5e <- GSEA(
  geneList = gene_list_entrez,
  TERM2GENE = selected_term2gene,
  minGSSize = 1,
  maxGSSize = Inf,
  pvalueCutoff = 1,
  eps = 0,
  verbose = FALSE
)

write_csv(
  as.data.frame(gsea_figure5e@result),
  file.path(
    output_dir,
    "Figure5E_GSEA_pathway_results.csv"
  )
)

pathway_order <- c(
  "Cholesterol homeostasis",
  "IL6 JAK STAT3 signaling",
  "Mitochondrial translation",
  "Antimicrobial peptides"
)

missing_pathways <- setdiff(
  pathway_order,
  gsea_figure5e@result$Description
)

if (length(missing_pathways) > 0) {
  stop(
    "These pathways were not returned by GSEA: ",
    paste(missing_pathways, collapse = ", ")
  )
}


# ==== Plot and save directly ====
draw_figure5e <- function() {
  enrichplot::gseaplot2(
    gsea_figure5e,
    geneSetID = pathway_order,
    color = c(
      "#147DB3",
      "#C5252B",
      "#ED1E79",
      "#393394"
    ),
    base_size = 13,
    rel_heights = c(1.7, 0.45, 0.55),
    subplots = 1:3,
    pvalue_table = FALSE,
    ES_geom = "line"
  )
}

grDevices::cairo_pdf(
  filename = file.path(output_dir, "Figure5E_GSEA_curves.pdf"),
  width = 10,
  height = 7.5
)
draw_figure5e()
grDevices::dev.off()

grDevices::png(
  filename = file.path(output_dir, "Figure5E_GSEA_curves.png"),
  width = 10,
  height = 7.5,
  units = "in",
  res = 300,
  type = "cairo"
)
draw_figure5e()
grDevices::dev.off()

message(
  "Saved Figure 5E to: ",
  normalizePath(output_dir)
)

message(
  "PDF size: ",
  file.info(file.path(output_dir, "Figure5E_GSEA_curves.pdf"))$size,
  " bytes"
)

message(
  "PNG size: ",
  file.info(file.path(output_dir, "Figure5E_GSEA_curves.png"))$size,
  " bytes"
)
