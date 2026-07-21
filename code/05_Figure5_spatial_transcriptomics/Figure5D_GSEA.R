suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Mm.eg.db)
  library(msigdbr)
  library(AnnotationDbi)
})

dea_file <- "outputs/05_Figure5_spatial_transcriptomics/Figure5B_volcano/Figure5B_DEA_HSV1_Antigen_Pos_vs_Neg.csv"
output_dir <- "outputs/05_Figure5_spatial_transcriptomics/Figure5D_GSEA"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==== Load DEA Results ====
dea_result <- read_csv(dea_file, show_col_types = FALSE)

required_columns <- c("gene", "log2foldchange")

if (!all(required_columns %in% colnames(dea_result))) {
  stop(
    "DEA file must contain columns: ",
    paste(required_columns, collapse = ", ")
  )
}

ranking_df <- dea_result %>%
  filter(
    !is.na(gene),
    !is.na(log2foldchange)
  ) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(desc(log2foldchange))

geneList <- setNames(
  ranking_df$log2foldchange,
  ranking_df$gene
)

geneList <- sort(geneList, decreasing = TRUE)

# ==== Convert Symbols to ENTREZ IDs ====
gene_df_mapped <- data.frame(
  symbol = names(geneList),
  log2FC = as.numeric(geneList),
  stringsAsFactors = FALSE
) %>%
  left_join(
    AnnotationDbi::select(
      org.Mm.eg.db,
      keys = names(geneList),
      columns = "ENTREZID",
      keytype = "SYMBOL"
    ),
    by = c("symbol" = "SYMBOL")
  ) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

geneList_entrez <- setNames(
  gene_df_mapped$log2FC,
  gene_df_mapped$ENTREZID
)

geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

# ==== Run GSEA ====
# Use pvalueCutoff = 1 so curated pathways are not removed before plotting.

msig_hallmark <- msigdbr(
  species = "Mus musculus",
  collection = "H"
)

term2gene_hallmark <- msig_hallmark %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  filter(!is.na(ncbi_gene)) %>%
  mutate(ncbi_gene = as.character(ncbi_gene)) %>%
  distinct()

gsea_hallmark <- GSEA(
  geneList = geneList_entrez,
  TERM2GENE = term2gene_hallmark,
  pvalueCutoff = 1,
  verbose = FALSE
)

gsea_kegg <- gseKEGG(
  geneList = geneList_entrez,
  organism = "mmu",
  pvalueCutoff = 1,
  verbose = FALSE
)

gsea_reactome <- gsePathway(
  geneList = geneList_entrez,
  organism = "mouse",
  pvalueCutoff = 1,
  verbose = FALSE
)

gsea_gobp <- gseGO(
  geneList = geneList_entrez,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pvalueCutoff = 1,
  verbose = FALSE
)

# ==== Export Full Supplemental GSEA Tables ====
write_csv(
  as.data.frame(gsea_hallmark@result),
  file.path(output_dir, "Supplementary_GSEA_Hallmark.csv")
)

write_csv(
  as.data.frame(gsea_kegg@result),
  file.path(output_dir, "Supplementary_GSEA_KEGG.csv")
)

write_csv(
  as.data.frame(gsea_reactome@result),
  file.path(output_dir, "Supplementary_GSEA_Reactome.csv")
)

write_csv(
  as.data.frame(gsea_gobp@result),
  file.path(output_dir, "Supplementary_GSEA_GO_BP.csv")
)

# ==== Combine All Results ====
add_database_label <- function(gsea_object, database_name) {
  as.data.frame(gsea_object@result) %>%
    filter(!is.na(NES)) %>%
    mutate(
      sign = ifelse(NES > 0, "activated", "suppressed"),
      db = database_name
    )
}

dot_data <- bind_rows(
  add_database_label(gsea_hallmark, "Hallmark"),
  add_database_label(gsea_kegg, "KEGG"),
  add_database_label(gsea_reactome, "Reactome"),
  add_database_label(gsea_gobp, "GO:BP")
)

# ==== Curated Pathways Used in Figure 5D ====
selected_pathways <- c(
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "Cytosolic DNA-sensing pathway - Mus musculus (house mouse)",
  "Herpes simplex virus 1 infection - Mus musculus (house mouse)",
  "NOD-like receptor signaling pathway - Mus musculus (house mouse)",
  "cellular response to type II interferon",
  "Nonsense-Mediated Decay (NMD)",
  "L13a-mediated translational silencing of Ceruloplasmin expression",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "O-linked glycosylation of mucins",
  "Keratinization",
  "glycolipid biosynthetic process",
  "glycosphingolipid biosynthetic process",
  "Amine ligand-binding receptors",
  "Mitochondrial translation",
  "positive regulation of glial cell migration",
  "Phospholipid metabolism",
  "Antimicrobial peptides",
  "response to dopamine"
)

selected_dot_data <- dot_data %>%
  filter(Description %in% selected_pathways) %>%
  mutate(
    gene_count = lengths(strsplit(core_enrichment, "/")),
    GeneRatio = gene_count / setSize,
    display_label = recode(
      Description,
      "HALLMARK_IL6_JAK_STAT3_SIGNALING" =
        "IL6 JAK STAT3 signaling",
      "HALLMARK_CHOLESTEROL_HOMEOSTASIS" =
        "cholesterol homeostasis",
      "HALLMARK_INFLAMMATORY_RESPONSE" =
        "inflammatory response",
      "HALLMARK_MYC_TARGETS_V1" =
        "MYC targets V1",
      "HALLMARK_IL2_STAT5_SIGNALING" =
        "IL2 STAT5 signaling",
      "Herpes simplex virus 1 infection - Mus musculus (house mouse)" =
        "Herpes simplex virus 1 infection",
      "Cytosolic DNA-sensing pathway - Mus musculus (house mouse)" =
        "Cytosolic DNA-sensing pathway",
      "NOD-like receptor signaling pathway - Mus musculus (house mouse)" =
        "NOD-like receptor signaling",
      "cellular response to type II interferon" =
        "Response to type II interferon",
      "O-linked glycosylation of mucins" =
        "O-linked glycosylation\nof mucins",
      "positive regulation of glial cell migration" =
        "positive regulation of glial\ncell migration",
      "glycosphingolipid biosynthetic process" =
        "glycosphingolipid\nbiosynthesis",
      "glycolipid biosynthetic process" =
        "glycolipid biosynthesis",
      "Amine ligand-binding receptors" =
        "Amine ligand-binding\nreceptors",
      "Mitochondrial translation" =
        "Mitochondrial\ntranslation"
    ),
    sign = factor(
      sign,
      levels = c("activated", "suppressed")
    )
  )

missing_pathways <- setdiff(
  selected_pathways,
  selected_dot_data$Description
)

if (length(missing_pathways) > 0) {
  warning(
    "The following curated pathways were not found in the current GSEA results:\n",
    paste(missing_pathways, collapse = "\n")
  )
}



# ==== Plot ====
combined_gsea <- ggplot(
  selected_dot_data,
  aes(
    x = GeneRatio,
    y = reorder(display_label, GeneRatio)
  )
) +
  geom_point(
    aes(
      size = gene_count,
      color = p.adjust
    )
  ) +
  scale_color_gradient(
    low = "red",
    high = "blue",
    limits = c(0, 0.25),
    breaks = c(0.05, 0.10, 0.15, 0.20, 0.25),
    oob = scales::squish,
    name = "Adjusted p-value"
  ) +
  scale_size(
    name = "Gene Count",
    range = c(3, 10)
  ) +
  scale_x_continuous(
    limits = c(0, 0.8),
    breaks = c(0, 0.4, 0.8),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  facet_wrap(
    ~sign,
    scales = "free_y",
    nrow = 1
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill = "gray90",
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      color = "black",
      size = 9
    ),
    axis.text.y = element_text(
      color = "black",
      size = 9
    ),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(
      color = "black",
      fill = NA
    ),
    panel.spacing.x = grid::unit(1.5, "lines"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  guides(
    color = guide_colorbar(
      order = 1,
      barheight = grid::unit(3, "cm")
    ),
    size = guide_legend(order = 2)
  ) +
  labs(
    title = NULL,
    x = "Gene Ratio",
    y = NULL
  )

ggsave(
  file.path(output_dir, "Figure5D_GSEA.pdf"),
  combined_gsea,
  width = 10,
  height = 8,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "Figure5D_GSEA.png"),
  combined_gsea,
  width = 10,
  height = 8,
  dpi = 300
)

message(
  "Saved Figure 5D and supplemental GSEA tables to: ",
  normalizePath(output_dir)
)
