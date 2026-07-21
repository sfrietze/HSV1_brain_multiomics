#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

ereg <- read.delim(
  "data/04_Figure4_SCENIC_footprinting/processed/eRegulons_extended.tsv",
  stringsAsFactors = FALSE
)

dea <- read.csv(
  "data/05_Figure5_spatial_transcriptomics/processed/HSV_Inf_vs_Uninf_DEA.csv",
  stringsAsFactors = FALSE
)

dea <- dea %>%
  filter(
    !is.na(gene),
    !is.na(log2foldchange),
    !is.na(p.value)
  )

spatial_up <- dea %>%
  filter(
    p.value < 0.05,
    log2foldchange > 0.5
  ) %>%
  pull(gene) %>%
  unique()

all_genes <- unique(dea$gene)
tf_list <- unique(ereg$TF)

enrichment_up <- lapply(tf_list, function(tf) {

  tf_targets <- ereg %>%
    filter(TF == tf) %>%
    pull(Gene) %>%
    unique()

  a <- length(intersect(tf_targets, spatial_up))
  b <- length(tf_targets) - a
  c <- length(spatial_up) - a
  d <- length(setdiff(all_genes, union(tf_targets, spatial_up)))

  fisher <- fisher.test(
    matrix(c(a, b, c, d), nrow = 2)
  )

  data.frame(
    TF = tf,
    overlap = a,
    odds_ratio = unname(fisher$estimate),
    pval = fisher$p.value
  )
}) %>%
  bind_rows() %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  arrange(FDR)

outdir <- "outputs/05_Figure5_spatial_transcriptomics/Figure5G_regulon_enrichment"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  enrichment_up,
  file.path(outdir, "Figure5G_regulon_enrichment_results.csv"),
  row.names = FALSE
)

tf_order <- c(
  "Klf6",
  "Etv6",
  "Egr1",
  "Mbd2",
  "Stat1",
  "Irf1",
  "Stat3",
  "Cebpb",
  "Stat2"
)

plot_df <- enrichment_up %>%
  filter(TF %in% tf_order) %>%
  mutate(
    TF = factor(TF, levels = rev(tf_order)),
    logFDR = -log10(FDR)
  )

p <- ggplot(plot_df, aes(x = odds_ratio, y = TF)) +
  geom_point(aes(size = overlap, color = logFDR)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  scale_color_gradient(
    low = "#FDB863",
    high = "#B2182B"
  ) +
  scale_size_continuous(
    range = c(3, 10),
    breaks = c(20, 40, 60, 80)
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Odds Ratio (Enrichment Strength)",
    y = NULL,
    color = expression(-log[10](FDR)),
    size = "Overlap Genes"
  ) +
  theme(
    axis.text.y = element_text(face = "plain"),
    legend.position = "right"
  )

ggsave(
  file.path(outdir, "Figure5G_regulon_enrichment.pdf"),
  p,
  width = 6,
  height = 5.5,
  device = cairo_pdf
)

print(
  plot_df %>%
    select(TF, overlap, odds_ratio, pval, FDR)
)
