# ============================================================
#  RNA-seq Analysis — CTLA4 Blockade with Ipilimumab
#  Dataset: GSE114716 (GEO)
#  CD4+ T cells — Melanoma patients
#  Baseline vs. Post-Ipilimumab (3 doses)
#  4 patients × 2 timepoints = 8 samples
#
#  Master's Course: Transcriptomics of RNA-seq
#  Verified pipeline — run step by step in class
# ============================================================

# ============================================================
# 1. Install packages (run once)
# ============================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install from CRAN (Comprehensive R Archive Network)
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readxl")

# Install from Bioconductor
BiocManager::install("DESeq2")
BiocManager::install("GEOquery")

# ============================================================
# 2. Load packages
# ============================================================
library(DESeq2)
library(GEOquery)
library(readxl)
library(ggplot2)
library(ggrepel)

# ============================================================
# 3. Download GEO metadata and raw counts
# ============================================================
gse_id <- "GSE114716"

gse <- getGEO(gse_id, GSEMatrix = TRUE)
meta <- pData(gse[[1]])

getGEOSuppFiles(gse_id, makeDirectory = TRUE)

counts_raw <- read_excel(
  "GSE114716/GSE114716_raw.counts.hs.xlsx"
)

# ============================================================
# 4. Prepare count matrix
# ============================================================
counts <- as.matrix(counts_raw[, -1])
rownames(counts) <- make.unique(as.character(counts_raw[[1]]))

storage.mode(counts) <- "numeric"
counts <- round(counts)
storage.mode(counts) <- "integer"

colnames(counts) <- as.character(names(counts_raw)[-1])

meta$description <- as.character(meta$description)
meta <- meta[match(colnames(counts), meta$description), ]

stopifnot(
  !anyNA(meta$description),
  identical(colnames(counts), meta$description)
)

# ============================================================
# 5. Define conditions and paired patients
# ============================================================
meta$condition <- ifelse(
  grepl("baseline", meta$title, ignore.case = TRUE),
  "Baseline", "Post_Ipi"
)

meta$condition <- factor(
  meta$condition,
  levels = c("Baseline", "Post_Ipi")
)

meta$patient <- factor(
  ceiling(as.numeric(meta$description) / 2)
)

rownames(meta) <- meta$description

stopifnot(
  all(table(meta$patient, meta$condition) == 1)
)

# ============================================================
# 6. Build DESeq2 dataset
# ============================================================
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ patient + condition
)

dds <- dds[rowSums(counts(dds) >= 10) >= 4, ]
dds <- DESeq(dds)

# ============================================================
# 7. Extract results
# ============================================================
res <- results(
  dds,
  contrast = c("condition", "Post_Ipi", "Baseline"),
  alpha    = 0.05
)

res_df        <- as.data.frame(res)
res_df$gene   <- rownames(res_df)

plot_df        <- res_df[!is.na(res_df$padj), ]
plot_df$group  <- "Not significant"

plot_df$group[
  plot_df$padj < 0.05 & plot_df$log2FoldChange >= 1
] <- "Upregulated"

plot_df$group[
  plot_df$padj < 0.05 & plot_df$log2FoldChange <= -1
] <- "Downregulated"

# --- Full DEG counts (source of truth) ---
# These numbers reflect the real DEG list used in GO analysis
n_up   <- sum(plot_df$group == "Upregulated")
n_down <- sum(plot_df$group == "Downregulated")

cat("=== DEG Summary (full list) ===\n")
cat("Upregulated genes   (padj < 0.05 & log2FC >= 1):", n_up,   "\n")
cat("Downregulated genes (padj < 0.05 & log2FC <= -1):", n_down, "\n")
cat("Total significant DEGs:", n_up + n_down, "\n\n")

# ============================================================
# 8. Select top 15 up/down for volcano plot labels only
# ============================================================
# Note: these are used ONLY for labeling — not for GO analysis
up_label   <- head(plot_df[plot_df$group == "Upregulated",  ][order(plot_df[plot_df$group == "Upregulated",  ]$padj), ], 15)
down_label <- head(plot_df[plot_df$group == "Downregulated",][order(plot_df[plot_df$group == "Downregulated",]$padj), ], 15)

top_genes          <- rbind(up_label, down_label)
top_genes$negLogP  <- -log10(pmax(top_genes$padj, 1e-300))

cat("Genes labeled in volcano plot:", nrow(top_genes),
    "(top 15 up + top 15 down by padj — labels only)\n\n")

print(top_genes[, c("gene", "log2FoldChange", "padj")])

# ============================================================
# 9. Volcano plot - Visualization
# ============================================================
plot_df$negLogP <- -log10(pmax(plot_df$padj, 1e-300))

ggplot(plot_df, aes(x = log2FoldChange, y = negLogP)) +
  geom_point(aes(color = group), alpha = 0.7, size = 1.5) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(
    "Not significant" = "grey70",
    "Upregulated"     = "red",
    "Downregulated"   = "blue"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title    = "RNA-seq: Post-Ipilimumab vs Baseline",
    subtitle = "GSE114716 | DESeq2 | FDR < 0.05 and |log2FC| >= 1",
    x        = "Log2 fold change",
    y        = "-Log10 adjusted p-value",
    color    = "Expression"
  ) +
  theme_classic()

ggsave("Volcano_DESeq2_GSE114716.pdf", width = 10, height = 8)

# ============================================================
# 10. Save all DESeq2 results
# ============================================================
write.csv(
  res_df[order(res_df$padj), ],
  "DESeq2_results_GSE114716.csv",
  row.names = FALSE
)

# End of script. Part 1.

# ============================================================
#  Functional Enrichment — GO Analysis + GSEA ranked list
#  Continuation of the DESeq2 pipeline (GSE114716)
# ============================================================

# Install from Bioconductor
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")

# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# ============================================================
# 11. Gene Ontology (GO) Enrichment Analysis
# ============================================================
# Using the FULL DEG lists from plot_df (not the top 15 label subset)

up_genes   <- plot_df[plot_df$group == "Upregulated",   ]
down_genes <- plot_df[plot_df$group == "Downregulated",  ]

cat("=== GO input — full DEG lists ===\n")
cat("Upregulated genes for GO:", nrow(up_genes), "\n")
cat("Downregulated genes for GO:", nrow(down_genes), "\n\n")

# Convert gene symbols to Entrez IDs
up_entrez <- bitr(
  up_genes$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

down_entrez <- bitr(
  down_genes$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

cat("Upregulated mapped to Entrez:", nrow(up_entrez), "\n")
cat("Downregulated mapped to Entrez:", nrow(down_entrez), "\n\n")

# Run GO enrichment
go_up <- enrichGO(
  gene          = up_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_down <- enrichGO(
  gene          = down_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cat("=== GO Upregulated Summary ===\n")
print(head(as.data.frame(go_up)[, c("Description", "pvalue", "p.adjust")], 10))

cat("\n=== GO Downregulated Summary ===\n")
print(head(as.data.frame(go_down)[, c("Description", "pvalue", "p.adjust")], 10))

# Dotplots — top 5 terms per ontology
dotplot(go_up, showCategory = 5, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  scale_color_gradient(low = "red", high = "#ffcccc") +
  labs(
    title    = "GO Enrichment — Upregulated genes",
    subtitle = "GSE114716 | Top 5 terms per ontology | BH-adjusted p < 0.05"
  )

ggsave("GO_dotplot_upregulated_GSE114716.pdf", width = 8, height = 5)

dotplot(go_down, showCategory = 5, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  scale_color_gradient(low = "blue", high = "#cce0ff") +
  labs(
    title    = "GO Enrichment — Downregulated genes",
    subtitle = "GSE114716 | Top 5 terms per ontology | BH-adjusted p < 0.05"
  )

ggsave("GO_dotplot_downregulated_GSE114716.pdf", width = 8, height = 5)

# Save GO results
write.csv(as.data.frame(go_up),   "GO_upregulated_results_GSE114716.csv",   row.names = FALSE)
write.csv(as.data.frame(go_down), "GO_downregulated_results_GSE114716.csv",  row.names = FALSE)

cat("GO results saved.\n")

# ============================================================
# 12. Ranked list for GSEA
# ============================================================
rnk_gsea <- data.frame(
  gene  = rownames(res_df),
  score = -log10(res_df$padj) * sign(res_df$log2FoldChange),
  stringsAsFactors = FALSE
)

rnk_gsea <- rnk_gsea[!is.na(rnk_gsea$score), ]

max_score <- max(rnk_gsea$score[is.finite(rnk_gsea$score)])
min_score <- min(rnk_gsea$score[is.finite(rnk_gsea$score)])

rnk_gsea$score[rnk_gsea$score ==  Inf] <- max_score
rnk_gsea$score[rnk_gsea$score == -Inf] <- min_score

rnk_gsea <- rnk_gsea[order(rnk_gsea$score, decreasing = TRUE), ]

cat("\nRanked list preview (top 5 and bottom 5):\n")
print(head(rnk_gsea, 5))
print(tail(rnk_gsea, 5))
cat("Total genes in ranked list:", nrow(rnk_gsea), "\n")

write.table(
  rnk_gsea,
  file      = "GSE114716_DESeq2.rnk",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep       = "\t"
)

cat("Ranked list saved to GSE114716_DESeq2.rnk\n")

# End of script. Part 2.
# ============================================================