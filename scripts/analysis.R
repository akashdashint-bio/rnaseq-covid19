# =============================================================================
# RNA-seq Differential Expression Analysis: COVID-19 vs Healthy Controls
# Dataset: GSE152641 (Thair et al., 2021, iScience)
# Author: Akash Dash (@AkashDash-bio)
# =============================================================================


# -----------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# -----------------------------------------------------------------------------

library(GEOquery)       # download data from NCBI GEO
library(DESeq2)         # differential expression analysis
library(ggplot2)        # plotting
library(ggrepel)        # non-overlapping volcano plot labels
library(pheatmap)       # heatmaps
library(org.Hs.eg.db)   # human gene annotation database
library(AnnotationDbi)  # gene ID conversion
library(clusterProfiler)# pathway enrichment analysis
library(R.utils)        # unzip .gz files


# -----------------------------------------------------------------------------
# 2. DOWNLOAD AND LOAD DATA
# -----------------------------------------------------------------------------

# Download metadata from GEO
gse <- getGEO("GSE152641", GSEMatrix = TRUE)
metadata_raw <- pData(gse[[1]])

# Download supplementary count matrix
getGEOSuppFiles("GSE152641")

# Unzip count matrix
gunzip("data/GSE152641_Inflammatix_COVID19_counts_entrez.csv.gz", remove = FALSE)

# Load count matrix (genes x samples)
counts <- read.csv("data/GSE152641_Inflammatix_COVID19_counts_entrez.csv",
                   row.names = 1)

# Load raw metadata
metadata_raw <- read.csv("data/metadata_raw.csv", row.names = 1)

# Inspect dimensions
dim(counts)        # expected: ~20000 genes x 86 samples
dim(metadata_raw)  # expected: 86 samples x 45 metadata columns


# -----------------------------------------------------------------------------
# 3. BUILD CLEAN METADATA
# -----------------------------------------------------------------------------

# Check available condition labels
unique(metadata_raw$characteristics_ch1)

# Build clean metadata dataframe with two conditions: COVID19 and control
metadata <- data.frame(
    sample_id = metadata_raw$title,
    condition = ifelse(
        metadata_raw$characteristics_ch1 == "disease: COVID19",
        "COVID19", "control"
    ),
    row.names = metadata_raw$title
)

# Verify sample counts per condition
table(metadata$condition)
# Expected: control = 24, COVID19 = 62

# Confirm sample names match between counts and metadata (must return TRUE)
all(colnames(counts) == rownames(metadata))


# -----------------------------------------------------------------------------
# 4. BUILD DESeq2 OBJECT AND FILTER
# -----------------------------------------------------------------------------

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = ~ condition   # compare by condition (COVID19 vs control)
)

# Filter out genes with very low counts across all samples (noise reduction)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dim(dds)  # expected: ~19936 genes x 86 samples


# -----------------------------------------------------------------------------
# 5. QUALITY CONTROL — PCA PLOT
# -----------------------------------------------------------------------------

# Variance-stabilising transformation for visualisation (blind = TRUE for QC)
vsd <- vst(dds, blind = TRUE)

# PCA plot — samples should cluster by condition if data is clean
pca_plot <- plotPCA(vsd, intgroup = "condition") +
    theme_bw() +
    labs(title = "PCA — COVID-19 vs Healthy Control",
         subtitle = "Whole blood transcriptomics (GSE152641)")

print(pca_plot)

# Save PCA plot
ggsave("results/figures/pca_plot.png", pca_plot, width = 8, height = 6, dpi = 300)


# -----------------------------------------------------------------------------
# 6. DIFFERENTIAL EXPRESSION ANALYSIS
# -----------------------------------------------------------------------------

# Run full DESeq2 pipeline:
# Step 1 — estimates size factors (normalises for sequencing depth)
# Step 2 — estimates dispersions (models gene-level variability)
# Step 3 — fits negative binomial model and tests each gene
dds <- DESeq(dds)

# Extract results: COVID19 vs control
res <- results(dds,
    contrast      = c("condition", "COVID19", "control"),
    alpha         = 0.05   # significance threshold
)

# Summary of results
summary(res)
# Expected: ~6294 upregulated, ~2900 downregulated

# Convert to dataframe and sort by adjusted p-value
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

# Filter significant DEGs (padj < 0.05, |log2FC| > 1)
sig_genes <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
nrow(sig_genes)  # total significant DEGs

# Save full results table
write.csv(res_df, "results/tables/DESeq2_results.csv")


# -----------------------------------------------------------------------------
# 7. VOLCANO PLOT
# -----------------------------------------------------------------------------

# Classify each gene by expression status
res_df$status <- "NS"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "Up"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

# Label top 15 most significant genes
res_df$label <- ""
res_df$label[1:15] <- res_df$gene[1:15]

# Plot
volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
    geom_point(alpha = 0.4, size = 1.2) +
    scale_color_manual(values = c("Up" = "#D85A30", "Down" = "#378ADD", "NS" = "#888780")) +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    theme_bw() +
    labs(
        title    = "COVID-19 vs Healthy Control",
        subtitle = "Whole blood transcriptomics (GSE152641)",
        x        = "Log2 Fold Change",
        y        = "-Log10 Adjusted P-value",
        color    = "Expression"
    )

print(volcano)

# Save volcano plot
ggsave("results/figures/volcano_plot.png", volcano, width = 8, height = 6, dpi = 300)


# -----------------------------------------------------------------------------
# 8. GENE ANNOTATION — CONVERT ENTREZ IDs TO GENE SYMBOLS
# -----------------------------------------------------------------------------

# Convert top 20 Entrez IDs to readable gene symbols
top_genes <- rownames(res_df)[1:20]

symbols <- mapIds(org.Hs.eg.db,
    keys      = top_genes,
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
)

# Print annotated top genes
data.frame(EntrezID = top_genes, Symbol = symbols)


# -----------------------------------------------------------------------------
# 9. PATHWAY ENRICHMENT ANALYSIS
# -----------------------------------------------------------------------------

# Extract significantly upregulated genes
sig_up <- rownames(subset(res_df, padj < 0.05 & log2FoldChange > 1))

# GO biological process enrichment
go_results <- enrichGO(
    gene          = sig_up,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",        # Biological Process
    pAdjustMethod = "BH",        # Benjamini-Hochberg multiple testing correction
    qvalueCutoff  = 0.05
)

# Plot top 15 enriched pathways
enrichment_plot <- dotplot(go_results,
    showCategory = 15,
    title        = "Top enriched pathways — COVID19 upregulated genes"
)

print(enrichment_plot)

# Save enrichment plot
ggsave("results/figures/pathway_enrichment.png", enrichment_plot,
       width = 8, height = 7, dpi = 300)

# Save enrichment results table
write.csv(as.data.frame(go_results), "results/tables/GO_enrichment.csv")


# =============================================================================
# END OF ANALYSIS
# Key findings:
#   - 9,168 significant DEGs (6,294 up, 2,900 down)
#   - Top gene: IFI27 (interferon response)
#   - Antiviral genes: APOBEC3A, APOBEC3A_B
#   - Cell cycle disruption: RRM2, CDC6, CLSPN, MCM10, POLQ
#   - Pathway enrichment: olfactory/sensory perception (anosmia signature),
#     antimicrobial and humoral immune response
# =============================================================================
