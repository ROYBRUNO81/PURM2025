library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(tibble)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)

# Read saved object
combined <- readRDS("combined_mouse_colon_final.rds")

combined$condition <- ifelse(
  startsWith(colnames(combined), "HC_"),
  "HC",
  "DSS"
)

table(combined$condition)

# Extract the numeric suffix (1,2,3) and combine with condition
combined$mouse <- with(
  data.frame(barcode = colnames(combined), condition = combined$condition),
  paste0(
    condition,
    sub(".*-(\\d+)$", "\\1", barcode)
  )
)

# Quick sanity check
table(combined$mouse)

# Aggregate raw counts by cell type × mouse sample
agg <- AggregateExpression(
  object        = combined,
  assays        = "RNA",
  slot          = "counts",
  group.by      = c("cell_type", "mouse"),
  return.seurat = FALSE
)
counts_pbulk <- agg$RNA

# Build sample metadata (colData) for DESeq2
coldata <- do.call(rbind, lapply(colnames(counts_pbulk), function(col) {
  parts <- strsplit(col, "_")[[1]]
  data.frame(
    cell_type = parts[1],
    mouse     = parts[2],
    condition = ifelse(grepl("^HC", parts[2]), "HC", "DSS"),
    row.names = col,
    stringsAsFactors = FALSE
  )
}))

# Ensure the rows line up
stopifnot(nrow(coldata) == ncol(counts_pbulk))
coldata$condition <- factor(coldata$condition, levels = c("HC", "DSS"))

# Round counts to integers
counts_pbulk_int <- round(counts_pbulk)
stopifnot(all(counts_pbulk_int == floor(counts_pbulk_int)))

# Initialize a list to hold each cell‐type’s results
all_res <- list()

# Loop over cell types, run DESeq2, and store results with cell_type column
for (ct in unique(coldata$cell_type)) {
  # Subset counts & metadata
  keep      <- coldata$cell_type == ct
  ct_counts <- counts_pbulk_int[, keep, drop = FALSE]
  ct_meta   <- coldata[keep, , drop = FALSE]
  
  # Build DESeqDataSet and run DE
  dds <- DESeqDataSetFromMatrix(
    countData = ct_counts,
    colData   = ct_meta,
    design    = ~ condition
  )
  dds <- DESeq(dds)
  
  # Extract and shrink results
  res <- results(dds, contrast = c("condition", "DSS", "HC"))
  res_shrunk <- lfcShrink(dds,
                          contrast = c("condition", "DSS", "HC"),
                          type     = "ashr")
  
  # Convert to data.frame and add gene + cell_type columns
  df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("gene") %>%
    mutate(cell_type = ct)
  
  # Store
  all_res[[ct]] <- df
}

# Bind all into one big data.frame
final_res <- do.call(rbind, all_res)

head(final_res)
rownames(final_res) <- NULL
head(final_res)

# Write to a single CSV without rownames
write.csv(final_res,
          file      = "DE_pseudobulk_all_celltypes_DSS_vs_HC.csv",
          row.names = FALSE)

# Filter for significance and effect size
sig_genes <- final_res %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%   # padj < 0.05 & |LFC| > 1
  arrange(cell_type, desc(abs(log2FoldChange)))      # rank by magnitude of change

sig_genes <- as_tibble(sig_genes)

#compute a log‐normalized pseudobulk matrix:
pb_data <- AggregateExpression(
  combined,
  group.by      = c("cell_type","mouse"),
  assay         = "RNA",
  slot          = "data",
  return.seurat = FALSE
)$RNA
pb_log <- log2(pb_data + 1)

# Make output directories
dir.create("MAplots",      showWarnings = FALSE)
dir.create("VolcanoPlots", showWarnings = FALSE)
dir.create("Heatmaps",     showWarnings = FALSE)

# Loop over each cell type
for (ct in unique(coldata$cell_type)) {
  # 1) Subset the pseudobulk counts & metadata
  keep      <- coldata$cell_type == ct
  ct_counts <- counts_pbulk_int[, keep, drop = FALSE]
  ct_meta   <- coldata[keep, , drop = FALSE]
  
  # 2) Build & run DESeq2 for this cell type
  dds <- DESeqDataSetFromMatrix(
    countData = ct_counts,
    colData   = ct_meta,
    design    = ~ condition      # model: ~ HC vs DSS
  )
  dds   <- DESeq(dds)            # fit negative-binomial GLM
  res   <- results(dds,          # raw DE results
                   contrast = c("condition","DSS","HC"))
  
  # 3) MA Plot via plotMA()
  #    • main: plot title
  #    • ylim: y-axis limits for log2 fold changes
  #    • alpha: significance cutoff (default α = 0.1)
  png(file.path("MAplots", paste0("plotMA_", ct, ".png")),
      width = 700, height = 500)
  plotMA(res,
         main = paste("MA Plot —", ct),
         ylim = c(-3, 3),
         alpha = 0.05)            # color points with padj < 0.05
  dev.off()
  
  # 4) Volcano Plot via EnhancedVolcano()
  #    • lab: row labels (genes)
  #    • x/y: which columns to use for axes
  #    • pCutoff: padj threshold
  #    • FCcutoff: |log2FC| threshold
  #    • pointSize/textSize: plotting accents
  # Capture the EnhancedVolcano plot in an object
  volcano_plot <- EnhancedVolcano(res,
                                  lab            = rownames(res),
                                  x              = 'log2FoldChange',
                                  y              = 'padj',
                                  pCutoff        = 0.05,
                                  FCcutoff       = 1,
                                  title          = paste("Volcano —", ct),
                                  subtitle       = "",
                                  pointSize      = 2.5,
                                  labSize        = 3.0,
                                  legendPosition = 'right',
                                  col            = c('grey30','forestgreen','royalblue','red2'))
  
  # Then save it
  ggsave(filename = file.path("VolcanoPlots", paste0("Volcano_", ct, ".png")),
         plot     = volcano_plot,
         width    = 12, height = 10, dpi = 300)
  
  # 5) Heatmap of Top 10 DE Genes (by |log2FC|):
  top10 <- as.data.frame(res) %>%
    na.omit() %>%
    tibble::rownames_to_column("gene") %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice_head(n = 10) %>%
    pull(gene)
  
  # Subset log-normalized pseudobulk matrix
  mat  <- pb_log[top10, grep(paste0("^", ct, "_"), colnames(pb_log))]
  mat_z <- t(scale(t(mat)))  # z-score per gene
  
  # Plot & save
  pheatmap(mat_z,
           cluster_cols  = FALSE,
           show_colnames = FALSE,
           main          = paste(ct, "Top 10 DE Genes"),
           filename      = file.path("Heatmaps", paste0("Heatmap_", ct, ".png")),
           width         = 6, height = 5)
}

# Saving seurat object
saveRDS(
  combined,
  file = "combined_mouse_colon_final.rds"
)



