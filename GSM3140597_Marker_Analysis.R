# Load required library
library(Seurat)
library(patchwork)
library(tibble)    # for rownames_to_column()
library(tidyr)
library(ggplot2)
library(scales)
library(dplyr)

# Read saved object and file
combined <- readRDS("combined_mouse_colon_final.rds")
markers <- read.csv("Mouse_Colon_cluster_markers.csv")

# Example for IL33
FeaturePlot(combined, features="Il33") + ggtitle("IL33 expression (cluster annotation)")

# Example for Sox6
FeaturePlot(combined, features="Sox6") + ggtitle("Sox6 expression (cluster annotation)")

# Define a named list of canonical markers per cell type
celltype_markers <- list(
  Epithelial        = c("Cdh1", "Epcam", "Alpi", "Muc2", "Vil1", "Lyz1", "Lgr5", "Olfm4", "Ccnd1"),
  PanImmune         = c("Ptprc"),
  BCell             = c("Ms4a1", "Cd19"),
  PlasmaCell        = c("Mzb1", "Ccr10", "Sdc1"),
  TCell             = c("Cd3e"),
  CD4_TCell         = c("Cd4"),
  CD8_TCell         = c("Cd8a", "Cd8b"),
  NKCell            = c("Nkg7", "Ncam1"),
  Myeloid           = c("Cd68", "Cd74"),
  Monocyte_Mac      = c("Cd14"),
  Dendritic         = c("Itgax"),
  MastCell          = c("Kit", "Il1rl1"),
  Treg              = c("Foxp3"),
  Endothelial       = c("Pecam1"),
  Lymphatic_Endo    = c("Lyve1"),
  Glial             = c("Plp1", "S100b"),
  Pericyte          = c("Rgs5", "Pdgfrb", "Adipoq"),
  SMC               = c("Myh11", "Actg2", "Myocd", "Des", "Acta2"),
  Fibroblast        = c("Pdgfra", "Sparc", "Dcn", "Lum", "Col1a1", "Col14a1", "Fgfr2", "Col3a1",
                        "Col4a5", "Col4a6", "Bmp5", "Bmp4", "Grem1", "Vcam1", "Ogn", "Mgp", "Sfrp2", "C3", "Dpt")
)

# Prepare to collect missing markers
missing_markers <- list()

# Output directory for plots
outdir <- "UMAP_celltype_markers"
if (!dir.exists(outdir)) dir.create(outdir)

# Loop over each cell‐type and its markers
for (ct in names(celltype_markers)) {
  genes <- celltype_markers[[ct]]
  
  # Check which markers are actually present in your data
  present <- genes[genes %in% rownames(combined)]
  missing <- setdiff(genes, present)
  missing_markers[[ct]] <- missing
  if (length(present) == 0) next
  
  # For each marker, make a FeaturePlot
  plots <- lapply(present, function(gene) {
    FeaturePlot(
      combined,
      features = gene,
      reduction = "umap",
      pt.size = 0.5
    ) + 
      ggtitle(paste(ct, "–", gene)) +
      theme(legend.position = "right")
  })
  
  # Combine into one panel and save
  combined_plot <- wrap_plots(plots, ncol = 2)
  ggsave(
    filename = file.path(outdir, paste0(ct, "_markers_umap.png")),
    plot = combined_plot,
    width = 8, height = 4 * ceiling(length(plots)/2)
  )
}

# Inspecting missing markers
missing_markers_df <- data.frame(
  cell_type = rep(names(missing_markers), lengths(missing_markers)),
  gene      = unlist(missing_markers),
  row.names = NULL
)
print(missing_markers_df)

# Search Ignoring Case
grep("cdh1", rownames(combined), value = TRUE, ignore.case = TRUE)

# Optionally save
write.csv(missing_markers_df, "missing_celltype_markers.csv", row.names = FALSE)

# Combined FeaturePlot: overlay all markers for each cell type in one UMAP,
# coloring each marker’s positive cells with a distinct hue.

# Create an output directory
dir.create("UMAP_combined_features", showWarnings = FALSE)

# Make sure your UMAP has been run:
if (!"umap" %in% names(combined@reductions)) {
  combined <- RunUMAP(combined, dims = 1:8)
}

# Extract UMAP embeddings once, as a data.frame
umap_df <- Embeddings(combined, "umap") %>%
  as.data.frame() %>%
  # Explicitly rename the first two columns to UMAP_1 and UMAP_2
  setNames(c("UMAP_1", "UMAP_2")) %>%
  rownames_to_column("cell")

# Loop per cell type
for (ct in names(celltype_markers)) {
  genes <- intersect(celltype_markers[[ct]], rownames(combined))
  if (length(genes) == 0) next
  
  # Pull expression for those genes
  expr_df <- GetAssayData(combined, slot = "data")[genes, , drop = FALSE] %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("cell")
  
  # Combine with UMAP coords
  df_long <- umap_df %>%
    inner_join(expr_df, by = "cell") %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr") %>%
    filter(expr > 0)
  
  # Assign each gene a unique color
  pal <- hue_pal()(length(genes))
  names(pal) <- genes
  
  p <- ggplot(df_long, aes(x = UMAP_1, y = UMAP_2, color = gene)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = pal) +
    theme_classic() +
    ggtitle(sprintf("%s — combined markers", ct)) +
    theme(legend.title = element_blank())
  
  ggsave(
    file.path("UMAP_combined_features", paste0(ct, "_combined_featureplot.png")),
    p, width = 6, height = 5
  )
}

# Individual Violin plots per gene per cell type (no points, just violins + median)
dir.create("Violin_per_marker", showWarnings = FALSE)

for (ct in names(celltype_markers)) {
  genes <- intersect(celltype_markers[[ct]], rownames(combined))
  if (length(genes) < 1) next
  
  p <- VlnPlot(
    combined,
    features = genes,
    group.by = "seurat_clusters",
    pt.size  = 0,
    combine  = TRUE
  ) +
    ggtitle(paste0(ct, " marker violins")) +
    theme_classic()
  
  ggsave(
    filename = file.path("Violin_per_marker", paste0(ct, "_violin_markers.png")),
    plot = p,
    width = 8, height = 4
  )
}

# Violin of summed expression per cell type: for each cell compute the sum of its markers,
#    then violin by cluster.

dir.create("Violin_sum_markers", showWarnings = FALSE)

# Add summed expression metadata for each cell type
for (ct in names(celltype_markers)) {
  genes <- intersect(celltype_markers[[ct]], rownames(combined))
  if (length(genes) < 1) next
  # Sum the (normalized) data slot
  combined[[paste0(ct, "_sum_expr")]] <- Matrix::colSums(
    GetAssayData(combined, slot="data")[genes, , drop=FALSE]
  )
}

# Now plot one violin per cell type sum
for (ct in names(celltype_markers)) {
  sum_col <- paste0(ct, "_sum_expr")
  if (!(sum_col %in% colnames(combined@meta.data))) next
  
  p <- VlnPlot(
    combined,
    features = sum_col,
    group.by = "seurat_clusters",
    pt.size  = 0
  ) +
    ggtitle(paste0(ct, " summed expression")) +
    theme_classic()
  
  ggsave(
    filename = file.path("Violin_sum_markers", paste0(ct, "_violin_sum.png")),
    plot = p,
    width = 8, height = 4
  )
}

# Plot UMAP, coloring by cluster
DimPlot(combined, reduction = "umap", label = TRUE)

# Create a mapping from old cluster IDs to new cell‐type names
new_labels <- c(
  `0` = "Fibroblasts",
  `1` = "Fibroblasts",
  `2` = "Fibroblasts",
  `3` = "Fibroblasts",
  `4` = "Fibroblasts",
  `5` = "Fibroblasts",
  `6` = "Fibroblasts",
  `7` = "Mast cells",
  `8` = "SMC",
  `9` = "Lymphatic endothelial",
  `10` = "Vascular endothelial",
  `11` = "SMC",
  `12` = "Pericytes"
)

# Apply the mapping to the Seurat object
Idents(combined) <- combined$seurat_clusters
combined <- RenameIdents(combined, new_labels)

# store it in metadata
combined$cell_type <- Idents(combined)

# Plot UMAP, grouping by the new cell_type names
cell_type_plot <- DimPlot(
  combined,
  reduction = "umap",
  group.by  = "cell_type",
  label     = TRUE,
  label.size= 5
) + NoLegend()

ggsave("UMAP_cell_types.png", plot = cell_type_plot, width = 6, height = 5)

saveRDS(
  combined,
  file = "combined_mouse_colon_final.rds"
)

