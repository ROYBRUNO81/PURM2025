# Load required libraries
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(scDblFinder)
library(tidyr)
library(tibble)    # for rownames_to_column()
library(scales)
library(ggnewscale)
library(DESeq2)

# Read in the count matrices
hc_counts <- read.table(
  "GSE114374_Mouse_HC_expression_matrix.txt.gz",
  header      = TRUE, row.names = 1,
  sep         = "\t", comment.char = "",
  check.names = FALSE
)
dss_counts <- read.table(
  "GSE114374_Mouse_DSS_expression_matrix.txt.gz",
  header      = TRUE, row.names = 1,
  sep         = "\t", comment.char = "",
  check.names = FALSE
)

# Create Seurat objects
hc  <- CreateSeuratObject(counts = hc_counts,  project = "HC",  min.cells = 3, min.features = 200)
dss <- CreateSeuratObject(counts = dss_counts, project = "DSS", min.cells = 3, min.features = 200)

################################################################################
#                              PROCESS HC SAMPLE                              #
################################################################################

## QC & Filtering for HC
hc[["percent.mt"]] <- PercentageFeatureSet(hc, pattern = "^mt-")

# Violin plot
VlnPlot(hc,
        features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol     = 3) +
  ggtitle("HC: QC metrics")

# Scatter plots
FeatureScatter(hc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("HC: nCount_RNA vs percent.mt")
FeatureScatter(hc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("HC: nCount_RNA vs nFeature_RNA")

# Filter cells
hc <- subset(hc, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)

## Remove zero‐only genes
hc[["RNA"]] <- JoinLayers(hc[["RNA"]])
hc_counts_mat <- hc@assays$RNA@layers[["counts"]]
hc_keep_genes <- Matrix::rowSums(hc_counts_mat > 0) > 0
hc <- subset(hc, features = rownames(hc_counts_mat)[hc_keep_genes])

## Pre‐clustering (for doublet detection)
hc <- NormalizeData(hc)
hc <- FindVariableFeatures(hc)
hc <- ScaleData(hc)
hc <- RunPCA(hc)
ElbowPlot(hc) + ggtitle("HC: PCA Elbow")

hc <- FindNeighbors(hc, dims = 1:9)
hc <- FindClusters(hc, resolution = 0.5)

## Doublet detection & filtering
hc_sce <- scDblFinder(as.SingleCellExperiment(hc), clusters = hc$seurat_clusters)
hc$doublet_class <- hc_sce$scDblFinder.class

DimPlot(hc, group.by = "doublet_class") +
  ggtitle("HC: Doublet calls")

hc <- subset(hc, subset = doublet_class == "singlet")

## Full normalization → PCA → UMAP
hc <- NormalizeData(hc, normalization.method = "LogNormalize", scale.factor = 1e4)
hc <- FindVariableFeatures(hc, selection.method = "vst", nfeatures = 2000)

# Variable feature plot
top10_hc <- head(VariableFeatures(hc), 10)
p_hc_vf <- VariableFeaturePlot(hc) + ggtitle("HC: Variable features")
LabelPoints(plot = p_hc_vf, points = top10_hc)

hc <- ScaleData(hc, features = rownames(hc))
hc <- RunPCA(hc, features = VariableFeatures(hc))

print(hc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(hc, dims = 1:2, reduction = "pca") + ggtitle("HC: PCA loadings")
DimHeatmap(hc, dims = 1, cells = 500, balanced = TRUE) + ggtitle("HC: PC1 heatmap")

ElbowPlot(hc) + ggtitle("HC: PCA Elbow (post-scaling)")

hc <- FindNeighbors(hc, dims = 1:9)
hc <- FindClusters(hc, resolution = 0.5)
hc <- RunUMAP(hc, dims = 1:8)

DimPlot(hc, reduction = "umap", label = TRUE) +
  ggtitle("HC: UMAP clusters")

## Marker detection & export
hc_markers <- FindAllMarkers(
  hc,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)
write.csv(hc_markers,
          file      = "Mouse_colon/Seperate_Analysis/HC_cluster_markers.csv",
          row.names = FALSE)

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

# Define blue→green→yellow→red palette
heatCols <- c("blue", "green", "yellow", "red")

# extract UMAP coords and make a data.frame once
hc_df <- Embeddings(hc, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell") 

# and pull out normalized expression slot
expr_mat_hc <- GetAssayData(hc, slot="data")

# Base output directory
outdir2 <- "Mouse_colon/Seperate_Analysis/UMAP_feature_plots_hc"
if (!dir.exists(outdir2)) dir.create(outdir2)

# Loop over each cell‐type and its markers
for (ct in names(celltype_markers)) {
  # restrict to genes actually present in HC
  genes <- intersect(celltype_markers[[ct]], rownames(expr_mat_hc))
  if (length(genes) == 0) next
  
  for (gene in genes) {
    # build a df with UMAP coords + this gene's expression
    df <- hc_df %>%
      mutate(expr = expr_mat_hc[gene, cell])
    
    # feature plot
    p <- ggplot(df, aes(x = umap_1, y = umap_2, color = expr)) +
      geom_point(size = 0.5) +
      scale_color_gradientn(
        colours = heatCols,
        name    = gene,
        limits  = range(df$expr, na.rm = TRUE)
      ) +
      labs(
        title = paste0("HC: ", ct, " — ", gene),
        x     = "UMAP 1",
        y     = "UMAP 2"
      ) +
      theme_classic() +
      theme(
        plot.title      = element_text(hjust = 0.5),
        legend.position = "right"
      )
    
    # save
    ggsave(
      filename = file.path(outdir2, paste0(ct, "_", gene, "_feature_plot.png")),
      plot     = p,
      width    = 6,
      height   = 5,
      dpi      = 300
    )
  }
}

# Create a mapping from old cluster IDs to new cell‐type names
new_labels_hc <- c(
  `0` = "Fibroblasts",
  `1` = "Fibroblasts",
  `2` = "Fibroblasts",
  `3` = "Fibroblasts",
  `4` = "Fibroblasts",
  `5` = "Fibroblasts",
  `6` = "SMC",
  `7` = "Fibroblasts",
  `8` = "Vascular endothelial",
  `9` = "Lymphatic endothelial",
  `10` = "Pericytes",
  `11` = "Fibroblasts"
)

# Apply the mapping to the Seurat object
Idents(hc) <- hc$seurat_clusters
hc <- RenameIdents(hc, new_labels_hc)

# store it in metadata
hc$cell_type <- Idents(hc)

# Plot UMAP, grouping by the new cell_type names
cell_type_plot <- DimPlot(
  hc,
  reduction = "umap",
  group.by  = "cell_type",
  label     = TRUE,
  label.size= 5
) + NoLegend()

ggsave("UMAP_cell_types_hc.png", plot = cell_type_plot, width = 6, height = 5)

saveRDS(hc, file = "Seperate_Analysis/HC_processed.rds")

################################################################################
#                             PROCESS DSS SAMPLE                             #
################################################################################

## QC & Filtering for DSS
dss[["percent.mt"]] <- PercentageFeatureSet(dss, pattern = "^mt-")

# Violin plot
VlnPlot(dss,
        features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol     = 3) +
  ggtitle("DSS: QC metrics")

# Scatter plots
FeatureScatter(dss, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("DSS: nCount_RNA vs percent.mt")
FeatureScatter(dss, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("DSS: nCount_RNA vs nFeature_RNA")

# Filter cells
dss <- subset(dss, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)

## Remove zero‐only genes
dss[["RNA"]] <- JoinLayers(dss[["RNA"]])
dss_counts_mat <- dss@assays$RNA@layers[["counts"]]
dss_keep_genes <- Matrix::rowSums(dss_counts_mat > 0) > 0
dss <- subset(dss, features = rownames(dss_counts_mat)[dss_keep_genes])

## Pre‐clustering (for doublet detection)
dss <- NormalizeData(dss)
dss <- FindVariableFeatures(dss)
dss <- ScaleData(dss)
dss <- RunPCA(dss)
ElbowPlot(dss) + ggtitle("DSS: PCA Elbow")

dss <- FindNeighbors(dss, dims = 1:8)
dss <- FindClusters(dss, resolution = 0.5)

## Doublet detection & filtering
dss_sce <- scDblFinder(as.SingleCellExperiment(dss), clusters = dss$seurat_clusters)
dss$doublet_class <- dss_sce$scDblFinder.class

DimPlot(dss, group.by = "doublet_class") +
  ggtitle("DSS: Doublet calls")

dss <- subset(dss, subset = doublet_class == "singlet")

## Full normalization → PCA → UMAP
dss <- NormalizeData(dss, normalization.method = "LogNormalize", scale.factor = 1e4)
dss <- FindVariableFeatures(dss, selection.method = "vst", nfeatures = 2000)

# Variable feature plot
top10_dss <- head(VariableFeatures(dss), 10)
p_dss_vf <- VariableFeaturePlot(dss) + ggtitle("DSS: Variable features")
LabelPoints(plot = p_dss_vf, points = top10_dss)

dss <- ScaleData(dss, features = rownames(dss))
dss <- RunPCA(dss, features = VariableFeatures(dss))

print(dss[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(dss, dims = 1:2, reduction = "pca") + ggtitle("DSS: PCA loadings")
DimHeatmap(dss, dims = 1, cells = 500, balanced = TRUE) + ggtitle("DSS: PC1 heatmap")

ElbowPlot(dss) + ggtitle("DSS: PCA Elbow (post-scaling)")

dss <- FindNeighbors(dss, dims = 1:8)
dss <- FindClusters(dss, resolution = 0.5)
dss <- RunUMAP(dss, dims = 1:8)

DimPlot(dss, reduction = "umap", label = TRUE) +
  ggtitle("DSS: UMAP clusters")

## Marker detection & export
dss_markers <- FindAllMarkers(
  dss,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)
write.csv(dss_markers,
          file      = "Mouse_colon/Seperate_Analysis/DSS_cluster_markers.csv",
          row.names = FALSE)

# extract UMAP coords and make a data.frame once
dss_df <- Embeddings(dss, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell") 

# and pull out normalized expression slot
expr_mat_dss <- GetAssayData(dss, slot="data")

# Base output directory
outdir2 <- "Mouse_colon/Seperate_Analysis/UMAP_feature_plots_dss"
if (!dir.exists(outdir2)) dir.create(outdir2)

# Loop over each cell‐type and its markers
for (ct in names(celltype_markers)) {
  # restrict to genes actually present in HC
  genes <- intersect(celltype_markers[[ct]], rownames(expr_mat_dss))
  if (length(genes) == 0) next
  
  for (gene in genes) {
    # build a df with UMAP coords + this gene's expression
    df <- dss_df %>%
      mutate(expr = expr_mat_dss[gene, cell])
    
    # feature plot
    p <- ggplot(df, aes(x = umap_1, y = umap_2, color = expr)) +
      geom_point(size = 0.5) +
      scale_color_gradientn(
        colours = heatCols,
        name    = gene,
        limits  = range(df$expr, na.rm = TRUE)
      ) +
      labs(
        title = paste0("DSS: ", ct, " — ", gene),
        x     = "UMAP 1",
        y     = "UMAP 2"
      ) +
      theme_classic() +
      theme(
        plot.title      = element_text(hjust = 0.5),
        legend.position = "right"
      )
    
    # save
    ggsave(
      filename = file.path(outdir2, paste0(ct, "_", gene, "_feature_plot.png")),
      plot     = p,
      width    = 6,
      height   = 5,
      dpi      = 300
    )
  }
}

# Plot UMAP, coloring by cluster
DimPlot(combined, reduction = "umap", label = TRUE)

# Reclustering cluster 6
# 1. Subset to only cluster 6 cells
dss_sub <- subset(dss, idents = 6)

# 2. Pre-process the subset exactly as before
#    (you can tweak parameters like variable feature count or PCA dims if needed)
dss_sub <- NormalizeData(dss_sub, normalization.method = "LogNormalize", scale.factor = 1e4)
dss_sub <- FindVariableFeatures(dss_sub, selection.method = "vst", nfeatures = 2000)
dss_sub <- ScaleData(dss_sub, features = rownames(dss_sub))
dss_sub <- RunPCA(dss_sub, features = VariableFeatures(dss_sub))

# 3. Choose PCs for subclustering (inspect ElbowPlot if you want)
ElbowPlot(dss_sub) + ggtitle("DSS cluster 6: PCA elbow")

# 4. Build neighborhood graph and find subclusters
dss_sub <- FindNeighbors(dss_sub, dims = 1:9)
dss_sub <- FindClusters(dss_sub, resolution = 0.5)  
# ↑ you may need a higher res (e.g. 0.8–1.2) to separate SMC vs pericytes

# 5. Run UMAP for visualization
dss_sub <- RunUMAP(dss_sub, dims = 1:9)

# 6. Plot the subclusters
DimPlot(dss_sub, reduction = "umap", label = TRUE) +
  ggtitle("DSS cluster 6: subclustering")

# 7. Check marker genes to annotate the new subclusters
sub_markers <- FindAllMarkers(dss_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(sub_markers)

# extract UMAP coords and make a data.frame once
dss_df_sub <- Embeddings(dss_sub, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell") 

# and pull out normalized expression slot
expr_mat_dss <- GetAssayData(dss_sub, slot="data")

# Base output directory
outdir2 <- "Seperate_Analysis/UMAP_feature_plots_dss_Cl6"
if (!dir.exists(outdir2)) dir.create(outdir2)

# Loop over each cell‐type and its markers
for (ct in names(celltype_markers)) {
  # restrict to genes actually present in HC
  genes <- intersect(celltype_markers[[ct]], rownames(expr_mat_dss))
  if (length(genes) == 0) next
  
  for (gene in genes) {
    # build a df with UMAP coords + this gene's expression
    df <- dss_df_sub %>%
      mutate(expr = expr_mat_dss[gene, cell])
    
    # feature plot
    p <- ggplot(df, aes(x = umap_1, y = umap_2, color = expr)) +
      geom_point(size = 0.5) +
      scale_color_gradientn(
        colours = heatCols,
        name    = gene,
        limits  = range(df$expr, na.rm = TRUE)
      ) +
      labs(
        title = paste0("DSS: ", ct, " — ", gene),
        x     = "UMAP 1",
        y     = "UMAP 2"
      ) +
      theme_classic() +
      theme(
        plot.title      = element_text(hjust = 0.5),
        legend.position = "right"
      )
    
    # save
    ggsave(
      filename = file.path(outdir2, paste0(ct, "_", gene, "_feature_plot.png")),
      plot     = p,
      width    = 6,
      height   = 5,
      dpi      = 300
    )
  }
}

new_labels_dss_Cl6 <- c(
  `0` = "SMC",
  `1` = "SMC",
  `2` = "Pericytes",
  `3` = "Fibroblasts"
)


# Create a mapping from old cluster IDs to new cell‐type names
new_labels_dss <- c(
  `0` = "Fibroblasts",
  `1` = "Fibroblasts",
  `2` = "Fibroblasts",
  `3` = "Fibroblasts",
  `4` = "Fibroblasts",
  `5` = "Fibroblasts",
  `6` = "NA",
  `7` = "Fibroblasts",
  `8` = "Lymphatic endothelial",
  `9` = "Vascular endothelial"
)

# Make sure Idents are set to the numeric clusters
Idents(dss) <- dss$seurat_clusters

# Rename the main DSS clusters using your combined-derived mapping
dss <- RenameIdents(dss, new_labels_dss)

# 3. Store those as metadata
dss$cell_type <- Idents(dss)

# 4. Now handle the refined sub-clustering of cluster 6:
#    - First, get the cells in cluster “6” before relabeling
cells6 <- rownames(dss@meta.data)[ dss@meta.data$seurat_clusters == 6 ]

#    - Subset out the original cluster-6 cells into dss_sub (you already have this)
#      and relabel them with new_labels_dss_Cl6
Idents(dss_sub) <- dss_sub$seurat_clusters
dss_sub <- RenameIdents(dss_sub, new_labels_dss_Cl6)
dss_sub$subcluster6 <- Idents(dss_sub)

#    - Now inject those sublabels back into the full dss metadata
dss$subcluster6 <- NA
dss$subcluster6[cells6] <- as.character(dss_sub$subcluster6)

# 5. Plot and save the UMAP colored by the “global” cell types
p1 <- DimPlot(
  dss,
  reduction  = "umap",
  group.by   = "cell_type",
  label      = TRUE,
  label.size = 5
) +
  ggtitle("DSS: Cell types (combined labels)") +
  NoLegend()

ggsave(
  filename = "Seperate_Analysis/DSS_UMAP_combined_cell_types.png",
  plot     = p1,
  width    = 6,
  height   = 5,
  dpi      = 300
)

# 6. Plot and save the UMAP colored by the refined cluster-6 sublabels
p2 <- DimPlot(
  dss,
  reduction  = "umap",
  group.by   = "subcluster6",
  label      = TRUE,
  label.size = 5
) +
  ggtitle("DSS: Subclusters within original cluster 6") +
  NoLegend()

ggsave(
  filename = "Seperate_Analysis/DSS_UMAP_cluster6_sublabels.png",
  plot     = p2,
  width    = 6,
  height   = 5,
  dpi      = 300
)

saveRDS(dss, file = "Seperate_Analysis/DSS_processed.rds")

################################################################################
#########  Compare Combined and Seperate Analysis Cluster Annotations  #########
################################################################################

library(mclust)
combined <- readRDS("combined_mouse_colon_final.rds")

# 1) Transfer HC labels onto combined (HC as reference)
anchors_hc <- FindTransferAnchors(
  reference          = hc,
  query              = combined,
  dims               = 1:30
)
pred_hc <- TransferData(
  anchorset = anchors_hc,
  refdata   = hc$cell_type,   # your HC labels
  dims      = 1:30
)
# pred_hc is a DataFrame with column "predicted.id"

# 2) Transfer DSS labels onto combined (DSS as reference)
anchors_dss <- FindTransferAnchors(
  reference          = dss,
  query              = combined,
  dims               = 1:30
)
pred_dss <- TransferData(
  anchorset = anchors_dss,
  refdata   = dss$cell_type,  # your DSS labels
  dims      = 1:30
)
# pred_dss has its own "predicted.id"

# 3) Build confusion tables *without* touching combined@meta.data
orig_clusters <- combined$seurat_clusters

hc_labels  <- pred_hc$predicted.id
dss_labels <- pred_dss$predicted.id

# 1) Subset to HC cells only
hc_only <- subset(combined, subset = condition == "HC")

# 2) Pull out the original clusters and the HC-derived predictions
orig_clusters_hc <- hc_only$seurat_clusters
hc_labels_hc     <- pred_hc[colnames(hc_only), "predicted.id"]

# 3) Table just HC cells
table(Combined = orig_clusters_hc,
      HC       = hc_labels_hc)

# Subset to DSS cells
dss_only <- subset(combined, subset = condition == "DSS")

# Get vectors for DSS
orig_clusters_dss <- dss_only$seurat_clusters
dss_labels_dss     <- pred_dss[colnames(dss_only), "predicted.id"]

# Table just DSS cells
table(Combined = orig_clusters_dss,
      DSS      = dss_labels_dss)

# 4) Compute ARI to summarize agreement
ari_hc  <- adjustedRandIndex(orig_clusters_hc, hc_labels_hc)
ari_dss <- adjustedRandIndex(orig_clusters_dss, dss_labels_dss)

cat(sprintf("\nARI (combined ↔ HC):  %.3f\n", ari_hc))
cat(sprintf("ARI (combined ↔ DSS): %.3f\n", ari_dss))

# Assuming hc_table and dss_table hold the two confusion matrices:
hc_table  <- table(Combined = orig_clusters_hc, HC  = hc_labels_hc)
dss_table <- table(Combined = orig_clusters_dss, DSS = dss_labels_dss)

# Convert to data.frame for saving
hc_df  <- as.data.frame(hc_table)
dss_df <- as.data.frame(dss_table)

# Create output directory if needed
outdir <- "Separate_Analysis/Confusion_Tables"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Write CSV files
write.csv(hc_df,
          file      = file.path(outdir, "HC_confusion_table.csv"),
          row.names = FALSE)

write.csv(dss_df,
          file      = file.path(outdir, "DSS_confusion_table.csv"),
          row.names = FALSE)




