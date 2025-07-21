# Load required library
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(scDblFinder)
library(tidyr)

# Read in the count matrices
hc_counts <- read.table(
  "GSE114374_Mouse_HC_expression_matrix.txt.gz",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  comment.char = "",
  check.names = FALSE
)

dss_counts <- read.table(
  "GSE114374_Mouse_DSS_expression_matrix.txt.gz",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  comment.char = "",
  check.names = FALSE
)

# Create Seurat objects, adding a 'condition' field
hc <- CreateSeuratObject(
  counts = hc_counts,
  project = "Colon",
  min.cells = 3,
  min.features = 200
)
hc

dss <- CreateSeuratObject(
  counts = dss_counts,
  project = "Colon",
  min.cells = 3,
  min.features = 200
)
dss

# Merge the two objects into one, here prefixing cell barcodes to keep them unique
combined <- merge(
  x = hc,
  y = dss,
  add.cell.ids = c("HC", "DSS"),
  project = "Mouse_Colon"
)

# Check that metadata carried over
table(combined$orig.ident)

# QC Metrics & Filtering

# Calculate percent.mt
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")

# Violin plots of QC metrics
p_violin <- VlnPlot(
  combined,
  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
  ncol = 3
)
p_violin

# Scatter plots: nCount_RNA vs percent.mt, and vs nFeature_RNA
plot1 <- FeatureScatter(combined, feature1="nCount_RNA", feature2="percent.mt")
plot1

plot2 <- FeatureScatter(combined, feature1="nCount_RNA", feature2="nFeature_RNA")
plot2

# Filtering
combined <- subset(
  combined,
  subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5
)

# Collapse all split counts and data layers into one of each
combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

# Extract the raw counts from the “counts” layer
counts_mat <- combined@assays$RNA@layers[["counts"]]

# Identify genes detected in at least one cell
keep_genes <- Matrix::rowSums(counts_mat > 0) > 0

# Subset the object to drop those all‐zero genes
combined <- subset(
  combined,
  features = rownames(counts_mat)[keep_genes]
)

# Preliminary Clustering (for scDblFinder)
combined <- NormalizeData(combined)             
combined <- FindVariableFeatures(combined) 
combined <- ScaleData(combined)   
combined <- RunPCA(combined)
ElbowPlot(combined)
combined <- FindNeighbors(combined, dims = 1:8) 
combined <- FindClusters(combined, resolution = 0.5)

# Doublet Detection
sce <- scDblFinder(
  as.SingleCellExperiment(combined),
  clusters = combined$seurat_clusters
)
combined$doublet_class <- sce$scDblFinder.class

# Visualize doublet calls
DimPlot(combined, group.by = "doublet_class")

# Keep only singlets
combined <- subset(combined, subset = doublet_class == "singlet")
combined

# Normalization & Variable Feature Selection

# Normalize
combined <- NormalizeData(combined, normalization.method="LogNormalize", scale.factor=10000)

# Variable features
combined <- FindVariableFeatures(combined, selection.method="vst", nfeatures=2000)

# Plot top 10 variable features
top10 <- head(VariableFeatures(combined), 10)
p_vf <- VariableFeaturePlot(combined)
p_lab <- LabelPoints(plot=p_vf, points=top10, repel=TRUE)
p_vf + p_lab

# Scale all genes so that mean = 0 and variance = 1 (for PCA)
all_genes <- rownames(combined)
DefaultAssay(combined) <- "RNA"
combined <- ScaleData(
  object   = combined,
  features = all_genes
)

# Run PCA on the scaled data
combined <- RunPCA(
  combined,
  features = VariableFeatures(combined)
)

# Examine and visualize PCA results a few different ways
print(combined[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(combined, dims = 1:2, reduction = "pca")

DimHeatmap(combined, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(combined, dims = 1:15, cells = 500, balanced = TRUE)

# Visualize variance explained per PC and choose cutoff
ElbowPlot(combined)

# Build shared nearest neighbor graph based on chosen PCs
combined <- FindNeighbors(combined, dims = 1:8)

# Identify clusters with resolution = 0.5
combined <- FindClusters(combined, resolution = 0.5)

# Run UMAP for 2D embedding
combined <- RunUMAP(combined, dims = 1:8)

# Plot UMAP, coloring by cluster
DimPlot(combined, reduction = "umap", label = TRUE)

################################################################################
#################### Marker Identification & CSV Export ########################
################################################################################

# 1. Unbiased cluster marker detection
# Finds genes that are upregulated in each cluster compared to all other cells
# only.pos = TRUE    : only returns genes positively enriched in each cluster
# min.pct = 0.25     : gene must be detected in at least 25% of cells in the cluster
# logfc.threshold = 0.25 : at least 0.25 log2-fold change

markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(
  markers,
  file      = "Mouse_Colon_cluster_markers.csv",
  row.names = FALSE
)

head(markers)

# Saving seurat object
saveRDS(
  combined,
  file = "combined_mouse_colon_final.rds"
)
