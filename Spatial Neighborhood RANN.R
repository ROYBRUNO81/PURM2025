########################################
# Install (one time) & load packages
########################################
# CRAN
install.packages(c("RANN", "uwot", "ggplot2", "dplyr", "tidyr", "Matrix", "pals"))

# Bioconductor (for reading .h5ad)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("zellkonverter", "SingleCellExperiment"), force = TRUE)

library(zellkonverter)
library(SingleCellExperiment)
library(RANN)
library(uwot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(pals)

set.seed(42)

########################################
# Load your AnnData (.h5ad) -> SCE
########################################
# Make sure this path points to your file:
h5ad_path <- "adata_hc.h5ad"

sce <- readH5AD(h5ad_path)   # returns a SingleCellExperiment
names(colData(sce))
cd  <- as.data.frame(colData(sce))

# Required columns in colData:
# x, y, Slice_ID, Tier1, Tier3
stopifnot(all(c("x","y","Slice_ID","Tier1","Tier3") %in% names(cd)))

# Clean labels (avoid NA issues)
cd$Tier1 <- as.character(cd$Tier1); cd$Tier1[is.na(cd$Tier1)] <- "Unknown"; cd$Tier1 <- factor(cd$Tier1)
cd$Tier3 <- as.character(cd$Tier3); cd$Tier3[is.na(cd$Tier3)] <- "Unknown"; cd$Tier3 <- factor(cd$Tier3)

t1_levels <- levels(cd$Tier1)
t3_levels <- levels(cd$Tier3)
n_cells   <- nrow(cd)

########################################
# 2) Parameters for RANN neighborhoods
########################################
# Radius in the same units as your x/y. Tune to your tissue scale.
r_radius <- 50

# Sampling size for estimating k (OSTA trick)
sample_n <- 1000

########################################
# 3) Helper to estimate k for radius search (OSTA)
########################################
estimate_k_radius <- function(coords, r, sample_n = 1000) {
  n <- nrow(coords)
  if (n < 2) return(1)
  q  <- sample(seq_len(n), min(sample_n, n))
  # Use a generous k during test; RANN will zero out neighbors beyond r
  test <- RANN::nn2(data = coords,
                    query = coords[q, , drop = FALSE],
                    k = ceiling(n/2),
                    searchtype = "radius",
                    radius = r)
  k_hat <- max(rowSums(test$nn.idx > 0))
  # OSTA recommends using ~2x the observed max to be safe
  k     <- max(10, min(n - 1, ceiling(2 * k_hat)))
  return(k)
}

########################################
# 4) Build composition vectors (Tier1 + Tier3)
########################################
# Initialize output matrix: rows=cells, cols=Tier1+Tier3
comp <- matrix(0, nrow = n_cells, ncol = length(t1_levels) + length(t3_levels))
rownames(comp) <- colnames(sce)
colnames(comp) <- c(paste0("T1_", t1_levels), paste0("T3_", t3_levels))

count_levels_vec <- function(idx, lev, nb) {
  # idx are integer codes (1..nb), tabulate to length nb
  out <- tabulate(idx, nbins = nb)
  as.numeric(out)
}

# Do it per slice (neighbors should be within a physical section)
slice_ids <- unique(cd$Slice_ID)
for (sid in slice_ids) {
  idx_slice <- which(cd$Slice_ID == sid)
  coords    <- as.matrix(cd[idx_slice, c("x","y")])
  if (nrow(coords) < 2) next
  
  # 4a) Estimate k via OSTA’s approach, then do radius search
  k_eff <- estimate_k_radius(coords, r = r_radius, sample_n = sample_n)
  
  nn <- RANN::nn2(
    data       = coords,
    query      = coords,
    k          = k_eff,
    searchtype = "radius",
    radius     = r_radius
  )
  nn_idx <- nn$nn.idx  # matrix [n_slice x k_eff], zeros where > r
  
  # Factor codes for fast tabulation
  t1_codes <- as.integer(factor(cd$Tier1[idx_slice], levels = t1_levels))  # 1..|T1|
  t3_codes <- as.integer(factor(cd$Tier3[idx_slice], levels = t3_levels))  # 1..|T3|
  
  # 4b) For each cell, tabulate neighbor Tier1/Tier3 within radius
  mat_slice <- vapply(seq_len(nrow(nn_idx)), function(i) {
    nbr <- nn_idx[i, ]
    nbr <- nbr[nbr > 0 & nbr != i]   # drop zeros (outside r) and self
    if (length(nbr) == 0L) {
      return(numeric(length(t1_levels) + length(t3_levels)))
    } else {
      t1_ct <- count_levels_vec(t1_codes[nbr], t1_levels, length(t1_levels))
      t3_ct <- count_levels_vec(t3_codes[nbr], t3_levels, length(t3_levels))
      vec   <- c(t1_ct, t3_ct)
      # convert to fractions (recommended)
      vec / length(nbr)
    }
  }, numeric(length(t1_levels) + length(t3_levels)))
  
  comp[idx_slice, ] <- t(mat_slice)
}

########################################
# 5) UMAP of neighborhood compositions + k-means contexts
########################################
# Scale, PCA (optional), then UMAP on compositions
comp_scaled <- scale(comp)
pca         <- prcomp(comp_scaled, center = FALSE, scale. = FALSE)
npc         <- min(30, ncol(pca$x))

um <- uwot::umap(
  X           = pca$x[, seq_len(npc), drop = FALSE],
  n_neighbors = 15,
  min_dist    = 0.3,
  verbose     = TRUE
)

# Define neighborhood labels via k-means on the same space (tune 'centers')
km  <- kmeans(pca$x[, seq_len(npc), drop = FALSE], centers = 25, nstart = 20)
ctx <- factor(km$cluster)

# Attach to colData for convenience
cd$RANN_ctx <- ctx

########################################
# 6) Plots
########################################
# 6a) UMAP of neighborhoods
umap_df <- data.frame(UMAP1 = um[,1], UMAP2 = um[,2], RANN_ctx = ctx)

ggplot(umap_df, aes(UMAP1, UMAP2, color = RANN_ctx)) +
  geom_point(size = 0.25, alpha = 0.8) +
  scale_color_manual(values = pals::glasbey(nlevels(umap_df$RANN_ctx))) +
  theme_void() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  ggtitle("Neighborhood UMAP (RANN radius; Tier1 + Tier3)")

# 6b) Spatial plot for one slice (pick an ID present in your data)
slice_example <- slice_ids[1]
df_slice <- cd %>%
  mutate(cell = colnames(sce)) %>%
  filter(Slice_ID == slice_example)

ggplot(df_slice, aes(x = x, y = y, color = RANN_ctx)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = pals::glasbey(nlevels(df_slice$RANN_ctx))) +
  coord_equal(expand = FALSE) +
  theme_void() +
  # flip y if your coordinates grow downward; otherwise remove:
  scale_y_reverse() +
  ggtitle(paste("Spatial neighborhoods — slice", slice_example))

# Saving Object

# Store UMAP and composition in reducedDims -> will appear in AnnData .obsm
reducedDims(sce)$X_umap_rann <- um
reducedDims(sce)$X_rann_comp <- comp   # optional but handy to keep

# Store labels in colData -> will appear in AnnData .obs
colData(sce)$RANN_ctx <- as.factor(ctx)

# (Optional) keep your parameters in metadata -> will appear in AnnData .uns
metadata(sce)$rann_params <- list(radius = r_radius, sampled_per_slice = sample_n)

# Write a new H5AD without touching your original file
writeH5AD(sce, "adata_hc_rann.h5ad")

