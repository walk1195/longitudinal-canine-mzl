#############################################################################################
############################################################################################
### ORv-01 R Spatial Analysis: ISCHIA
### Samples: ORv-01 Visium v1 Spatial Sample (D10)
### Author: Grace Walker
### Date: Feb 20, 2026
#############################################################################################
#############################################################################################

# Set up environment
# -------------------------------------------------------------------------------------------

# Load libs
library(Seurat)
library(ISCHIA)
library(tidyverse)
library(patchwork) # For plot wrapping
library(matrixStats) # Because we had to edit a function

# Dirs
objDir <- ""
resDir <- ""

####################################################################################
# Standard clustering
####################################################################################

# Read in raw seurat obj
spatial <- readRDS(glue::glue("{objDir}/orv01_spatial_clustered.rds"))

# Set Idents
spatial@meta.data['sample'] <- "ORv-01"
Idents(spatial) <- "sample"

# Visualize data
## Total counts
p1 <- VlnPlot(spatial, features = "total_counts", pt.size = 0.1) + NoLegend()
p2 <- SpatialFeaturePlot(spatial, features = "total_counts") + theme(legend.position = "right")
wrap_plots(p1, p2)
ggsave(glue::glue("{resDir}qc_total_counts.png"), dpi=400, height=5, width=10)

## Unique genes
p1 <- VlnPlot(spatial, features = "n_genes_by_counts", layer='counts',pt.size = 0.1) + NoLegend()
p2 <- SpatialFeaturePlot(spatial, features = "n_genes_by_counts") + theme(legend.position = "right")
wrap_plots(p1, p2)
ggsave(glue::glue("{resDir}qc_n_genes_by_counts.png"), dpi=400, height=5, width=10)

# Standard normalization
spatial <- NormalizeData(spatial, assay = "Spatial")
spatial <- ScaleData(spatial, assay = "Spatial")

# Check normalized expression
SpatialFeaturePlot(spatial, features = c("CD79B", "GSN"), pt.size.factor = 1.3)

# Dim reduction and clustering
spatial <- spatial %>%
  FindVariableFeatures() %>%
  RunPCA() %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(reduction = "pca", dims = 1:30)

# Visualize
p1 <- DimPlot(spatial, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(spatial, label = TRUE, label.size = 3)
p1 + p2
ggsave(glue::glue("{resDir}standard_seurat_clustering_res0.8.png"), dpi=400, height=5, width=10)


# Separate by group
SpatialDimPlot(spatial,
               cells.highlight = CellsByIdentities(object = spatial, idents = c(0, 6, 4, 1, 5, 3)),
               facet.highlight = TRUE,
               ncol = 3)
ggsave(glue::glue("{resDir}seurat_individual_clusters_0.8.png"), dpi=400, height=5, width=10)

# Interactive plotting if desired
# SpatialDimPlot(spatial, interactive = TRUE)
# SpatialFeaturePlot(spatial, features = "CD79A", interactive = TRUE)

# Plotting individual cell types
# Vector of cell types
cell_predictions <- colnames(spatial@meta.data)[grep('cells', colnames(spatial@meta.data))]

# Loop over cell types
for (cell_type in cell_predictions) {
  SpatialFeaturePlot(spatial, features = cell_type, pt.size.factor = 1.3, crop = TRUE)
  # Fix backslash bug
  if (grepl("/", cell_type)) {
    cell_type <- gsub("/", "_", cell_type)
  }
  # Save image
  ggsave(glue::glue("{resDir}cell_predictions_{cell_type}.png"), dpi=400, height=5, width=7)
}

# Save object
saveRDS(spatial, glue::glue("{objDir}/orv01_spatial_clustered.rds"))

####################################################################################
# ISCHIA
####################################################################################

# Read seurat obj
spatial <- readRDS(glue::glue("{objDir}orv01_spatial_clustered.rds"))

# Extract deconv cell prediction matrix
deconv.mat <- as.matrix(spatial@meta.data[,cell_predictions])

# Spot-based normalization
deconv.norm <- deconv.mat / rowSums(deconv.mat)

# Prep object for plotting normalized labels
deconv.norm.labelled <- deconv.norm
colnames(deconv.norm.labelled) <- paste0(colnames(deconv.norm), "_norm") 
spatial_normalized <- spatial
spatial_normalized <- AddMetaData(spatial_normalized, metadata = deconv.norm.labelled)

# Replot
for (cell_type in cell_predictions) {
  SpatialFeaturePlot(spatial_normalized, features = paste0(cell_type, "_norm"), pt.size.factor = 1.3, crop = TRUE)
  # Fix backslash bug
  if (grepl("/", cell_type)) {
    cell_type <- gsub("/", "_", cell_type)
  }
  # Save image
  ggsave(glue::glue("{resDir}cell_predictions_{cell_type}_normalized.png"), dpi=400, height=5, width=7)
}

# Decide k value
Composition.cluster.k(deconv.norm, 20)

# Run composition clustering of the deconvoluted spatial spots
spatial <- Composition.cluster(spatial, deconv.norm, 5)

# Inspect number of clusters
table(spatial$CompositionCluster_CC)

# Plot
SpatialDimPlot(spatial, group.by = c("CompositionCluster_CC")) + 
  scale_fill_manual(values = c("dodgerblue2", "royalblue4", "orange","mediumvioletred","seagreen1")) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  labs(fill = "Composition Cluster") +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))
ggsave(glue::glue("{resDir}ISCHIA/composition_cluster_spatial_plot.png"), dpi=400, height=5, width=7)


####################################################################################
# Gotta do a quick function fix
f <- ISCHIA::Composition_cluster_enrichedCelltypes
print(f, max = 50)

Composition_cluster_enrichedCelltypes2 <- function (Spatial.object, COI, Celltype_deconvolved_probs) {
  COI.spots <- rownames(Spatial.object@meta.data[which(Spatial.object@meta.data$CompositionCluster_CC == 
                                                         COI), ])
  COI.topic.probs <- Celltype_deconvolved_probs[COI.spots, 
  ]
  COI.topic.probs.medians <- sort(colMedians(COI.topic.probs), 
                                  decreasing = T)
  COI.leading.topics <- names(COI.topic.probs.medians[1:2])
  COI.topic.probs.melt <- reshape2::melt(COI.topic.probs) ### Perform the fix here
  colnames(COI.topic.probs.melt) <- c("spot_id", "Celltype", 
                                      "prob")
  Topic.prob.all.df <- COI.topic.probs.melt
  Topic.prob.all.df$Topic <- factor(Topic.prob.all.df$Celltype)
  p <- ggplot(Topic.prob.all.df, aes(x = Topic, y = prob, fill = Topic)) + 
    geom_boxplot(show.legend = FALSE) + 
    labs(x = "Cell type", y = "Predicted Proportion (Norm)") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, face = "bold", colour = "black"), 
                                                           axis.text.y = element_text(hjust = 1, size = 17, face = "bold", 
                                                                                      colour = "black"), panel.border = element_blank(), 
                                                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                           axis.line = element_line(colour = "black")) + ggtitle(paste("Composition Cluster", 
                                                                                                                       COI, sep = " "))
  return(p)
}
####################################################################################

# Plot cell type proportions of each cluster
for (cluster in unique(spatial$CompositionCluster_CC)) {
  Composition_cluster_enrichedCelltypes2(spatial, cluster, deconv.norm)
  ggsave(glue::glue("{resDir}ISCHIA/{cluster}_proportions_boxplot.png"), dpi=400, height=5, width=7)
}

# Plot UMAP
spatial.umap <- Composition_cluster_umap(spatial, deconv.norm)
#> Plotting scatterpies for 2185 pixels with 20 cell-types...this could take a while if the dataset is large.
spatial.umap$umap.cluster.gg
ggsave(glue::glue("{resDir}ISCHIA/composition_cluster_umap.png"), dpi=400, height=5, width=7)

# Deconvoluted umap
spatial.umap$umap.deconv.gg
ggsave(glue::glue("{resDir}ISCHIA/composition_cluster_umap_deconvolved.png"), dpi=400, height=5, width=7)


# Cell co-occurrence analysis
### CC1
CC1.celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=spatial,
                                                     deconv.prob.mat=deconv.norm,
                                                     COI="CC1", prob.th= 0.05, Condition=unique(spatial$orig.ident))
plot.celltype.cooccurence(CC1.celltype.cooccur)
ggsave(glue::glue("{resDir}ISCHIA/coocurrence_CC1.png"), dpi=400, height=6, width=7)

### CC2
CC2.celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=spatial,
                                                     deconv.prob.mat=deconv.norm,
                                                     COI="CC2", prob.th= 0.05, Condition=unique(spatial$orig.ident))
plot.celltype.cooccurence(CC2.celltype.cooccur)
ggsave(glue::glue("{resDir}ISCHIA/coocurrence_CC2.png"), dpi=400, height=6, width=7)

### CC3
CC3.celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=spatial,
                                                     deconv.prob.mat=deconv.norm,
                                                     COI="CC3", prob.th= 0.05, Condition=unique(spatial$orig.ident))
plot.celltype.cooccurence(CC3.celltype.cooccur)
ggsave(glue::glue("{resDir}ISCHIA/coocurrence_CC3.png"), dpi=400, height=6, width=7)

### CC4
CC4.celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=spatial,
                                                     deconv.prob.mat=deconv.norm,
                                                     COI="CC4", prob.th= 0.05, Condition=unique(spatial$orig.ident))
plot.celltype.cooccurence(CC4.celltype.cooccur)
ggsave(glue::glue("{resDir}ISCHIA/coocurrence_CC4.png"), dpi=400, height=6, width=7)

### CC5
CC5.celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=spatial,
                                                     deconv.prob.mat=deconv.norm,
                                                     COI="CC5", prob.th= 0.05, Condition=unique(spatial$orig.ident))
plot.celltype.cooccurence(CC5.celltype.cooccur)
ggsave(glue::glue("{resDir}ISCHIA/coocurrence_CC5.png"), dpi=400, height=6, width=7)


# Horizontal bar plot of cell type proportions for CCs
cluster_assignments <- data.frame(spot = rownames(deconv.norm), cluster = spatial$CompositionCluster_CC)

# Combine cluster info with spot proportions
deconv.df <- deconv.norm %>%
  as.data.frame() %>%
  mutate(spot = rownames(deconv.norm)) %>%
  left_join(cluster_assignments, by = "spot")

# Compute average
avg.comp <- deconv.df %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric), mean))  # average per cell type

avg.comp.long <- avg.comp %>%
  pivot_longer(
    cols = -cluster,
    names_to = "cell_type",
    values_to = "avg_proportion"
  )

ggplot(avg.comp.long, aes(x = cluster, y = avg_proportion, fill = cell_type)) +
  geom_col(width = 0.6,color = "black", size = 0.2) +
  theme_classic(base_size = 14) +
  coord_flip() +              # make it horizontal
  labs(x = "Composition Cluster", y = "Avg Cell Proportion", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(glue::glue("{resDir}ISCHIA/composition_cluster_stacked_barchart.png"), dpi=400, height=5, width=7)

# -------------------------------------------------------------------------------------------
# Rerunning for k = 4 (4 clusters)
# -------------------------------------------------------------------------------------------

# Plot spatially
SpatialDimPlot(spatial, group.by = c("CompositionCluster_CC")) + 
  scale_fill_manual(values = c("dodgerblue2", "orange","mediumvioletred","seagreen1")) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  labs(fill = "Composition Cluster") +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))
ggsave(glue::glue("{resDir}ISCHIA/composition_cluster_spatial_plot_4.png"), dpi=400, height=5, width=7)


# Plot cell type proportions of each cluster
for (cluster in unique(spatial$CompositionCluster_CC)) {
  Composition_cluster_enrichedCelltypes2(spatial, cluster, deconv.norm)
  ggsave(glue::glue("{resDir}ISCHIA/{cluster}_proportions_boxplot_4.png"), dpi=400, height=5, width=7)
}

# Plot UMAP
spatial.umap <- Composition_cluster_umap(spatial, deconv.norm)
#> Plotting scatterpies for 2185 pixels with 20 cell-types...this could take a while if the dataset is large.
spatial.umap$umap.cluster.gg
ggsave(glue::glue("{resDir}ISCHIA/composition_cluster_umap_4.png"), dpi=400, height=5, width=7)

# Deconvoluted umap
spatial.umap$umap.deconv.gg
ggsave(glue::glue("{resDir}ISCHIA/composition_cluster_umap_deconvolved_4.png"), dpi=400, height=5, width=7)


# Horizontal bar plot of cell type proportions for CCs
cluster_assignments <- data.frame(spot = rownames(deconv.norm), cluster = spatial$CompositionCluster_CC)

# Combine cluster info with spot proportions
deconv.df <- deconv.norm %>%
  as.data.frame() %>%
  mutate(spot = rownames(deconv.norm)) %>%
  left_join(cluster_assignments, by = "spot")

# Compute average
avg.comp <- deconv.df %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric), mean))  # average per cell type

avg.comp.long <- avg.comp %>%
  pivot_longer(
    cols = -cluster,
    names_to = "cell_type",
    values_to = "avg_proportion"
  )

ggplot(avg.comp.long, aes(x = cluster, y = avg_proportion, fill = cell_type)) +
  geom_col(width = 0.6,color = "black", size = 0.2) +
  theme_classic(base_size = 14) +
  coord_flip() +              # make it horizontal
  labs(x = "Composition Cluster", y = "Avg Cell Proportion", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(glue::glue("{resDir}ISCHIA/composition_cluster_stacked_barchart_4.png"), dpi=400, height=5, width=7)







