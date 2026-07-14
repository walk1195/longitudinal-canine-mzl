#############################################################################################
############################################################################################
### Processing of ORv-01 Longitudinal scRNAseq: Doublet Filtering, SCT Normalization, & Clustering
### Samples: ORv-01 Longitudinal scRNAseq (n=7)
### Author: Grace Walker
### Date: April 23, 2025
#############################################################################################

#############################################################################################
## SCTransform: Normalization / Scaling / FindVariableFeatures -- makes it so that high abundance genes are not dominating your downstream analysis
# Uses a Pearson's residuals transformation to account more for the distribution of expression
# Gives greater weight to genes with specific cell type expression than genes with broad expression
# Better at removing technical effects
#############################################################################################

# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(tidyverse) # dplyr and ggplot
library(Seurat) # for seurat analysis
library(patchwork) # for plotting orientation
library(ggplot2) # for visualization
library(scales) # for plotting
library(glmGamPoi) # for improved sctransform
library(DoubletFinder) # for detecting doublets
library(presto) # for faster FindMarkers 
library(clustree) # for optimal clustering

# Params
options(future.globals.maxSize = 4 * 1024^3)  # Set to 4 GB
gc() # Free up mem

# Dirs
objDir <- ""
figDir <- ""

# -------------------------------------------------------------------------------------------

# Read in data
s1 <- readRDS(objDir, "merged_filtered.rds")

# Resort samples chronologically
s1@meta.data$orig.ident <- factor(s1@meta.data$orig.ident, levels = unique(s1@meta.data$orig.ident))

#############################################################################################
# Inspecting the effect of regressing out MT expression
#############################################################################################

###### (1) Including MT ######
include_mt <- SCTransform(s1, verbose = TRUE) # Include MT percent expression
saveRDS(include_mt, paste0(objDir,'merged_processed_withMT.rds')) # Save progress

# Simple processing to visualize
include_mt <- RunPCA(include_mt, assay="SCT") %>%
  FindNeighbors(dims = 1:30) %>% # Default to 30 PCs for initial vis
  FindClusters(resolution=0.1) %>% # Initial vis
  RunUMAP(dims = 1:30)

# Visualize
DimHeatmap(include_mt, dims = 1:9, cells = 500, balanced = TRUE,ncol=3) # Top 10 PCs
ggsave(paste(figDir, 'top10pcs_withMT.png', sep=""), width = 8, height = 5, dpi=400)
DimHeatmap(include_mt, dims = 20:30, cells = 500, balanced = TRUE,ncol=3) # 20-30 PCs

ElbowPlot(include_mt, ndims = 40) # Elbow plot
ggsave(paste(figDir, 'elbow_plot_merged_withMT.png', sep=""), width = 8, height = 5, dpi=400)

DimPlot(include_mt, reduction = "umap", group.by = c("orig.ident", "seurat_clusters","sample_type"),
        alpha=0.4, ncol=3) # UMAP
ggsave(paste(figDir, 'umaps_withMT.png', sep=""), width = 12, height = 4, dpi=400)

FeaturePlot(include_mt, features='percent.mt')
ggsave(paste(figDir, 'umap_MTdistribution_withMT.png', sep=""), width = 5, height = 5, dpi=400)

VlnPlot(include_mt, features = "percent.mt", split.by = "seurat_clusters") # MT percent distribution across clusters
ggsave(paste(figDir, 'MTbycluster_withMT.png', sep=""), width = 8, height = 5, dpi=400)

# Save progress
saveRDS(include_mt, paste0(objDir,'merged_processed_withMT.rds'))

###### (2) Regressing out MT #######
regress_mt <- SCTransform(s1, vars.to.regress = "percent.mt", verbose = TRUE) # Regressing out MT % expression
saveRDS(regress_mt, paste0(objDir,'merged_processed_regressMT.rds')) # Save progress

# Simple processing to visualize
include_mt <- RunPCA(regress_mt, assay="SCT") %>%
  FindNeighbors(dims = 1:30) %>% # Default to 30 PCs for initial vis
  FindClusters(resolution=0.1) %>% # Initial vis
  RunUMAP(dims = 1:30)

# Visualize
ElbowPlot(regress_mt, ndims = 40) # Elbow plot
ggsave(paste(figDir, 'elbow_plot_merged_regressMT.png', sep=""), width = 8, height = 5, dpi=400)

# PC heatmaps
png(filename = paste(figDir, 'top10pcs_regressMT.png', sep=""), width = 8, height = 5, units = "in", res = 400)
DimHeatmap(regress_mt, dims = 1:9, cells = 500, balanced = TRUE,ncol=3) # Top 10 PCs
dev.off()

png(filename = paste(figDir, '10_20pcs_regressMT.png', sep=""), width = 8, height = 5, units = "in", res = 400)
DimHeatmap(regress_mt, dims = 10:20, cells = 500, balanced = TRUE,ncol=3) # 20-30 PCs
dev.off()

# Resort clusters
regress_mt$seurat_clusters <- factor(regress_mt$seurat_clusters, levels = rev(sort(unique(regress_mt$seurat_clusters))))
p1 <- DimPlot(regress_mt, reduction = "umap", group.by = "orig.ident", alpha = 0.4)
p2 <- DimPlot(regress_mt, reduction = "umap", group.by = "seurat_clusters", alpha = 0.4) 
p3 <- DimPlot(regress_mt, reduction = "umap", group.by = "sample_type", order = c('FNA', 'resect'), alpha = 0.4)
p1 + p2 + p3
ggsave(paste(figDir, 'umaps_regressMT.png', sep=""), width = 12, height = 4, dpi=400)

FeaturePlot(regress_mt, features='percent.mt')
ggsave(paste(figDir, 'umap_MTdistribution_regressMT.png', sep=""), width = 5, height = 5, dpi=400)

VlnPlot(regress_mt, features = "percent.mt", split.by = "seurat_clusters") # MT percent distribution across clusters
ggsave(paste(figDir, 'MTbycluster_regressMT.png', sep=""), width = 8, height = 5, dpi=400)

# UMAP across timepoints
DimPlot(regress_mt, reduction = "umap", group.by = "seurat_clusters", split.by = "orig.ident", alpha = 0.4, ncol=3) 
ggsave(paste(figDir, 'cluster_across_timepoints.png', sep=""), width = 8, height = 5, dpi=400)

# Sample contribution to each cluster
dittoBarPlot(var = "orig.ident", object = regress_mt, group.by = "seurat_clusters", x.labels.rotate=FALSE)
ggsave(paste(figDir, 'sample_contribution_by_cluster.png', sep=""), width = 8, height = 5, dpi=400)

# Save progress
saveRDS(regress_mt, paste0(objDir,'merged_processed_regressMT.rds')) 


#############################################################################################
# Running Doublet Detection
#############################################################################################

# Set annotations (clusters)
s1@meta.data$seurat_clusters <- s1@meta.data$cell_type
annotations <- s1@meta.data$seurat_clusters

# Split object
samp_split <- SplitObject(s1, split.by = "orig.ident")

## Custom DoubletFinder function
# Returns a dataframe with the cell IDs and a column with either 'Singlet' or 'Doublet'
run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL){
  # Print sample number
  print(paste0("Sample ", unique(seu_sample_subset[['orig.ident']]), '...........')) 
  
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # Set sample and PCs
  sample <- seu_sample_subset
  min_pc <- 25
  
  # pK identification (no ground-truth) 
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = TRUE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  ## Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # Run DoubletFinder
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj, # number of expected real doublets
                          sct=TRUE) # sct assay
  # Change name of metadata column with Singlet/Doublet information
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
  
  # Get res
  double_finder_res <- sample@meta.data['doublet_finder'] # get the metadata column with singlet, doublet info
  double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
  return(double_finder_res)
}


### Run DoubletFinder on all samples
samp_split_doubletfinder <- lapply(samp_split, run_doubletfinder_custom)

# Post process
doublet_results <- data.frame(bind_rows(samp_split_doubletfinder)) # merge results to a single dataframe
rownames(doublet_results) <- doublet_results$row_names # assign cell IDs to row names to ensure match
doublet_results$row_names <- NULL
head(doublet_results)

# Add to obj
s1 <- AddMetaData(s1, doublet_results, col.name = 'doublet_finder')

# Summary table
# Get doublets per sample
doublets_summary <- s1@meta.data %>% 
  group_by(orig.ident, doublet_finder) %>% 
  summarise(total_count = n(),.groups = 'drop') %>% as.data.frame() %>% ungroup() %>%
  group_by(orig.ident) %>%
  mutate(countT = sum(total_count)) %>%
  group_by(doublet_finder, .add = TRUE) %>%
  mutate(percent = paste0(round(100 * total_count/countT, 2),'%')) %>%
  dplyr::select(-countT)
write.table(doublets_summary, file=paste0(figDir, 'doubletfinder_results.txt'), quote = FALSE, row.names = FALSE, sep = '\t')

# Visualize results
VlnPlot(s1, group.by = 'orig.ident', split.by = "doublet_finder",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave(paste0(figDir, 'doublet_qc_plots.png'), width = 8, height = 5, dpi=400)

DimPlot(s1, reduction = "umap", group.by = c("orig.ident", "seurat_clusters", 'doublet_finder'), alpha=0.4, ncol=3) # UMAP
ggsave(paste0(figDir, 'doublet_umap_distribution.png'), width = 12, height = 5, dpi=400)

# Save
saveRDS(s1, paste0(objDir, 'merged_clustered_object.rds'))

#############################################################################################
# Optimizing number of PCs and cluster resolution
#############################################################################################

figDir = ""

# Run PCA
s1 <- RunPCA(s1, assay="SCT")
 
# Check elbow plot 
ElbowPlot(s1, ndims = 40)
ggsave(paste(figDir, 'elbow_plot.png', sep=""), width=8, height=5, dpi=400)

# Testing different # of PCs
pcs <- c(15, 20, 25, 30, 40, 45)
for (n in pcs) {
  s1 <- FindNeighbors(s1, dims = 1:n)
  s1 <- RunUMAP(s1, dims = 1:n)
  DimPlot(s1, reduction = "umap", group.by = "orig.ident", alpha = 0.4)
  ggsave(paste0(figDir, 'pc_testing/umap_', n, '_pcs.png'), width = 6, height = 5, dpi=400)
}

# Set optimal PCs
s1 <- FindNeighbors(s1, dims = 1:25)
s1 <- RunUMAP(s1, dims = 1:25)

# Generate clusters at different resolutions
resolutions = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
for (res in resolutions){
  # Cluster
  s1 <- FindClusters(s1, resolution=res, cluster.name=paste0('res',res)) # Setting smallest resolution
  # Plot UMAP
  DimPlot(s1, reduction = "umap", group.by = c("orig.ident", paste0('res',res)), alpha = 0.4)
  ggsave(paste0(figDir, 'umaps/umap', res, '.png'), width = 8, height = 5, dpi=400)
  
}

# Find top DEGs for clusters @ ideal resolution
cluster_markers <- FindAllMarkers(s1, group.by = 'res0.05', only.pos = TRUE)
write.csv(cluster_markers, paste0(figDir,'cluster_markers.csv')) # Save unfiltered cluster markers

# Save object
saveRDS(s1, paste0(objDir, 'merged_clustered_object.rds'))

# Examine cluster stability with clustree
clusterings <- s1@meta.data
clustree(clusterings, prefix = "res")
ggsave(paste(figDir, 'clustree_diagram.png', sep=""), width = 12, height = 7, dpi=400)


########################################################################
#           QC Plots -- Post Clustering
########################################################################

DimPlot(s1, group.by="orig.ident")
ggsave(paste0(figDir, 'umap_by_sample.png'), width = 6, height = 5, dpi=400)

FeaturePlot(s1, features='nFeature_RNA')
ggsave(paste0(figDir, 'umap_nFeature.png'), width = 6, height = 5, dpi=400)

FeaturePlot(s1, features='nCount_RNA')
ggsave(paste0(figDir, 'umap_nCount.png'), width = 6, height = 5, dpi=400)

FeaturePlot(s1, features='percent.mt')
ggsave(paste0(figDir, 'umap_MT.png'), width = 6, height = 5, dpi=400)

DimPlot(s1, group.by = "res0.05", split.by = "orig.ident", alpha = 0.4, ncol=3) 
ggsave(paste0(figDir, 'clusters_across_timepoints.png'), width = 8, height = 5, dpi=400)

#############################################################################################
# Plotting Top Cluster DEGs
#############################################################################################

figDir = ""
res = ""

# Cluster markers
cluster_markers <- read.csv(paste0(figDir,'cluster_markers.csv'))

# Filter for sig genes (p_val_adj < 0.05)
significant <- cluster_markers %>%
  filter(p_val_adj<0.05)

# Get top 10 for each cluster
top10 <- significant %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()

# Save
write.table(significant, file=paste0(figDir, 'cluster_markers_all_', res, ".tsv"), quote=FALSE, row.names=TRUE)
write.table(top10, file=paste0(figDir, 'cluster_markers_top10_', res, ".tsv"), quote=FALSE, row.names=TRUE)

# Heatmap plot
png(filename = paste0(figDir,'marker_heatmap_', res, '.png'), width = 8, height = 6, units = "in", res = 400)
p <- DoHeatmap(s1, features = top10$gene, draw.lines=TRUE, label=FALSE)
p + scale_fill_gradient(colors = c("cornflowerblue", "ivory1", "darkorange"))
dev.off()


########################################################################
#           Marker Gene Plots -- Post Clustering
########################################################################
figDir <- ""

bcell_markers <- c('CD79A', 'CD79B', 'PAX5', 'CD19', "CD20")
FeaturePlot(s1, features=bcell_markers)
ggsave(paste0(figDir, 'umap_bcell_markers.png'), width = 6, height = 5, dpi=400)

tcell_markers <- c('CD3E', 'CD3D', 'CD4', 'CD8A')
FeaturePlot(s1, features=tcell_markers)
ggsave(paste0(figDir, 'umap_tcell_markers.png'), width = 6, height = 5, dpi=400)

myeloid <- c('CD163', 'CD68', 'S100A8', 'CSF1R')
FeaturePlot(s1, features=myeloid)
ggsave(paste0(figDir, 'umap_myeloid_markers.png'), width = 6, height = 5, dpi=400)

prolif <- c('MKI67', 'TOP2A', 'CENPF','BUB1')
FeaturePlot(s1, features=prolif)
ggsave(paste0(figDir, 'umap_prolif_markers.png'), width = 6, height = 5, dpi=400)

ribo <- c('RPS19', 'RPL32-2', 'RPS26', 'RPS8')
ggsave(paste0(figDir, 'umap_ribo_markers.png'), width = 6, height = 5, dpi=400)



########################################################################
# SingleR Automated Annotation
########################################################################

figDir <- ""

# Import reference datasets
ref1 <- celldex::HumanPrimaryCellAtlasData() # Human primary scRNAseq atlas 
ref2 <- celldex::BlueprintEncodeData() # Human bulk RNAseq
ref3 <- celldex::ImmGenData() # Mouse immune bulk RNAseq

# Convert seurat object to sce object
DefaultAssay(s1) <- "SCT"
integrated_sce <- as.SingleCellExperiment(s1, assay = "SCT")

# Ensure log normalized counts
integrated_sce <- logNormCounts(integrated_sce)

# Subset to common genes
common1 <- intersect(rownames(integrated_sce),rownames(ref1))
common2 <- intersect(rownames(integrated_sce),rownames(ref2)) 
common3 <- intersect(rownames(integrated_sce),rownames(ref3))  

# Check common genes
length(common1) 
length(common2)
length(common3)

## Inspect gene names of each ##
head(rownames(ref3))
head(rownames(integrated_sce))

## Convert mouse gene names to uppercase
rownames(ref3) <- toupper(rownames(ref3))

# Subset
integrated_sce1 <- integrated_sce[common1,]
integrated_sce2 <- integrated_sce[common2,]
integrated_sce3 <- integrated_sce[common3,]

# Run SingleR with different reference datasets -- outputs dataframes
hpca_results_main <- SingleR(test=integrated_sce1, ref=ref1, labels=ref1$label.main, assay.type.ref = "logcounts")
hpca_results_fine <- SingleR(test=integrated_sce1, ref=ref1, labels=ref1$label.fine, assay.type.ref = "logcounts")
blue_results <- SingleR(test=integrated_sce2, ref=ref2, labels=ref2$label.main, assay.type.ref = "logcounts")
immgen_results <- SingleR(test=integrated_sce3, ref=ref3, labels=ref3$label.main, assay.type.ref = "logcounts")

# Save results as csv
write.table(hpca_results_main, paste0(figDir, "SingleR_Results_HPCA_Main.tsv"), sep = "\t", quote = FALSE)
write.table(hpca_results_fine, paste0(figDir, "SingleR_Results_HPCA_Fine.tsv"), sep = "\t", quote = FALSE)
write.table(blue_results, paste0(figDir, "SingleR_Results_BlueEncode.tsv"), sep = "\t", quote = FALSE)
write.table(immgen_results, paste0(figDir, "SingleR_Results_ImmGen.tsv"), sep = "\t", quote = FALSE)


### Plot Results
# Plot score heatmaps
png(paste0(figDir, "SingleR_Heatmap_CellAssignmentScores_HPCA_Main.png"),height=600,width=900)
plotScoreHeatmap(hpca_results_main)
dev.off()

png(paste0(figDir, "SingleR_Heatmap_CellAssignmentScores_HPCA_Fine.png"),height=800,width=1100)
plotScoreHeatmap(hpca_results_fine)
dev.off()

png(paste0(figDir, "SingleR_Heatmap_CellAssignmentScores_BlueEncode.png"),height=600,width=900)
plotScoreHeatmap(blue_results)
dev.off()

png(paste0(figDir, "SingleR_Heatmap_CellAssignmentScores_ImmGen.png"),height=600,width=900)
plotScoreHeatmap(immgen_results)
dev.off()

# Plot delta distributions
png(paste0(figDir, "delta_distributions_hpca_main.png"),height=1000,width=800)
plotDeltaDistribution(hpca_results_main, ncol = 4, dots.on.top = FALSE)
dev.off()

png(paste0(figDir, "delta_distributions_hpca_fine.png"),height=1000,width=800)
plotDeltaDistribution(hpca_results_fine, ncol = 4, dots.on.top = FALSE)
dev.off()

png(paste0(figDir, "delta_distributions_hpca_fine.png"),height=1000,width=800)
plotDeltaDistribution(blue_results, ncol = 4, dots.on.top = FALSE)
dev.off()

png(paste0(figDir, "delta_distributions_immgen.png"),height=1000,width=800)
plotDeltaDistribution(immgen_results, ncol = 4, dots.on.top = FALSE)
dev.off()


### Merge results with object
rownames(hpca_results_main)[1:5] # Checking cell IDs
s1 <- AddMetaData(s1, hpca_results_main$labels, col.name = 'hpca_main_raw')
s1 <- AddMetaData(s1, hpca_results_main$pruned.labels, col.name = 'hpca_main_pruned')

rownames(hpca_results_fine)[1:5] # Checking cell IDs
s1 <- AddMetaData(s1, hpca_results_fine$labels, col.name = 'hpca_fine_raw')
s1 <- AddMetaData(s1, hpca_results_fine$pruned.labels, col.name = 'hpca_fine_pruned')

rownames(blue_results)[1:5] # Checking cell IDs
s1 <- AddMetaData(s1, blue_results$labels, col.name = 'blue_raw')
s1 <- AddMetaData(s1, blue_results$pruned.labels, col.name = 'blue_pruned')

rownames(immgen_results)[1:5] # Checking cell IDs
s1 <- AddMetaData(s1, immgen_results$labels, col.name = 'immgen_raw')
s1 <- AddMetaData(s1, immgen_results$pruned.labels, col.name = 'immgen_pruned')


# UMAPs - HPCA main
Idents(s1) <- s1@meta.data$hpca_fine_pruned
DimPlot(s1, label = F , repel = T, label.size = 3)
ggsave(paste0(figDir, 'hpca_main_umap.png'), height=5, width=8, dpi=400)

# B cells
bcells <- c("B_cell:Memory", "B_cell:Naive", "B_cell", "B_cell:Germinal_center")
s1$bcells <- ifelse(s1$hpca_fine_pruned %in% bcells, 
                         s1$hpca_fine_pruned, 
                         "Other")
DimPlot(s1, group.by = "bcells") + ggtitle('B Cells')
ggsave(paste0(figDir, 'umap_bcells.png'), height=6, width=8, dpi=400)

# Immune cells
immune <- c("B_cell:Memory", "B_cell:Naive", "B_cell", "B_cell:Germinal_center", "T_cell:CD8+", "T_cell:CD4+", "T_cell:CD4+_Naive", "T_cell:Treg:Naive", "NK_cell")
s1$immune <- ifelse(s1$hpca_fine_pruned %in% immune, 
                            s1$hpca_fine_pruned, 
                            "Other")
DimPlot(s1, reduction='umap.harmony', group.by = "immune") + ggtitle('Immune Cells')
ggsave(paste0(figDir, 'umap_immunecells.png'), height=6, width=8, dpi=400)


# UMAPs - BlueEncode
Idents(s1) <- "blue_pruned"
DimPlot(s1, reduction='umap.harmony')
ggsave(paste0(figDir, 'umap_blueencode.png'), height=7, width=8, dpi=400)

# UMAPs - ImmGen
Idents(s1) <- "immgen_pruned"
DimPlot(s1, reduction='umap.harmony')
ggsave(paste0(figDir, 'umap_immgen.png'), height=5, width=8, dpi=400)


