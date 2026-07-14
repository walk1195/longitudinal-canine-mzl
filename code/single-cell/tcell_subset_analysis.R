#############################################################################################
############################################################################################
### ORv01 Longitudinal scRNAseq Analysis: Subclustering of T cell Population
### Samples: ORv-01 longitudinal scRNAseq (n=7), T cells only
### Author: Grace Walker
### Date: Sep 11, 2025
#############################################################################################
#############################################################################################

# Set up environment
# -------------------------------------------------------------------------------------------

# Load libs
library(tidyverse) # dplyr and ggplot
library(Seurat) # for seurat analysis
library(patchwork) # for plotting orientation
library(ggplot2) # for visualization
library(scran) # singler normalization
library(EnhancedVolcano) # for dge visualization
library(clustree) # for cluster optimization
library(celldex) # for ref datasets
library(SingleR) # for automated cell typing
library(dittoSeq) # for stacked bar plots
library(jsonlite) # for json handling
library(ggpubr) # for boxplots
library(ggalluvial)  # for alluvial plots
library(clusterProfiler) # for GO analysis
library(RColorBrewer) # for a colourful plot
library(pheatmap) # for heatmaps
library(biomaRt) # for gene symbol conversions
library(UCell) # for gene signature scoring
library(enrichplot) # for pathway results visualizations
library(fgsea) # for gene set enrichment analysis
library(msigdbr) # for test gene sets
library(presto) # for gene set enrichment analysis
library(scales)

# Params
options(future.globals.maxSize = 4 * 1024^3)  # Set to 4 GB
gc() # free up memory

# Set volcano plot theme
theme_set(theme_classic(base_size = 10) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(0.5), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(0.6), color = 'black'),
              plot.title = element_text(hjust = 0.5)))

# Dirs
figDir <- ""
objDir <- ""
resDir <- ""

# -------------------------------------------------------------------------------------------

# Read in data
s1 <- readRDS(paste0(objDir, 'processed_object_final.rds'))
s1.tcells <- readRDS(paste0(objDir, 'tcell_subset_object_final.rds'))


##############################################################################################
# Subset to T cell cluster
##############################################################################################
# Subset
s1.tcells <- subset(s1, subset=cell_type=='T cells')

##############################################################################################
# Initial Plots
##############################################################################################
Idents(s1.tcells) <- s1.tcells$cell_type_final
s1.tcells$orig.ident <- factor(s1.tcells$orig.ident, levels = c('PreTx', 'D3FNA', 'D7FNA', 'D10FNA', 'D10resect', 'D28FNA', 'D49FNA'))


# Cluster contributions by sample
dittoBarPlot(var = "orig.ident", object = s1.tcells, group.by = "cell_type_final", x.labels.rotate=FALSE)
ggsave(paste0(figDir, 'cluster_contribution_by_sample.png'), width = 8, height = 5, dpi=400)

# Sample contributions by cluster
dittoBarPlot(var = "cell_type_final", object = s1.tcells, group.by = "orig.ident", x.labels.rotate=FALSE)
ggsave(paste0(figDir, 'sample_contribution_by_cluster.png'), width = 8, height = 5, dpi=400)

# Cell type by sample
p <- dittoBarPlot(var = "cell_type_final", object = s1.tcells, group.by = "orig.ident") +
  scale_x_discrete(limits = c("PreTx", "D3FNA", "D7FNA", "D10FNA", "D10resect", "D28FNA", "D49FNA"))
ggsave(paste0(figDir, 'cell_type_by_sample_dittobar.png'), plot=p, width = 8, height = 5, dpi=400)


### Plot T cell counts over time
s1 <- s1.tcells
# Format data
meta_data <- s1@meta.data %>%
  as_tibble() %>%
  mutate(orig.ident = as.character(orig.ident)) # Convert factor to character
tcell_counts <- meta_data %>%
  filter(cell_type == "T cells") %>%  
  dplyr::count(orig.ident, name = "tcell_count") # Specify dplyr count function

# Plot
tcell_counts$orig.ident <- factor(tcell_counts$orig.ident, levels = unique(s1$orig.ident)) # chronological order
ggplot(tcell_counts, aes(x = orig.ident, y = tcell_count, fill = orig.ident)) +
  geom_bar(stat = "identity", color='black') +
  labs(title = "T Cell Count By Sample", x = "Sample", y = "Cell Count") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5, face="bold"))
ggsave(paste0(figDir, 't_cells_counts_by_sample.png'), dpi=400, height=5, width=7)


### Plot percent of total cells (normalized cell count)
total_counts <- meta_data %>%
  dplyr::count(orig.ident, name = "total_cells")
# Merge total counts with t cell counts
tcell_counts <- left_join(tcell_counts, total_counts, by = "orig.ident") %>%
  mutate(pct_tcells = (tcell_count / total_cells) * 100)  # Calculate percentage
tcell_counts$orig.ident <- factor(tcell_counts$orig.ident, levels = unique(s1$orig.ident))
# Plot
ggplot(tcell_counts, aes(x = orig.ident, y = pct_tcells, fill = orig.ident)) +
  geom_bar(stat = "identity", color='black') +
  labs(title = "T Cell Proportion By Sample", x = "Sample", y = "% of total cells") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5, face="bold"))
ggsave(paste0(figDir, 't_cell_proportion_by_sample.png'), dpi=400, height=5, width=7)


##############################################################################################
# Normalization & Clustering
##############################################################################################
figDir <- ""

s1.tcells <- SCTransform(s1.tcells, verbose = TRUE) # Normalize
s1.tcells <- RunPCA(s1.tcells, assay="SCT") 
ElbowPlot(s1.tcells, ndims = 40) # Elbow plot
s1.tcells <- FindNeighbors(s1.tcells, dims = 1:30)
s1.tcells <- RunUMAP(s1.tcells, dims = 1:30) # only do this once
DimPlot(s1.tcells, reduction = "umap", group.by = c("orig.ident", 'cell_type'), alpha = 0.4)
ggsave(paste0(figDir, 'subset_umap_unlabelled.png'), dpi=400, height=5, width=9)

# Cluster
resolutions = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
for (res in resolutions){
  # Cluster
  s1.tcells <- FindClusters(s1.tcells, resolution=res, cluster.name=paste0('res',res)) # Setting smallest resolution
  # Plot UMAP
  #s1.tcells <- RunUMAP(s1.tcells, dims = 1:25) # only do this once
  DimPlot(s1.tcells, reduction = "umap", group.by = paste0('res',res), alpha = 0.4)
  ggsave(paste0(figDir, 'umaps/clusters_', res, '.png'), width = 7, height = 5, dpi=400)
}
# Clustree to optimize
clusterings <- s1.tcells@meta.data
clustree(clusterings, prefix = "res")
ggsave(paste(figDir, 'clustree_diagram_tcell_subset.png', sep=""), width = 7, height = 9, dpi=400)

### Plot nFeatures distribution
FeaturePlot(s1.tcells, features="nFeature_RNA", alpha = 0.4)
ggsave(paste(figDir, 'feature_distribution.png', sep=""), width = 7, height = 9, dpi=400)

# Save object
saveRDS(s1.tcells, paste0(objDir, '	tcell_subset_object_final.rds'))


##############################################################################################
# Inspecting Markers
##############################################################################################
figDir <- paste0(figDir, 'markers/')

rownames(s1.tcells)[grep('ITGAX', rownames(s1.tcells))] # To check for a marker

# Markers
tcell_markers <- c('CD3E', 'CD4', 'CD8A', 'CD8B')
nk_markers <- c('GZMB', 'GZMA', 'KLRB1', 'NCAM1', 'NCAM2') # KLRB1 = CD161
treg_markers <- c('FOXP3', 'IL2RA')
ifng <- c("IFNGR1", "IFNG", "IFNGR2")
th17_markers <- c("IL17A", "IL17F", "IL22", "IL23R", "IL1R1")
bcell_markers <- c('CD19', 'CD79A', 'CD79B', 'PAX5', 'MS4A1')


# Feature plots
FeaturePlot(s1.tcells, features=tcell_markers)
ggsave(paste0(figDir, 'tcell_markers.png'), width = 7, height = 5, dpi=400)
FeaturePlot(s1.tcells, features=nk_markers)
ggsave(paste0(figDir, 'nk_markers.png'), width = 7, height = 5, dpi=400)
FeaturePlot(s1.tcells, features=treg_markers)
ggsave(paste0(figDir, 'treg_markers.png'), width = 7, height = 5, dpi=400)
FeaturePlot(s1.tcells, features=ifng)
ggsave(paste0(figDir, 'ifng_markers.png'), width = 7, height = 5, dpi=400)
FeaturePlot(s1.tcells, features=th17_markers)
ggsave(paste0(figDir, 'th17_markers.png'), width = 7, height = 5, dpi=400)
FeaturePlot(s1.tcells, features=bcell_markers)
ggsave(paste0(figDir, 'bcell_markers.png'), width = 7, height = 5, dpi=400)

# Individual plots
FeaturePlot(s1.tcells, features='KLRB1') + labs(title='KLRB1 (CD161)')
ggsave(paste0(figDir, 'cd161_nk_marker.png'), width = 7, height = 5, dpi=400)
FeaturePlot(s1.tcells, features='IL2RA') + labs(title='IL2RA (CD25)') 
ggsave(paste0(figDir, 'cd25_treg_marker.png'), width = 7, height = 5, dpi=400)


# Cluster markers
figDir <- ""

resolutions <- c('res0.06')

for (res in resolutions) {
  # Set Idents
  Idents(s1.tcells) <- s1.tcells@meta.data[[res]]
  # Prep SCT
  s1.tcells <- PrepSCTFindMarkers(s1.tcells,verbose=T)
  # Get markers
  markers <- FindAllMarkers(s1.tcells, only.pos = TRUE, recorrect_umi=FALSE)

  # Filter for significant genes (p_val_adj < 0.05)
  significant <- markers %>%
    filter(p_val_adj<0.05)
  
  # Heatmap plot
  png(filename = paste0(figDir, 'marker_heatmap_', res, '.png'), width = 8, height = 6, units = "in", res = 400)
  top10 <- significant %>%
    group_by(cluster) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # Plot
  p <- DoHeatmap(subset(s1.tcells,  downsample = 6000), features = top10$gene, draw.lines=TRUE, label=FALSE)
  p + scale_fill_gradientn(colors = c("cornflowerblue", "ivory1", "darkorange"))
  dev.off()
  
  # Save
  write.table(significant, file=paste0(figDir, 'cluster_markers_all_', res, ".tsv"), quote=FALSE, row.names=TRUE)
  write.table(top10, file=paste0(figDir, 'cluster_markers_top10_', res, ".tsv"), quote=FALSE, row.names=TRUE)
  
  
  # Plot over time
  s1.tcells@meta.data$orig.ident <- factor(s1.tcells@meta.data$orig.ident, levels = c('PreTx', 'D3FNA', 'D7FNA', 'D10FNA', 'D10resect', 'D28FNA', 'D49FNA'))
  
  # Plot
  DimPlot(s1.tcells, group.by = res, split.by = "orig.ident", alpha = 0.4, ncol=3) +
    ggtitle('Cell type by sample')
  ggsave(paste0(figDir, 'tcell_clusters_across_timepoints_', res, '.png'), width = 6, height = 5, dpi=400)
}

# Assign annotations
### FINAL SUBSET -- res 0.06
cluster_annotations <- c("0" = "CD4+ T cells", 
                         "1" = "NK cells", 
                         "2" = "Naive B cells", 
                         "3" = "Naive T cells", 
                         "4" = "Proliferating T cells",
                         "5" = "CD8+ T cells",
                         "6" = "IL23R+ T cells")

# Map annotations to clusters
s1.tcells@meta.data$cell_type_final <- cluster_annotations[as.character(s1.tcells@meta.data$res0.06)]

### Plotting
# Count the number of cells per sample in each cluster
cluster_counts <- s1.tcells@meta.data %>%
  group_by(cell_type, orig.ident) %>%
  summarise(cell_count = n(), .groups = "drop")

# Compute the proportion of each sample within the cluster
cluster_counts <- cluster_counts %>%
  group_by(cell_type) %>%
  mutate(proportion = cell_count / sum(cell_count))

# Plot
ggplot(cluster_counts, aes(x = orig.ident, y = proportion, fill = orig.ident)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~cell_type, scales = "free_y", labeller = labeller(cell_type = function(x) {
    paste0(x, " (n=", cluster_total_cells$total_cells[cluster_total_cells$cell_type == x], ")")
  })) +
  labs(title = "Proportion of Cells per Sample in Each Cluster",
       x = "Sample (orig.ident)",
       y = "Proportion of Total Cells in Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(figDir, 'tcell_types_in_each_sample.png', sep=""), width = 8, height = 5, dpi=400)

# Additional plots
DimPlot(s1.tcells, group.by='cell_type_final')
ggsave(paste(figDir, 'umap_cell_types.png', sep=""), width = 6, height = 5, dpi=400)

DimPlot(s1.tcells, group.by='res0.06')
ggsave(paste(figDir, 'umap_clusters.png', sep=""), width = 6, height = 5, dpi=400)

FeaturePlot(s1.tcells, features='nCount_RNA')
ggsave(paste(figDir, 'count_distribution.png', sep=""), width = 7, height = 5, dpi=400)

FeaturePlot(s1.tcells, features='nFeature_RNA')
ggsave(paste(figDir, 'feature_distribution.png', sep=""), width = 7, height = 5, dpi=400)

FeaturePlot(s1.tcells, features='percent.mt')
ggsave(paste(figDir, 'mt_distribution.png', sep=""), width = 7, height = 5, dpi=400)


##############################################################################################
# SingleR Automated Annotation
##############################################################################################

# Import reference datasets
ref1 <- celldex::HumanPrimaryCellAtlasData() # Human primary scRNAseq atlas 
ref2 <- celldex::ImmGenData() # Mouse immune bulk RNAseq

# Convert seurat object to sce object
DefaultAssay(s1.tcells) <- "SCT"
integrated_sce <- as.SingleCellExperiment(s1.tcells, assay = "SCT")
# Ensure log normalized counts
integrated_sce <- logNormCounts(integrated_sce)

# Subset to common genes
common1 <- intersect(rownames(integrated_sce),rownames(ref1))
common2 <- intersect(rownames(integrated_sce),rownames(ref2)) 

## Inspect gene names of each ##
head(rownames(ref2))
head(rownames(integrated_sce))
## Convert mouse gene names to uppercase
rownames(ref2) <- toupper(rownames(ref2))

# Check common genes
length(common1) 
length(common2)

# Subset
integrated_sce1 <- integrated_sce[common1,]
integrated_sce2 <- integrated_sce[common2,]


# Run SingleR with different reference datasets -- outputs dataframes
hpca_results_main <- SingleR(test=integrated_sce1, ref=ref1, labels=ref1$label.main, assay.type.ref = "logcounts")
hpca_results_fine <- SingleR(test=integrated_sce1, ref=ref1, labels=ref1$label.fine, assay.type.ref = "logcounts")
immgen_results <- SingleR(test=integrated_sce2, ref=ref2, labels=ref2$label.main, assay.type.ref = "logcounts")

# Save results as csv
write.table(hpca_results_main, paste0(figDir, "SingleR_Results_HPCA_Main.tsv"), sep = "\t", quote = FALSE)
write.table(hpca_results_fine, paste0(figDir, "SingleR_Results_HPCA_Fine.tsv"), sep = "\t", quote = FALSE)
write.table(immgen_results, paste0(figDir, "SingleR_Results_ImmGen.tsv"), sep = "\t", quote = FALSE)


### Plotting Results ###
# Plot score heatmaps
png(paste0(figDir,"SingleR_Heatmap_CellAssignmentScores_HPCA_Main.png"),height=600,width=900)
plotScoreHeatmap(hpca_results_main, legend_labels=c(round(min(hpca_results_main$scores), 2), round(max(hpca_results_main$scores),2)))
dev.off()

png(paste0(figDir,"SingleR_Heatmap_CellAssignmentScores_HPCA_Fine.png"),height=800,width=1100)
plotScoreHeatmap(hpca_results_fine, legend_labels=c(round(min(hpca_results_fine$scores), 2), round(max(hpca_results_fine$scores),2)))
dev.off()

png(paste0(figDir,"SingleR_Heatmap_CellAssignmentScores_ImmGen.png"),height=600,width=900)
plotScoreHeatmap(immgen_results, legend_labels=c(round(min(immgen_results$scores), 2), round(max(immgen_results$scores),2)))
dev.off()

# Plot delta distributions
png("delta_distributions_hpca_main.png",height=1000,width=800)
plotDeltaDistribution(hpca_results_main, ncol = 4, dots.on.top = FALSE)
dev.off()

png("delta_distributions_hpca_fine.png",height=1000,width=800)
plotDeltaDistribution(hpca_results_fine, ncol = 4, dots.on.top = FALSE)
dev.off()

png("delta_distributions_immgen.png",height=1000,width=800)
plotDeltaDistribution(immgen_results, ncol = 4, dots.on.top = FALSE)
dev.off()


### Merge results with object ###
rownames(hpca_results_main)[1:5] # Checking cell IDs
s1.tcells <- AddMetaData(s1.tcells, hpca_results_main$labels, col.name = 'hpca_main_raw')
s1.tcells <- AddMetaData(s1.tcells, hpca_results_main$pruned.labels, col.name = 'hpca_main_pruned')

rownames(hpca_results_fine)[1:5] # Checking cell IDs
s1.tcells <- AddMetaData(s1.tcells, hpca_results_fine$labels, col.name = 'hpca_fine_raw')
s1.tcells <- AddMetaData(s1.tcells, hpca_results_fine$pruned.labels, col.name = 'hpca_fine_pruned')

rownames(immgen_results)[1:5] # Checking cell IDs
s1.tcells <- AddMetaData(s1.tcells, immgen_results$labels, col.name = 'immgen_raw')
s1.tcells <- AddMetaData(s1.tcells, immgen_results$pruned.labels, col.name = 'immgen_pruned')


# Plot - HPCA main
Idents(s1.tcells) <- s1.tcells@meta.data$hpca_main_pruned
DimPlot(s1.tcells, label = F , repel = T, label.size = 3, cols = palette30)
ggsave(paste0(figDir, 'hpca_main_umap.png'), height=5, width=7, dpi=400)


# Plot - HPCA fine
Idents(s1.tcells) <- s1.tcells@meta.data$hpca_fine_pruned

celltypes <- levels(s1.tcells)
palette30 <- Polychrome::createPalette(length(celltypes), seedcolors = c("#000000", "#FF0000"))
names(palette30) <- celltypes

DimPlot(s1.tcells, label = F , repel = T, label.size = 3, cols=palette30)
ggsave(paste0(figDir, 'hpca_fine_umap.png'), height=5, width=12, dpi=400)

# Just immune cells
immune <- c("Pro-Myelocyte", "HSC_-G-CSF", "CMP", "BM", "Endothelial_cells:blood_vessel", "Endothelial_cells:lymphatic", "GMP", "NA", "Platelets")
s1.tcells$immune <- ifelse(s1.tcells$hpca_fine_pruned %in% immune, 
                           "Other", 
                           s1.tcells$hpca_fine_pruned)
DimPlot(s1.tcells, group.by = "immune") + ggtitle('Immune Cells')
ggsave(paste0(figDir, 'hpca_fine_umap_immune_only.png'), height=6, width=12, dpi=400)


# Plot - ImmGen
Idents(s1.tcells) <- "immgen_pruned"
DimPlot(s1.tcells)
ggsave(paste0(figDir, 'umap_immgen.png'), height=5, width=8, dpi=400)


# Save object
saveRDS(s1.tcells, paste0(objDir, 'tcell_object_final.rds'))

##############################################################################################
# Cluster predicted cell type histograms 
##############################################################################################

# For each cluster (there are 6), we want a histogram of the the SingleR HPCA fine cell type predictions
for (cluster in unique(s1.tcells$res0.06)) {
  
  # Subset cells for the current cluster
  cells_in_cluster <- s1.tcells@meta.data %>%
    filter(res0.06 == cluster)
  # Ensure that hpca_fine_pruned is a factor
  cells_in_cluster$hpca_fine_pruned <- as.factor(cells_in_cluster$hpca_fine_pruned)
  # Use table to count cell types
  celltype_counts <- table(cells_in_cluster$hpca_fine_pruned) %>%
    as.data.frame() %>%
    arrange(desc(Freq))  # Sort by frequency
  # Filter out cell types with fewer than 10 cells
  celltype_counts <- celltype_counts %>%
    filter(Freq >= 10)
  # Reorder factor levels by frequency
  cells_in_cluster$hpca_fine_pruned <- factor(
    cells_in_cluster$hpca_fine_pruned,
    levels = celltype_counts$Var1)
  
  ggplot(cells_in_cluster, aes(x = hpca_fine_pruned)) +
    geom_bar(fill = "steelblue", color = "black") +
    theme_classic() +
    labs(title = paste("Cluster", cluster, 'SingleR Predictions (HPCA Fine)'),
         x = "Cell Type", y = "Cell Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(figDir, 'histograms/cluster', cluster, '_hpca_predictions.png'), width=9, height=5, dpi=400)
}

##############################################################################################
# DGEA 
##############################################################################################
### We are interested in defining differential expression between T cells of all pairwise timepoints

# Add a col grouping the T cells together by sample
s1.tcells$celltype.sample <- paste(s1.tcells$cell_type, s1.tcells$orig.ident, sep = "_")
# Set idents
Idents(s1.tcells) <- s1.tcells$celltype.sample

# Subset 1: 500 cells per dataset
num_cells_per_sample <- 500 # No. cells
set.seed(163)  # For reproducibility
cells_to_keep <- unlist(lapply(unique(s1.tcells$orig.ident), function(sample_id) {
  # Sample 500 cells from each sample
  sample(Cells(s1.tcells)[s1.tcells$orig.ident == sample_id], num_cells_per_sample)
}))
s1.tcells.subset <- subset(s1.tcells, cells = cells_to_keep)


# Subset 2: Remove cells annotated as B cells
annotation_cols <- c("immgen_pruned", "hpca_fine_pruned", "hpca_main_pruned")
# Remove if all 3 cols contain B cell annotation
cells_to_remove <- apply(s1.tcells@meta.data[, annotation_cols], 1, function(x) {
  all(grepl("B[ _]cell", x, ignore.case = TRUE))  # Matches "B cell" OR "B_cell"
})
# Subset object
s1.tcells.cleaned <- subset(s1.tcells, cells = colnames(s1.tcells)[!cells_to_remove])


# Subset 3: Immgen predicted T cells only (based on SingleR annotations of entire dataset)
s1.tcells.singler <- subset(s1, subset=immgen_pruned=='T cells')

### Now we have 3 objects:
######### s1.tcells.subset -> 500 cells per sample, cluster 2 only (manual 'T cells')
######### s1.tcells.cleaned -> not subset, but removed cells SingleR predicted as B cell
######### s1.tcells.singler -> 500 cells per sample, SingleR Immgen 'T cells' only

# Group T cells together by sample
s1.tcells.singler$celltype.sample <- paste(s1.tcells.singler$immgen_pruned, s1.tcells.singler$orig.ident, sep = "_")
# Set idents
Idents(s1.tcells.singler) <- s1.tcells.singler$celltype.sample

# Prep SCT - skip for subsets (minimum UMI unchanged)
s1.tcells.subset <- PrepSCTFindMarkers(s1.tcells.subset,verbose=T) # General T cells

# Set the correct path
figDir <- ""

# Initialize sample vector - general T cells
sample_ids <- unique(s1.tcells.subset$celltype.sample) # Subset
sample_ids <- unique(s1.tcells.cleaned$celltype.sample) # No subset but cleaned
sample_ids <- unique(s1.tcells.singler$celltype.sample) # No subset but cleaned


# Set global obj
#curr_object <- s1.tcells.subset
#curr_object <- s1.tcells.cleaned
curr_object <- s1.tcells.singler

#### Run pairwise comparisons - General T cells
while (length(sample_ids) > 0) {
  
  # Get current sample (entire ID)
  curr_sample1 <- sample_ids[1]
  # Remove from list
  sample_ids <- sample_ids[-1]
  
  # Get curr sample (just name)
  curr_sample_name <- strsplit(curr_sample1, "_")[[1]][2]
  
  # Start of loop
  for (sample in sample_ids) {
    print(paste('Getting markers for', curr_sample1, 'vs', sample))
    
    # Get sample name
    split_str <- strsplit(sample, "_")
    timepoint <- split_str[[1]][2]
    
    # Get markers
    markers <- FindMarkers(object = curr_object, 
                           ident.1 = curr_sample1, 
                           ident.2 = sample,
                           recorrect_umi=FALSE)
    
    # Sort by log2FC
    markers_ordered <- markers[order(-markers$avg_log2FC), ]
    significant_genes <- markers_ordered %>%
      filter(p_val_adj < 0.05)
  
    # Save filtered gene list
    write.table(significant_genes, file=paste0(figDir, curr_sample_name, '_', timepoint, '_significant_genes.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
    # Save unfiltered gene list
    write.table(markers, file=paste0(figDir, curr_sample_name, '_', timepoint, '_unfiltered_genes.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
    
    
    # Plot
    # Create gene name col
    markers$genes <- rownames(markers)
    # Create col for up/down regulated genes
    markers$diffexpressed <- "NO"
    markers$diffexpressed[markers$avg_log2FC > 0.6 & markers$p_val_adj < 0.05] <- "UP"
    markers$diffexpressed[markers$avg_log2FC < -0.6 & markers$p_val_adj < 0.05] <- "DOWN"
    # Create label col for top 30
    markers$top30 <- ifelse(markers$genes %in% head(markers[order(markers$p_val_adj), "genes"], 30) & markers$diffexpressed != "NO", markers$genes, NA)
    # Add genes with largest log2fc
    top_log2fc <- c(rownames(tail(significant_genes, n=10)), rownames(head(significant_genes, n=10)))
    markers$top30 <- ifelse(is.na(markers$top30) & markers$genes %in% top_log2fc & markers$diffexpressed != "NO", markers$genes, markers$top30)
    
  
    # Volcano plot
    markers$p_val_adj[markers$p_val_adj == 0] <- 1e-300
    y_max <- max(-log10(markers$p_val_adj), na.rm = TRUE) + 10  # Get max y val with buffer of 10
    
    ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = top30)) +
      geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
      geom_point(size = 0.5) +
      coord_cartesian(ylim = c(0, y_max), xlim = c(-7, 7)) +
      scale_color_manual(values = c("blue", "grey", "red"),
                         labels = c("Downregulated", "Not significant", "Upregulated")) +
      scale_x_continuous(breaks = seq(-10, 10, 2)) +
      labs(color = 'Severe') +
      geom_text_repel(size=3, max.overlaps = Inf) +
      ggtitle(paste(curr_sample1, "vs", timepoint))
    ggsave(paste0(figDir, curr_sample_name, '_', timepoint, '_volcano.png'), width=7, height=5, dpi=400)
  }
}

### We also are interested in running DGEA for T cell subtypes between timepoints
# This object contains T cell subtype labels
s1.tcells <- PrepSCTFindMarkers(s1.tcells,verbose=T)

# Initialize sample vector - subtypes
sample_ids <- unique(s1.tcells$celltype.sample)

curr_object <- s1.tcells

figDir <- ""

# Group T cells by subtype
s1.tcells$celltype.sample <- paste(s1.tcells$cell_type_final, s1.tcells$orig.ident, sep = "_")
# Set idents
Idents(s1.tcells) <- s1.tcells$celltype.sample

# Prep
s1.tcells <- PrepSCTFindMarkers(s1.tcells,verbose=T)

# Set sample vectors
sample_ids <- c("PreTx", "D3FNA", "D7FNA", "D10FNA", "D10resect", 'D28FNA', 'D49FNA')
subtype_sample_ids <- unique(s1.tcells$celltype.sample)

### Run pairwise comparisons - Subtype specific
while (length(sample_ids) > 0) {
  
  # Get current sample (entire ID)
  curr_sample1 <- sample_ids[1]
  # Remove from list
  sample_ids <- sample_ids[-1]
  
  # Get curr sample (just name)
  curr_sample_name <- curr_sample1
  
  # Start of loop
  for (sample in sample_ids) {
    
    # Make results dir
    resDir <- figDir
    subresDir <- paste0(curr_sample1, '_', sample)

    if (file.exists(subresDir)){
      setwd(file.path(resDir, subresDir))
    } else {
      dir.create(file.path(resDir, subresDir)) # Create if doesn't exist
      setwd(file.path(resDir, subresDir)) # Set as current wd
    }
    
    for (subtype in subtype_sample_ids) {

      # Get ids for subtypes
      curr_sample1_subtype <- paste0(subtype, '_', curr_sample1)
      sample_subtype <- paste0(subtype, '_', sample)
      
      
      print(paste('Getting markers for', curr_sample1, 'vs', sample, subtype))

      # Get markers
      markers <- FindMarkers(object = curr_object, 
                             ident.1 = curr_sample1_subtype, 
                             ident.2 = sample_subtype,
                             recorrect_umi=FALSE)
      
      # Sort by log2FC
      markers_ordered <- markers[order(-markers$avg_log2FC), ]
      significant_genes <- markers_ordered %>%
        filter(p_val_adj < 0.05)
      
      # Save filtered gene list
      if (grepl("\\+", subtype)) {
        subtype_fixed <- gsub("\\+", "", subtype)
      } else {
        subtype_fixed <- subtype
      }
      
      write.table(significant_genes, file=paste0(subtype_fixed, '_significant_genes.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
      # Save unfiltered gene list
      write.table(markers, file=paste0(subtype_fixed, '_unfiltered_genes.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
      
      
      # Plot
      # Create gene name col
      markers$genes <- rownames(markers)
      # Create col for up/down regulated genes
      markers$diffexpressed <- "NO"
      markers$diffexpressed[markers$avg_log2FC > 0.6 & markers$p_val_adj < 0.05] <- "UP"
      markers$diffexpressed[markers$avg_log2FC < -0.6 & markers$p_val_adj < 0.05] <- "DOWN"
      # Create label col for top 30
      markers$top30 <- ifelse(markers$genes %in% head(markers[order(markers$p_val_adj), "genes"], 30) & markers$diffexpressed != "NO", markers$genes, NA)
      # Add genes with largest log2fc
      top_log2fc <- c(rownames(tail(significant_genes, n=10)), rownames(head(significant_genes, n=10)))
      markers$top30 <- ifelse(is.na(markers$top30) & markers$genes %in% top_log2fc & markers$diffexpressed != "NO", markers$genes, markers$top30)
      
      
      # Volcano plot
      markers$p_val_adj[markers$p_val_adj == 0] <- 1e-300
      y_max <- max(-log10(markers$p_val_adj), na.rm = TRUE) + 10  # Get max y val with buffer of 10
      
      ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = top30)) +
        geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
        geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
        geom_point(size = 0.5) +
        coord_cartesian(ylim = c(0, y_max), xlim = c(-7, 7)) +
        scale_color_manual(values = c("blue", "grey", "red"),
                           labels = c("Downregulated", "Not significant", "Upregulated")) +
        scale_x_continuous(breaks = seq(-10, 10, 2)) +
        labs(color = 'Severe') +
        geom_text_repel(size=3, max.overlaps = Inf) +
        ggtitle(paste(subtype, curr_sample1, "vs", sample))
      ggsave(paste0(subtype_fixed, '_volcano.png'), width=7, height=5, dpi=400)
    }
  }
}

# Save object
saveRDS(s1.tcells, file=paste0(objDir, 'tcell_subsets/tcell_object_final.rds'))


##############################################################################################
# Gene set enrichment
##############################################################################################

# Proliferation (G2M) genes from Seurat
prolif_genes <- read.table("g2m_genes.txt", header = FALSE, sep = "\t")
prolif_genes <- prolif_genes[["V1"]]

# Compute scores
s1.tcells <- AddModuleScore(object = s1.tcells, features = list(prolif_genes), name = "prolif")


# Visualize
VlnPlot(s1.tcells, features = "prolif1", group.by = "orig.ident")

# Individual genes 
VlnPlot(s1.tcells, features = "GZMA", group.by = "orig.ident")
VlnPlot(s1.tcells, features = "NFATC1", group.by = "orig.ident")
VlnPlot(s1.tcells, features = "NFKB1", group.by = "orig.ident")
VlnPlot(s1.tcells, features = "IFNG", group.by = "orig.ident")

##############################################################################################
# Individual Gene Plots
##############################################################################################

figDir <- ""

# Genes of interest
gene_list <- c('CD4', 'NKG7', 'TOP2A', 'MKI67', 'IL23R', 'CENPF', 'IFNG', 'CD8A', 'CD8B', 'GZMA',
               'RPL37', 'CD3E', 'CD19')

# Check for a gene
rownames(s1.tcells)[grep('CD184', rownames(s1.tcells))] # To check for a marker

# Plot expression of each gene
for (gene in gene_list) {
  # Extract expression data from each timepoint
  expr_data <- FetchData(s1.tcells, vars = c(gene, "orig.ident"))  # SCT normalized data
  # Convert to factor
  expr_data$orig.ident <- factor(expr_data$orig.ident, levels = sort(unique(expr_data$orig.ident)))
  
  # Get cell counts 
  cell_counts <- expr_data %>%
    group_by(orig.ident) %>%
    summarise(n_cells = n(), .groups = "drop")
  
  # Run anova on expression
  anova_result <- aov(get(gene) ~ orig.ident, data = expr_data)
  # Check p-value
  #summary(anova_result) 
  
  # Plot
  ggplot(expr_data %>% mutate(orig.ident = factor(orig.ident, levels = unique(s1.tcells$orig.ident))),
         aes(x = orig.ident, y = get(gene), fill = orig.ident)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    #geom_jitter(width = .01, alpha = 0.7, size = 0.5) +
    geom_point(position = position_nudge(x = 0), alpha = 0.7, size = 1.5) +  # No horizontal shift
    stat_compare_means(method = "anova", label.y = max(expr_data[[gene]]) * 1.1) +  # ANOVA p-value
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste('Expression across samples: ', gene),
         y = "Expression (SCT norm)",
         x = "Timepoint") +
    expand_limits(y = -0.05) +
    geom_text(data = cell_counts, aes(x = orig.ident, y = max(expr_data[[gene]], na.rm = TRUE) * 1.05, label = n_cells), 
              size = 4, color='grey50') + 
    theme(legend.position = "none", text = element_text(size = 14), axis.title.x = element_text(margin = margin(t = 10)))
  ggsave(paste0(figDir, gene, '.png'), dpi=400, height=5, width=7)
}


########################################################################################################################
#                                       Pathway Enrichment with ClusterProfiler
########################################################################################################################

degs_path <- ""
out_path <- ""

# Set sample pairs -----------------------------------------

# PreTx v all
sample_list1 <- c('PreTx')
sample_list2 <- c('D28FNA', 'D7FNA', 'D10resect', 'D28FNA', 'D49FNA')

# Other pairs of interest
sample_list1 <- c('D10FNA')
sample_list2 <- c('D28FNA')

------------------------------------------------------------
resDir <- ""

for (sample1 in sample_list1) {
  for(sample2 in sample_list2) {
    # Set results directory
    subresDir <- paste0(sample1, '_', sample2)
    
    if (file.exists(subresDir)){
      setwd(file.path(resDir, subresDir))
    } else {
      dir.create(file.path(resDir, subresDir)) # Create if doesn't exist
      setwd(file.path(resDir, subresDir)) # Set as current wd
    }
    
    print(paste0('Running Pathway Analysis for ', sample1, ' vs ', sample2, '....'))

    # Load in df of DEGs results for curr pairwise grouop
    df <- read.table(paste0(degs_path, sample1, '_', sample2, '_significant_genes.tsv'), sep='\t')
    
    # Get only upregulated genes
    #upregulated <- df[df$avg_log2FC > 0, ]
    # Save
    #write.table(upregulated, file=paste0(out_path, '_upregulated_only.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
    # Override
    #df <- upregulated
    
    # Initialize cutoffs and output file
    log2FC_lower <- -0.05
    log2FC_upper <- 0.05
    pval_lim <- 0.001
    
    sink('deg_filtering.txt')
    
    # Set log2FC cutoff until < 500 DEGs in list
    while (nrow(df) > 500) {

      # Print initial DEG amount and cutoffs
      cat(paste0('Initial DEG count:', nrow(df),'\n'))
      cat(paste0('Filtering at (', log2FC_lower, ', ', log2FC_upper, ') log2FC cutoffs... ','\n'))
      
      # Filter
      df <- df %>%
        filter(avg_log2FC < log2FC_lower | avg_log2FC > log2FC_upper) %>%
        filter(p_val_adj < pval_lim)
      
      # Print resulting DEG count
      cat(paste0('DEG count: ', nrow(df),'\n'))
      
      # Update cutoffs
      log2FC_lower <- -(abs(log2FC_lower) + 0.1)
      log2FC_upper <- log2FC_upper + 0.1
    }
    # End output redirection
    sink()
    
    sink('parameters.txt')
    cat(paste0(sample1, ' v ', sample2, ' : ', nrow(df), ' DEGs being used for ORA','\n'))
    sink()


    ### Format input df
    # Annotate according to differential expression
    df <- df %>% mutate(diffexpressed = case_when(
      avg_log2FC > 0 & p_val_adj < 0.05 ~ 'Upregulated',
      avg_log2FC < 0 & p_val_adj < 0.05 ~ 'Downregulated',
    ))

    
    # Split into 2 dfs: upregulated and downregulated
    deg_results_list <- split(df, df$diffexpressed)
    
    # Get deg names for up and downregulated dfs
    genes_in_data_up <-rownames(deg_results_list$Upregulated)
    genes_in_data_down <- rownames(deg_results_list$Downregulated)
    
    
    # Convert query gene symbols to ENTREZ IDs
    library("org.Cf.eg.db")
    entrez_genes_up <- bitr(genes_in_data_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Cf.eg.db")
    entrez_genes_down <- bitr(genes_in_data_down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Cf.eg.db")
    
    # Convert background gene symbols to ENTREZ IDs
    universe_genes <- bitr(rownames(s1.tcells), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Cf.eg.db")

    ################# Run ClusterProfiler #################
    # Settings
    padj_cutoff <- 0.25 # p-adjusted threshold, used to filter out pathways
    GS_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways
    
    # Write params to outfile
    sink("parameters.txt", append = TRUE)
    cat(paste0('padj cutoff: ', padj_cutoff, '\n'))
    cat(paste0('GS cutoff: ', GS_cutoff, '\n'))
    sink()
    
    # Run enrichKEGG
    # Downregulated
    results_kegg_down <- enrichKEGG(gene = entrez_genes_down$ENTREZID, # applies the enricher function to each gene 
                                              universe = universe_genes$ENTREZID, # set background genes
                                              minGSSize = GS_cutoff, # min genes required for a pathway
                                              pvalueCutoff = padj_cutoff, # pval cutoff for enriched pathways
                                              organism = 'cfa') # cfa for canine

    # Upregulated
    results_kegg_up <-  tryCatch({enrichKEGG(gene = entrez_genes_up$ENTREZID, 
                                    universe = universe_genes$ENTREZID, 
                                    minGSSize = GS_cutoff, 
                                    pvalueCutoff = padj_cutoff, 
                                    organism = 'cfa') 
    }, error = function(e) {
      message("KEGG Upregulated enrichment failed: ", e$message)
      return(NULL)
    })
      
    # Plot results -- downreg
    if (!is.null(results_kegg_down)) {
      # Barplot
      p <- barplot(results_kegg_down, showCategory = 15) + ggtitle(paste0(sample1, ' v ', sample2, ' (KEGG Downregulated)'))
      ggsave("downreg_kegg_barplot.png", plot = p, width = 8, height = 7, dpi=400)
      
      # Dotplot
      p <- dotplot(results_kegg_down, showCategory = 15) + ggtitle(paste0(sample1, ' v ', sample2, ' (KEGG Downregulated)'))
      ggsave("downreg_kegg_dotplot.png", plot = p, width = 8, height = 7, dpi=400)
    }
    
    # Plot results -- upreg
    if (!is.null(results_kegg_up)) {
      # Barplot
      p <- barplot(results_kegg_up, showCategory = 15) + ggtitle(paste0(sample1, ' v ', sample2, ' (KEGG Upregulated)'))
      ggsave("upreg_kegg_barplot.png", plot = p, width = 8, height = 7, dpi=400)
      
      # Dot plot
      p <- dotplot(results_kegg_up, showCategory = 15) + ggtitle(paste0(sample1, ' v ', sample2, ' (KEGG Upregulated)'))
      ggsave("upreg_kegg_dotplot.png", plot = p, width = 8, height = 7, dpi=400)
    }
    

    # Run enrichGO 
    # Downregulated
    results_go_down <- enrichGO(gene = entrez_genes_down$ENTREZID,
                                           universe = universe_genes$ENTREZID, 
                                           ont='ALL',
                                           minGSSize = GS_cutoff, 
                                           pvalueCutoff = padj_cutoff, 
                                           OrgDb = 'org.Cf.eg.db')
    # Upregulated
    results_go_up <-  tryCatch({enrichGO(gene = entrez_genes_up$ENTREZID,
                                universe = universe_genes$ENTREZID, 
                                ont='ALL',
                                minGSSize = GS_cutoff, 
                                pvalueCutoff = padj_cutoff, 
                                OrgDb = 'org.Cf.eg.db')
    }, error = function(e) {
      message("GO Upregulated enrichment failed: ", e$message)
      return(NULL)
    })
    
  
    # Plot results
    if (!is.null(results_go_down)) {
      # Barplot
      p <- barplot(results_go_down, showCategory = 12, font.size=13) + ggtitle(paste0(sample1, ' v ', sample2, ' (GO Downregulated)'))
      ggsave("downreg_go_barplot.png", plot = p, width = 8, height = 8, dpi=400)
      
      # Dotplot
      p <- dotplot(results_go_down, showCategory = 12, font.size=13) + ggtitle(paste0(sample1, ' v ', sample2, ' (GO Downregulated)'))
      ggsave("downreg_go_dotplot.png", plot = p, width = 9, height = 8, dpi=400)
    }
    
    # Plot results
    if (!is.null(results_go_up)) {
      # Barplot
      p <- barplot(results_go_up, showCategory = 20) + ggtitle(paste0(sample1, ' v ', sample2, ' (GO Upregulated)'))
      ggsave("upreg_go_barplot.png", plot = p, width = 8, height = 7, dpi=400)
      
      # Dotplot
      p <- dotplot(results_go_up, showCategory = 20) + ggtitle(paste0(sample1, ' v ', sample2, ' (GO Upregulated)'))
      ggsave("upreg_go_dotplot.png", plot = p, width = 8, height = 7, dpi=400)
    }
    
    # Save results
    # Format into condensed df
    df_kegg_down <- if (!is.null(results_kegg_down)) as.data.frame(results_kegg_down) else NULL
    df_kegg_up <- if (!is.null(results_kegg_up)) as.data.frame(results_kegg_up) else NULL
    df_kegg_na_uni <- if (!is.null(results_kegg_down_na_uni)) as.data.frame(results_kegg_down_na_uni) else NULL
    df_go_down <- if (!is.null(results_go_down)) as.data.frame(results_go_down) else NULL
    df_go_up <- if (!is.null(results_go_up)) as.data.frame(results_go_up) else NULL
    
    # Save
    if (!dir.exists("results")) {
      dir.create("results")
    }
    write.table(df_kegg_down, 'results/kegg_results_downreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_kegg_up, 'results/kegg_results_upreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_go_down, 'results/go_results_downreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_go_up, 'results/go_results_upreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_go_up, 'results/go_results_upreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    
  }
} 


########################################################################################################################
#                                       Updating main object with filtered T cells
########################################################################################################################
figDir = ""

# Checking if populations match
original_tcells <- s1[,s1$cell_type == 'T cells']
all(Cells(original_tcells) %in% Cells(s1.tcells))

# Extract subtype annotations
subtype_annotations <- s1.tcells$cell_type_final
names(subtype_annotations) <- Cells(s1.tcells)

# Update cell type final annotations
s1$cell_type_final <- s1$cell_type  
s1$cell_type_final[names(subtype_annotations)] <- as.character(subtype_annotations)
unique(s1$cell_type_final) # Check

# Plot 
DimPlot(s1, group.by='cell_type_final')
ggsave(paste0(figDir, 'umap_fine_annotations.png'), dpi=400, height=5, width=7)

# Save
saveRDS(s1, paste0(objDir, 'processed_object_final.rds')) # with final T cell annotations





