#############################################################################################
############################################################################################
### Comparison of ORv-01 MZL Longitudinal Samples to Adjacent Healthy Spleen / WBCs
### Samples: ORv-01 longitudinal scRNAseq (n=7), DHSA_2101 (spleen), WBCs (n=3)
### Author: Grace Walker
### Date: July 10, 2025
#############################################################################################
#############################################################################################

# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(ggplot2) # for plotting
library(tidyverse) # for data manipulation
library(Seurat) # for Seurat processing
library(patchwork) # for plotting
library(scales)
library(clustree) # for clustering
library(celldex) # for reference datasets
library(SingleR) # for automated annotation
library(scran) # for normalization
library(EnhancedVolcano) # for dge visualization
library(ggpubr) # for boxplots

# Params
options(future.globals.maxSize = 4 * 1024^3)  # Set to 4 GB
gc() # free up memory

# Set dirs
objDir <- ""
spleenDir <- "" 
wbcDir <- ""

# -------------------------------------------------------------------------------------------

# Read in data
s1 <- readRDS(paste0(objDir, 'processed_object_final.rds'))
s2 <- readRDS(paste0(objDir, 'spleen/spleen_processed.rds'))
s3 <- readRDS(paste0(objDir, 'wbc/wbc_merged_processed.rds'))

#############################################################################################
# Prep data
#############################################################################################

### Spleen dataset
spleen <- Read10X(data.dir = spleenDir)
s2 <- CreateSeuratObject(counts = spleen, project = 'spleen', min.cells = 3, min.features = 200)


### WBC datasets
wbc_datasets <- c('cellranger_Modiano_061_K9_WBC_Dog_1_YP_GEX',
                  'cellranger_Modiano_064_K9_WBC_141517_GEX_FL',
                  'cellranger_Modiano_064_K9_WBC_1832517_GEX_FL')
objects <- c()
for (i in 1:length(wbc_datasets)) {
  umgc_id <- wbc_datasets[[i]]
  full_dir <- paste0(wbcDir, umgc_id, '/outs/filtered_feature_bc_matrix/')
  counts <- Read10X(data.dir = full_dir)
  sample_id <- paste0('WBC0',i)
  s3 <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200)
  objects[[sample_id]] <- s3
}

########################################################################
#           Preprocessing & QC of spleen and WBC datasets
########################################################################
### Spleen 
res_dir <- ""

# Calculate percent mt
s2[["percent.mt"]] <- PercentageFeatureSet(s2, pattern = "^MT-")

# QC plots
VlnPlot(s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), layer="counts", ncol = 3)
ggsave('QC/violin_plots.png', dpi=400, width=9, height=6)

# Scatter plots
plot1 <- FeatureScatter(s2, feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle('nCount_RNA x percent.mt')
plot2 <- FeatureScatter(s2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle('nCount_RNA x nFeature_RNA')
plot3 <- FeatureScatter(s2, feature1 = "percent.mt", feature2 ="nFeature_RNA") + ggtitle('percent.mt x nFeature_RNA')
plot1 + plot2 + plot3
ggsave('QC/scatter_plots.png', dpi=400, width=15, height=6)

# Density plots
# Thresholds
mt = 10
feature_min = 300
count_min = 400
# Plot
p1 <- s2@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() + 
  scale_x_log10() + 
  geom_vline(xintercept = feature_min,color="red",linetype="dotted") +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none") +
  ggtitle("nFeature")
p2 <- s2@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() + 
  scale_x_log10() + 
  geom_vline(xintercept = count_min,color="red",linetype="dotted") +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none") +
  ggtitle("nCount")
p3 <- s2@meta.data %>% 
  ggplot(aes(x=percent.mt,fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10(labels = label_comma()) + 
  geom_vline(xintercept = mt,color="red",linetype="dotted") +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("MT Count")
p1 + p2 + p3
ggsave('QC/density_plots.png', dpi=400, width=15, height=6)

# Filter cells
prefilter_count = length(Cells(s2))
s2@meta.data$keep <- with(s2@meta.data, ifelse(nFeature_RNA > feature_min & nCount_RNA > count_min & percent.mt < mt, TRUE, FALSE))
s2 <- subset(s2, subset = keep == TRUE)
# Count
postfilter_count = length(Cells(s2))
total = prefilter_count - postfilter_count
print(paste("Cells removed:", total))

## Save the filtered object
saveRDS(s2, paste0(objDir, "/spleen/spleen_processed.rds"))


### WBCs
res_dir <- ""
for (i in 1:length(objects)) {
  # Get object
  s2 <- objects[[i]]
  sample_id <- names(objects)[i]
  # Set curr res dir
  curr_dir <- paste0(res_dir,sample_id, '/')
  # Create subdir
  dir.create(paste0(curr_dir, 'QC'))
             
  
  # Calculate percent mt
  s2[["percent.mt"]] <- PercentageFeatureSet(s2, pattern = "^MT-")
  
  # QC plots
  VlnPlot(s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), layer="counts", ncol = 3)
  ggsave(paste0(curr_dir,'QC/violin_plots.png'), dpi=400, width=9, height=6)
  
  # Scatter plots
  plot1 <- FeatureScatter(s2, feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle('nCount_RNA x percent.mt')
  plot2 <- FeatureScatter(s2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle('nCount_RNA x nFeature_RNA')
  plot3 <- FeatureScatter(s2, feature1 = "percent.mt", feature2 ="nFeature_RNA") + ggtitle('percent.mt x nFeature_RNA')
  plot1 + plot2 + plot3
  ggsave(paste0(curr_dir,'QC/scatter_plots.png'), dpi=400, width=15, height=6)
  
  # Density plots
  # Thresholds
  mt = 6
  feature_min = 350
  count_min = 400
  # Plot
  p1 <- s2@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() + 
    scale_x_log10() + 
    geom_vline(xintercept = feature_min,color="red",linetype="dotted") +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(legend.position = "none") +
    ggtitle("nFeature")
  p2 <- s2@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() + 
    scale_x_log10() + 
    geom_vline(xintercept = count_min,color="red",linetype="dotted") +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(legend.position = "none") +
    ggtitle("nCount")
  p3 <- s2@meta.data %>% 
    ggplot(aes(x=percent.mt,fill=orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10(labels = label_comma()) + 
    geom_vline(xintercept = mt,color="red",linetype="dotted") +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("MT Count")
  p1 + p2 + p3
  ggsave(paste0(curr_dir,'QC/density_plots.png'), dpi=400, width=15, height=6)
}

# Filter each dataset
s2 <- objects[[1]]
s2 <- objects[[2]]
s2 <- objects[[3]]

prefilter_count = length(Cells(s2))
s2@meta.data$keep <- with(s2@meta.data, ifelse(nFeature_RNA > feature_min & nCount_RNA > count_min & percent.mt < mt, TRUE, FALSE))
s2 <- subset(s2, subset = keep == TRUE)
# Count
postfilter_count = length(Cells(s2))
total = prefilter_count - postfilter_count
print(paste("Cells removed:", total))

# Save
saveRDS(s2, paste0(objDir, "/wbc/WBC01_processed.rds"))
saveRDS(s2, paste0(objDir, "/wbc/WBC02_processed.rds"))
saveRDS(s2, paste0(objDir, "/wbc/WBC03_processed.rds"))


# Read back in and merge
wbc1 <- readRDS(paste0(objDir, 'wbc/WBC01_processed.rds'))
wbc2 <- readRDS(paste0(objDir, 'wbc/WBC02_processed.rds'))
wbc3 <- readRDS(paste0(objDir, 'wbc/WBC03_processed.rds'))

# Merge
s3 <- merge(wbc1, y = c(wbc2, wbc3), add.cell.ids = c("wbc1", "wbc2", "wbc3"), project = "WBCs")
# Save
saveRDS(s3, paste0(objDir, 'wbc/wbc_merged_processed.rds'))

########################################################################
#           Normalization & Clustering
########################################################################
subres_dir <- 'WBC03'
s2 <- readRDS(paste0(objDir, 'wbc/', subres_dir, '_processed.rds'))

dir.create(paste0(res_dir, subres_dir, '/clustering'))
setwd(paste0(res_dir, subres_dir, '/clustering'))

# SCT normalize
s2 <- SCTransform(s2, verbose = TRUE)
# PCA
s2 <- RunPCA(s2, assay="SCT") 
# Elbow plot
ElbowPlot(s2, ndims = 40) # Elbow plot
ggsave('elbow_plot.png', width = 8, height = 5, dpi=400)
# Neighbors
s2 <- FindNeighbors(s2, dims = 1:30)
# Clustering
resolutions = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)

for (res in resolutions){
  # Cluster
  s2 <- FindClusters(s2, resolution=res, cluster.name=paste0('res',res)) # Setting smallest resolution
  # Plot UMAP
  #s2 <- RunUMAP(s2, dims = 1:30) # only do this once
  DimPlot(s2, reduction = "umap", group.by = paste0('res',res), alpha = 0.7)
  ggsave(paste0('umaps/clusters_', res, '.png'), width = 6, height = 5, dpi=400)
}


# Clustree to choose resolution
clusterings <- s2@meta.data
clustree(clusterings, prefix = "res")
ggsave(paste('clustree_diagram.png', sep=""), width = 7, height = 9, dpi=400)

# Save progress
saveRDS(s2, paste0(objDir, 'wbc/', subres_dir, '_processed.rds'))

########################################################################
#           SingleR Automated Annotation
########################################################################
res_dir < "" # for spleen
res_dir <- paste0(resDir, subres_dir, "/singler/") # for wbcs

# Import reference datasets
ref1 <- celldex::HumanPrimaryCellAtlasData() # Human primary scRNAseq atlas 
ref2 <- celldex::BlueprintEncodeData() # Human bulk RNAseq
ref3 <- celldex::ImmGenData() # Mouse immune bulk RNAseq

# Get normalized counts
# Convert seurat object to sce object
DefaultAssay(s2) <- "SCT"
integrated_sce <- as.SingleCellExperiment(s2,assay = "SCT")
# Ensure log normalized counts
integrated_sce <- logNormCounts(integrated_sce)

## Convert mouse gene names to uppercase
rownames(ref3) <- toupper(rownames(ref3))

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
write.table(hpca_results_main, paste0(res_dir, "SingleR_Results_HPCA_Main.tsv"), sep = "\t", quote = FALSE)
write.table(hpca_results_fine, paste0(res_dir, "SingleR_Results_HPCA_Fine.tsv"), sep = "\t", quote = FALSE)
write.table(blue_results, paste0(res_dir, "SingleR_Results_BlueEncode.tsv"), sep = "\t", quote = FALSE)
write.table(immgen_results, paste0(res_dir, "SingleR_Results_ImmGen.tsv"), sep = "\t", quote = FALSE)


### Plotting Results ###
# Plot score heatmaps
png(paste0(res_dir,"SingleR_Heatmap_CellAssignmentScores_HPCA_Main.png"),height=2600,width=4000,res=400)
plotScoreHeatmap(hpca_results_main)
dev.off()

png(paste0(res_dir,"SingleR_Heatmap_CellAssignmentScores_HPCA_Fine.png"),height=2600,width=5000, res=400)
plotScoreHeatmap(hpca_results_fine)
dev.off()

png(paste0(res_dir,"SingleR_Heatmap_CellAssignmentScores_BlueEncode.png"),height=2600,width=4000, res=400)
plotScoreHeatmap(blue_results)
dev.off()

png(paste0(res_dir,"SingleR_Heatmap_CellAssignmentScores_ImmGen.png"),height=2600,width=4000, res=400)
plotScoreHeatmap(immgen_results)
dev.off()


### Merge results with object ###
rownames(hpca_results_main)[1:5] # Checking cell IDs
s2 <- AddMetaData(s2, hpca_results_main$labels, col.name = 'hpca_main_raw')
s2 <- AddMetaData(s2, hpca_results_main$pruned.labels, col.name = 'hpca_main_pruned')

s2 <- AddMetaData(s2, hpca_results_fine$labels, col.name = 'hpca_fine_raw')
s2 <- AddMetaData(s2, hpca_results_fine$pruned.labels, col.name = 'hpca_fine_pruned')

s2 <- AddMetaData(s2, blue_results$labels, col.name = 'blue_raw')
s2 <- AddMetaData(s2, blue_results$pruned.labels, col.name = 'blue_pruned')

s2 <- AddMetaData(s2, immgen_results$labels, col.name = 'immgen_raw')
s2 <- AddMetaData(s2, immgen_results$pruned.labels, col.name = 'immgen_pruned')


# Plot - HPCA main
Idents(s2) <- s2@meta.data$hpca_main_pruned
DimPlot(s2, label = F , repel = T, label.size = 3) + ggtitle('HPCA Main Annotations')
ggsave(paste0(res_dir, 'hpca_main_umap.png'), height=6, width=9, dpi=400) 

# Plot - HPCA fine
Idents(s2) <- s2@meta.data$hpca_fine_pruned
DimPlot(s2, label = F , repel = T, label.size = 3) + ggtitle('HPCA Fine Annotations')
ggsave(paste0(res_dir, 'hpca_fine_umap.png'), height=6, width=20, dpi=400) 

# B cells
bcells <- c("B_cell:Memory", "B_cell:Naive", "B_cell", "B_cell:Germinal_center")
cleaned$bcells <- ifelse(cleaned$hpca_fine_pruned %in% bcells, 
                         cleaned$hpca_fine_pruned, 
                         "Other")
DimPlot(cleaned, group.by = "bcells") + ggtitle('B Cells')
ggsave(paste0(fig_path, 'umap_bcells.png'), height=6, width=8, dpi=400)

# Immune cells
immune <- c("B_cell:Memory", "B_cell:Naive", "B_cell", "B_cell:Germinal_center", "T_cell:CD8+", "T_cell:CD4+", "T_cell:CD4+_Naive", "T_cell:Treg:Naive", "NK_cell")
integrated$immune <- ifelse(integrated$hpca_fine_pruned %in% immune, 
                            integrated$hpca_fine_pruned, 
                            "Other")
DimPlot(integrated, reduction='umap.harmony', group.by = "immune") + ggtitle('Immune Cells')
ggsave(paste0(fig_path, 'umap_immunecells.png'), height=6, width=8, dpi=400)


# Plot - BlueEncode
Idents(s2) <- "blue_pruned"
DimPlot(s2, reduction='umap') + ggtitle('BlueEncode Annotations')
ggsave(paste0(res_dir, 'umap_blueencode.png'), height=7, width=20, dpi=400)

# Plot - ImmGen
Idents(s2) <- "immgen_pruned"
DimPlot(s2, reduction='umap') + ggtitle('ImmGen Annotations')
ggsave(paste0(res_dir, 'umap_immgen.png'), height=5, width=8, dpi=400)

# Save progress
saveRDS(s2, paste0(objDir, "/spleen/spleen_processed.rds"))
saveRDS(s2, paste0(objDir, 'wbc/', subres_dir, '_processed.rds'))


########################################################################
#           Cluster DEGs
########################################################################
res_dir <- "" # for spleen
res_dir <- paste0(resDir, subres_dir, "/cluster_markers/") # for wbcs

resolutions <- c('res0.1', 'res0.09') # spleen
resolutions <- c('res0.04') # wbc1
resolutions <- c('res0.04') # wbc2
resolutions <- c('res0.02') # wbc3


# Plotting top DEGs per cluster
for (res in resolutions) {
  # Set Idents
  Idents(s2) <- s2@meta.data[[res]]
  # Prep SCT
  s2 <- PrepSCTFindMarkers(s2,verbose=T)
  # Get markers
  markers <- FindAllMarkers(s2, only.pos = TRUE, recorrect_umi=FALSE)
  
  # Filter for significant genes (p_val_adj < 0.05)
  significant <- markers %>%
    filter(p_val_adj<0.05)
  
  # Get top 10 for each cluster
  top10 <- significant %>%
    group_by(cluster) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # Save
  write.table(significant, file=paste0(res_dir, 'cluster_markers_all_', res, ".tsv"), quote=FALSE, row.names=TRUE)
  write.table(top10, file=paste0(res_dir, 'cluster_markers_top10_', res, ".tsv"), quote=FALSE, row.names=TRUE)
  
  # Heatmap plot
  png(filename = paste0(res_dir, 'marker_heatmap_', res, '.png'), width = 8, height = 6, units = "in", res = 400)
  p <- DoHeatmap(s2, features = top10$gene, draw.lines=TRUE, label=FALSE)
  p + scale_fill_gradient(colors = c("cornflowerblue", "ivory1", "darkorange"))
  dev.off()
}

# Save 
saveRDS(s2, paste0(objDir, 'wbc/', subres_dir, '_processed.rds'))

########################################################################
#           Known markers
########################################################################
res_dir <- "" # for spleen
res_dir <- paste0(resDir, subres_dir, "/cluster_markers/markers/") # for wbcs

clusters <- 'res0.1'
clusters <- 'res0.04'


# Immune markers
bcell <- c('CD79A', 'CD79B', 'CD19', 'CD22')
tcell <- c('CD3E', 'CD4', 'CD8A', 'CD8B')
macro <- c('CD68', 'CD14', 'CD163', 'CD11B')
plasma <- c('CD138', 'CD38', 'CD27', 'IGA')

# Plot expression
FeaturePlot(s2, features=bcell)
ggsave(paste(res_dir, 'bcell_markers.png', sep=""), width = 8, height = 5, dpi=400)
FeaturePlot(s2, features=tcell)
ggsave(paste(res_dir, 'tcell_markers.png', sep=""), width = 8, height = 5, dpi=400)
FeaturePlot(s2, features=macro)
ggsave(paste(res_dir, 'macro_markers.png', sep=""), width = 8, height = 5, dpi=400)
FeaturePlot(s2, features=mzb)
ggsave(paste(res_dir, 'mzb_marker.png', sep=""), width = 8, height = 5, dpi=400)
FeaturePlot(s2, features=plasma)
ggsave(paste(res_dir, 'plasma_marker.png', sep=""), width = 8, height = 5, dpi=400)


########################################################################
#           Differential Expression Analysis of B cells
########################################################################
# Create col for broad cell type labels
s1$cell_type_broad <- dplyr::case_when(
  s1$cell_type %in% c("+Ribo/-Prolif B cells", "-Ribo B cells", "+Prolif B cells") ~ "B cells",
  s1$cell_type == "T cells" ~ "T cells",
  s1$cell_type == "Immune cells" ~ "Immune cells",
  TRUE ~ "Other")

# Subsetting objects
# Subset ORv01 obj to 500 PreTx and D10res B cells
pretx <- subset(s1, subset=orig.ident=='PreTx')
d10res <- subset(s1, subset=orig.ident=='D10resect')

# Subset to 500 pretx
pretx.bcells <- subset(pretx, subset=cell_type_broad=="B cells")
set.seed(67)  # For reproducibility
sampled_cells <- unlist(lapply(unique(pretx.bcells$cell_type), function(subtype) {
  cells <- Cells(pretx.bcells)[pretx.bcells$cell_type == subtype]
  sample(cells, size = min(200, length(cells)))}))
ORv01_PreTx <- subset(pretx.bcells, cells = sampled_cells)

# Subset to 500 d10res
d10res.bcells <- subset(d10res, subset=cell_type_broad=="B cells")
set.seed(123)  # For reproducibility
sampled_cells <- unlist(lapply(unique(d10res.bcells$cell_type), function(subtype) {
  cells <- Cells(d10res.bcells)[d10res.bcells$cell_type == subtype]
  sample(cells, size = min(200, length(cells)))}))
ORv01_D10res <- subset(d10res.bcells, cells = sampled_cells)

# Subset spleen dataset to only B cells
DHSA_104 <- subset(s2, subset=res0.1==1)

# Subset WBCs to only 500 B cells
Idents(wbc1) <- 'res0.04'
all.cells <- WhichCells(wbc1, idents = 2)  # Assuming Idents(wbc1) is set to res0.04
# Grab 500 
set.seed(123)  # for reproducibility
sampled.cells <- sample(all.cells, 600)
# Subset
WBC_01 <- subset(wbc1, cells = sampled.cells)


# Check length of all objects
print(length(Cells(WBC_01)))
print(length(Cells(DHSA_104)))
print(length(Cells(ORv01_PreTx)))
print(length(Cells(ORv01_D10res)))


##################### Grab raw counts of B cells from each
DHSA_104$condition <- "DHSA_104"
ORv01_D10res$condition <- "ORv01_D10res"
ORv01_PreTx$condition <- "ORv01_PreTx"
WBC_01$condition <- "WBC_01"


# Merge datasets
merged <- merge(
  x = DHSA_104,
  y = list(ORv01_PreTx, ORv01_D10res, WBC_01),
  add.cell.ids = c("DHSA", "PreTx", "D10res", "WBC_01"))

# Fix NAs for WBC and spleen
merged$cell_type[is.na(merged$cell_type) & merged$condition == 'WBC_01'] <- "WBC_01 B cells" 
merged$cell_type[is.na(merged$cell_type) & merged$condition == 'DHSA_104'] <- "DHSA_104 B cells"

# Clean up object
DefaultAssay(merged) <- "RNA"
merged[["SCT"]] <- NULL
merged <- JoinLayers(merged)


# Renormalize
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)


# Varibale Features
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

# Plot top 10
top10 <- head(VariableFeatures(merged), 10)
plot1 <- VariableFeaturePlot(merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scale
all.genes <- rownames(merged)
merged <- ScaleData(merged, features = all.genes)

# Reprocess
merged <- RunPCA(merged)

# Visualize PCs
VizDimLoadings(merged, dims = 1:2, reduction = "pca")

ElbowPlot(merged)

merged <- FindNeighbors(merged, dims = 1:20)

merged <- RunUMAP(merged, dims = 1:20) # only do this once

# Plot
DimPlot(merged, reduction = "umap", group.by = c('condition', 'cell_type'), alpha = 0.8)
ggsave(paste0('merged_umaps/bcell_umap.png'), width = 12, height = 5, dpi=400)


# Save merged object
saveRDS(merged, paste0(objDir, "/spleen/spleen_merged_processed.rds"))


##################### Pairwise comparisons 
res_dir <- ""
samples <- unique(merged$condition)

while (length(samples) > 0) {
  
  # Get current sample (entire ID)
  curr_sample <- samples[1]
  # Remove from list
  samples <- samples[-1]
  
  # Loop over each sample vs spleen dataset
  for (sample in samples) {
    
    # Make results dir
    resDir <- res_dir
    subresDir <- paste0(curr_sample, '_', sample)
    
    if (file.exists(subresDir)){
      setwd(file.path(resDir, subresDir))
    } else {
      dir.create(file.path(resDir, subresDir)) # Create if doesn't exist
      setwd(file.path(resDir, subresDir)) # Set as current wd
    }
    
    print(paste('Getting markers for', curr_sample, 'vs', sample))
    
    # Get markers
    markers <- FindMarkers(object = merged, 
                           ident.1 = curr_sample, 
                           ident.2 = sample,
                           recorrect_umi=FALSE)
    
    # Sort by log2FC
    markers_ordered <- markers[order(-markers$avg_log2FC), ]
    significant_genes <- markers_ordered %>%
      filter(p_val_adj < 0.05)
    
    # Save filtered gene list
    write.table(significant_genes, file=paste0(curr_sample, '_', sample, '_sig_genes.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
    # Save unfiltered gene list
    write.table(markers, file=paste0(curr_sample, '_', sample, '_unfiltered_genes.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
    
    # Volcano plot
    auto_genes <- head(rownames(markers)[abs(markers$avg_log2FC) > 0.6 & markers$p_val_adj < 0.05], n=15)
    custom_genes <- c('KLF2') # Edit if genes of interest not getting labelled
    genes_to_plot <- c(auto_genes,custom_genes)

    # Plot
    EnhancedVolcano(
      markers,
      lab = rownames(markers),
      x = 'avg_log2FC',
      y = 'p_val_adj',
      selectLab = genes_to_plot,
      #xlim = c(-5, 6.5),
      #ylim = c(0, 40),
      pCutoff = 0.01,
      FCcutoff = 0.6,
      pointSize = 3,
      labSize = 4,
      axisLabSize = 13,
      #colAlpha = 1,
      legendLabels = c('NS', 'log_2FC', 'pval_adj', 'pval_adj and log2FC'),
      #legendPosition = 'right',
      #legendLabSize = 15,
      title = paste0(curr_sample, "vs", sample, 'B cells'),
      subtitle='',
      caption = paste0('n = ', nrow(markers), ' genes'),
      captionLabSize = 12,
      max.overlaps = Inf,
      boxedLabels = TRUE,
      drawConnectors = TRUE
    ) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0(res_dir, curr_sample, '_', sample, '_bcell_volcano.png'), width=7, height=7, dpi=400)
  }
}


#####################  Plot expression of specific genes
res_dir <- ""

genes <- c('NFKB1', 'NFKBIA', 'KLF2', 'MAP3K1', 'RIPOR2', 'BACH2')

# Boxplots
for (gene in genes) {
  # Extract expression data from each timepoint
  expr_data <- FetchData(merged, vars = c(gene, "condition"), layer='scale.data')  # Raw data
  # Convert to factor
  expr_data$condition <- factor(expr_data$condition, levels = sort(unique(expr_data$condition)))
  
  # Get cell counts 
  cell_counts <- expr_data %>%
    group_by(condition) %>%
    summarise(n_cells = n(), .groups = "drop")
  
  # Run anova on expression
  anova_result <- aov(get(gene) ~ condition, data = expr_data)
  
  # Plot
  ggplot(expr_data %>% mutate(condition = factor(condition, levels = unique(merged$condition))),
         aes(x = condition, y = get(gene), fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    #geom_jitter(width = .01, alpha = 0.7, size = 0.5) +
    geom_point(position = position_nudge(x = 0), alpha = 0.7, size = 1.5) +  # No horizontal shift
    stat_compare_means(method = "anova", label.y = max(expr_data[[gene]]) * 1.1) +  # ANOVA p-value
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = gene, y = "Expression", x = "Dataset") +
    expand_limits(y = -0.05) +
    geom_text(data = cell_counts, aes(x = condition, y = max(expr_data[[gene]], na.rm = TRUE) * 1.05, label = n_cells), 
              size = 4, color='grey50') + 
    theme(legend.position = "none", text = element_text(size = 14), axis.title.x = element_text(margin = margin(t = 10)))
  ggsave(paste0(res_dir, gene, '.png'), dpi=400, height=5, width=7)
}


#############################################################################################
# Comparing each PreTx B-cell subpopulation to WBCs
### Date: February 27, 2026
############################################################################################

# -------------------------------------------------------------------------------------------
# Function to sample 500 cells from a cell type
sample_cells <- function(group_name) {
  cells <- WhichCells(pretx, expression = cell_type == group_name)
  sample(cells, 600)
}
# -------------------------------------------------------------------------------------------

# Read in objects
s1 <- readRDS(paste0(objDir,'processed_object_final.rds'))
wbc1 <- readRDS(paste0(objDir, 'wbcs/WBC01_processed.rds')) # Cleanest WBC dataset

# Subset to prex only
pretx <- subset(s1, subset=orig.ident=='PreTx')

# Define group names
groups <- c("+Ribo/-Prolif B cells", "-Ribo B cells", "+Prolif B cells")

# Sample random subset of each B cell subpopulation
cells_prolif <- sample_cells("+Prolif B cells")
cells_normal <- sample_cells("+Ribo/-Prolif B cells")
cells_noribo <- sample_cells("-Ribo B cells")

PreTx_Prolif <- subset(pretx, cells = cells_prolif)
PreTx_Normal <- subset(pretx, cells = cells_normal)
PreTx_NoRibo <- subset(pretx, cells = cells_noribo)

# Subset WBCs to only 500 B cells (called by Immgen)
Idents(wbc1) <- 'immgen_pruned'
all.cells <- WhichCells(wbc1, idents = 'B cells')  # Assuming Idents(wbc1) is set to res0.04
# Grab 500 
sampled.cells <- sample(all.cells, 600)
# Subset
WBC_01 <- subset(wbc1, cells = sampled.cells)

# Check length of all objects
print(length(Cells(WBC_01)))
print(length(Cells(PreTx_Prolif)))
print(length(Cells(PreTx_Normal)))
print(length(Cells(PreTx_NoRibo)))


# Grab raw counts of B cells from each
WBC_01$condition <- "WBC_01"
PreTx_Prolif$condition <- "+Prolif B cells"
PreTx_Normal$condition <- "+Ribo_-Prolif B cells"
PreTx_NoRibo$condition <- "-Ribo B cells"

# Merge datasets
merged <- merge(
  x = WBC_01,
  y = list(PreTx_Prolif, PreTx_Normal, PreTx_NoRibo),
  add.cell.ids = c("WBC_01", "PreTx_Prolif", "PreTx_Normal", "PreTx_NoRibo"))

# Fix NAs for WBC
merged$cell_type[is.na(merged$cell_type) & merged$condition == 'WBC_01'] <- "WBC_01_Bcells" 

# Clean up object
DefaultAssay(merged) <- "RNA"
merged[["SCT"]] <- NULL
merged <- JoinLayers(merged)

# Renormalize
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)

# Varibale Features
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

# Scale
all.genes <- rownames(merged)
merged <- ScaleData(merged, features = all.genes)

# Reprocess
merged <- RunPCA(merged)

# Visualize PCs
VizDimLoadings(merged, dims = 1:2, reduction = "pca")

# Elbow plot
ElbowPlot(merged)

# Neighbors
merged <- FindNeighbors(merged, dims = 1:15)

# UMAP
merged <- RunUMAP(merged, dims = 1:15) # only do this once

# Plot
DimPlot(merged, reduction = "umap", group.by = c('condition', 'cell_type'), alpha = 0.8)
ggsave(paste0(resDir,'bcells_wbc_umap.png'), width = 12, height = 5, dpi=400)

# Save merged object
saveRDS(merged, paste0(objDir, "/wbcs/bcells_and_wbcs.rds"))

###############################################################
# Run DGEA
###############################################################

# Read in object if needed
merged <- readRDS(paste0(objDir, "/wbcs/bcells_and_wbcs.rds"))

# Create new dir
dir.create(glue::glue("{resDir}clustering"))

# Cluster
resolutions = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
for (res in resolutions){
  # Cluster
  merged <- FindClusters(merged, resolution=res, cluster.name=paste0('res',res)) # Setting smallest resolution
  # Plot UMAP
  DimPlot(merged, reduction = "umap", group.by = paste0('res',res), alpha = 0.7)
  ggsave(paste0(resDir,'clustering/clusters_', res, '.png'), width = 6, height = 5, dpi=400)
}


# Plotting
DimPlot(merged, group.by=c('res0.06', 'condition'))
ggsave(paste0(resDir,'clustering/clusters_and_condition.png'), width = 10, height = 5, dpi=400)


# DGEA by clustering vs condition
samples <- c("+Prolif B cells", "+Ribo_-Prolif B cells", "-Ribo B cells") 
clusters <- c("1", "3", "0") 

# -------------------------------------------------------------------------------------------
# By Bcell group
# -------------------------------------------------------------------------------------------
# Set idents
Idents(merged) <- merged$condition

samples <- c("+Prolif B cells", "+Ribo_-Prolif B cells", "-Ribo B cells") 
while (length(samples) > 0) {
  
  # Get current sample (entire ID)
  curr_sample <- samples[1]
  # Remove from list
  samples <- samples[-1]
  
  # Loop over each sample vs spleen dataset
  # for (sample in samples) {
  
  # Make results dir
  
  print(paste('Running DGEA for', curr_sample, 'vs', 'WBC_01 (Healthy B cells)'))
  
  # Get markers
  markers <- FindMarkers(object = merged, 
                         ident.1 = curr_sample, 
                         ident.2 = 'WBC_01',
                         recorrect_umi=FALSE)
  
  # Sort by log2FC
  markers_ordered <- markers[order(-markers$avg_log2FC), ]
  significant_genes <- markers_ordered %>%
    filter(p_val_adj < 0.05)
  
  # Save filtered gene list
  write.table(significant_genes, file=paste0(resDir,curr_sample, '_WBCs_sig_genes.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
  # Save unfiltered gene list
  write.table(markers, file=paste0(resDir,curr_sample, '_WBCs_unfiltered_genes.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
  
  # Plot
  EnhancedVolcano(
    markers,
    lab = rownames(markers),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    #    selectLab = genes_to_plot,
    #xlim = c(-5, 6.5),
    #ylim = c(0, 40),
    pCutoff = 0.05,
    FCcutoff = 0.6,
    pointSize = 3,
    labSize = 3.5,
    axisLabSize = 13,
    #colAlpha = 1,
    legendLabels = c('NS', 'log_2FC', 'pval_adj', 'pval_adj and log2FC'),
    #legendPosition = 'right',
    #legendLabSize = 15,
    title = paste(curr_sample, "vs", 'Healthy B cells'),
    subtitle='',
    caption = paste0('n = ', nrow(markers), ' genes'),
    captionLabSize = 12,
    max.overlaps = 30,
    #  boxedLabels = TRUE,
    #  drawConnectors = TRUE
  ) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(resDir, curr_sample, '_volcano.png'), width=6, height=6, dpi=400)
}


