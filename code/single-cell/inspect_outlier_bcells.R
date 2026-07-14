#############################################################################################
############################################################################################
### Inspecting Outlier B cell Population in ORv01 Longitudinal Data
### Samples: ORv01 longitudinal scRNAseq (n=7)
### Author: Grace Walker
### Date: April 29, 2025
#############################################################################################
#############################################################################################

# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(tidyverse) # dplyr and ggplot
library(Seurat) # for seurat analysis
library(ggplot2) # for visualization
library(patchwork) # for plotting
library(CellChat) # for cell communication analysis

# Params
options(future.globals.maxSize = 50 * 1024^3)  # Set to 5 GB
gc() # free up memory
options(stringsAsFactors = FALSE) # required for cellchat

# Dirs
fig_path_tcells <- ""
fig_path <- ""
obj_path <- ""

# Load objects
s1 <- readRDS(paste0(obj_path, 'clustered_object_cleaned_d49.rds'))
tcell_subset <- readRDS(paste0(obj_path, 'tcell_subsets/tcell_subset_object_clustered_2pass.rds'))


############################################################################## 
# Preprocessing
############################################################################## 

############# Merge T cell subtype annotations with large dataset ############# 
# Check for annotations
unique(tcell_subset@meta.data$cell_type)
unique(s1@meta.data$cell_type)

# Get overlapping cells
common_cells <- intersect(Cells(s1), Cells(tcell_subset))

# Extract subtype annotations
subtype_annotations <- tcell_subset$cell_type
names(subtype_annotations) <- Cells(tcell_subset)

# Update main seurat obj
s1$cell_type_final <- s1$cell_type  
s1$cell_type_final[names(subtype_annotations)] <- as.character(subtype_annotations)
unique(s1$cell_type_final) # Check

# Update B cells as outlier B cells
s1$cell_type_final[s1$cell_type_final == 'B cells'] <- 'Outlier B cells'

# Plot to check
DimPlot(s1, group.by='cell_type_final')
ggsave(paste(fig_path, 'umap_outlier_bcells.png', sep=""), width = 7, height = 5, dpi=400)

# Save object
saveRDS(s1, file=paste0(obj_path, 'processed_object_final.rds'))


# Plot original broad labels
DimPlot(s1, group.by='cell_type')
ggsave(paste(fig_path, 'umap_broad_annotations.png', sep=""), width = 7, height = 5, dpi=400)

# Plot only the outlier B cells
celltypes <- unique(s1$cell_type_final)
color_vector <- setNames(rep("lightgray", length(celltypes)), celltypes)
color_vector["Outlier B cells"] <- "dodgerblue"
DimPlot(s1, group.by = "cell_type_final", cols = color_vector)


# Plot histogram of nFeature and nCount
outlier_meta <- s1@meta.data[s1$cell_type_final == "B cells", ]
# nFeature
ggplot(outlier_meta, aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 100, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "nFeatureRNA Distribution", x = "Features", y = "Count")
# nCount
ggplot(outlier_meta, aes(x = nCount_RNA)) +
  geom_histogram(binwidth = 100, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "nFeatureRNA Distribution", x = "Features", y = "Count")


# Get top markers for outlier B cells when compared against all cell types (B and T cells)
Idents(s1) <- s1@meta.data$cell_type_final
# Prep SCT
s1 <- PrepSCTFindMarkers(s1,verbose=T)
# Diffexp
markers_outliers <- FindMarkers(s1, ident.1='Outlier B cells', ident.2=c('-Ribo B cells', '+Prolif B cells', '+Ribo/-Prolif B cells'), only.pos = TRUE, recorrect_umi=FALSE)
# Filter
markers_outliers_filtered <- markers_outliers %>%
  dplyr::filter(avg_log2FC > 1) %>% # filter out less variable genes
  slice_head(n = 100) %>% # top 100
  ungroup() 
# Save
write.table(markers_outliers_filtered, file=paste0(fig_path,'outlier_bcell_vs_other_bcell_subtypes_markers_top100.tsv'), sep='\t', row.names = TRUE, quote = FALSE)


# Get top markers for all cell types
markers_all <- FindAllMarkers(s1, only.pos = TRUE, recorrect_umi=FALSE)
# Filter
markers_all_filtered_30 <- markers_all %>%
  group_by(cluster) %>% # group by cluster
  dplyr::filter(avg_log2FC > 1) %>% # filter out less variable genes
  slice_head(n = 30) %>% # top 20
  ungroup()
markers_all_filtered_50 <- markers_all %>%
  group_by(cluster) %>% # group by cluster
  dplyr::filter(avg_log2FC > 1) %>% # filter out less variable genes
  slice_head(n = 50) %>% # top 20
  ungroup()
# Save
write.table(markers_all_filtered_50, file=paste0(fig_path,'one_vs_all_others_markers_top50.tsv'), sep='\t', row.names = TRUE, quote = FALSE)

# Plot heatmap
png(filename = paste0(fig_path,'one_vs_all_heatmap.png'), width = 8, height = 12, units = "in", res = 400)
dev.off()

# Get only outlier B cells
outlier_markers <- markers_all_filtered_50[markers_all_filtered_50$cluster == "Outlier B cells", ]
write.table(outlier_markers, file=paste0(fig_path,'outlier_bcell_vs_all_others_markers_top50.tsv'), sep='\t', row.names = TRUE, quote = FALSE)


# Get top markers for the outlier B cell cluster when compared against only other T cell subtypes
Idents(tcell_subset) <- tcell_subset@meta.data$cell_type
# Prep SCT
tcell_subset <- PrepSCTFindMarkers(tcell_subset,verbose=T)
# Diffexp
markers <- FindAllMarkers(tcell_subset, only.pos = TRUE, recorrect_umi=FALSE)

# Filter to top 30 markers
markers_filtered <- markers %>%
  group_by(cluster) %>% # group by cluster
  dplyr::filter(avg_log2FC > 1) %>% # filter out less variable genes
  slice_head(n = 100) %>% # top 20
  ungroup() 
# Save results df
write.table(markers_filtered, file=paste0(fig_path, 'tcell_subtype_cluster_markers.tsv'), sep='\t', row.names = TRUE, quote = FALSE)

# Filter to only outlier B cell markers
outlier_markers <- markers_filtered[markers_filtered$cluster == "B cells", ]
# Save
write.table(outlier_markers, file=paste0(fig_path, 'outlier_bcell_vs_tcell_subtypes_markers_top100.tsv'), sep='\t', row.names = TRUE, quote = FALSE)

