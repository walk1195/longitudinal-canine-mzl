#############################################################################################
############################################################################################
### Initial Preprocessing & QC of ORv-01 Longitudinal scRNAseq Data
### Samples: ORv-01 Longitudinal scRNAseq (n=7)
### Author: Grace Walker
### Date: February 28, 2025
#############################################################################################
#############################################################################################

# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(ggplot2)
library(tidyverse)
library(Seurat)
library(patchwork)
library(scales)

# Set dirs
dataDir = ""
resDir = ""
outDir = ""

#############################################################################################
# Generating Seurat objects from CellRanger output (per sample)
#############################################################################################

# Read in sample names
samples <- readLines("samples.txt")

# Iterate over samples
sample_objects <- list() # Empty list to store objects
for (sample_id in samples) {
  # Read in filtered data
  print(paste('Reading in ', sample_id, sep=""))
  # Get full dir
  directory = paste(dataDir, sample_id,  "/outs/filtered_feature_bc_matrix/", sep="")
  data <- Read10X(data.dir = directory)
  
  # Generate Seurat object
  batch_key = substr(sample_id, 17, nchar(sample_id)-7) # Grab only unique ID
  seurat_obj <- CreateSeuratObject(counts = data, project = batch_key, min.cells = 3, min.features = 200)
  
  # Append curr object to list
  sample_objects[[batch_key]] <- seurat_obj

  # Save object
  print(paste("Saving ", sample_id, sep=""))
  saveRDS(seurat_obj, file = paste0(outDir, batch_key, ".rds"))
}

#############################################################################################
# Merge samples
##############################################################################################

sample_names = c("PreTx", "D3FNA", "D7FNA", "D10FNA", "D10resect", "D28FNA", "D49FNA")

sample_objects <- list()
for (i in sample_names) {
  file_path = paste0(outDir, i, '.rds')
  curr_object = readRDS(file_path)
  sample_objects[[i]] <- curr_object
}

# Fix D3 name because it was mislabelled in sequencing
names(sample_objects)[names(sample_objects) == "D3"] <- "D3FNA"

# Merge
merged <- merge(sample_objects[[1]], y = sample_objects[2:length(sample_objects)], 
             add.cell.ids = names(sample_objects), project="ORv-01")

# Add sample type condition
merged$sample_type <- ifelse(str_detect(merged@meta.data$orig.ident, "^FNA"),
                        "FNA","resect")
# Correcting PreTx col
merged$sample_type <- ifelse(str_detect(merged@meta.data$orig.ident, "PreTx"),
                             "FNA", merged$sample_type)

# Save
saveRDS(merged, paste0(outDir,"/samples_merged.rds"))


#############################################################################################
# Visualize Initial Cell Counts per sample
##############################################################################################
# Preserve sample order
merged@meta.data$orig.ident <- factor(merged@meta.data$orig.ident, levels = unique(merged@meta.data$orig.ident))

# Plot
merged@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar(color="black") +
  stat_count(geom = "text", colour = "black", size = 3.5, 
             aes(label = after_stat(count)),
             position=position_stack(vjust=1.05))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Sample") +
  labs(x = "Sample", fill = "Sample")
ggsave(paste0(figDir,"cells_per_sample.png"), width = 8, height = 5, dpi=400)


#############################################################################################
# Individual Sample Filtering
##############################################################################################

# Read in single dataset
sample = "D49FNA"
s1 = readRDS(glue::glue("{outDir}/{sample}.rds"))

# Set curr res dir
curr_dir <- paste0(figDir,'Preprocessing/',sample, '/')
dir.create(curr_dir)

# Recalculate percent MT
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^MT-")

# Visualize filter thresholds
# Violin plots
VlnPlot(s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), layer="counts", ncol = 3)

# Scatter plots
plot1 <- FeatureScatter(s1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(s1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(s1, feature1 = "percent.mt", feature2 ="nFeature_RNA")
plot1 + plot2 + plot3

# Density plots
# Set thresholds based on visual inspection and QC plots
mt = 5
feature_min = 350
count_min = 650

p1 <- s1@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() + 
  scale_x_log10() + 
  geom_vline(xintercept = feature_min,color="red",linetype="dotted") +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none") +
  ggtitle("nFeature")
p2 <- s1@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() + 
  scale_x_log10() + 
  geom_vline(xintercept = count_min,color="red",linetype="dotted") +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none") +
  ggtitle("nCount")
p3 <- s1@meta.data %>% 
  ggplot(aes(x=percent.mt,fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10(labels = label_comma()) + 
  geom_vline(xintercept = mt,color="red",linetype="dotted") +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("MT Count")
p1 + p2 + p3
ggsave(paste0(figDir, sample, "_densities.png"), width = 10, height = 4, dpi=400)

# Scatter (all metrics + thresholds)
ggplot(s1@meta.data) +
  geom_point(aes(x=nCount_RNA,y=nFeature_RNA,fill=percent.mt > 5),shape=21,alpha=0.4) + 
  theme_classic() +
  scale_x_log10()+
  scale_y_log10()+
  facet_grid(.~orig.ident) +
  geom_vline(xintercept = count_min,color="red",linetype="dotted")+
  geom_hline(yintercept=feature_min,color="red", linetype="dotted")+
  scale_fill_manual(values=c("FALSE"="lightblue", "TRUE"="purple"))  # Customize colors
ggsave(paste0(figDir, sample, "_thresholds.png"), width = 10, height = 4, dpi=400)


# Apply QC filtering
prefilter_count = length(Cells(merged))
merged@meta.data$keep <- with(merged@meta.data, ifelse(orig.ident == sample & 
                                                         nFeature_RNA > feature_min & nCount_RNA > count_min & percent.mt < mt,
                                                       TRUE, FALSE))
merged_filtered <- subset(merged, subset = keep == TRUE | orig.ident != sample)
postfilter_count = length(Cells(merged_filtered))
total = prefilter_count - postfilter_count
print(paste("Cells removed:", total))

# Save filtered object
saveRDS(merged_filtered, paste0(outDir,"/merged_filtered.rds"))

