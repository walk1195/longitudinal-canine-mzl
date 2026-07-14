#############################################################################################
############################################################################################
### Converting Seurat to adata (for interoperability with Python-based analysis)
### Author: Grace Walker
### Date: April 21, 2026
#############################################################################################
#############################################################################################

# Libs
library(Seurat)
library(rhdf5)
library(Matrix)

# Set dirs
objDir = ""
outDir = ""

# Read in data
s1 <- readRDS(paste0(objDir, 'processed_object_final.rds'))

# Format for saving 
DefaultAssay(object = s1) <- "RNA" # Set defaults assay
s1[["SCT"]] <- NULL # Remove SCT assay

# Merge raw counts layers
s1 <- JoinLayers(s1)

# Normalize 
s1 <- NormalizeData(s1, normalization.method = "LogNormalize")
s1 <- ScaleData(s1)

# Write outfiles
writeMM(s1[["RNA"]]$counts, file=paste0(outDir,'ORv01_counts_matrix.mtx'))
write.csv((s1@meta.data), file = paste0(outDir,"ORv01_obs.csv"))
write.csv(Features(s1, assay = "RNA"), file = paste0(outDir,"ORv01_gene_names.csv"))


