#############################################################################################
############################################################################################
### ORv01 Spatial Analysis: Converting ST adata object to Seurat
### Samples: ORv01 D10 Visium v1 Sample
### Author: Grace Walker
### Date: Feb 19, 2026
#############################################################################################
#############################################################################################

# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(reticulate) # For python
library(SeuratDisk) # Requires hdf5 module to be loaded
library(anndata) # For building Seurat object manually
library(Seurat)
library(scater)
library(ISCHIA)
library(dplyr)

# Dirs
dataDir <- ''
objDir <- ''
resDir <- ''

####################################################################################
# Build Seurat object from scratch
####################################################################################

# Build seurat obj
adata <- ad$read_h5ad(glue::glue("{objDir}/orv01_spatial_clustered.h5ad")) # Requires anndata library

# Read in raw 10x data
spatial <- Load10X_Spatial(
  dataDir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "orv01",
  bin.size = NULL,
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL)

# Filter based on the spots in adata
filtered_barcodes <- py_to_r(adata$obs_names) # Extract spatial coords

# Check
head(filtered_barcodes)
head(colnames(spatial))

# Subset
cells_to_keep <- intersect(filtered_barcodes, colnames(spatial))
spatial <- subset(spatial, cells = cells_to_keep)

# Transfer metadata
spatial@meta.data <- py_to_r(adata$obs)


# Function to clean object up
clean_visium_object <- function(obj, assay = "Spatial") {
  stopifnot(assay %in% names(obj@assays))
  
  keep_cells <- rownames(obj@meta.data)
  
  # Extract counts matrix (this is the source of truth)
  counts <- GetAssayData(obj, assay = assay, layer = "counts")[, keep_cells, drop = FALSE]
  
  # Create a NEW Seurat object
  new_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = obj@meta.data,
    assay = assay
  )
  
  # Copy images (safe in VisiumV2)
  new_obj@images <- obj@images
  
  # Copy misc / tools if needed
  new_obj@misc <- obj@misc
  
  # Filter images 
  cells.keep <- Cells(obj)
  
  for (img in names(obj@images)) {
    obj@images[[img]] <- subset(
      obj@images[[img]],
      cells = cells.keep
    )
  }
  
  return(new_obj)
}

# Subset
spatial <- clean_visium_object(spatial)

# Check
length(Cells(spatial))
length(rownames(spatial@meta.data))

# Save
saveRDS(spatial, glue::glue("{objDir}/orv01_spatial_clustered.rds"))
