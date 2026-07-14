#############################################################################################
############################################################################################
### Running InferCNV for CNV detection of ORv01 Longitudinal scRNAseq Data: BATCH
### Samples: ORv01 Longitudinal scRNAseq (n=7)
### Author: Grace Walker
### Date: January 22, 2026
#############################################################################################
#############################################################################################

# -------------------------------------------------------------------------------------------
# The following environment set up is required to enable successful import of inferCNV :

# module load R/4.4.0-openblas-rocky8
# module load JAGS/4.3.0-gcc7.2.0

# This is called in the batch script, to be submitted to SLURM
# -------------------------------------------------------------------------------------------

# Load libs
library(tidyverse) # dplyr and ggplot
library(Seurat) # for seurat analysis
library(patchwork) # for plotting orientation
library(ggplot2) # for visualization
library(scales) # for plotting
library(presto) # for faster FindMarkers 
library(dittoSeq) # for dittoplot
library(infercnv) # for cnv detection
library(optparse) # for batching samples

# Params
options(future.globals.maxSize = 4 * 1024^3)  # Set to 4 GB
gc() # free up memory

options(scipen = 100) # ensure no scientific notation, required for inferCNV
options("preferRaster" = FALSE) # no raster for plotting heatmaps

# -------------------------------------------------------------------------------------------

##### Get current sample ID
option_list = list(make_option(c("-s", "--sampleID"),
                               type = "character",
                               default = NULL,
                               help = "Sample ID of current sample",
                               metavar = "character"))

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);
curr_sample = opt$sampleID

#### Set file paths
objDir <- "" # Must be filled in for script to run
outDir <- glue::glue("/infercnv/batch_results/",curr_sample, "/")

# Create output dir and set as wd
dir.create(file.path(outDir), showWarnings = FALSE) # Suppress warning if dir already exists
setwd(outDir)

#### Read in large seurat object
s1 <- readRDS(paste0(objDir, 'processed_object_final.rds')) # with D49

############################################################################################################

# Subset to current sample ID
s1.subset <- subset(s1, subset = orig.ident == curr_sample)

# Get gene expr matrix
# Join layers
s1.subset <- JoinLayers(s1.subset, assay = "RNA")

# Extract raw expr data
exp.rawdata <- GetAssayData(s1.subset, assay = "RNA", layer = "counts")

# Print checks
message("Processing sample: ", curr_sample)
message("Cell Type Distribution:", table(s1.subset$cell_type))


#### Create annotation file -- col1=cell_names, col2=cell_type labels
Idents(s1.subset) <- s1.subset$cell_type
annotations <- as.matrix(Idents(s1.subset))
colnames(annotations) <- NULL


#### Create object
# First grab list of all chrs
gene_positions <- read.table(paste0(gene_pos_path,"gene_positions.txt"), sep='\t')
chrs <- unique(gene_positions$V2)

# Subset to only chrs we care about
chrs_to_exclude <- chrs[grepl('chrUn_*',chrs)]


#### Create InferCNV object
message('Generating infercnv object .....')
infercnv_obj_pretx = CreateInfercnvObject(raw_counts_matrix=exp.rawdata, # Raw counts
                                    annotations_file=annotations, # Cell type annotations
                                    delim="\t",
                                    gene_order_file=paste0(gene_pos_path,"gene_positions.txt"), # Gene positions
                                    ref_group_names=c("T cells", "Myeloid cells"),
                                    chr_exclude = chrs_to_exclude) # Customized chr list excluding certain chrs


#### Run infercnv 
message('Running infercnv ......')

### On the cutoff parameter:
# The cutoff value determines which genes will be used for the infercnv analysis.
# Genes with a mean number of counts across cells will be excluded.
# For smart-seq (full-length transcript sequencing, typically using cell plate assays rather than droplets), a value of 1 works well.
# For 10x (and potentially other 3'-end sequencing and droplet assays, where the count matrix tends to be more sparse), a value of 0.1 is found to generally work well.


infercnv_obj = infercnv::run(infercnv_obj_pretx,
                             cutoff=0.1, # Recommended for 10x Genomics data
                             out_dir=outDir, # set output dir
                             cluster_by_groups=TRUE, # clustering by cell types
                             denoise=TRUE,
                             HMM=TRUE,
                             BayesMaxPNormal=0.3, # More stringent posterior threshold
                             num_threads=128, # Optimize computational load
                             leiden_resolution=0.0002, # Only for rerunning the subclustering
                            # up_to_step=15
                             ) 
message(glue::glue('InferCNV complete. Results folder can be found at ', outDir))


### Results:

# infercnv.preliminary.png : the preliminary inferCNV view (prior to denoising or HMM prediction)
# infercnv.png : the final heatmap generated by inferCNV with denoising methods applied.
# infercnv.references.txt : the 'normal' cell matrix data values.
# infercnv.observations.txt : the tumor cell matrix data values
# infercnv.observation_groupings.txt : group memberships for the tumor cells as clustered.
# infercnv.observations_dendrogram.txt : the newick formatted dendrogram for the tumor cells that matches the heatmap.


##### Update seurat obj with inferCNV results
s1.subset <- add_to_seurat(seurat_obj = s1.subset,
                           infercnv_output_path = outDir)

# Save subsetted object with cnv details
obj_name = glue::glue(curr_sample,'_infercnv_object.rds')
saveRDS(s1.subset, paste0(outDir, obj_name))


#### Plot results
message("Plotting InferCNV results....")

# Set output dir
outDir <- glue::glue(outDir, 'umaps/')
dir.create(file.path(outDir))
setwd(outDir)

# Plotting cell types first
DimPlot(s1.subset, group.by='cell_type') + ggtitle(paste0('Cell Types (',curr_sample, ')'))
ggsave('umap_subset_cell_types.png', dpi=400, height=5, width=6)


# Plotting cell types first
DimPlot(s1.subset, group.by=c('has_cnv_chr13', 'has_dupli_chr13', 'has_loss_chr13'))
ggsave('umap_chr13_cnvs.png', dpi=400, height=5, width=12)

FeaturePlot(s1.subset, features=c('proportion_scaled_dupli_chr13'))
ggsave('featureplot_chr13_scaled_dupli.png', dpi=400, height=5, width=6)


# What cells have the most CNVs?
DimPlot(s1.subset, group.by=c('infercnv_subcluster')) + ggtitle("InferCNV subclusters")
ggsave('umap_subclusters.png', dpi=400, height=5, width=6)


message('InferCNV complete.')


