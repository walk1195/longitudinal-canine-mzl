#############################################################################################
############################################################################################
### DGEA and Pathway Analysis of ORv-01 Longitudinal scRNAseq Data
### Samples: ORv-01 longitudinal scRNAseq (n=7)
### Author: Grace Walker
### Date: March 30, 2026
#############################################################################################
#############################################################################################

# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(tidyverse) # dplyr and ggplot
library(Seurat) # for seurat analysis
library(patchwork) # for plotting orientation
library(scales) # for plotting
library(presto) # for faster FindMarkers 
library(clustree) # for optimal clustering
library(harmony) # for integration
library(dittoSeq) # for dittoplot
library(celldex) # for reference datasets
library(SingleR) # for automated annotation
library(scran) # for normalization
library(ggpubr) # for stats tests 
library(EnhancedVolcano) # for dge visualization
library(clusterProfiler) # for pathway anlaysis

# Params
options(future.globals.maxSize = 4 * 1024^3)  # Set to 4 GB
gc() # free up memory

# Set dirs
objDir = ""

# -------------------------------------------------------------------------------------------

# Read in data
s1 <- readRDS(paste0(objDir, 'processed_object_final.rds'))
s1.tcells <- readRDS(paste0(objDir, '/tcell_subsets/tcell_subset_object_final.rds'))

#############################################################################################
# Differential Gene Expression Analysis
#############################################################################################

### B cells
resDir <- ""

# Add col of broad cell type
cell_types <- c(
  '0' = "B cells",
  '1' = "B cells",
  '2' = "T cells",
  '3' = "B cells",
  '4' = "Myeloid cells",
  '5' = "B cells")

# Assign new metadata column with mapped values
s1$cell_type_broad <- plyr::mapvalues(
  x = as.character(s1$res0.05),
  from = names(cell_types),
  to = cell_types
)

# Add a col grouping by cell type and sample
s1$celltype.sample <- paste(s1$cell_type_broad, s1$orig.ident, sep = "_")
Idents(s1) <- s1$celltype.sample

# Subset method 1: 500 per Bcell type, per sample (1500 cells per sample total)
# bcells_only <- subset(s1, subset=s1$cell_type_broad %in% c('B cells')) # Get only B cells
# samples <- unique(s1$orig.ident)
# bcell_types <- unique(bcells_only$cell_type_final)
# selected_cells <- c()
# set.seed(123)
# for (sample in samples) {
#   for (bcell_type in bcell_types) {
#     cells <- WhichCells(bcells_only, expression = (orig.ident == !!sample & cell_type_final == !!bcell_type))
#     sampled_cells <- if (length(cells) >= 500) sample(cells, 500) else cells
#     selected_cells <- c(selected_cells, sampled_cells)
#   }
# }
# s1.subset <- subset(bcells_only, cells = selected_cells)

# Subset method 2: 500 Bcells per sample total, regardless of subtype
bcells_only <- subset(s1, subset=s1$cell_type_broad %in% c('B cells')) # Get only B cells

# Subset to 500 cells per sample
num_cells_per_sample <- 500 # No. cells
set.seed(123)  # For reproducibility
cells_to_keep <- unlist(lapply(unique(bcells_only$celltype.sample), function(sample_id) {
  # Sample 500 cells from each sample
  sample(Cells(bcells_only)[bcells_only$celltype.sample == sample_id], num_cells_per_sample)
}))
s1.subset <- subset(bcells_only, cells = cells_to_keep)

# Check cell type distributions
table(s1.subset$celltype.sample)

# Prep SCT 
s1.subset <- PrepSCTFindMarkers(s1.subset,verbose=T)

# Initialize sample vector - all B cells
sample_ids <- unique(s1.subset$celltype.sample) # To get all

#### Run pairwise comparisons 
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
    markers <- FindMarkers(object = s1.subset, 
                           ident.1 = curr_sample1, 
                           ident.2 = sample,
                           recorrect_umi=FALSE)
    
    # Sort by log2FC
    markers_ordered <- markers[order(-markers$avg_log2FC), ]
    sig_genes <- markers_ordered %>%
      filter(p_val_adj < 0.05)
    
    # Get upreg genes only
    upreg_genes <- sig_genes %>%
      filter(avg_log2FC>0.6)
    
    # Save sig genes
    write.table(sig_genes, file=paste0(resDir, curr_sample_name, '_', timepoint, '_genes_significant.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
    # Save total gene list
    write.table(markers, file=paste0(resDir, curr_sample_name, '_', timepoint, '_genes_all.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
    # Save upreg genes
    write.table(upreg_genes, file=paste0(resDir, curr_sample_name, '_', timepoint, '_genes_upreg.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
    
    
    # Plot 
    # Get auto plotted genes (highest sig) to cross check
    # auto_genes <- head(rownames(markers)[abs(markers$avg_log2FC) > 0.8 & markers$p_val_adj < 0.04], n=20)
    # # Genes to plot
    # custom_genes <- c('TOP2A', 'CENPF', 'MKI67', 'PCLAF', 'PFN1',
    #                   'LOC119863879', 'LOC100687749', 'RPL7', 'RPS6', 'TPT1', 'LTB', 'LYZ', 'CD74', 'MS4A1', 'RBM3', 'ELAVL1', 'H5C9', 'ABCA1') #KLRF1, NKG7
    # 
    # genes_to_plot <- c(custom_genes,auto_genes)
    
    # Plot
    EnhancedVolcano(
      markers,
      lab = rownames(markers),
      x = 'avg_log2FC',
      y = 'p_val_adj',
      #selectLab = custom_genes,
      pCutoff = 0.05,
      FCcutoff = 0.6,
      pointSize = 3,
      labSize = 4,
      axisLabSize = 13,
      legendLabels = c('NS', 'log_2FC', 'pval_adj', 'pval_adj and log2FC'),
      title = paste('B cells', curr_sample_name, 'vs', timepoint),
      subtitle='',
      caption = paste0('n = ', nrow(markers), ' genes'),
      captionLabSize = 12,
      max.overlaps = Inf,
      #boxedLabels = TRUE,
      #drawConnectors = TRUE
    ) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0(resDir, curr_sample_name, '_', timepoint, '_bcells_volcano.png'), width=7, height=7, dpi=400)
  }
}


### T cells
resDir <- ""

# Add a col grouping by cell type and sample
s1.tcells$celltype.sample <- paste(s1.tcells$cell_type, s1.tcells$orig.ident, sep = "_")
Idents(s1.tcells) <- s1.tcells$celltype.sample

# Subset to 500 cells per sample
num_cells_per_sample <- 500 # No. cells
set.seed(72)  # For reproducibility
cells_to_keep <- unlist(lapply(unique(s1.tcells$orig.ident), function(sample_id) {
  # Sample 500 cells from each sample
  sample(Cells(s1.tcells)[s1.tcells$orig.ident == sample_id], num_cells_per_sample)
}))
s1.tcells.subset <- subset(s1.tcells, cells = cells_to_keep)

# Prep SCT 
s1.tcells.subset <- PrepSCTFindMarkers(s1.tcells.subset,verbose=T)

# Initialize sample vector
sample_ids <- unique(s1.tcells.subset$celltype.sample) # To get all

#### Run pairwise comparisons 
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
    markers <- FindMarkers(object = s1.tcells, 
                           ident.1 = curr_sample1, 
                           ident.2 = sample,
                           recorrect_umi=FALSE)
    
    # Sort by log2FC
    markers_ordered <- markers[order(-markers$avg_log2FC), ]
    sig_genes <- markers_ordered %>%
      filter(p_val_adj < 0.05)
    
    # Get upreg genes only
    upreg_genes <- sig_genes %>%
      filter(avg_log2FC>0.6)
    
    # Save sig genes
    write.table(sig_genes, file=paste0(resDir, curr_sample_name, '_', timepoint, '_genes_significant.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
    # Save total gene list
    write.table(markers, file=paste0(resDir, curr_sample_name, '_', timepoint, '_genes_all.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
    # Save upreg genes
    write.table(upreg_genes, file=paste0(resDir, curr_sample_name, '_', timepoint, '_genes_upreg.tsv'), sep='\t', row.names = TRUE, quote = FALSE)
    
    
    # Plot 
    # Get auto plotted genes (highest sig) to cross check
    # auto_genes <- head(rownames(markers)[abs(markers$avg_log2FC) > 0.8 & markers$p_val_adj < 0.04], n=20)
    # # Genes to plot
    # custom_genes <- c('TOP2A', 'CENPF', 'MKI67', 'PCLAF', 'PFN1',
    #                   'LOC119863879', 'LOC100687749', 'RPL7', 'RPS6', 'TPT1', 'LTB', 'LYZ', 'CD74', 'MS4A1', 'RBM3', 'ELAVL1', 'H5C9', 'ABCA1') #KLRF1, NKG7
    # 
    # genes_to_plot <- c(custom_genes,auto_genes)
    
    # Plot
    EnhancedVolcano(
      markers,
      lab = rownames(markers),
      x = 'avg_log2FC',
      y = 'p_val_adj',
      #selectLab = custom_genes,
      pCutoff = 0.05,
      FCcutoff = 0.6,
      pointSize = 3,
      labSize = 4,
      axisLabSize = 13,
      legendLabels = c('NS', 'log_2FC', 'pval_adj', 'pval_adj and log2FC'),
      title = paste('B cells', curr_sample_name, 'vs', timepoint),
      subtitle='',
      caption = paste0('n = ', nrow(markers), ' genes'),
      captionLabSize = 12,
      max.overlaps = Inf,
      #boxedLabels = TRUE,
      #drawConnectors = TRUE
    ) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0(resDir, curr_sample_name, '_', timepoint, '_tcells_volcano.png'), width=7, height=7, dpi=400)
  }
}


################################################################################################################
# Pathway Enrichment with ClusterProfiler
################################################################################################################

# B cells
degsDir <- ""
resDir <- ""

# T cells
degsDir <- ""
resDir <- ""


# Set sample pairs -----------------------------------------

sample_list2 <- c() # Input sample comparisons we want to perform
sample_list2 <- c()

# --------------------------------------- B cells ---------------------------------------
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
    
    # Load in df of sig DEGs for curr pairwise group
    df <- read.table(paste0(degsDir, sample1, '_', sample2, '_genes_significant.tsv'), sep='\t')
    
    # Initialize cutoffs and output file
    log2FC_lower <- -0.6
    log2FC_upper <- 0.6
    
    # Print initial DEG count
    cat(paste0('Initial DEG count:', nrow(df),'\n'))
    
    
    sink('deg_filtering.txt')
    
    # Set log2FC cutoff until < 500 DEGs in list
    
    while (nrow(df) > 500) {
      
      # Print initial DEG amount and cutoffs
      cat(paste0('Initial DEG count:', nrow(df),'\n'))
      cat(paste0('Filtering at (', log2FC_lower, ', ', log2FC_upper, ') log2FC cutoffs... ','\n'))
      
      # Filter
      df <- df %>%
        filter(avg_log2FC < log2FC_lower | avg_log2FC > log2FC_upper)
      
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
    universe_genes <- bitr(rownames(s1), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Cf.eg.db")
    
    ################# Run ClusterProfiler #################
    # Settings
    padj_cutoff <- 0.15 # p-adjusted threshold, used to filter out pathways
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
    #df_kegg_na_uni <- if (!is.null(results_kegg_down_na_uni)) as.data.frame(results_kegg_down_na_uni) else NULL
    df_go_down <- if (!is.null(results_go_down)) as.data.frame(results_go_down) else NULL
    df_go_up <- if (!is.null(results_go_up)) as.data.frame(results_go_up) else NULL
    
    # Save
    if (!dir.exists("Results")) {
      dir.create("Results")
    }
    write.table(df_kegg_down, 'Results/kegg_results_downreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_kegg_up, 'Results/kegg_results_upreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_go_down, 'Results/go_results_downreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_go_up, 'Results/go_results_upreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_go_up, 'Results/go_results_upreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    
  }
}

# --------------------------------------- T cells ---------------------------------------
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
    
    # Load in df of sig DEGs for curr pairwise group
    df <- read.table(paste0(degsDir, sample1, '_', sample2, '_genes_significant.tsv'), sep='\t')
    
    # Initialize cutoffs and output file
    log2FC_lower <- -0.6
    log2FC_upper <- 0.6
    
    # Print initial DEG count
    cat(paste0('Initial DEG count:', nrow(df),'\n'))
    
    
    sink('deg_filtering.txt')
    
    # Set log2FC cutoff until < 500 DEGs in list
    
    while (nrow(df) > 500) {
      
      # Print initial DEG amount and cutoffs
      cat(paste0('Initial DEG count:', nrow(df),'\n'))
      cat(paste0('Filtering at (', log2FC_lower, ', ', log2FC_upper, ') log2FC cutoffs... ','\n'))
      
      # Filter
      df <- df %>%
        filter(avg_log2FC < log2FC_lower | avg_log2FC > log2FC_upper)
      
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
    padj_cutoff <- 0.15 # p-adjusted threshold, used to filter out pathways
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
    #df_kegg_na_uni <- if (!is.null(results_kegg_down_na_uni)) as.data.frame(results_kegg_down_na_uni) else NULL
    df_go_down <- if (!is.null(results_go_down)) as.data.frame(results_go_down) else NULL
    df_go_up <- if (!is.null(results_go_up)) as.data.frame(results_go_up) else NULL
    
    # Save
    if (!dir.exists("Results")) {
      dir.create("Results")
    }
    write.table(df_kegg_down, 'Results/kegg_results_downreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_kegg_up, 'Results/kegg_results_upreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_go_down, 'Results/go_results_downreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_go_up, 'Results/go_results_upreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    write.table(df_go_up, 'Results/go_results_upreg.tsv', row.names = FALSE, sep='\t', quote=FALSE)
    
  }
} 


################################################################################################################
# Replotting results to fix figure dims
################################################################################################################

### Read in results for pair of interest
# PreTx v D10res T cells
resDir <- ""
results_go_up <-read.table(resDir, header=TRUE, sep='\t')
sample1 <- "PreTx"
sample2 <- "D10res"

# First make a new object
res <- new("enrichResult", result = results_go_up)

# Plot
resDir <- ""
p <- dotplot(res, showCategory = 15, font.size = 13) + ggtitle(paste0(sample1, ' v ', sample2, ' (GO Upregulated)'))
ggsave(glue::glue("{resDir}PreTx_D10res_upreg_go_dotplot.png"), plot = p, width = 8, height = 8, dpi=400)


# D28FNA v D49FNA B cells
resDir <- ""
sample1 <- "D28FNA"
sample2 <- "D49FNA"

# First make a new object
results_df <-read.table(resDir, header=TRUE, sep='\t')
res <- new("enrichResult", result = results_df)

# Plot
resDir <- ""
p <- dotplot(res, showCategory = 15, font.size = 13) + ggtitle(paste0(sample1, ' v ', sample2, ' (GO Upregulated)'))
ggsave(glue::glue("{resDir}D28FNA_D49FNA_upreg_go_dotplot.png"), plot = p, width = 8, height = 8, dpi=400)



# D3FNA v D49FNA B cells
resDir <- ""
sample1 <- "D3FNA"
sample2 <- "D49FNA"

# First make a new object
results_df <-read.table(resDir, header=TRUE, sep='\t')
res <- new("enrichResult", result = results_df)

# Plot
resDir <- ""
p <- dotplot(res, showCategory = 15, font.size = 13) + ggtitle(paste0(sample1, ' v ', sample2, ' (GO Upregulated)'))
ggsave(glue::glue("{resDir}D3FNA_D49FNA_upreg_go_dotplot.png"), plot = p, width = 8, height = 8, dpi=400)



