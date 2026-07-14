#############################################################################################
############################################################################################
### Running FeatureCounts to Quantify Presence of VSV Transcripts in ORv01 Longitudinal scRNAseq Data
### Samples: ORv01 Longitudinal scRNAseq (n=7)
### Author: Grace Walker
### Date: April 9, 2025
#############################################################################################
#############################################################################################

# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(Rsubread)


# Params
options(future.globals.maxSize = 4 * 1024^3)  # Set to 4 GB
gc() # Free up mem

# Dirs
annoDir <- ""
outDir <- ""
bamDir <- ""

# -------------------------------------------------------------------------------------------

# Read in BAM files
samples <- list.files(pattern="*.bam$", recursive=FALSE)


# Iterate over each file
for (curr_sample in samples) {
	sample_name <- sub(".bam", "", curr_sample)

	# Summarize paired-end reads and counting fragments (instead of reads):
	results <- featureCounts(files=curr_sample, isPairedEnd=TRUE, annot.ext=annoDir,
	isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id")

	# Get results
	raw_counts <- results$counts
	results_df <- results$annotation # df of read counts 
	counts_summary <- results$stat # df of successful counts results

	# Append counts to results df
	results_df$counts <- raw_counts[,1]

	# Save output as 2 tab delimited files
	write.table(results_df, file=paste0(outDir, sample_name, '_counts.txt'), sep='\t', row.names=FALSE)
	write.table(counts_summary, file=paste0(outDir, sample_name, '_summary.txt'), sep='\t', row.names=FALSE)

	# Also output save as rds
	saveRDS(results_df, file=paste0(outDir, sample_name, '_counts.rds'))
}


