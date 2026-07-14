#!/usr/bin/env python

### ALL FILE PATHS NEED TO BE UPDATED BEFORE SCRIPT CAN RUN ###

# Catch any errors
import sys
import os

log_file = 'c2l_part2.log'
error_file = 'c2l_part2_errors.err'

sys.stdout = open(log_file, 'w') # Logs go to directory where script is
sys.stderr = open(error_file, 'w')

# Packages
import scanpy as sc
import numpy as np
import cell2location as c2l
import matplotlib
import scipy
from matplotlib import rcParams
import matplotlib.pyplot as plt
import scipy.sparse

# Settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor="white")


# Figure settings
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
rcParams['figure.figsize'] = 6,6

# File paths
data_dir = ''

# Set up results directory and paths to results folders
results = ''

pt1_run_name = f'{results}/Reference_signatures'
pt2_run_name = f'{results}/Spatial_mapping'

# Spatial
adata_st = sc.read_h5ad(f'{data_dir}orv01_spatial_clustered.h5ad')

# Get rid of MT genes in spatial
adata_st.var["MT_gene"] = [gene.startswith("MT") for gene in adata_st.var_names]
adata_st.obs['sample'] = list(adata_st.uns['spatial'].keys())[0]

# Remove genes
adata_st.obsm["MT"] = adata_st[:, adata_st.var["MT_gene"].values].X.toarray()
adata_st = adata_st[:, ~adata_st.var["MT_gene"].values]

# Extracting reference cell types signatures as a pd.DataFrame (inf_aver)

# For spatial mapping we just need the estimated expression of every gene in every cell type
# Read in the reference to get gene signatures
adata_file = f"{pt1_run_name}/reference_signatures.h5ad"
adata_ref_fil = sc.read_h5ad(adata_file)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref_fil.varm.keys():
    inf_aver = adata_ref_fil.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref_fil.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref_fil.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref_fil.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref_fil.uns['mod']['factor_names']
inf_aver.iloc[0:5]

# RUN MODEL


# Find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_st.var_names, inf_aver.index)
adata_st = adata_st[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()


# Setting raw counts
adata_st.X = adata_st.layers['counts'].copy()

# Prepare anndata for cell2location model
c2l.models.Cell2location.setup_anndata(adata=adata_st) # no batch key

# create model
mod = c2l.models.Cell2location(
    adata_st, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=8, # standard Visium data (spot-based)
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)

print('Training the model....')
mod.train(max_epochs=15000, # Testing first at 15000 
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          accelerator='gpu'
         )
print('Model training completed. Saving model.')

# Exporting model to results file
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_st = mod.export_posterior(
    adata_st, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
)

# Save model
mod.save(f"{pt2_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{pt2_run_name}/c2l_label_predications.h5ad"
adata_st.write(adata_file)

print('Done! Script ran successfully.')