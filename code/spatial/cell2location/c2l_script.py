#!/usr/bin/env python

# Catch any errors
import sys
import os

log_file = 'c2l_errors.log'
sys.stdout = open(log_file, 'w') # Logs go to directory where script is
sys.stderr = open(log_file, 'w')

# Packages
import scanpy as sc
# import squidpy as sq
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
sc.settings.figdir = '' # Update with desired output dir for figures

# Figure settings
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
rcParams['figure.figsize'] = 6,6

# Set up results directory and paths to results folders
results = '' # Update with full path to store results

pt1_run_name = f'{results}/Reference_signatures'

# Read in reference object
## Reference has already been filtered and formatted for the model
adata_ref_fil = sc.read_h5ad('Reference_c2l_filtered.h5ad') # Update with full path to object

# (1) Generating reference gene signatures

# Prepare anndata for the regression model
c2l.models.RegressionModel.setup_anndata(adata=adata_ref_fil,
                        # batch_key='Source', # Only a single sample in reference
                        # cell type, covariate used for constructing signatures
                        labels_key='cell_type_final',
                        # no covariates to add
                       )

# Create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref_fil)

# view anndata_setup as a sanity check
mod.view_anndata_setup()

print('Training the model....')
# Training the model
mod.train(max_epochs=250)
print('Model training completed. Saving model.')

# Exporting the model to results file
adata_ref_fil = mod.export_posterior(
    adata_ref_fil, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)

# Save model
mod.save(f"{pt1_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{pt1_run_name}/reference_signatures.h5ad"
adata_ref_fil.write(adata_file)

print('Script ran successfully.')