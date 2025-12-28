from behav3d.utils import load_behav3d_metadata, check_behav3d_metadata

import pandas as pd
from pathlib import Path

import pandas as pd
import numpy as np

import umap

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random

from sklearn.cluster import KMeans, HDBSCAN
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold


from pathlib import Path

import seaborn as sns
# import hdbscan

from behav3d.utils.rolling_classification import *
from behav3d.utils.rolling_classification.features import *
from behav3d.utils.rolling_classification.filtering import *
from behav3d.utils.rolling_classification.clustering import *
from behav3d.utils.rolling_classification.plotting import *
from behav3d.utils.rolling_classification.videos import *

import numpy as np
#%matplotlib inline

seed = 123
random.seed(seed)
np.random.seed(seed)

ssd_dir = r"F:/"
ssd_dir = Path(ssd_dir)
outfolder = Path(ssd_dir, r"BHVD_BEHAV3D\BEHAV3D_python\rolling_classification")
adata_full = sc.read_h5ad(Path(outfolder,"adata_full.h5ad"))

min_length = 100
max_length = 100
adata_filt = filter_and_truncate_tracks(
    adata_full,
    groupby_cols=["sample_name", "TrackID"],   # <-- set to your actual column
    time_col="position_t",      # <-- set to your actual column (e.g., "t", "time", "frame"); or None
    min_length=min_length,
    max_length=max_length
)

state_col = 'ClusterID'
test, blocks = extract_descibing_track_state_features(
    adata_filt,
    state_col=state_col
)

"""
Add L2 normalization??
"""
test_norm = l2_normalize_features_blocks(test, blocks=blocks)

to_run = test_norm
to_run = run_pca(
    to_run,
    ncomps=len(to_run.var_names),
    pca_var=0.95, 
    seed=None
    )


# to_run = run_leiden_clustering(
#     to_run, 
#     n_neighbors=50,
#     resolution=0.2, 
#     metric="cosine",
#     method="umap",
#     use_rep="X_pca",
#     key_added="ClusterID",
#     random_state=seed
#     )


to_run = run_leiden_clustering(
    to_run, 
    n_neighbors=30,
    resolution=0.4, 
    # metric="cosine",
    method="umap",
    use_rep="X",
    key_added="ClusterID",
    random_state=seed
    )

sc.tl.umap(
    to_run,
    min_dist=0.1,
    random_state=seed,
)

sc.pl.umap(to_run, color="ClusterID", size=4, alpha=0.5)


sc.tl.dendrogram(to_run, groupby="ClusterID")
sc.pl.heatmap(
    to_run,
    var_names=to_run.var_names,
    # var_group_labels=list(feature_dict.keys()),
    # var_group_rotation=0,
    groupby="ClusterID",
    standard_scale="var",
    figsize=(25, 20),
    swap_axes=True,
    dendrogram=True,
    show_gene_labels=True,
)

sc.pl.matrixplot(
    to_run,
    var_names=to_run.var_names,    
    groupby="ClusterID",
    standard_scale="var",
    figsize=(20, 20),
    swap_axes=True,
    dendrogram=True
)

plot_exemplar_tracks_by_cluster(
    adata_filt,
    to_run,
    n_per_cluster=5,
    state_key="ClusterID",
)

plot_exemplar_tracks_by_cluster(
    adata_filt,
    to_run,
    n_per_cluster=5,
    state_key="ClusterID_filt5",
)
