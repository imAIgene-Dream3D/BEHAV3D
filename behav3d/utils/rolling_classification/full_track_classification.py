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

min_length = 100
max_length = 100
adata_filt = filter_and_truncate_tracks(
    adata_full,
    groupby_cols=["sample_name", "TrackID"],   # <-- set to your actual column
    time_col="position_t",      # <-- set to your actual column (e.g., "t", "time", "frame"); or None
    min_length=min_length,
    max_length=max_length
)


test,blocks = extract_descibing_track_state_features(
    adata_filt
)

test = run_leiden_clustering(
    test, 
    n_neighbors=50,
    resolution=0.2, 
    metric="cosine",
    method="umap",
    use_rep="X",
    key_added="ClusterID",
    random_state=seed
    )

sc.tl.umap(
    test,
    min_dist=0.1,
    random_state=seed,
)

sc.pl.umap(test, color="ClusterID", size=6, alpha=0.5)

sc.tl.dendrogram(test, groupby="ClusterID")
sc.pl.heatmap(
    test,
    var_names=test.var_names,
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
    test,
    var_names=test.var_names,    
    groupby="ClusterID",
    standard_scale="var",
    figsize=(20, 20),
    swap_axes=True,
    dendrogram=True
)


plot_exemplar_tracks_by_cluster(
    adata_filt,
    test,
    n_per_cluster=20
)
