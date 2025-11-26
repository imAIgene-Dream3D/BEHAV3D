

"""
#############################

TSFRESH

#############################
"""


from test.tsfresh import extract_features, select_features
from tsfresh.utilities.dataframe_functions import impute

id_columns = ["sample_name", "TrackID"]
# features=[
#     # "elongation",
#     # "sphericity",
#     "percentage_dead_mask",
#     # "nr_dead_mask_pixels",
#     "organoid_contact",
#     "tcell_contact",
#     # "displacement",
#     # "mean_square_displacement",
#     "speed",
#     # # "dead",
#     # "active_tcell_contact",
#     # "position_t"
# ]

df_tracks_tsfresh = df_tracks.copy()

min_track_length=100
max_track_length=100

if min_track_length is not None:
    df_tracks_tsfresh=df_tracks_tsfresh.groupby(id_columns).filter(lambda group: len(group) >= min_track_length).reset_index(drop=True)
if max_track_length is not None:
    df_tracks_tsfresh=df_tracks_tsfresh.groupby(id_columns).apply(lambda group: group.iloc[:max_track_length]).reset_index(drop=True)
    
    
df_tracks_tsfresh["composite_id"] = (
    df_tracks_tsfresh["sample_name"] + "--" + df_tracks_tsfresh["TrackID"].astype(str)
)

df_tracks_tsfresh = df_tracks_tsfresh[features+["composite_id","position_t"]]
df_tracks_tsfresh["organoid_contact_pixels"] = df_tracks_tsfresh["organoid_contact_pixels"].astype(np.int16)
df_tracks_tsfresh["tcell_contact_pixels"] = df_tracks_tsfresh["tcell_contact_pixels"].astype(np.int16)

extracted_features = extract_features(
    df_tracks_tsfresh, 
    column_id="composite_id", 
    column_sort="position_t",
    impute_function=impute,
    n_jobs=0
    )

X_reduced, dropped, report = drop_highly_correlated_features(
    df=extracted_features,
    feature_cols=extracted_features.columns.tolist(),
    threshold=0.95
)
print(f"Dropped {len(dropped)} features.")
len(X_reduced)  # see which were kept vs dropped and why

df_tsfresh = X_reduced.reset_index()
df_tsfresh[["sample_name", "TrackID"]] = df_tsfresh["index"].str.split("--", expand=True)

X_scaled = StandardScaler().fit_transform(X_reduced)

pca_model = PCA(n_components=0.95, random_state=seed)
X_pca = pca_model.fit_transform(X_scaled)

# df_pca = pd.DataFrame(X_pca[:, :2], columns=["PC1", "PC2"], index=df_windows_descriptive.index)
# df_analysis = pd.concat([df_windows_descriptive, df_pca], axis=1)
df_tsfresh["PC1"] = X_pca[:, 0]
df_tsfresh["PC2"] = X_pca[:, 1]

"""
UMAP on raw features
"""
umap_model = umap.UMAP(
            n_components=2, 
            n_neighbors=100, 
            min_dist=0.1, 
            # init="random", 
            metric = "cosine",
            random_state=seed,
            )
umap_embedding = umap_model.fit_transform(X_pca)
df_tsfresh["UMAP1"] = umap_embedding[:,0]
df_tsfresh["UMAP2"] = umap_embedding[:,1]

hdbscan_clusterer=HDBSCAN(
    min_cluster_size=100,
    min_samples=50,
    metric="euclidean",
    alpha=1.0,
    cluster_selection_method="eom",
    cluster_selection_epsilon=0.0,
    allow_single_cluster=False,
    leaf_size=40,
    algorithm="auto",                # sklearn’s default
    n_jobs=None,
    copy=False
    )
cluster_labels = hdbscan_clusterer.fit_predict(umap_embedding)
df_tsfresh["cluster_label_hdbscan"] = cluster_labels
df_tsfresh["ClusterID"] = cluster_labels

n_clusters = 8  # you can tune this
kmeans_clusterer = KMeans(
    n_clusters=n_clusters, 
    random_state=seed, 
    n_init="auto"
    )
cluster_labels = kmeans_clusterer.fit_predict(umap_embedding)
df_tsfresh["cluster_label_kmeans"] = cluster_labels

# ---------- 1) PCA & UMAP quick looks ----------
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
sns.scatterplot(data=df_tsfresh, x="PC1", y="PC2", s=8, alpha=0.6, edgecolor=None)
plt.title("PCA (PC1 vs PC2)")

plt.subplot(1,2,2)
sns.scatterplot(data=df_tsfresh, x="UMAP1", y="UMAP2", s=8, alpha=0.6, edgecolor=None)
plt.title("UMAP (2D)")
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,5))
sns.scatterplot(
    data=df_tsfresh, x="UMAP1", y="UMAP2",
    hue="cluster_label_kmeans", palette="tab20",
    s=10, alpha=0.8, edgecolor=None
)
plt.title("UMAP colored by KMeans clusters")
plt.show()

plt.figure(figsize=(6,5))
sns.scatterplot(
    data=df_tsfresh, x="UMAP1", y="UMAP2",
    hue="cluster_label_hdbscan", palette="tab20",
    s=10, alpha=0.8, edgecolor=None
)
plt.title("UMAP colored by HDBSCAN clusters (−1 = noise)")
plt.show()

test = [c for c in X_reduced.columns if c.endswith("mean")]
plot_umap_feature_grid(df_tsfresh, feature_cols=test)








### Z-scale certain features
scaler = StandardScaler()
# scaler = RobustScaler()
df_tracks_tsfresh[features] = pd.DataFrame(scaler.fit_transform(pd.DataFrame(df_tracks_tsfresh[features])), columns=features)


extracted_features = extract_features(
    df_tracks_tsfresh, 
    column_id="composite_id", 
    column_sort="position_t",
    impute_function=impute,
    n_jobs=0
    )



extracted_features_scaled = extracted_features.copy()
scaler = StandardScaler()
extracted_features_scaled = scaler.fit_transform(extracted_features)
extracted_features_scaled = pd.DataFrame(extracted_features_scaled, 
                                         index=extracted_features.index, 
                                         columns=extracted_features.columns)
# Start with your extracted tsfresh features
features = extracted_features_scaled.copy()  # just to be safe

# 1. Drop columns that are all NaN
features_clean = features.dropna(axis=1, how='all')

# 2. Drop constant columns (no variance)
features_clean = features_clean.loc[:, features_clean.nunique() > 1]

# 3. Drop highly correlated columns (e.g., correlation > 0.95)
corr_matrix = features_clean.corr().abs()
upper_triangle = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
to_drop = [col for col in upper_triangle.columns if any(upper_triangle[col] > 0.95)]
features_clean = features_clean.drop(columns=to_drop)

# # 4. (Optional) Standardize features (for later clustering)
# scaler = StandardScaler()
# features_scaled = scaler.fit_transform(features_clean)

# # 5. Keep as DataFrame with original index and column names
# features_clean_scaled = pd.DataFrame(features_scaled, 
#                                      index=features_clean.index, 
#                                      columns=features_clean.columns)

# extracted_features.index.name = 'composite_id'
# features_clean = features_clean.reset_index()
# df_extracted_features.index.name = 'index'
# summarized_features = df_tracks_tsfresh.groupby(["composite_id"])[features].mean().reset_index()
# summarized_features = extracted_features.groupby(["composite_id"]).mean().reset_index()
df_extracted_features = features_clean.copy()
df_extracted_features = df_extracted_features.reset_index()

df_extracted_features[["sample_name", "TrackID"]] = df_extracted_features["composite_id"].str.split("--", expand=True)
df_extracted_features["sample_name"] = df_extracted_features["sample_name"].astype(str)
df_extracted_features["TrackID"] = df_extracted_features["TrackID"].astype(int)

umap_features = df_extracted_features.drop(columns=["composite_id", "sample_name", "TrackID"])
pca = PCA(n_components=20)
X_pca = pca.fit_transform(umap_features)

umap_model = umap.UMAP(
        n_components=2, 
        n_neighbors=30, 
        min_dist=0.1, 
        init="random", 
        random_state=seed,
        # metric="precomputed", 
        )

# umap_embedding = umap_model.fit_transform(umap_features)
umap_embedding = umap_model.fit_transform(X_pca)
# df_extracted_features = df_extracted_features.reset_index().rename(columns={"index": "composite_id"})
df_extracted_features[["sample_name","TrackID"]] = df_extracted_features["composite_id"].str.split("--", expand=True)
df_extracted_features["sample_name"] = df_extracted_features["sample_name"].astype(str)
df_extracted_features["TrackID"] = df_extracted_features["TrackID"].astype(int)

umap_embedding = pd.DataFrame(umap_embedding, columns=['UMAP1', 'UMAP2'])

df_umap = umap_embedding
df_umap["sample_name"] = df_extracted_features["sample_name"]
df_umap["TrackID"] = df_extracted_features["TrackID"]
df_umap = pd.merge(df_umap, df_extracted_features, how="left", on=["TrackID", "sample_name"])

sns.scatterplot(
            data=df_umap,
            x="UMAP1",
            y="UMAP2",
            # hue=feature,
            s=40,
            alpha=0.95,
            palette="viridis",
            # legend=False
        )

summarized_features = df_tracks_tsfresh.groupby(["composite_id"])[features].mean().reset_index()

df_umap_plot = pd.merge(df_umap, summarized_features, how="left", on=["composite_id"])
for feature in features:
    sns.scatterplot(
            data=df_umap_plot,
            x="UMAP1",
            y="UMAP2",
            hue=feature,
            s=40,
            alpha=0.95,
            palette="viridis",
            # legend=False
        )
    plt.show()


# --- Clustering ---
# Option A: KMeans
# n_clusters = 4  # you can tune this
# kmeans = KMeans(n_clusters=n_clusters, random_state=seed, n_init="auto")
# cluster_labels = kmeans.fit_predict(X_pca)
# df_analysis["cluster_label_kmeans"] = cluster_labels


# # Option B: (alternative) HDBSCAN (for variable density clusters)
# clusterer = hdbscan.HDBSCAN(min_cluster_size=50, metric='euclidean')
# cluster_labels = clusterer.fit_predict(X_pca)
# df_analysis["cluster_label_hdbscan"] = cluster_labels

# --- UMAP ---
# umap_model = umap.UMAP(
#             n_components=2, 
#             n_neighbors=100, 
#             min_dist=0.1, 
#             # init="random", 
#             metric = "cosine",
#             random_state=seed,
#             )

# umap_embedding = umap_model.fit_transform(X_pca)
# df_analysis["UMAP1"] = umap_embedding[:,0]
# df_analysis["UMAP2"] = umap_embedding[:,1]
    