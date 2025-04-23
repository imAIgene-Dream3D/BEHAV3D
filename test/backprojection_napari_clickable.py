### !!!! COloring the labels based on the features in the table does not work in napari for 3D rendering
### Go 
import numpy as np
import pandas as pd
import zarr
import napari
from qtpy.QtWidgets import QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem, QComboBox
from napari.layers import Labels
from napari.viewer import Viewer
import matplotlib.pyplot as plt

# --- File paths ---
seg_path = 'D:/BHVD_BEHAV3D/BEHAV3D_python/test/images/SMI_JM1_Exp004_Img14/SMI_JM1_Exp004_Img14_organoid_tracked.zarr'
csv_path = 'D:/BHVD_BEHAV3D/BEHAV3D_python/test/analysis/organoid/track_features/BEHAV3D_organoid_combined_track_features.csv'

# --- Load segmentation from Zarr (T, C, Z, Y, X) ---
seg_zarr = zarr.open(seg_path, mode='r')
seg = np.array(seg_zarr)  # shape: (T, C, Z, Y, X)
print(f"Segmentation shape: {seg.shape}")

# --- Load track table ---
df = pd.read_csv(csv_path)

# Filter out rows with NaNs in t or TrackID (just in case)
# df = df.dropna(subset=['position_t', 'TrackID'])

# Ensure integers
df['position_t'] = df['position_t'].astype(int)
df['TrackID'] = df['TrackID'].astype(int)

# --- Create unique label ID for each (t, TrackID) ---

df['label_id'] = df['position_t'] * 10000 + df['TrackID']

# --- Relabel segmentation ---
seg_relabel = np.zeros_like(seg, dtype=np.int32)


# for t, stack_img in enumerate(seg):
seg.shape
for _, row in df.iterrows():
    t = row['position_t']
    track_id = row['TrackID']
    new_label = row['label_id']

    mask = seg[t] == track_id
    seg_relabel[t][mask] = new_label
    
# np.unique(seg_relabel[1])
# df[df["label_id"]==2990008]
# --- Extract features to attach to labels ---
feature_cols = [col for col in df.columns if col not in [
    'position_t', 
    'SegmentID',
    'position_x',
    'position_y',
    'position_z',
    'pixel_position_x',
    'pixel_position_y',
    'pixel_position_z',
    'relative_time',
    'time',
    'distance_unit',
    'time_unit',
    'interpolated',
    'touching_tcells',
    'touching_organoids',
     'sample_name',
     'label_id'
    ]
]
features_df = df.set_index('label_id')[feature_cols]
features_df.index = features_df.index.astype(np.int32)
# features_df.to_dict()
if 0 not in features_df.index:
    dummy_row = pd.DataFrame({col: [np.nan] for col in features_df.columns}, index=[0])
    features_df = pd.concat([dummy_row, features_df])

# print(features_df.loc[0])
# --- Load into napari ---


class FeatureTableWidget(QWidget):
    def __init__(self, layer: Labels, features_df: pd.DataFrame):
        super().__init__()
        self.layer = layer
        self.features_df = features_df
        self.table = QTableWidget()
        layout = QVBoxLayout()
        layout.addWidget(self.table)
        self.setLayout(layout)

        # Connect label selection event
        self.layer.mouse_drag_callbacks.append(self.update_table_from_click)

    def update_table_from_click(self, layer, event):
        position = tuple(int(round(p)) for p in event.position)
        # coords = tuple(int(round(p)) for p in position[::-1])  # reverse to (X, Y, Z...) order
        print(f"Clicked at: {position}")

        # Get the full coord based on number of dimensions
        label_val = layer.data[position]
        print(f"Clicked on label: {label_val}")
        if label_val in self.features_df.index:
            row_data = self.features_df.loc[label_val]
            self.populate_table(row_data)
        else:
            print("Label not found in feature table.")
            self.table.clearContents()
            self.table.setRowCount(0)
            self.table.setColumnCount(0)

    def populate_table(self, row_data):
        if isinstance(row_data, pd.Series):
            self.table.setRowCount(len(row_data))
            self.table.setColumnCount(2)
            self.table.setHorizontalHeaderLabels(['Feature', 'Value'])
            for i, (key, value) in enumerate(row_data.items()):
                self.table.setItem(i, 0, QTableWidgetItem(str(key)))
                self.table.setItem(i, 1, QTableWidgetItem(str(value)))
            
                
viewer = napari.Viewer()

layer = viewer.add_labels(
    seg_relabel,
    name="Tracked Labels",
    properties=features_df.reset_index(),
    features=features_df.reset_index(),
    # color={col: 'viridis' for col in features_df.columns},  # allows dropdown
)

# Create a ComboBox for feature selection
feature_combo = QComboBox()
feature_combo.addItems(features_df.columns.tolist())  # Add feature columns to the combo box

def update_colors(selected_feature):
    label_ids = features_df.index.values
    df_selected_feature = features_df[selected_feature]
    df_selected_feature = df_selected_feature.dropna()
    
    if selected_feature in ['TrackID', 'label_id']:
        track_ids = df_selected_feature.values
        unique_track_ids = np.unique(track_ids)

        # Create a color per unique TrackID
        cmap = plt.get_cmap('tab20')
        rgba = cmap(np.linspace(0, 1, len(unique_track_ids)))
        track_to_color = {track: tuple(color[:3]) for track, color in zip(unique_track_ids, rgba)}

        # Map each label_id to its track's color
        color_map = {
            int(label_id): track_to_color[track_id]
            for label_id, track_id in zip(label_ids, track_ids)
        }
    else:
        values = df_selected_feature.values.astype(np.float32)
        normed = (values - values.min()) / (values.max() - values.min() + 1e-6)  # avoid divide by 0
        cmap = plt.get_cmap('viridis')
        rgba = cmap(normed)
        color_map = {int(label): tuple(color[:3]) for label, color in zip(label_ids, rgba)}

    # Add background color (label 0)
    # Add background color for 0 label
    color_map[0] = (0, 0, 0)
    layer._set_colormap(color_map)
    # Update the layer
    layer.color = color_map
    layer.color_mode = "direct"  # this must be set for the color dict to apply

    # Trigger redraw
    layer.refresh()
    layer.update()
    layer._update_thumbnail()
    layer.visible = False
    layer.visible = True
    layer.data = layer.data.copy()

    print("TrackID color mapping applied.")

def on_change():
    update_colors(feature_combo.currentText())
# Connect the combo box to the function
feature_combo.currentTextChanged.connect(on_change)

# Create a QWidget for the drop-down menu and layout
menu_widget = QWidget()
layout = QVBoxLayout()
layout.addWidget(feature_combo)
menu_widget.setLayout(layout)

# Add the widget to the viewer's GUI
viewer.window.add_dock_widget(menu_widget, area='right')
table_widget = FeatureTableWidget(layer, features_df)
viewer.window.add_dock_widget(table_widget, area='right', name='Track Features')
napari.run()


# viewer.window.add_dock_widget(layer.label_properties, area='right')
# @layer.mouse_move_callbacks.append
# def show_tooltip(layer, event):
#     try:
#         # Get mouse position in data coordinates
#         position = tuple(int(round(p)) for p in event.position)
#         # Safety check for shape
#         if all(0 <= p < s for p, s in zip(position, layer.data.shape)):
#             label_val = layer.data[position]
#             if label_val in features_df.index:
#                 row = features_df.loc[label_val]
#                 print(f"Hovering on label {label_val}: {row.to_dict()}")
#         else:
#             print("Mouse outside data bounds.")
#     except Exception as e:
#         print(f"Tooltip error: {e}")

# viewer = napari.Viewer()

# layer = viewer.add_labels(
#     seg,
# )
# layer = viewer.add_image(
#     raw,
# )
# napari.run()

