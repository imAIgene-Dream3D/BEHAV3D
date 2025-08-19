# behav3d_widgets.py
import ipywidgets as widgets
from IPython.display import display, clear_output
from behav3d.utils import load_behav3d_metadata, check_behav3d_metadata
import builtins


##Define widgets
#Widnget for output dir and metadata
output_dir_widget = widgets.Text(
    value=r"D:/SURF_2/Scripts/Python/BEHAV3D3_0/output",
    description='Output Dir:', style={'description_width': '70px'},
    layout=widgets.Layout(width='100%')
)

metadata_path_widget = widgets.Text(
    value=r"D:/SURF_2/Scripts/Python/BEHAV3D3_0/metadata.csv",
    description='Metadata:',
    layout=widgets.Layout(width='100%')
)

load_button1 = widgets.Button(description="Load Metadata", button_style='success')
output_area1 = widgets.Output()
output_dir =None
metadata_csv_path = None
metadata = None

def on_load_clicked1(b):
    global_vars = globals()  # get the global namespace of this module
    # or better, assign to the notebook's global namespace explicitly:
    import __main__
    global_vars = __main__.__dict__

    output_dir_val = output_dir_widget.value
    metadata_csv_path_val = metadata_path_widget.value

    metadata_val = load_behav3d_metadata(metadata_csv_path_val)

    with output_area1:
        clear_output()
        print(f"Loading metadata from:\nOutput Dir: {output_dir_val}\nMetadata CSV: {metadata_csv_path_val}")
        display(metadata_val)
        print("Checking metadata...")
        check_behav3d_metadata(metadata_val)
        print("✅ Metadata checks passed!")

    # Update notebook globals directly:
    global_vars['output_dir'] = output_dir_val
    global_vars['metadata_csv_path'] = metadata_csv_path_val
    global_vars['metadata'] = metadata_val


#widget for dimensions
dim_order_options = [
    ('T, C, Z, Y, X', ('T', 'C', 'Z', 'Y', 'X')),
    ('T, Z, C, Y, X', ('T', 'Z', 'C', 'Y', 'X')),
    ('Z, C, T, Y, X', ('Z', 'C', 'T', 'Y', 'X')),
    ('Z, T, C, Y, X', ('Z', 'T', 'C', 'Y', 'X')),
    ('C, Z, T, Y, X', ('C', 'Z', 'T', 'Y', 'X')),
    ('C, T, Z, Y, X', ('C', 'T', 'Z', 'Y', 'X')),
]

# Dropdown widget with labels and values
manual_dim_order_widget = widgets.Dropdown(
    options=dim_order_options,
    value=('T', 'C', 'Z', 'Y', 'X'),  # default
    description='Select Order:',
    style={'description_width': 'initial'},
    layout=widgets.Layout(width='300px')
)

# Variable to store the selected order
manual_dim_order = manual_dim_order_widget.value

# Handler to update variable on change
def on_dim_order_change(change):
    new_val = change['new']

    import __main__
    global_vars = __main__.__dict__

    # Update local variable (optional)
    global manual_dim_order
    manual_dim_order = new_val

    # Update notebook global
    global_vars['manual_dim_order'] = new_val

    print(f"Selected manual_dim_order: {new_val}")


# Widgets for parameters
exp_duration_widget = widgets.IntText(
    value=1000,  # or your default experiment duration
    description='Exp Duration:',
    style={'description_width': '110px'},
    layout=widgets.Layout(width='300px')
)

min_track_length_widget = widgets.IntSlider(
    value=100, min=0, max=1000, step=1,
    description='Min Track Length:',
    style={'description_width': '130px'},
    layout=widgets.Layout(width='300px')
)

max_track_length_widget = widgets.IntSlider(
    value=500, min=0, max=1000, step=1,
    description='Max Track Length:',
    style={'description_width': '130px'},
    layout=widgets.Layout(width='300px')
)

set_params_button = widgets.Button(description="Set Parameters", button_style='info')
output_area_params = widgets.Output()

def on_set_params_clicked(b):
    exp_duration_val = exp_duration_widget.value
    tcell_min_track_length_val = min_track_length_widget.value
    tcell_max_track_length_val = max_track_length_widget.value

    with output_area_params:
        clear_output()
        print(f"Parameters set:")
        print(f"  Experiment Duration: {exp_duration_val}")
        print(f"  Min Track Length: {tcell_min_track_length_val}")
        print(f"  Max Track Length: {tcell_max_track_length_val}")

    import __main__
    global_vars = __main__.__dict__
    global_vars['exp_duration'] = exp_duration_val
    global_vars['tcell_min_track_length'] = tcell_min_track_length_val
    global_vars['tcell_max_track_length'] = tcell_max_track_length_val


## Widgets to select the parameters for classification
all_features = [
    'nr_pixels', 'volume', 'bbox_volume', 'extent', 'solidity', 'equivalent_diameter', 
    'major_axis_length', 'minor_axis_length', 'elongation', 'surface_area', 'sphericity',
    'convex_volume', 'orientation_vector', 'mean_intensity_ch1', 'mean_intensity_ch2', 
    'mean_dead_dye', 'mean_intensity_ch4', 'percentage_dead_mask', 'nr_dead_mask_pixels', 
    'increase_dead_mask', 'organoid_contact', 'organoid_contact_pixels', 'touching_organoids',
    'tcell_contact', 'tcell_contact_pixels', 'touching_tcells', 'interpolated', 'displacement',
    'cumulative_displacement', 'displacement_from_origin', 'mean_square_displacement', 'speed',
    'mean_speed', 'dead', 'active_tcell_contact', 'z_mean_square_displacement', 'z_speed', 
    'z_mean_dead_dye'
]

default_selected = {
    "z_mean_square_displacement", 
    "z_speed", 
    "z_mean_dead_dye", 
    "tcell_contact", 
    "organoid_contact"
}

checkboxes = []
for feature in all_features:
    cb = widgets.Checkbox(
        value=(feature in default_selected),
        description=feature,
        indent=False,
        layout=widgets.Layout(width='300px')
    )
    checkboxes.append(cb)

checkboxes_box = widgets.VBox(checkboxes, layout=widgets.Layout(
    overflow='auto',
    height='300px',
    border='1px solid gray',
    padding='5px'
))

save_button2 = widgets.Button(description="Save DTW Features", button_style='success', layout=widgets.Layout(width='200px'))
output_area2 = widgets.Output()

def on_save_clicked(b):
    selected_features = [cb.description for cb in checkboxes if cb.value]

    with output_area2:
        clear_output()
        print("Selected DTW features:")
        for f in selected_features:
            print(f"- {f}")

    import __main__
    global_vars = __main__.__dict__
    global_vars['dtw_features'] = selected_features

# Widgets for UMAP and clustering parameters
umap_distance_widget = widgets.FloatSlider(
    value=0.1, min=0.0, max=1.0, step=0.05, 
    style={'description_width': '64px'}, description='UMAP dist:'
)

umap_neighbors_widget = widgets.IntSlider(
    value=5, min=2, max=50, step=1, 
    style={'description_width': '101px'}, description='UMAP neighbors:'
)

num_clusters_widget = widgets.IntSlider(
    value=6, min=2, max=20, step=1, 
    style={'description_width': '80px'}, description='Num Clusters:'
)

save_button3 = widgets.Button(description="UMAP&Clustering Params", button_style='success',layout=widgets.Layout(width='250px'))
output_area3 = widgets.Output()

def on_save_clicked3(b):
    umap_dist_val = umap_distance_widget.value
    umap_neighbors_val = umap_neighbors_widget.value
    nr_clusters_val = num_clusters_widget.value

    with output_area3:
        clear_output()
        print(f"Saved parameters:\nUMAP distance: {umap_dist_val}\nUMAP neighbors: {umap_neighbors_val}\nNumber of clusters: {nr_clusters_val}")

    import __main__
    global_vars = __main__.__dict__
    global_vars['umap_minimal_distance'] = umap_dist_val
    global_vars['umap_n_neighbors'] = umap_neighbors_val
    global_vars['nr_of_clusters'] = nr_clusters_val


# Define placeholders
# Create widgets once, globally
sample_dropdown = widgets.Dropdown(description='Sample:', layout=widgets.Layout(width='50%'), style={'description_width': 'initial'})
select_button4 = widgets.Button(description="Select Sample", button_style='success', layout=widgets.Layout(width='30%'))
sample_output = widgets.Output()

sample_to_backproject = None  # will hold selected sample name

def setup_sample_selection_widgets(metadata):
    # Update dropdown options based on metadata
    sample_names = metadata["sample_name"].unique().tolist()
    sample_dropdown.options = sample_names
    if sample_names:
        sample_dropdown.value = sample_names[0]

def on_select_clicked4(b):
    import __main__
    global sample_to_backproject
    sample_to_backproject = sample_dropdown.value
    __main__.sample_to_backproject = sample_to_backproject  # push to notebook globals

    with sample_output:
        clear_output()
        print(f"✅ Sample selected for backprojection: {sample_to_backproject}")

# Attach callback once
select_button4.on_click(on_select_clicked4)
# --------------------

__all__ = [
    "output_dir_widget", "metadata_path_widget", "load_button1", "output_area1",
    "manual_dim_order_widget", "manual_dim_order",
    "output_dir", "metadata_csv_path", "metadata",
    "exp_duration_widget", "min_track_length_widget", "max_track_length_widget",
    "set_params_button", "output_area_params",
    "checkboxes", "checkboxes_box", "save_button2", "output_area2",
    "umap_distance_widget", "umap_neighbors_widget", "num_clusters_widget",
    "save_button3", "output_area3","on_load_clicked1",
    "on_set_params_clicked",
    "on_save_clicked",
    "on_save_clicked3",
    "on_dim_order_change",
    "on_select_clicked4",
    "sample_dropdown",
    "select_button4",
    "sample_output",
    "sample_to_backproject",
    "setup_sample_selection_widgets"
]