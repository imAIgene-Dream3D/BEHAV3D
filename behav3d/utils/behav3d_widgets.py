# behav3d_widgets.py
import ipywidgets as widgets
from ipyfilechooser import FileChooser
from IPython.display import display, clear_output
from behav3d.utils import load_behav3d_metadata, check_behav3d_metadata
from behav3d.preprocessing import convert_input_files_to_zarr
import builtins
from pathlib import Path
from traitlets import Any, Unicode, Bool
import os
##Define widgets
class PathPicker(widgets.VBox):
    """
    Displayable widget for picking a file or directory.
    - mode: 'file' or 'dir'
    - start_dir: starting directory for the chooser
    - default: default filename (file mode) or starting path (dir mode)
    - filter_pattern: e.g. '*.csv' (file mode only)
    """
    def __init__(
        self,
        mode='file',
        start_dir='.',
        default='',
        description='Path:',
        button_description='Browse…',
        placeholder='Type a path or click Browse…',
        filter_pattern=None,
        description_width='90px',
        width='100%',
    ):
        assert mode in ('file', 'dir')
        self._mode = mode

        # Text box (the only persistent UI)
        self.text = widgets.Textarea(
            value=str(default or ''),
            description=description,
            placeholder=placeholder,
            style={'description_width': description_width},
            layout=widgets.Layout(width=width, height='32px'),
        )

        # Browse button
        self.button = widgets.Button(
            description=button_description, icon='folder-open',
            tooltip=f"Choose a {'folder' if mode=='dir' else 'file'}"
        )

        # Resolve starting directory & filename
        start_dir = start_dir or '.'
        filename = ''
        if mode == 'file' and default:
            p = Path(default)
            if p.is_absolute():
                start_dir = str(p.parent if p.parent.exists() else start_dir)
                filename = p.name
            else:
                filename = str(default)
        elif mode == 'dir' and default and os.path.isdir(default):
            start_dir = str(default)

        # FileChooser (hidden until button click)
        self.fc = FileChooser(path=start_dir, filename=filename)
        self.fc.title = f"<b>Select {'directory' if mode=='dir' else 'file'}</b>"
        self.fc.show_only_dirs = (mode == 'dir')
        self.fc.use_dir_icons = True
        if filter_pattern and mode == 'file':
            self.fc.filter_pattern = filter_pattern

        # Container: row with Text + Button (chooser injected below when needed)
        self._row = widgets.HBox([self.text, self.button])
        super().__init__([self._row])

        # Wire events
        self.button.on_click(self._toggle_chooser)
        self.fc.register_callback(self._on_select)
        self.text.observe(self._on_text_change, names='value')

    # ---- public API ----
    @property
    def value(self) -> str:
        """Current selected path (string)."""
        return self.text.value

    @value.setter
    def value(self, path: str):
        """Set the path programmatically (updates chooser location if possible)."""
        self.text.value = str(path or '')

    @property
    def mode(self) -> str:
        return self._mode

    # ---- internal helpers ----
    def _toggle_chooser(self, _=None):
        if self.fc in self.children:
            # hide
            self.children = tuple(c for c in self.children if c is not self.fc)
        else:
            # show
            self.children = tuple(list(self.children) + [self.fc])

    def _on_select(self, chooser):
        # Fill text and hide chooser
        if self._mode == 'dir':
            self.text.value = chooser.selected_path
        else:
            self.text.value = chooser.selected
        if self.fc in self.children:
            self.children = tuple(c for c in self.children if c is not self.fc)

    def _on_text_change(self, change):
        # Keep chooser in sync when user types/pastes a path
        newv = (change.get('new') or '').strip()
        if not newv:
            return
        try:
            if self._mode == 'dir':
                if os.path.isdir(newv):
                    self.fc.reset(path=newv)
            else:
                # file mode: split into directory + filename
                parent = os.path.dirname(newv) or '.'
                fname = os.path.basename(newv)
                if os.path.isdir(parent):
                    if fname:
                        self.fc.reset(path=parent, filename=fname)
                    else:
                        self.fc.reset(path=parent)
        except Exception:
            # Don't crash UI if reset fails; ignore silently
            pass
class MetadataLoader(widgets.VBox):
    """
    A VBox widget that wraps a file PathPicker and a 'Load' button.
    After clicking the button, `value` holds the loaded pandas DataFrame,
    and `path` holds the CSV path. Works without globals; read .value later.

    Parameters
    ----------
    file_picker : a widget with a `.value` path (e.g., your PathPicker)
        If None, you can pass one later via set_file_picker().
    button_description : str
        Button label text.
    """
    _busy = Bool(default_value=False).tag(sync=False)
    _handler_bound = Bool(default_value=False).tag(sync=False)
    
    def __init__(
        self,
        metadata_path_picker,
        output_dir_picker,
        button_description: str = "Load metadata",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.file_picker = metadata_path_picker
        self.output_dir_picker = output_dir_picker
        
        self.button = widgets.Button(description=button_description, button_style='success')
        self.out = widgets.Output()

        # self.button.on_click(self._on_click)
        self._click_handler = self._on_click
        self.button.on_click(self._click_handler)
        # Build UI (file_picker may be None initially)
        children = []
        if self.output_dir_picker is not None:
            children.append(self.output_dir_picker)
        if self.file_picker is not None:
            children.append(self.file_picker)  
        children += [self.button, self.out]
        self.children = tuple(children)

   
    def set_file_picker(self, file_picker: widgets.Widget):
        """Attach/replace the file picker widget (must have a `.value` path)."""
        self.file_picker = file_picker
        # Rebuild children with the new picker at the top
        # self.children = (self.file_picker, self.button, self.out)

    def load(self, path = None):
        """Programmatic load (same as clicking the button)."""
        if path is None:
            if self.file_picker is None:
                raise ValueError("No file_picker attached and no path provided.")
            path = self.file_picker.value

        with self.out:
            clear_output(wait=True)
                
            if not path or not str(path).lower().endswith(".csv"):
                print("Please choose a .csv file.")
                return

            self.metadata = load_behav3d_metadata(path)
            self.output_dir = self.output_dir_picker.value
            self.metadata_csv_path = str(path)
            
            check_behav3d_metadata(self.metadata)

            print("✅ Checks passed!")
            display(self.metadata)
           
            
            # import __main__
            # global_vars = __main__.__dict__
            # global_vars['output_dir'] = self.output_dir_picker.value
            # global_vars['metadata_csv_path'] = self.path
            # global_vars['metadata'] = self.metadata
   
    # Button handler
    def _on_click(self, _):
        if self._busy:
            return  # re-entrancy guard prevents double execution
        self._busy = True
        try:
            self.load()
        except Exception as e:
            with self.out:
                print(f"❌ Error: {e}")
        finally:
            self._busy = False

def convert_zarr_button(
        MetadataLoader
    ):
    btn = widgets.Button(
    description="Convert to Zarr",
    button_style="success",  # 'success', 'info', 'warning', 'danger' or ''
    icon="cogs"
    )
    out = widgets.Output()

    def _run_conversion(_):
        with out:
            out.clear_output()
            print("Starting conversion…")
            btn.disabled = True
            try:
                result = convert_input_files_to_zarr(
                    metadata=MetadataLoader.metadata,
                    output_dir=MetadataLoader.output_dir
                    )  # add args here if needed
                print("Done ✅")
            except Exception as e:
                import traceback; traceback.print_exc()
            finally:
                btn.disabled = False

    btn.on_click(_run_conversion)
    display(widgets.VBox([btn, out]))

def set_dim_order(
        MetadataLoader
    ):
    manual_dim_order_widget = widgets.Dropdown(
        options=dim_order_options,
        value=('T', 'C', 'Z', 'Y', 'X'),  # default
        description='Select Order:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='300px')
    )
    out = widgets.Output()

    def _run_conversion(_):
        with out:
            out.clear_output()
            print("Starting conversion…")
            btn.disabled = True
            try:
                result = convert_input_files_to_zarr(
                    metadata=MetadataLoader.metadata,
                    output_dir=MetadataLoader.output_dir
                    )  # add args here if needed
                print("Done ✅")
            except Exception as e:
                import traceback; traceback.print_exc()
            finally:
                btn.disabled = False

    btn.on_click(_run_conversion)
    display(widgets.VBox([btn, out]))
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

# __all__ = [
#     "output_dir_widget", "metadata_path_widget", "load_button1", "output_area1",
#     "manual_dim_order_widget", "manual_dim_order",
#     "output_dir", "metadata_csv_path", "metadata",
#     "exp_duration_widget", "min_track_length_widget", "max_track_length_widget",
#     "set_params_button", "output_area_params",
#     "checkboxes", "checkboxes_box", "save_button2", "output_area2",
#     "umap_distance_widget", "umap_neighbors_widget", "num_clusters_widget",
#     "save_button3", "output_area3","on_load_clicked1",
#     "on_set_params_clicked",
#     "on_save_clicked",
#     "on_save_clicked3",
#     "on_dim_order_change",
#     "on_select_clicked4",
#     "sample_dropdown",
#     "select_button4",
#     "sample_output",
#     "sample_to_backproject",
#     "setup_sample_selection_widgets"
# ]