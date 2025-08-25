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
from behav3d.utils.fileio import load_image, get_image_shape, get_image_dimension_order
import asyncio
import pandas as pd
from behav3d.preprocessing.segmentation.napari_pixelclassifier import train_pixel_classifier, run_pixel_classifier_segmentation
import traceback
from behav3d.preprocessing.tracking import visualize_tracks
import json
from copy import deepcopy

# ===============================
# JSON CONFIG (loader + defaults)
# ===============================
CONFIG_PATH = Path("behav3d_config.json")  # adjust if you prefer a different location

_DEFAULT_CONFIG = {
    "paths": {
        "metadata_csv": r"/Volumes/T7_Sam//BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata.csv",
        "output_dir": r"/Volumes/T7_Sam//BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/"
    },
    "dim_order": {
        "default_apply_all": "TCZYX"
    },
    "pixel_classifier": {
        "examples_per_sample": 3,
        "sample_specific_classifier": False,
        "workers": 8,
        "organoid_edt_threshold": 12.0,
        "use_all_timepoints": True,
        "tp_start": 0,
        "tp_end": 0
    },
    "tracking": {
        "method": "lap",
        "overwrite": False,
        "lap": {
            "track_cost_px": 45,
            "gap_close_cost_px": 60,
            "gap_close_max_frames": 3,
            "merging_cost_px": 0,
            "splitting_cost_px": 0
        },
        "trackpy": {
            "search_range_px": 31,
            "memory_frames": 2,
            "adaptive_stop": 10.0,
            "adaptive_step": 0.95
        }
    },
    "tracking_visualization": {
        "sample_name": None,
        "use_range": False,
        "start_t": 0,
        "end_t": 0,
        "channel_colors": ["cyan", "yellow", "red", "green", "magenta", "blue"]
    }
}

def _deep_merge(a: dict, b: dict) -> dict:
    out = deepcopy(a)
    for k, v in (b or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out

def _load_config(path: Path = CONFIG_PATH) -> dict:
    cfg = deepcopy(_DEFAULT_CONFIG)
    if path.exists():
        try:
            user = json.loads(path.read_text())
            cfg = _deep_merge(cfg, user)
        except Exception:
            # ignore malformed files; keep defaults
            pass
    return cfg

def _cfg_get(cfg: dict, dotted_key: str, default=None):
    cur = cfg
    for part in dotted_key.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur

CFG = _load_config()

def _mk_timepoint_range(use_all: bool, start: int, end: int):
    if use_all:
        return None
    return (int(start), int(end))
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

        # Prefer JSON defaults if caller didn't provide one
        if not default:
            if mode == 'file':
                default = _cfg_get(CFG, "paths.metadata_csv", "") or ""
            else:  # dir
                default = _cfg_get(CFG, "paths.output_dir", "") or ""

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
        self.metadata = None
        self.output_dir = None
        self.metadata_csv_path = None
            
        self.file_picker = metadata_path_picker
        self.output_dir_picker = output_dir_picker

        # initialize pickers with JSON defaults if empty
        try:
            if hasattr(self.file_picker, "value") and not self.file_picker.value:
                self.file_picker.value = _cfg_get(CFG, "paths.metadata_csv", "") or ""
        except Exception:
            pass
        try:
            if hasattr(self.output_dir_picker, "value") and not self.output_dir_picker.value:
                self.output_dir_picker.value = _cfg_get(CFG, "paths.output_dir", "") or ""
        except Exception:
            pass
        
        self.button = widgets.Button(description=button_description, button_style='success')
        self.out = widgets.Output()

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
            self.output_dir = self.output_dir_picker.value if self.output_dir_picker is not None else ""
            self.metadata_csv_path = str(path)
            
            check_behav3d_metadata(self.metadata)
            print("✅ Checks passed!")
            display(self.metadata)
   
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

def convert_zarr_button(metadata_loader, dim_order_widget):
    btn = widgets.Button(description="Convert to Zarr", button_style="success", icon="cogs")
    out = widgets.Output()

    def _run_conversion(_):
        with out:
            out.clear_output()
            print("Starting conversion…")
            btn.disabled = True
            # make sure current dim-order choices are written first (optional but handy)
            try:
                dim_order_widget.write_dimorder_to_metadata()
            except Exception:
                pass

            # ensure `result` is always defined even if conversion errors
            result = metadata_loader.metadata
            try:
                result = convert_input_files_to_zarr(
                    metadata=metadata_loader.metadata,
                    output_dir=metadata_loader.output_dir
                )
            except Exception:
                import traceback; traceback.print_exc()
            finally:
                metadata_loader.metadata = result
                metadata_loader.metadata.to_csv(metadata_loader.metadata_csv_path, index=False)
                print("Done ✅")
                btn.disabled = False

    btn.on_click(_run_conversion)
    # ⬇️ return the widget; do not display here
    return widgets.VBox([btn, out])
class DimOrderTable:
    DIM_ORDER_OPTIONS = [
        "TCZYX",
        "TZCYX",
        "ZCTYX",
        "ZTCYX",
        "CZTYX",
        "CTZYX",
    ]
    DEFAULT_ORDER = "TCZYX"

    def __init__(
        self,
        metadata_loader,
        sample_col="sample_name",
        path_col="raw_image_path",
        widths=("40%","30%","30%"),
        dim_col="dimension_order",
        auto_write=False
    ):
        self.metadata_loader = metadata_loader
        self.sample_col = sample_col
        self.path_col = path_col
        self.dim_col = dim_col
        self.auto_write = auto_write

        self._width_sample, self._width_shape, self._width_order = widths

        # State
        self._rows = []
        self._selections = {}

        # Precompute allowed tuples
        self._allowed_orders = self.DIM_ORDER_OPTIONS

        # UI skeleton
        self._status = widgets.HTML("<b>Waiting for user to load metadata…</b>")
        self._refresh_btn = widgets.Button(description="Refresh", tooltip="Build/Update table from metadata")
        self._refresh_btn.on_click(self._on_refresh)

        _apply_default = str(_cfg_get(CFG, "dim_order.default_apply_all", self.DEFAULT_ORDER))
        if _apply_default not in self.DIM_ORDER_OPTIONS:
            _apply_default = self.DEFAULT_ORDER

        self._apply_all_dd = widgets.Dropdown(
            options=self.DIM_ORDER_OPTIONS,
            value=_apply_default,
            description="Apply to all:",
            style={'description_width': 'initial'},
            layout=widgets.Layout(width="350px", margin="6px 8px 6px 0")
        )
        self._apply_all_btn = widgets.Button(description="Apply", layout=widgets.Layout(width="100px"))
        self._apply_all_btn.on_click(lambda _: self.set_all(self._apply_all_dd.value))

        self._header = widgets.HBox([
            widgets.Label("Sample", layout=widgets.Layout(width=self._width_sample, overflow="hidden")),
            widgets.Label("Image shape", layout=widgets.Layout(width=self._width_shape, overflow="hidden")),
            widgets.Label("Dim order", layout=widgets.Layout(width=self._width_order, overflow="hidden")),
        ])
        self._table_body = widgets.VBox()
        self._out = widgets.Output()

        self.widget = widgets.VBox([
            self._status,
            widgets.HBox([self._refresh_btn]),
            self._header,
            self._table_body,
            widgets.HBox([self._apply_all_dd, self._apply_all_btn]),
            self._out
        ])

        # NEW: auto-load immediately if metadata already exists
        self._maybe_autoload()

    # NEW: helper that only shows “waiting” if metadata is None; otherwise loads
    def _maybe_autoload(self):
        df = getattr(self.metadata_loader, "metadata", None)
        if df is None:
            self._status.value = "<b>Waiting for user to load metadata…</b>"
            return
        if not isinstance(df, pd.DataFrame):
            self._status.value = "<b>metadata_loader.metadata must be a pandas DataFrame.</b>"
            return
        missing = [c for c in (self.sample_col, self.path_col) if c not in df.columns]
        if missing:
            self._status.value = f"<b>Metadata missing columns: {missing}</b>"
            return
        self._build_rows_from(df)
        self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"

    # UPDATED: only show “waiting” when metadata is None; otherwise load/refresh
    def _on_refresh(self, _btn):
        df = getattr(self.metadata_loader, "metadata", None)
        if df is None:
            self._status.value = "<b>Waiting for user to load metadata…</b>"
            return
        if not isinstance(df, pd.DataFrame):
            self._status.value = "<b>metadata_loader.metadata must be a pandas DataFrame.</b>"
            return

        missing = [c for c in (self.sample_col, self.path_col) if c not in df.columns]
        if missing:
            self._status.value = f"<b>Metadata missing columns: {missing}</b>"
            return

        self._build_rows_from(df)
        self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"
        
    def display(self):
        display(self.widget)

    def get_selections(self) -> dict:
        """Return dict: sample_name -> tuple like ('T','C','Z','Y','X')"""
        return dict(self._selections)

    def get_selections_str(self) -> dict:
        """Return dict: sample_name -> 'TCZYX'"""
        return {s: self._order_to_str(o) for s, o in self._selections.items()}

    def write_dimorder_to_metadata(self, col=None):
        col = col or self.dim_col
        df = getattr(self.metadata_loader, "metadata", None)
        if not isinstance(df, pd.DataFrame):
            raise ValueError("metadata_loader.metadata is not a DataFrame yet.")
        if df.empty:
            raise ValueError("metadata_loader.metadata is empty.")
        if self.sample_col not in df.columns:
            raise ValueError(f"metadata is missing '{self.sample_col}' column.")
        # Map selections by sample name
        str_map = self.get_selections_str()
        self.metadata_loader.metadata[col] = df[self.sample_col].map(str_map)
        return self.metadata_loader.metadata

    def set_all(self, order_tuple):
        """Programmatically set all rows to the given order tuple or string."""
        if isinstance(order_tuple, str):
            order_tuple = tuple(order_tuple)
        for row in self._rows:
            row['dd'].value = order_tuple  # triggers observers

    def refresh_shapes(self):
        """Recompute shape labels (e.g., after mounting a drive)."""
        for row in self._rows:
            row['shape_label'].value = self._probe_image_shape(row['path'])

    def to_dataframe(self):
        data = []
        for row in self._rows:
            tup = row["dd"].value
            data.append({
                "sample_name": row["sample"],
                "raw_image_path": row["path"],
                "image_shape": row["shape_label"].value,
                "dim_order": tup,
                self.dim_col: self._order_to_str(tup),
            })
        return pd.DataFrame(data)

    def _on_refresh(self, _btn):
        df = getattr(self.metadata_loader, "metadata", None)
        if not (isinstance(df, pd.DataFrame) and not df.empty):
            self._status.value = "<b>Waiting for user to load metadata…</b>"
            return
        missing = [c for c in (self.sample_col, self.path_col) if c not in df.columns]
        if missing:
            self._status.value = f"<b>Metadata missing columns: {missing}</b>"
            return

        self._build_rows_from(df)
        self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"

    def _build_rows_from(self, df: pd.DataFrame):
        # Reset state
        self._rows.clear()
        self._selections.clear()
        row_boxes = []

        # Optional: prefill from existing string column if present
        prefill = None
        if self.dim_col in df.columns:
            prefill = df[self.dim_col].astype(str).to_dict()  # index -> string

        for idx, r in df.iterrows():
            sample = str(r[self.sample_col])
            path = r[self.path_col]

            probed_shape = self._probe_image_shape(path)
            sample_lbl = widgets.Label(sample, layout=widgets.Layout(width=self._width_sample, overflow="hidden"))
            shape_lbl  = widgets.Label(probed_shape, layout=widgets.Layout(width=self._width_shape, overflow="hidden"))

            # Determine default dropdown value
            default_val = get_image_dimension_order(path)
            if default_val is None:
                default_val = self.DEFAULT_ORDER
            if prefill is not None:
                s = prefill.get(idx, None)
                if isinstance(s, str):
                    tu = tuple(s.upper())
                    if tu in self._allowed_orders:
                        default_val = tu

            dd = widgets.Dropdown(
                options=self.DIM_ORDER_OPTIONS,
                value=default_val,
                layout=widgets.Layout(width=self._width_order),
            )

            # keep live selection
            self._selections[sample] = dd.value
            dd.observe(lambda ch, s=sample: self._on_dd_changed(s, ch), names='value')

            self._rows.append({"sample": sample, "path": path, "shape_label": shape_lbl, "dd": dd})
            row_boxes.append(widgets.HBox([sample_lbl, shape_lbl, dd]))

        # Swap UI body
        self._table_body.children = tuple(row_boxes)

    def _on_dd_changed(self, sample_name, change):
        if change.get('name') == 'value':
            new_tuple = change['new']
            self._selections[sample_name] = new_tuple
            if self.auto_write:
                df = getattr(self.metadata_loader, "metadata", None)
                if isinstance(df, pd.DataFrame) and self.sample_col in df.columns:
                    mask = df[self.sample_col].astype(str) == str(sample_name)
                    df.loc[mask, self.dim_col] = self._order_to_str(new_tuple)

    def _probe_image_shape(self, path):
        path = Path(path)
        if not path.exists():
            return "⛔ missing"
        try:
            shape = get_image_shape(path)
            return "×".join(map(str, shape)) if shape else "unknown"
        except Exception as e:
            return f"⚠️ {type(e).__name__}"

    @staticmethod
    def _order_to_str(order_tuple):
        return "".join(order_tuple)
class PixelClassifierPanel:
    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        pc = _cfg_get(CFG, "pixel_classifier", {})

        # -------- Train controls --------
        self.examples_per_sample = widgets.IntSlider(
            description="Examples per sample",
            value=int(pc.get("examples_per_sample", 3)),
            min=1, max=50, step=1, continuous_update=False
        )
        self.sample_specific_classifier = widgets.Checkbox(
            description="Sample-specific classifier",
            value=bool(pc.get("sample_specific_classifier", False))
        )
        self.n_workers = widgets.IntSlider(
            description="Workers",
            value=int(pc.get("workers", (os.cpu_count() or 8))),
            min=1, max=max(8, (os.cpu_count() or 8)),
            step=1, continuous_update=False
        )

        # -------- Apply controls --------
        self.organoid_edt_threshold = widgets.FloatSlider(
            description="Organoid EDT thr",
            value=float(pc.get("organoid_edt_threshold", 12.0)),
            min=0.0, max=50.0, step=0.5, readout_format=".1f",
            continuous_update=False,
            style={'description_width': '160px'}
        )
        self.use_all_timepoints = widgets.Checkbox(
            description="Process ALL timepoints",
            value=bool(pc.get("use_all_timepoints", True))
        )
        self.tp_start = widgets.IntText(
            description="Start t",
            value=int(pc.get("tp_start", 0))
        )
        self.tp_end   = widgets.IntText(
            description="End t",
            value=int(pc.get("tp_end", 0))
        )

        # Show/hide start/end boxes based on checkbox
        self.use_all_timepoints.observe(self._toggle_timepoint_inputs, names='value')
        self._toggle_timepoint_inputs()  # initialize

        # Buttons + log
        self.btn_train = widgets.Button(description="Train", button_style="primary")
        self.btn_apply = widgets.Button(description="Apply", button_style="success")
        self.btn_train.on_click(self._on_train_clicked)
        self.btn_apply.on_click(self._on_apply_clicked)
        self.out = widgets.Output()

        # Layout
        train_box = widgets.VBox([
            widgets.HTML("<b>Train pixel classifier</b>"),
            widgets.HBox([self.examples_per_sample, self.n_workers]),
            self.sample_specific_classifier,
            self.btn_train,
        ])

        self.tp_row = widgets.HBox([self.use_all_timepoints, self.tp_start, self.tp_end])

        apply_box = widgets.VBox([
            widgets.HTML("<b>Apply segmentation</b>"),
            self.organoid_edt_threshold,
            self.tp_row,
            self.btn_apply,
        ])

        self.ui = widgets.VBox([train_box, widgets.HTML("<hr>"), apply_box, widgets.HTML("<hr>"), self.out])

    def _toggle_timepoint_inputs(self, change=None):
        show = not self.use_all_timepoints.value
        disp = None if show else 'none'
        self.tp_start.layout.display = disp
        self.tp_end.layout.display   = disp
        self.tp_start.disabled = not show
        self.tp_end.disabled   = not show

    def display(self):
        display(self.ui)

    def _lock(self, state: bool):
        for w in [
            self.btn_train, self.btn_apply,
            self.examples_per_sample, self.sample_specific_classifier, self.n_workers,
            self.organoid_edt_threshold, self.use_all_timepoints, self.tp_start, self.tp_end,
        ]:
            w.disabled = state

    # ---- callbacks ----
    def _on_train_clicked(self, _):
        self._lock(True)
        with self.out:
            self.out.clear_output()
            try:
                odir = Path(self.metadata_loader.output_dir).expanduser()
                odir.mkdir(parents=True, exist_ok=True)

                print("▶️ Training pixel classifier…", flush=True)
                print(f"  output_dir={odir}")
                print(f"  examples_per_sample={self.examples_per_sample.value}")
                print(f"  sample_specific_classifier={self.sample_specific_classifier.value}")
                print(f"  n_workers={self.n_workers.value}")

                train_pixel_classifier(
                    output_dir=str(odir),
                    metadata=self.metadata_loader.metadata,
                    examples_per_sample=int(self.examples_per_sample.value),
                    sample_specific_classifier=bool(self.sample_specific_classifier.value),
                    n_workers=int(self.n_workers.value),
                )
                print("✅ Training finished.", flush=True)
            except Exception:
                import traceback; traceback.print_exc()
            finally:
                self._lock(False)

    def _on_apply_clicked(self, _):
        self._lock(True)
        with self.out:
            self.out.clear_output()
            try:
                odir = Path(self.metadata_loader.output_dir).expanduser()
                odir.mkdir(parents=True, exist_ok=True)

                tpr = _mk_timepoint_range(
                    use_all=bool(self.use_all_timepoints.value),
                    start=int(self.tp_start.value),
                    end=int(self.tp_end.value),
                )

                print("▶️ Applying pixel classifier segmentation…", flush=True)
                print(f"  output_dir={odir}")
                print(f"  organoid_edt_threshold={self.organoid_edt_threshold.value}")
                print(f"  timepoint_range={tpr}", flush=True)

                new_md = run_pixel_classifier_segmentation(
                    output_dir=str(odir),
                    metadata=self.metadata_loader.metadata,
                    organoid_edt_threshold=float(self.organoid_edt_threshold.value),
                    timepoint_range=tpr
                )

                # Update loader + ALWAYS save CSV
                self.metadata_loader.metadata = new_md
                csv_path = Path(self.metadata_loader.metadata_csv_path).expanduser()
                print(f"💾 Saving metadata to: {csv_path}", flush=True)
                new_md.to_csv(csv_path, sep=",", index=False)

                print("✅ Apply finished.", flush=True)
            except Exception:
                import traceback; traceback.print_exc()
            finally:
                self._lock(False)
class TrackingPanel:
    """
    Minimal UI for tracking (LAP, TrackPy, or Propagation).
    - Uses metadata_loader.output_dir and metadata_loader.metadata_csv_path directly
    - Updates metadata_loader.metadata and always saves CSV
    - cell_type is supplied in __init__
    """
    def __init__(self, metadata_loader, cell_type="tcell"):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        tcfg = _cfg_get(CFG, "tracking", {})

        # Method selector
        self.method = widgets.Dropdown(
            description="Tracking method",
            options=[("LAP (laptrack)", "lap"),
                     ("TrackPy", "trackpy"),
                     ("Propagation", "prop")],
            value=str(tcfg.get("method", "lap")),
            style={'description_width': '160px'}
        )

        # Shared controls
        self.overwrite = widgets.Checkbox(
            description="Overwrite existing",
            value=bool(tcfg.get("overwrite", False))
        )

        # LAP params (distances in pixels; squared when calling)
        lap = tcfg.get("lap", {})
        self.track_cost_dist = widgets.IntSlider(
            description="Track cost (px)",
            value=int(lap.get("track_cost_px", 45)),
            min=1, max=200, step=1,
            continuous_update=False, style={'description_width':'160px'}
        )
        self.gap_cost_dist = widgets.IntSlider(
            description="Gap close cost (px)",
            value=int(lap.get("gap_close_cost_px", 60)),
            min=1, max=300, step=1,
            continuous_update=False, style={'description_width':'160px'}
        )
        self.gap_max_frames = widgets.IntSlider(
            description="Gap close max frames",
            value=int(lap.get("gap_close_max_frames", 3)),
            min=0, max=20, step=1,
            continuous_update=False, style={'description_width':'160px'}
        )
        self.merging_cost_dist = widgets.IntText(
            description="Merging cost (px)",
            value=int(lap.get("merging_cost_px", 0)),
            style={'description_width':'160px'}
        )
        self.splitting_cost_dist = widgets.IntText(
            description="Splitting cost (px)",
            value=int(lap.get("splitting_cost_px", 0)),
            style={'description_width':'160px'}
        )
        self.lap_params = widgets.VBox([
            widgets.HTML("<b>LAP parameters</b>"),
            self.track_cost_dist, self.gap_cost_dist, self.gap_max_frames,
            self.merging_cost_dist, self.splitting_cost_dist
        ])

        # TrackPy params
        tpy = tcfg.get("trackpy", {})
        self.tp_search_range = widgets.IntSlider(
            description="Search range (px)",
            value=int(tpy.get("search_range_px", 31)),
            min=1, max=200, step=1,
            continuous_update=False, style={'description_width':'160px'}
        )
        self.tp_memory = widgets.IntSlider(
            description="Memory (frames)",
            value=int(tpy.get("memory_frames", 2)),
            min=0, max=20, step=1,
            continuous_update=False, style={'description_width':'160px'}
        )
        self.tp_adaptive_stop = widgets.FloatSlider(
            description="Adaptive stop",
            value=float(tpy.get("adaptive_stop", 10.0)),
            min=0.0, max=50.0, step=0.5,
            readout_format=".1f", continuous_update=False, style={'description_width':'160px'}
        )
        self.tp_adaptive_step = widgets.FloatSlider(
            description="Adaptive step",
            value=float(tpy.get("adaptive_step", 0.95)),
            min=0.1, max=1.0, step=0.01,
            readout_format=".2f", continuous_update=False, style={'description_width':'160px'}
        )
        self.tp_params = widgets.VBox([
            widgets.HTML("<b>TrackPy parameters</b>"),
            self.tp_search_range, self.tp_memory, self.tp_adaptive_stop, self.tp_adaptive_step
        ])

        # Propagation params (none needed; just an info box)
        self.prop_params = widgets.VBox([
            widgets.HTML("<b>Propagation tracking</b>"),
            widgets.HTML("<i>No tunable parameters.</i>")
        ])

        # Swap param groups on method change
        self.param_container = widgets.VBox()
        self.method.observe(self._toggle_param_groups, names='value')
        self._toggle_param_groups()

        # Run + log
        self.btn_run = widgets.Button(description="Run tracking", button_style="success")
        self.btn_run.on_click(self._on_run_clicked)
        self.out = widgets.Output()

        # Layout
        self.ui = widgets.VBox([
            widgets.HTML("<b>Tracking</b>"),
            self.method,
            widgets.HBox([self.overwrite]),
            self.param_container,
            self.btn_run,
            widgets.HTML("<hr>"),
            self.out
        ])

    def display(self):
        display(self.ui)

    def _toggle_param_groups(self, change=None):
        if self.method.value == "lap":
            self.param_container.children = (self.lap_params,)
        elif self.method.value == "trackpy":
            self.param_container.children = (self.tp_params,)
        else:  # "prop"
            self.param_container.children = (self.prop_params,)

    def _lock(self, state: bool):
        for w in [
            self.method, self.overwrite,
            self.track_cost_dist, self.gap_cost_dist, self.gap_max_frames,
            self.merging_cost_dist, self.splitting_cost_dist,
            self.tp_search_range, self.tp_memory, self.tp_adaptive_stop, self.tp_adaptive_step,
            self.btn_run
        ]:
            if hasattr(w, "disabled"):
                w.disabled = state

    def _on_run_clicked(self, _):
        self._lock(True)
        with self.out:
            self.out.clear_output()
            try:
                out_dir = Path(self.metadata_loader.output_dir).expanduser()
                out_dir.mkdir(parents=True, exist_ok=True)

                csv_path = Path(self.metadata_loader.metadata_csv_path).expanduser()

                method = self.method.value
                if method == "lap":
                    # Square distances for laptracking (function expects px^2)
                    tc = int(self.track_cost_dist.value) ** 2
                    gc = int(self.gap_cost_dist.value) ** 2
                    mc = int(self.merging_cost_dist.value)
                    sc = int(self.splitting_cost_dist.value)
                    merging  = (mc ** 2) if mc > 0 else False
                    splitting = (sc ** 2) if sc > 0 else False

                    print("▶️ LAP tracking…", flush=True)
                    print(f"  track_cost_cutoff={tc}, gap_closing_cost_cutoff={gc}")
                    print(f"  gap_closing_max_frame_count={self.gap_max_frames.value}")
                    print(f"  merging_cost_cutoff={merging}, splitting_cost_cutoff={splitting}", flush=True)

                    from __main__ import run_tcell_laptracking
                    new_md = run_tcell_laptracking(
                        metadata=self.metadata_loader.metadata,
                        output_dir=str(out_dir),
                        track_cost_cutoff=tc,
                        gap_closing_cost_cutoff=gc,
                        gap_closing_max_frame_count=int(self.gap_max_frames.value),
                        merging_cost_cutoff=merging,
                        splitting_cost_cutoff=splitting,
                        cell_type=self.cell_type,
                        overwrite=bool(self.overwrite.value),
                    )

                elif method == "trackpy":
                    print("▶️ TrackPy tracking…", flush=True)
                    print(f"  search_range={self.tp_search_range.value}, memory={self.tp_memory.value}")
                    print(f"  adaptive_stop={self.tp_adaptive_stop.value}, adaptive_step={self.tp_adaptive_step.value}", flush=True)

                    from __main__ import run_tcell_trackpy_tracking
                    new_md = run_tcell_trackpy_tracking(
                        metadata=self.metadata_loader.metadata,
                        output_dir=str(out_dir),
                        overwrite=bool(self.overwrite.value),
                        cell_type=self.cell_type,
                        search_range=int(self.tp_search_range.value),
                        memory=int(self.tp_memory.value),
                        adaptive_stop=float(self.tp_adaptive_stop.value),
                        adaptive_step=float(self.tp_adaptive_step.value),
                    )

                else:  # "prop" propagation
                    print("▶️ Propagation tracking…", flush=True)
                    from __main__ import run_propagation_tracking
                    new_md = run_propagation_tracking(
                        metadata=self.metadata_loader.metadata,
                        output_dir=str(out_dir),
                        cell_type=self.cell_type,
                        overwrite=bool(self.overwrite.value),
                    )

                # Update loader + ALWAYS save CSV
                self.metadata_loader.metadata = new_md
                print(f"💾 Saving metadata to: {csv_path}", flush=True)
                new_md.to_csv(csv_path, sep=",", index=False)
                print("✅ Tracking finished.", flush=True)

            except Exception:
                import traceback; traceback.print_exc()
            finally:
                self._lock(False)
class TrackingVisualizationPanel:
    def __init__(
        self,
        metadata_loader,
        default_time_range=None,
        channel_colors=None,
    ):
        self.metadata_loader = metadata_loader
        # prefer JSON channel colors if not passed
        self.channel_colors = tuple(channel_colors or _cfg_get(
            CFG, "tracking_visualization.channel_colors",
            ["cyan", "yellow", "red", "green", "magenta", "blue"]
        ))

        self._viewer = None
        
        # --- UI: status + refresh
        self._status = widgets.HTML("<b>Waiting for user to load metadata…</b>")
        self._refresh_btn = widgets.Button(
            description="Refresh",
            tooltip="Build/Update selector from metadata_loader.metadata",
        )
        self._refresh_btn.on_click(self._on_refresh_clicked)

        # --- Main controls (disabled until metadata present)
        self.sample_dropdown = widgets.Dropdown(
            options=[],
            value=None,
            description="Sample:",
            layout=widgets.Layout(width="350px"),
            disabled=True,
        )

        # Tickbox to enable/disable time range (from JSON)
        self.use_range = widgets.Checkbox(
            description="Use custom time range",
            value=bool(_cfg_get(CFG, "tracking_visualization.use_range", False)),
            indent=False,
        )

        # Start/End boxes (defaults from constructor OR JSON)
        if isinstance(default_time_range, (tuple, list)) and len(default_time_range) == 2:
            _start_default, _end_default = map(int, default_time_range)
        else:
            _start_default = int(_cfg_get(CFG, "tracking_visualization.start_t", 0))
            _end_default   = int(_cfg_get(CFG, "tracking_visualization.end_t", 0))

        self.start_t = widgets.IntText(
            description="Start T:",
            value=_start_default,
            layout=widgets.Layout(width="180px")
        )
        self.end_t   = widgets.IntText(
            description="End T:",
            value=_end_default,
            layout=widgets.Layout(width="180px")
        )

        self.range_box = widgets.HBox([self.start_t, self.end_t])
        self.range_box.layout.display = "flex" if self.use_range.value else "none"

        def _toggle_range_visibility(change):
            self.range_box.layout.display = "flex" if change["new"] else "none"
        self.use_range.observe(_toggle_range_visibility, names="value")

        self.open_button = widgets.Button(
            description="Visualize segmentation and tracking results",
            button_style="primary",
            tooltip="Launch Napari for the selected sample",
            icon="eye",
            layout=widgets.Layout(width="300px"),
            disabled=True,
        )

        self.close_button = widgets.Button(
            description="Close viewer",
            button_style="danger",
            icon="stop",
            tooltip="Close the active Napari viewer",
            layout=widgets.Layout(width="200px", display="none"),  # hidden by default
        )
        
        self.msg = widgets.Output()

        # Events
        self.open_button.on_click(self._on_open_clicked)
        self.close_button.on_click(self._on_close_clicked)

        # Layout
        self._panel = widgets.VBox(
            [
                widgets.HBox([self._status, self._refresh_btn]),
                self.sample_dropdown,
                self.use_range,
                self.range_box,
                widgets.HBox([self.open_button, self.close_button]),
                self.msg,
            ]
        )

        # Try immediate build if metadata already present
        self._maybe_build_from_loader()

    # ---------------- Public API ----------------
    def display(self):
        display(self._panel)

    def get_selected_row(self) -> pd.Series:
        self._ensure_metadata_ready()
        name = self.sample_dropdown.value
        df = self.metadata_loader.metadata
        row = df[df["sample_name"].astype(str) == str(name)]
        if row.empty:
            raise KeyError(f"sample_name '{name}' not found in metadata.")
        return row.iloc[0]

    def open_viewer(self):
        row = self.get_selected_row()

        # only use range if tickbox checked
        if self.use_range.value:
            start_t, end_t = int(self.start_t.value), int(self.end_t.value)
            if end_t < start_t:
                start_t, end_t = end_t, start_t
            time_range = (start_t, end_t)
        else:
            time_range = None

        with self.msg:
            self.msg.clear_output()
            try:
                print(f"Opening viewer for: {row['sample_name']}")
                self._viewer = visualize_tracks(
                    metadata_row=row,
                    timepoint_range=time_range,
                    channel_colors=self.channel_colors,
                )
            except Exception as e:
                print(f"Error: {e}")
            self.close_button.layout.display = "inline-block" if self._viewer is not None else "none"

    # ---------------- Internal helpers ----------------
    def _ensure_metadata_ready(self):
        df = getattr(self.metadata_loader, "metadata", None)
        if not isinstance(df, pd.DataFrame) or df.empty:
            raise RuntimeError("Metadata not loaded yet. Click 'Refresh' once metadata_loader.metadata is set.")
   
    def _maybe_build_from_loader(self):
        df = getattr(self.metadata_loader, "metadata", None)
        if isinstance(df, pd.DataFrame) and not df.empty and ("sample_name" in df.columns):
            self._build_from_metadata(df)
            self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"
        else:
            self._status.value = "<b>Waiting for user to load metadata…</b>"
            self.sample_dropdown.disabled = True
            self.open_button.disabled = True

    def _build_from_metadata(self, df: pd.DataFrame):
        sample_names = df["sample_name"].astype(str).unique().tolist()
        if not sample_names:
            self._status.value = "<b>No sample_name values found in metadata.</b>"
            self.sample_dropdown.options = []
            self.sample_dropdown.value = None
            self.sample_dropdown.disabled = True
            self.open_button.disabled = True
            return

        desired = _cfg_get(CFG, "tracking_visualization.sample_name", None)
        self.sample_dropdown.options = sample_names
        self.sample_dropdown.value = desired if (desired in sample_names) else sample_names[0]
        self.sample_dropdown.disabled = False
        self.open_button.disabled = False

    # ---------------- Callbacks ----------------
    def _on_open_clicked(self, _):
        self.open_viewer()

    def _on_close_clicked(self, _):
        # Only relevant when we’re in non-blocking mode and hold a viewer
        with self.msg:
            try:
                if self._viewer is not None:
                    print("Closing viewer…")
                    self._viewer.close()
                    self._viewer = None
                else:
                    print("No active viewer to close (viewer was likely opened in blocking mode).")
            finally:
                # Hide the button regardless
                self.close_button.layout.display = "none"
                
    def _on_refresh_clicked(self, _):
        with self.msg:
            self.msg.clear_output()
            try:
                self._maybe_build_from_loader()
            except Exception as e:
                print(f"Refresh failed: {e}")
                
# ##Define widgets
# class PathPicker(widgets.VBox):
#     """
#     Displayable widget for picking a file or directory.
#     - mode: 'file' or 'dir'
#     - start_dir: starting directory for the chooser
#     - default: default filename (file mode) or starting path (dir mode)
#     - filter_pattern: e.g. '*.csv' (file mode only)
#     """
#     def __init__(
#         self,
#         mode='file',
#         start_dir='.',
#         default='',
#         description='Path:',
#         button_description='Browse…',
#         placeholder='Type a path or click Browse…',
#         filter_pattern=None,
#         description_width='90px',
#         width='100%',
#     ):
#         assert mode in ('file', 'dir')
#         self._mode = mode

#         # Text box (the only persistent UI)
#         self.text = widgets.Textarea(
#             value=str(default or ''),
#             description=description,
#             placeholder=placeholder,
#             style={'description_width': description_width},
#             layout=widgets.Layout(width=width, height='32px'),
#         )

#         # Browse button
#         self.button = widgets.Button(
#             description=button_description, icon='folder-open',
#             tooltip=f"Choose a {'folder' if mode=='dir' else 'file'}"
#         )

#         # Resolve starting directory & filename
#         start_dir = start_dir or '.'
#         filename = ''
#         if mode == 'file' and default:
#             p = Path(default)
#             if p.is_absolute():
#                 start_dir = str(p.parent if p.parent.exists() else start_dir)
#                 filename = p.name
#             else:
#                 filename = str(default)
#         elif mode == 'dir' and default and os.path.isdir(default):
#             start_dir = str(default)

#         # FileChooser (hidden until button click)
#         self.fc = FileChooser(path=start_dir, filename=filename)
#         self.fc.title = f"<b>Select {'directory' if mode=='dir' else 'file'}</b>"
#         self.fc.show_only_dirs = (mode == 'dir')
#         self.fc.use_dir_icons = True
#         if filter_pattern and mode == 'file':
#             self.fc.filter_pattern = filter_pattern

#         # Container: row with Text + Button (chooser injected below when needed)
#         self._row = widgets.HBox([self.text, self.button])
#         super().__init__([self._row])

#         # Wire events
#         self.button.on_click(self._toggle_chooser)
#         self.fc.register_callback(self._on_select)
#         self.text.observe(self._on_text_change, names='value')

#     # ---- public API ----
#     @property
#     def value(self) -> str:
#         """Current selected path (string)."""
#         return self.text.value

#     @value.setter
#     def value(self, path: str):
#         """Set the path programmatically (updates chooser location if possible)."""
#         self.text.value = str(path or '')

#     @property
#     def mode(self) -> str:
#         return self._mode

#     # ---- internal helpers ----
#     def _toggle_chooser(self, _=None):
#         if self.fc in self.children:
#             # hide
#             self.children = tuple(c for c in self.children if c is not self.fc)
#         else:
#             # show
#             self.children = tuple(list(self.children) + [self.fc])

#     def _on_select(self, chooser):
#         # Fill text and hide chooser
#         if self._mode == 'dir':
#             self.text.value = chooser.selected_path
#         else:
#             self.text.value = chooser.selected
#         if self.fc in self.children:
#             self.children = tuple(c for c in self.children if c is not self.fc)

#     def _on_text_change(self, change):
#         # Keep chooser in sync when user types/pastes a path
#         newv = (change.get('new') or '').strip()
#         if not newv:
#             return
#         try:
#             if self._mode == 'dir':
#                 if os.path.isdir(newv):
#                     self.fc.reset(path=newv)
#             else:
#                 # file mode: split into directory + filename
#                 parent = os.path.dirname(newv) or '.'
#                 fname = os.path.basename(newv)
#                 if os.path.isdir(parent):
#                     if fname:
#                         self.fc.reset(path=parent, filename=fname)
#                     else:
#                         self.fc.reset(path=parent)
#         except Exception:
#             # Don't crash UI if reset fails; ignore silently
#             pass
# class MetadataLoader(widgets.VBox):
#     """
#     A VBox widget that wraps a file PathPicker and a 'Load' button.
#     After clicking the button, `value` holds the loaded pandas DataFrame,
#     and `path` holds the CSV path. Works without globals; read .value later.

#     Parameters
#     ----------
#     file_picker : a widget with a `.value` path (e.g., your PathPicker)
#         If None, you can pass one later via set_file_picker().
#     button_description : str
#         Button label text.
#     """
#     _busy = Bool(default_value=False).tag(sync=False)
#     _handler_bound = Bool(default_value=False).tag(sync=False)
    
#     def __init__(
#         self,
#         metadata_path_picker,
#         output_dir_picker,
#         button_description: str = "Load metadata",
#         **kwargs,
#     ):
#         super().__init__(**kwargs)
#         self.metadata = None
#         self.output_dir = None
#         self.metadata_csv_path = None
            
#         self.file_picker = metadata_path_picker
#         self.output_dir_picker = output_dir_picker
        
#         self.button = widgets.Button(description=button_description, button_style='success')
#         self.out = widgets.Output()

#         # self.button.on_click(self._on_click)
#         self._click_handler = self._on_click
#         self.button.on_click(self._click_handler)
#         # Build UI (file_picker may be None initially)
#         children = []
#         if self.output_dir_picker is not None:
#             children.append(self.output_dir_picker)
#         if self.file_picker is not None:
#             children.append(self.file_picker)  
#         children += [self.button, self.out]
#         self.children = tuple(children)

   
#     def set_file_picker(self, file_picker: widgets.Widget):
#         """Attach/replace the file picker widget (must have a `.value` path)."""
#         self.file_picker = file_picker
#         # Rebuild children with the new picker at the top
#         # self.children = (self.file_picker, self.button, self.out)

#     def load(self, path = None):
#         """Programmatic load (same as clicking the button)."""
#         if path is None:
#             if self.file_picker is None:
#                 raise ValueError("No file_picker attached and no path provided.")
#             path = self.file_picker.value

#         with self.out:
#             clear_output(wait=True)
                
#             if not path or not str(path).lower().endswith(".csv"):
#                 print("Please choose a .csv file.")
#                 return

#             self.metadata = load_behav3d_metadata(path)
#             self.output_dir = self.output_dir_picker.value
#             self.metadata_csv_path = str(path)
            
#             check_behav3d_metadata(self.metadata)

#             print("✅ Checks passed!")
#             display(self.metadata)
   
#     # Button handler
#     def _on_click(self, _):
#         if self._busy:
#             return  # re-entrancy guard prevents double execution
#         self._busy = True
#         try:
#             self.load()
#         except Exception as e:
#             with self.out:
#                 print(f"❌ Error: {e}")
#         finally:
#             self._busy = False

# def convert_zarr_button(
#         metadata_loader,
#         dim_order_widget
#     ):
#     btn = widgets.Button(
#     description="Convert to Zarr",
#     button_style="success",  # 'success', 'info', 'warning', 'danger' or ''
#     icon="cogs"
#     )
#     out = widgets.Output()

#     def _run_conversion(_):
#         with out:
#             out.clear_output()
#             print("Starting conversion…")
#             btn.disabled = True
#             try:
#                 result = convert_input_files_to_zarr(
#                     metadata=metadata_loader.metadata,
#                     output_dir=metadata_loader.output_dir
#                     )  # add args here if needed 
#             except Exception as e:
#                 import traceback; traceback.print_exc()
#             finally:
#                 metadata_loader.metadata = result
#                 metadata_loader.metadata.to_csv(metadata_loader.metadata_csv_path, index=False)
#                 print("Done ✅")
#                 btn.disabled = False

#     btn.on_click(_run_conversion)
#     display(widgets.VBox([btn, out]))
# class DimOrderTable:
#     DIM_ORDER_OPTIONS = [
#         "TCZYX",
#         "TZCYX",
#         "ZCTYX",
#         "ZTCYX",
#         "CZTYX",
#         "CTZYX",
#     ]
#     DEFAULT_ORDER = "TCZYX"


#     def __init__(
#         self,
#         metadata_loader,
#         sample_col="sample_name",
#         path_col="raw_image_path",
#         widths=("40%","30%","30%"),
#         dim_col="dimension_order",                 # <-- column to write string like 'TCZYX'
#         auto_write=False                    # <-- if True, update metadata table on each change
#     ):
#         self.metadata_loader = metadata_loader
#         self.sample_col = sample_col
#         self.path_col = path_col
#         self.dim_col = dim_col
#         self.auto_write = auto_write

#         self._width_sample, self._width_shape, self._width_order = widths

#         # State
#         self._rows = []          # [{'sample','path','shape_label','dd'}, ...]
#         self._selections = {}    # sample_name -> dim order tuple

#         # Precompute allowed tuples for validation / prefill
#         self._allowed_orders = self.DIM_ORDER_OPTIONS

#         # UI skeleton
#         self._status = widgets.HTML("<b>Waiting for user to load metadata…</b>")
#         self._refresh_btn = widgets.Button(description="Refresh", tooltip="Build/Update table from metadata")
#         self._refresh_btn.on_click(self._on_refresh)

#         self._apply_all_dd = widgets.Dropdown(
#             options=self.DIM_ORDER_OPTIONS,
#             value=self.DEFAULT_ORDER,
#             description="Apply to all:",
#             style={'description_width': 'initial'},
#             layout=widgets.Layout(width="350px", margin="6px 8px 6px 0")
#         )
#         self._apply_all_btn = widgets.Button(description="Apply", layout=widgets.Layout(width="100px"))
#         self._apply_all_btn.on_click(lambda _: self.set_all(self._apply_all_dd.value))

#         self._header = widgets.HBox([
#             widgets.Label("Sample", layout=widgets.Layout(width=self._width_sample, overflow="hidden")),
#             widgets.Label("Image shape", layout=widgets.Layout(width=self._width_shape, overflow="hidden")),
#             widgets.Label("Dim order", layout=widgets.Layout(width=self._width_order, overflow="hidden")),
#         ])
#         self._table_body = widgets.VBox()
#         self._out = widgets.Output()

#         self.widget = widgets.VBox([
#             self._status,
#             widgets.HBox([self._refresh_btn]),
#             self._header,
#             self._table_body,
#             widgets.HBox([self._apply_all_dd, self._apply_all_btn]),
#             self._out
#         ])
        
#     def display(self):
#         display(self.widget)

#     def get_selections(self) -> dict:
#         """Return dict: sample_name -> 'TCZYX' """
#         return dict(self._selections)

#     def get_selections_str(self) -> dict:
#         """Return dict: sample_name -> 'TCZYX'"""
#         return {s: self._order_to_str(o) for s, o in self._selections.items()}

#     def write_dimorder_to_metadata(self, col=None):
#         """
#         Write string dimorder (e.g. 'TCZYX') into metadata_loader.metadata[col].
#         Returns the updated DataFrame.
#         """
#         col = col or self.dim_col
#         df = getattr(self.metadata_loader, "metadata", None)
#         if not isinstance(df, pd.DataFrame):
#             raise ValueError("metadata_loader.metadata is not a DataFrame yet.")
#         if df.empty:
#             raise ValueError("metadata_loader.metadata is empty.")
#         if self.sample_col not in df.columns:
#             raise ValueError(f"metadata is missing '{self.sample_col}' column.")
#         # Map selections by sample name
#         str_map = self.get_selections_str()
#         self.metadata_loader.metadata[col] = df[self.sample_col].map(str_map)
#         return self.metadata_loader.metadata

#     def set_all(self, order_tuple):
#         """Programmatically set all rows to the given order tuple."""
#         for row in self._rows:
#             row['dd'].value = order_tuple  # triggers observers -> updates _selections (+ auto_write)

#     def refresh_shapes(self):
#         """Recompute shape labels (e.g., after mounting a drive)."""
#         for row in self._rows:
#             row['shape_label'].value = self._probe_image_shape(row['path'])

#     def to_dataframe(self):
#         """
#         Return a pandas.DataFrame with:
#           sample_name, raw_image_path, image_shape, dim_order (tuple), dimorder (string)
#         """
#         data = []
#         for row in self._rows:
#             tup = row["dd"].value
#             data.append({
#                 "sample_name": row["sample"],
#                 "raw_image_path": row["path"],
#                 "image_shape": row["shape_label"].value,
#                 "dim_order": tup,
#                 self.dim_col: self._order_to_str(tup),
#             })
#         return pd.DataFrame(data)

#     def _on_refresh(self, _btn):
#         df = getattr(self.metadata_loader, "metadata", None)
#         if not (isinstance(df, pd.DataFrame) and not df.empty):
#             self._status.value = "<b>Waiting for user to load metadata…</b>"
#             return
#         missing = [c for c in (self.sample_col, self.path_col) if c not in df.columns]
#         if missing:
#             self._status.value = f"<b>Metadata missing columns: {missing}</b>"
#             return

#         self._build_rows_from(df)
#         self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"

#     def _build_rows_from(self, df: pd.DataFrame):
#         # Reset state
#         self._rows.clear()
#         self._selections.clear()
#         row_boxes = []

#         # Optional: prefill from existing string column if present
#         prefill = None
#         if self.dim_col in df.columns:
#             prefill = df[self.dim_col].astype(str).to_dict()  # index -> string

#         for idx, r in df.iterrows():
#             sample = str(r[self.sample_col])
#             path = r[self.path_col]

#             probed_shape = self._probe_image_shape(path)
#             sample_lbl = widgets.Label(sample, layout=widgets.Layout(width=self._width_sample, overflow="hidden"))
#             shape_lbl  = widgets.Label(probed_shape, layout=widgets.Layout(width=self._width_shape, overflow="hidden"))

#             # Determine default dropdown value (from prefill if valid)
#             default_val = get_image_dimension_order(path)
#             if default_val is None:
#                 default_val = self.DEFAULT_ORDER
#             if prefill is not None:
#                 s = prefill.get(idx, None)
#                 if isinstance(s, str):
#                     tu = tuple(s.upper())
#                     if tu in self._allowed_orders:
#                         default_val = tu

#             dd = widgets.Dropdown(
#                 options=self.DIM_ORDER_OPTIONS,
#                 value=default_val,
#                 layout=widgets.Layout(width=self._width_order),
#             )

#             # keep live selection
#             self._selections[sample] = dd.value
#             dd.observe(lambda ch, s=sample: self._on_dd_changed(s, ch), names='value')

#             self._rows.append({"sample": sample, "path": path, "shape_label": shape_lbl, "dd": dd})
#             row_boxes.append(widgets.HBox([sample_lbl, shape_lbl, dd]))

#         # Swap UI body
#         self._table_body.children = tuple(row_boxes)

#     def _on_dd_changed(self, sample_name, change):
#         if change.get('name') == 'value':
#             new_tuple = change['new']
#             self._selections[sample_name] = new_tuple
#             if self.auto_write:
#                 # live-write into metadata table (string form) for matching rows
#                 df = getattr(self.metadata_loader, "metadata", None)
#                 if isinstance(df, pd.DataFrame) and self.sample_col in df.columns:
#                     mask = df[self.sample_col].astype(str) == str(sample_name)
#                     df.loc[mask, self.dim_col] = self._order_to_str(new_tuple)

#     def _probe_image_shape(self, path):
#         path = Path(path)
#         if not path.exists():
#             return "⛔ missing"
#         try:
#             shape = get_image_shape(path)
#             return "×".join(map(str, shape)) if shape else "unknown"
#         except Exception as e:
#             return f"⚠️ {type(e).__name__}"

#     @staticmethod
#     def _order_to_str(order_tuple):
#         return "".join(order_tuple)

# def _mk_timepoint_range(use_all: bool, start: int, end: int):
#     if use_all:
#         return None
#     return (int(start), int(end))

# class PixelClassifierPanel:
#     def __init__(self, metadata_loader):
#         self.metadata_loader = metadata_loader

#         # -------- Train controls --------
#         self.examples_per_sample = widgets.IntSlider(
#             description="Examples per sample", value=3, min=1, max=50, step=1, continuous_update=False
#         )
#         self.sample_specific_classifier = widgets.Checkbox(
#             description="Sample-specific classifier", value=False
#         )
#         self.n_workers = widgets.IntSlider(
#             description="Workers", value=(os.cpu_count() or 8), min=1, max=max(8, (os.cpu_count() or 8)),
#             step=1, continuous_update=False
#         )

#         # -------- Apply controls --------
#         self.organoid_edt_threshold = widgets.FloatSlider(
#             description="Organoid EDT thr",
#             value=12.0, min=0.0, max=50.0, step=0.5, readout_format=".1f",
#             continuous_update=False,
#             style={'description_width': '160px'}  # wider label
#         )
#         self.use_all_timepoints = widgets.Checkbox(
#             description="Process ALL timepoints", value=True
#         )
#         self.tp_start = widgets.IntText(description="Start t", value=0)
#         self.tp_end   = widgets.IntText(description="End t", value=0)

#         # Show/hide start/end boxes based on checkbox
#         self.use_all_timepoints.observe(self._toggle_timepoint_inputs, names='value')
#         self._toggle_timepoint_inputs()  # initialize

#         # Buttons + log
#         self.btn_train = widgets.Button(description="Train", button_style="primary")
#         self.btn_apply = widgets.Button(description="Apply", button_style="success")
#         self.btn_train.on_click(self._on_train_clicked)
#         self.btn_apply.on_click(self._on_apply_clicked)
#         self.out = widgets.Output()

#         # Layout
#         train_box = widgets.VBox([
#             widgets.HTML("<b>Train pixel classifier</b>"),
#             widgets.HBox([self.examples_per_sample, self.n_workers]),
#             self.sample_specific_classifier,
#             self.btn_train,
#         ])

#         self.tp_row = widgets.HBox([self.use_all_timepoints, self.tp_start, self.tp_end])

#         apply_box = widgets.VBox([
#             widgets.HTML("<b>Apply segmentation</b>"),
#             self.organoid_edt_threshold,
#             self.tp_row,
#             self.btn_apply,
#         ])

#         self.ui = widgets.VBox([train_box, widgets.HTML("<hr>"), apply_box, widgets.HTML("<hr>"), self.out])

#     def _toggle_timepoint_inputs(self, change=None):
#         show = not self.use_all_timepoints.value
#         disp = None if show else 'none'
#         self.tp_start.layout.display = disp
#         self.tp_end.layout.display   = disp
#         self.tp_start.disabled = not show
#         self.tp_end.disabled   = not show

#     def display(self):
#         display(self.ui)

#     def _lock(self, state: bool):
#         for w in [
#             self.btn_train, self.btn_apply,
#             self.examples_per_sample, self.sample_specific_classifier, self.n_workers,
#             self.organoid_edt_threshold, self.use_all_timepoints, self.tp_start, self.tp_end,
#         ]:
#             w.disabled = state

#     # ---- callbacks ----
#     def _on_train_clicked(self, _):
#         self._lock(True)
#         with self.out:
#             self.out.clear_output()
#             try:
#                 odir = Path(self.metadata_loader.output_dir).expanduser()
#                 odir.mkdir(parents=True, exist_ok=True)

#                 print("▶️ Training pixel classifier…", flush=True)
#                 print(f"  output_dir={odir}")
#                 print(f"  examples_per_sample={self.examples_per_sample.value}")
#                 print(f"  sample_specific_classifier={self.sample_specific_classifier.value}")
#                 print(f"  n_workers={self.n_workers.value}")

#                 train_pixel_classifier(
#                     output_dir=str(odir),
#                     metadata=self.metadata_loader.metadata,
#                     examples_per_sample=int(self.examples_per_sample.value),
#                     sample_specific_classifier=bool(self.sample_specific_classifier.value),
#                     n_workers=int(self.n_workers.value),
#                 )
#                 print("✅ Training finished.", flush=True)
#             except Exception:
#                 traceback.print_exc()
#             finally:
#                 self._lock(False)

#     def _on_apply_clicked(self, _):
#         self._lock(True)
#         with self.out:
#             self.out.clear_output()
#             try:
#                 odir = Path(self.metadata_loader.output_dir).expanduser()
#                 odir.mkdir(parents=True, exist_ok=True)

#                 tpr = _mk_timepoint_range(
#                     use_all=bool(self.use_all_timepoints.value),
#                     start=int(self.tp_start.value),
#                     end=int(self.tp_end.value),
#                 )

#                 print("▶️ Applying pixel classifier segmentation…", flush=True)
#                 print(f"  output_dir={odir}")
#                 print(f"  organoid_edt_threshold={self.organoid_edt_threshold.value}")
#                 print(f"  timepoint_range={tpr}", flush=True)

#                 new_md = run_pixel_classifier_segmentation(
#                     output_dir=str(odir),
#                     metadata=self.metadata_loader.metadata,
#                     organoid_edt_threshold=float(self.organoid_edt_threshold.value),
#                     timepoint_range=tpr
#                 )

#                 # Update loader + ALWAYS save CSV
#                 self.metadata_loader.metadata = new_md
#                 csv_path = Path(self.metadata_loader.metadata_csv_path).expanduser()
#                 print(f"💾 Saving metadata to: {csv_path}", flush=True)
#                 new_md.to_csv(csv_path, sep=",", index=False)

#                 print("✅ Apply finished.", flush=True)
#             except Exception:
#                 traceback.print_exc()
#             finally:
#                 self._lock(False)

# class TrackingPanel:
#     """
#     Minimal UI for tracking (LAP, TrackPy, or Propagation).
#     - Uses metadata_loader.output_dir and metadata_loader.metadata_csv_path directly
#     - Updates metadata_loader.metadata and always saves CSV
#     - cell_type is supplied in __init__
#     """
#     def __init__(self, metadata_loader, cell_type="tcell"):
#         self.metadata_loader = metadata_loader
#         self.cell_type = str(cell_type).strip()

#         # Method selector (now includes Propagation)
#         self.method = widgets.Dropdown(
#             description="Tracking method",
#             options=[("LAP (laptrack)", "lap"),
#                      ("TrackPy", "trackpy"),
#                      ("Propagation", "prop")],
#             value="lap",
#             style={'description_width': '160px'}
#         )

#         # Shared controls
#         self.overwrite = widgets.Checkbox(description="Overwrite existing", value=False)

#         # LAP params (distances in pixels; squared when calling)
#         self.track_cost_dist = widgets.IntSlider(
#             description="Track cost (px)", value=45, min=1, max=200, step=1,
#             continuous_update=False, style={'description_width':'160px'}
#         )
#         self.gap_cost_dist = widgets.IntSlider(
#             description="Gap close cost (px)", value=60, min=1, max=300, step=1,
#             continuous_update=False, style={'description_width':'160px'}
#         )
#         self.gap_max_frames = widgets.IntSlider(
#             description="Gap close max frames", value=3, min=0, max=20, step=1,
#             continuous_update=False, style={'description_width':'160px'}
#         )
#         self.merging_cost_dist = widgets.IntText(
#             description="Merging cost (px)", value=0, style={'description_width':'160px'}
#         )
#         self.splitting_cost_dist = widgets.IntText(
#             description="Splitting cost (px)", value=0, style={'description_width':'160px'}
#         )
#         self.lap_params = widgets.VBox([
#             widgets.HTML("<b>LAP parameters</b>"),
#             self.track_cost_dist, self.gap_cost_dist, self.gap_max_frames,
#             self.merging_cost_dist, self.splitting_cost_dist
#         ])

#         # TrackPy params
#         self.tp_search_range = widgets.IntSlider(
#             description="Search range (px)", value=31, min=1, max=200, step=1,
#             continuous_update=False, style={'description_width':'160px'}
#         )
#         self.tp_memory = widgets.IntSlider(
#             description="Memory (frames)", value=2, min=0, max=20, step=1,
#             continuous_update=False, style={'description_width':'160px'}
#         )
#         self.tp_adaptive_stop = widgets.FloatSlider(
#             description="Adaptive stop", value=10.0, min=0.0, max=50.0, step=0.5,
#             readout_format=".1f", continuous_update=False, style={'description_width':'160px'}
#         )
#         self.tp_adaptive_step = widgets.FloatSlider(
#             description="Adaptive step", value=0.95, min=0.1, max=1.0, step=0.01,
#             readout_format=".2f", continuous_update=False, style={'description_width':'160px'}
#         )
#         self.tp_params = widgets.VBox([
#             widgets.HTML("<b>TrackPy parameters</b>"),
#             self.tp_search_range, self.tp_memory, self.tp_adaptive_stop, self.tp_adaptive_step
#         ])

#         # Propagation params (none needed; just an info box)
#         self.prop_params = widgets.VBox([
#             widgets.HTML("<b>Propagation tracking</b>"),
#             widgets.HTML("<i>No tunable parameters.</i>")
#         ])

#         # Swap param groups on method change
#         self.param_container = widgets.VBox()
#         self.method.observe(self._toggle_param_groups, names='value')
#         self._toggle_param_groups()

#         # Run + log
#         self.btn_run = widgets.Button(description="Run tracking", button_style="success")
#         self.btn_run.on_click(self._on_run_clicked)
#         self.out = widgets.Output()

#         # Layout
#         self.ui = widgets.VBox([
#             widgets.HTML("<b>Tracking</b>"),
#             self.method,
#             widgets.HBox([self.overwrite]),
#             self.param_container,
#             self.btn_run,
#             widgets.HTML("<hr>"),
#             self.out
#         ])

#     def display(self):
#         display(self.ui)

#     def _toggle_param_groups(self, change=None):
#         if self.method.value == "lap":
#             self.param_container.children = (self.lap_params,)
#         elif self.method.value == "trackpy":
#             self.param_container.children = (self.tp_params,)
#         else:  # "prop"
#             self.param_container.children = (self.prop_params,)

#     def _lock(self, state: bool):
#         for w in [
#             self.method, self.overwrite,
#             self.track_cost_dist, self.gap_cost_dist, self.gap_max_frames,
#             self.merging_cost_dist, self.splitting_cost_dist,
#             self.tp_search_range, self.tp_memory, self.tp_adaptive_stop, self.tp_adaptive_step,
#             self.btn_run
#         ]:
#             if hasattr(w, "disabled"):
#                 w.disabled = state

#     def _on_run_clicked(self, _):
#         self._lock(True)
#         with self.out:
#             self.out.clear_output()
#             try:
#                 out_dir = Path(self.metadata_loader.output_dir).expanduser()
#                 out_dir.mkdir(parents=True, exist_ok=True)

#                 csv_path = Path(self.metadata_loader.metadata_csv_path).expanduser()

#                 method = self.method.value
#                 if method == "lap":
#                     # Square distances for laptracking (function expects px^2)
#                     tc = int(self.track_cost_dist.value) ** 2
#                     gc = int(self.gap_cost_dist.value) ** 2
#                     mc = int(self.merging_cost_dist.value)
#                     sc = int(self.splitting_cost_dist.value)
#                     merging  = (mc ** 2) if mc > 0 else False
#                     splitting = (sc ** 2) if sc > 0 else False

#                     print("▶️ LAP tracking…", flush=True)
#                     print(f"  track_cost_cutoff={tc}, gap_closing_cost_cutoff={gc}")
#                     print(f"  gap_closing_max_frame_count={self.gap_max_frames.value}")
#                     print(f"  merging_cost_cutoff={merging}, splitting_cost_cutoff={splitting}", flush=True)

#                     from __main__ import run_tcell_laptracking
#                     new_md = run_tcell_laptracking(
#                         metadata=self.metadata_loader.metadata,
#                         output_dir=str(out_dir),
#                         track_cost_cutoff=tc,
#                         gap_closing_cost_cutoff=gc,
#                         gap_closing_max_frame_count=int(self.gap_max_frames.value),
#                         merging_cost_cutoff=merging,
#                         splitting_cost_cutoff=splitting,
#                         cell_type=self.cell_type,
#                         overwrite=bool(self.overwrite.value),
#                     )

#                 elif method == "trackpy":
#                     print("▶️ TrackPy tracking…", flush=True)
#                     print(f"  search_range={self.tp_search_range.value}, memory={self.tp_memory.value}")
#                     print(f"  adaptive_stop={self.tp_adaptive_stop.value}, adaptive_step={self.tp_adaptive_step.value}", flush=True)

#                     from __main__ import run_tcell_trackpy_tracking
#                     new_md = run_tcell_trackpy_tracking(
#                         metadata=self.metadata_loader.metadata,
#                         output_dir=str(out_dir),
#                         overwrite=bool(self.overwrite.value),
#                         cell_type=self.cell_type,
#                         search_range=int(self.tp_search_range.value),
#                         memory=int(self.tp_memory.value),
#                         adaptive_stop=float(self.tp_adaptive_stop.value),
#                         adaptive_step=float(self.tp_adaptive_step.value),
#                     )

#                 else:  # "prop" propagation
#                     print("▶️ Propagation tracking…", flush=True)
#                     from __main__ import run_propagation_tracking
#                     new_md = run_propagation_tracking(
#                         metadata=self.metadata_loader.metadata,
#                         output_dir=str(out_dir),
#                         cell_type=self.cell_type,      # e.g., "organoid" or "tcell"
#                         overwrite=bool(self.overwrite.value),
#                     )

#                 # Update loader + ALWAYS save CSV
#                 self.metadata_loader.metadata = new_md
#                 print(f"💾 Saving metadata to: {csv_path}", flush=True)
#                 new_md.to_csv(csv_path, sep=",", index=False)
#                 print("✅ Tracking finished.", flush=True)

#             except Exception:
#                 traceback.print_exc()
#             finally:
#                 self._lock(False)

# class TrackingVisualizationPanel:
#     def __init__(
#         self,
#         metadata_loader,
#         default_time_range=None,
#         channel_colors=("cyan", "yellow", "red", "green", "magenta", "blue"),
#     ):
#         self.metadata_loader = metadata_loader
#         self.channel_colors = channel_colors

#         self._viewer = None
        
#         # --- UI: status + refresh
#         self._status = widgets.HTML("<b>Waiting for user to load metadata…</b>")
#         self._refresh_btn = widgets.Button(
#             description="Refresh",
#             tooltip="Build/Update selector from metadata_loader.metadata",
#         )
#         self._refresh_btn.on_click(self._on_refresh_clicked)

#         # --- Main controls (disabled until metadata present)
#         self.sample_dropdown = widgets.Dropdown(
#             options=[],
#             value=None,
#             description="Sample:",
#             layout=widgets.Layout(width="350px"),
#             disabled=True,
#         )

#         # Tickbox to enable/disable time range
#         self.use_range = widgets.Checkbox(
#             description="Use custom time range",
#             value=False,
#             indent=False,
#         )

#         # Start/End boxes (wrapped in a container we can show/hide)
#         self.start_t = widgets.IntText(description="Start T:", layout=widgets.Layout(width="180px"))
#         self.end_t   = widgets.IntText(description="End T:",   layout=widgets.Layout(width="180px"))

#         if isinstance(default_time_range, (tuple, list)) and len(default_time_range) == 2:
#             self.start_t.value = int(default_time_range[0])
#             self.end_t.value   = int(default_time_range[1])
#             self.use_range.value = True
#         else:
#             self.start_t.value = 0
#             self.end_t.value   = 0
#             self.use_range.value = False

#         self.range_box = widgets.HBox([self.start_t, self.end_t])
#         self.range_box.layout.display = "flex" if self.use_range.value else "none"

#         def _toggle_range_visibility(change):
#             self.range_box.layout.display = "flex" if change["new"] else "none"
#         self.use_range.observe(_toggle_range_visibility, names="value")

#         self.open_button = widgets.Button(
#             description="Visualize segmentation and tracking results",
#             button_style="primary",
#             tooltip="Launch Napari for the selected sample",
#             icon="eye",
#             layout=widgets.Layout(width="300px"),
#             disabled=True,
#         )

#         self.close_button = widgets.Button(
#             description="Close viewer",
#             button_style="danger",
#             icon="stop",
#             tooltip="Close the active Napari viewer",
#             layout=widgets.Layout(width="200px", display="none"),  # hidden by default
#         )
        
#         self.msg = widgets.Output()

#         # Events
#         self.open_button.on_click(self._on_open_clicked)
#         self.close_button.on_click(self._on_close_clicked)

#         # Layout
#         self._panel = widgets.VBox(
#             [
#                 widgets.HBox([self._status, self._refresh_btn]),
#                 self.sample_dropdown,
#                 self.use_range,
#                 self.range_box,
#                 widgets.HBox([self.open_button, self.close_button]),
#                 self.msg,
#             ]
#         )

#         # Try immediate build if metadata already present
#         self._maybe_build_from_loader()

#     # ---------------- Public API ----------------
#     def display(self):
#         display(self._panel)

#     def get_selected_row(self) -> pd.Series:
#         self._ensure_metadata_ready()
#         name = self.sample_dropdown.value
#         df = self.metadata_loader.metadata
#         row = df[df["sample_name"].astype(str) == str(name)]
#         if row.empty:
#             raise KeyError(f"sample_name '{name}' not found in metadata.")
#         return row.iloc[0]

#     def open_viewer(self):
#         row = self.get_selected_row()

#         # only use range if tickbox checked
#         if self.use_range.value:
#             start_t, end_t = int(self.start_t.value), int(self.end_t.value)
#             if end_t < start_t:
#                 start_t, end_t = end_t, start_t
#             time_range = (start_t, end_t)
#         else:
#             time_range = None

#         with self.msg:
#             self.msg.clear_output()
#             try:
#                 print(f"Opening viewer for: {row['sample_name']}")
#                 self._viewer = visualize_tracks(
#                     metadata_row=row,
#                     timepoint_range=time_range,
#                     channel_colors=self.channel_colors,
#                 )
#             except Exception as e:
#                 print(f"Error: {e}")
#             self.close_button.layout.display = "inline-block" if self._viewer is not None else "none"

#     # ---------------- Internal helpers ----------------
#     def _ensure_metadata_ready(self):
#         df = getattr(self.metadata_loader, "metadata", None)
#         if not isinstance(df, pd.DataFrame) or df.empty:
#             raise RuntimeError("Metadata not loaded yet. Click 'Refresh' once metadata_loader.metadata is set.")
   
#     def _maybe_build_from_loader(self):
#         df = getattr(self.metadata_loader, "metadata", None)
#         if isinstance(df, pd.DataFrame) and not df.empty and ("sample_name" in df.columns):
#             self._build_from_metadata(df)
#             self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"
#         else:
#             self._status.value = "<b>Waiting for user to load metadata…</b>"
#             self.sample_dropdown.disabled = True
#             self.open_button.disabled = True

#     def _build_from_metadata(self, df: pd.DataFrame):
#         sample_names = df["sample_name"].astype(str).unique().tolist()
#         if not sample_names:
#             self._status.value = "<b>No sample_name values found in metadata.</b>"
#             self.sample_dropdown.options = []
#             self.sample_dropdown.value = None
#             self.sample_dropdown.disabled = True
#             self.open_button.disabled = True
#             return

#         self.sample_dropdown.options = sample_names
#         self.sample_dropdown.value = sample_names[0]
#         self.sample_dropdown.disabled = False
#         self.open_button.disabled = False

#     # ---------------- Callbacks ----------------
#     def _on_open_clicked(self, _):
#         self.open_viewer()

#     def _on_close_clicked(self, _):
#         # Only relevant when we’re in non-blocking mode and hold a viewer
#         with self.msg:
#             try:
#                 if self._viewer is not None:
#                     print("Closing viewer…")
#                     self._viewer.close()
#                     self._viewer = None
#                 else:
#                     print("No active viewer to close (viewer was likely opened in blocking mode).")
#             finally:
#                 # Hide the button regardless
#                 self.close_button.layout.display = "none"
                
#     def _on_refresh_clicked(self, _):
#         with self.msg:
#             self.msg.clear_output()
#             try:
#                 self._maybe_build_from_loader()
#             except Exception as e:
#                 print(f"Refresh failed: {e}")
           
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