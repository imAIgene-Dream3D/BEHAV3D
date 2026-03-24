import ipywidgets as widgets
from pathlib import Path
import pandas as pd
import traceback
import yaml
from copy import deepcopy
from .utils import (
    _cfg_get, 
    spinning_loader
)
from behav3d.preprocessing import convert_input_files_to_zarr
from behav3d.io.images import load_image, get_image_shape, get_image_dimension_order
#from behav3d.preprocessing.unmixing import visualize_unmix
#from behav3d.preprocessing.unmixing.signal_unmixing import signal_unmixing

def convert_zarr_button(metadata_loader, dim_order_widget):
    btn = widgets.Button(
        description="Convert to Zarr",
        button_style="success",
        icon="cogs",
        layout=widgets.Layout(
            width="fit-content",
            flex="0 0 auto"
        )
    )
    spinner = widgets.HTML(value=spinning_loader)
    spinner.layout.display = "none"
    out = widgets.Output()

    selector = widgets.RadioButtons(
        options=['All timepoints', 'Custom time range'],
        value='All timepoints',
        description='Temporal Range:',
        disabled=False,
        style={'description_width': 'initial'}
    )
    
    start_t = widgets.IntText(value=0, description='Start T:', layout=widgets.Layout(width='150px'))
    end_t = widgets.IntText(value=0, description='End T:', layout=widgets.Layout(width='150px'))
    
    range_box = widgets.HBox([start_t, end_t], layout=widgets.Layout(display='none', margin='0 0 10px 0'))

    def _on_selector_change(change):
        range_box.layout.display = 'flex' if change['new'] == 'Custom time range' else 'none'

    selector.observe(_on_selector_change, names='value')

    row = widgets.HBox([btn, spinner], layout=widgets.Layout(align_items="center", gap="8px"))

    def _run_conversion(_):
        with out:
            out.clear_output()
            btn.disabled = True
            spinner.layout.display = None

            try:
                dim_order_widget.write_dimorder_to_metadata()
            except Exception:
                pass

            t_start = None
            t_end = None
            if selector.value == 'Custom time range':
                if start_t.value < 0 or end_t.value < 0:
                    print("Error: Timepoints cannot be negative.")
                    btn.disabled = False
                    spinner.layout.display = "none"
                    return
                if start_t.value > end_t.value:
                    print(f"Error: Start T ({start_t.value}) must be less than or equal to End T ({end_t.value}).")
                    btn.disabled = False
                    spinner.layout.display = "none"
                    return
                
                t_start = start_t.value
                t_end = end_t.value
                print(f"Starting conversion for timepoints {t_start} to {t_end}…")
            else:
                print("Starting conversion for all timepoints…")

            result = metadata_loader.metadata
            try:
                result = convert_input_files_to_zarr(
                    metadata=metadata_loader.metadata,
                    output_dir=metadata_loader.output_dir,
                    t_start=t_start,
                    t_end=t_end
                )
            except Exception:
                traceback.print_exc()
            finally:
                metadata_loader.metadata = result
                try:
                    metadata_loader.metadata.to_csv(metadata_loader.metadata_csv_path, index=False)
                except Exception:
                    traceback.print_exc()
                print("Done ✅")
                btn.disabled = False
                spinner.layout.display = "none"

    btn.on_click(_run_conversion)
    return widgets.VBox([selector, range_box, row, out])

class DimOrderTable:
    DIM_ORDER_OPTIONS = [
        "TCZYX", "TZCYX", "ZCTYX", "ZTCYX", "CZTYX", "CTZYX",
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

        self._rows = []
        self._selections = {}

        self._status = widgets.HTML("<b>Waiting for user to load metadata…</b>")
        self._refresh_btn = widgets.Button(description="Refresh", tooltip="Build/Update table from metadata")
        self._refresh_btn.on_click(self._on_refresh)

        _apply_default = str(_cfg_get(self.metadata_loader.behav3d_parameters, "dim_order.default_apply_all", self.DEFAULT_ORDER))
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

        self._maybe_autoload()

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
        widgets.display(self.widget)

    def get_selections(self) -> dict:
        return dict(self._selections)

    def get_selections_str(self) -> dict:
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
        str_map = self.get_selections_str()
        self.metadata_loader.metadata[col] = df[self.sample_col].map(str_map)
        return self.metadata_loader.metadata

    def set_all(self, order_tuple):
        if isinstance(order_tuple, str):
            order_tuple = tuple(order_tuple)
        for row in self._rows:
            row['dd'].value = order_tuple

    def refresh_shapes(self):
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

    def _build_rows_from(self, df: pd.DataFrame):
        self._rows.clear()
        self._selections.clear()
        row_boxes = []

        prefill = None
        if self.dim_col in df.columns:
            prefill = df[self.dim_col].astype(str).to_dict()

        for idx, r in df.iterrows():
            sample = str(r[self.sample_col])
            path = r[self.path_col]

            probed_shape = self._probe_image_shape(path)
            sample_lbl = widgets.Label(sample, layout=widgets.Layout(width=self._width_sample, overflow="hidden"))
            shape_lbl  = widgets.Label(probed_shape, layout=widgets.Layout(width=self._width_shape, overflow="hidden"))

            default_val = get_image_dimension_order(path)
            if default_val is None:
                default_val = self.DEFAULT_ORDER
            if prefill is not None:
                s = prefill.get(idx, None)
                if isinstance(s, str):
                    tu = tuple(s.upper())
                    if tu in self.DIM_ORDER_OPTIONS:
                        default_val = tu

            dd = widgets.Dropdown(
                options=self.DIM_ORDER_OPTIONS,
                value=default_val,
                layout=widgets.Layout(width=self._width_order),
            )

            self._selections[sample] = dd.value
            dd.observe(lambda ch, s=sample: self._on_dd_changed(s, ch), names='value')

            self._rows.append({"sample": sample, "path": path, "shape_label": shape_lbl, "dd": dd})
            row_boxes.append(widgets.HBox([sample_lbl, shape_lbl, dd]))

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
    