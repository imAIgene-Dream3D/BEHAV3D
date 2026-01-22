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
from behav3d.preprocessing.unmixing import visualize_unmix
from behav3d.preprocessing.unmixing.signal_unmixing import signal_unmixing

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

    row = widgets.HBox([btn, spinner], layout=widgets.Layout(align_items="center", gap="8px"))

    def _run_conversion(_):
        with out:
            out.clear_output()
            print("Starting conversion…")
            btn.disabled = True
            spinner.layout.display = None

            try:
                dim_order_widget.write_dimorder_to_metadata()
            except Exception:
                pass

            result = metadata_loader.metadata
            try:
                result = convert_input_files_to_zarr(
                    metadata=metadata_loader.metadata,
                    output_dir=metadata_loader.output_dir
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
    return widgets.VBox([row, out])

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

class SignalUnmixingPanel:
    def __init__(
            self, 
            metadata_loader,
            default_time_range=None,
            channel_colors=None):
        
        self.metadata_loader = metadata_loader
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing", {})

        self.sink_channel = widgets.IntText(
            description="Sink channel",
            value=int(pc.get("sink_channel", 2))
        )
        self.source_channels = widgets.Text(
            value=pc.get("source_channels", "0,1"), 
            description="Source channels",
            style={"description_width": "initial"}, 
        )
        self.train_z = widgets.Text(
            value=pc.get("train_z", "4,5,6"), 
            description="Train Z", 
        )
        self.train_t = widgets.Text(
            value=pc.get("train_t", "1,3,6"), 
            description="Train T", 
        )
        self.bg_percentile = widgets.IntText(
            value=int(pc.get("bg_percentile", 1)), 
            description="BG percentile", 
        )
        self.median_size = widgets.IntText(
            value=int(pc.get("median_size", 3)), 
            description="Median size",
        )
        self.gaussian_sigma = widgets.IntText(
            value=int(pc.get("gaussian_sigma", 4)), 
            description="Gaussian sigma",
            style={"description_width": "initial"}, 
        )  
        self.use_all_timepoints = widgets.Checkbox(
            description="Process ALL timepoints",
            value=bool(pc.get("use_all_timepoints", True))
        )
        self.tp_n = widgets.IntText(
            description="Number of Timepoints", 
            value=int(pc.get("timepoints", 0)),
            style={"description_width": "initial"}
        )

        self.use_all_timepoints.observe(self._toggle_timepoint_inputs, names='value')
        self._toggle_timepoint_inputs()

        self.training_log = widgets.Checkbox(
            description="Show Training Log",
            value=bool(pc.get("training_log", True))
        )

        sum_scale_widgets = []
        rows = []
        file_list = self.metadata_loader.metadata.sample_name.to_list()

        header = widgets.HBox([
            widgets.Label(value="File Name", layout=widgets.Layout(width='200px')),
            widgets.Label(value="sum_scale (source1, source2)", layout=widgets.Layout(width='200px')),
            widgets.Label(value="")
        ])
        rows.append(header) 

        def apply_to_all_callback(btn):
            first_value = sum_scale_widgets[0]['input'].value
            for i, row in enumerate(sum_scale_widgets):
                if i == 0: continue
                row['input'].value = first_value

        sum_scale_config = pc.get("sum_scale", "20,30")

        for i, file in enumerate(file_list):
            file_label = widgets.Label(value=file, layout=widgets.Layout(width='200px'))
            if isinstance(sum_scale_config, dict):
                file_sum = sum_scale_config.get(file, [20, 30])
                value_str = ", ".join(str(x) for x in file_sum)
            elif isinstance(sum_scale_config, str):
                value_str = sum_scale_config
            else:
                value_str = "20,30"

            self.sum_scale_input = widgets.Text(
                value=value_str,
                description='',
                layout=widgets.Layout(width='150px')
            )
            self.sum_scale_input.style = {'description_width': 'initial'}

            if i == 0:
                self.apply_button = widgets.Button(description="Apply to all", layout=widgets.Layout(width='120px'))
                self.apply_button.on_click(apply_to_all_callback)
            else:
                self.apply_button = widgets.Label(value="")

            sum_scale_widgets.append({
                'label': file_label,
                'input': self.sum_scale_input,
                'button': self.apply_button
            })
            row = widgets.HBox([file_label, self.sum_scale_input, self.apply_button])
            rows.append(row)

        self.sum_scale_widgets = sum_scale_widgets

        self.btn_unmix = widgets.Button(
            description="Run signal unmixing",
            button_style="primary",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )
        self.spinner_unmix = widgets.HTML(value=spinning_loader)
        self.spinner_unmix.layout.display = "none"

        self.unmix_row = widgets.HBox(
            [self.btn_unmix, self.spinner_unmix],
            layout=widgets.Layout(align_items="center", gap="8px")
        )

        self.channel_colors = tuple(channel_colors or _cfg_get(
            self.metadata_loader.behav3d_parameters, "signal_unmixing.channel_colors",
            ["cyan", "yellow", "red", "green", "magenta", "blue"]
        ))

        self._viewer = None
        self._status = widgets.HTML("<b>Waiting for user to load metadata…</b>")
        self._refresh_btn = widgets.Button(
            description="Refresh",
            tooltip="Build/Update selector from metadata_loader.metadata",
        )
        self._refresh_btn.on_click(self._on_refresh_clicked)

        self.sample_dropdown = widgets.Dropdown(
            options=[],
            value=None,
            description="Sample:",
            layout=widgets.Layout(width="350px"),
            disabled=True,
        )

        self.use_range = widgets.Checkbox(
            description="Use custom time range",
            value=bool(_cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing.use_range", False)),
            indent=False,
        )

        if isinstance(default_time_range, (tuple, list)) and len(default_time_range) == 2:
            _start_default, _end_default = map(int, default_time_range)
        else:
            _start_default = int(_cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing.start_t", 0))
            _end_default   = int(_cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing.end_t", 0))

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
            description="Visualize Signal Unmixing results",
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
            layout=widgets.Layout(width="200px", display="none"),
        )

        self._maybe_build_from_loader()

        self.btn_save = widgets.Button(
            description="Save Signal Unmixing",
            button_style="success",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )
        self.use_unmix_path = widgets.Checkbox(
            description="Use Unmix as raw path",
            value=bool(_cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing.use_unmix_path", False)),
            indent=False,
        )
        self.spinner_save = widgets.HTML(value=spinning_loader)
        self.spinner_save.layout.display = "none"
        self.save_row = widgets.HBox(
            [self.btn_save, self.spinner_save],
            layout=widgets.Layout(align_items="center", gap="8px")
        )

        self.btn_unmix.on_click(self._on_btn_unmix_clicked)
        self.open_button.on_click(self._on_open_clicked)
        self.close_button.on_click(self._on_close_clicked)
        self.btn_save.on_click(self._on_save_clicked)

        self.out = widgets.Output()

        unmix_box = widgets.VBox([
            widgets.HTML("<b>Run Signal Unmixing</b>"),
            widgets.VBox(rows),
            widgets.HBox([self.sink_channel, self.source_channels]),
            widgets.HBox([self.train_z, self.train_t]),
            widgets.HBox([self.bg_percentile, self.median_size, self.gaussian_sigma]),
            widgets.HBox([self.training_log]),
            widgets.HBox([self.use_all_timepoints, self.tp_n]),
            self.unmix_row,
        ])

        visualize_box = widgets.VBox([
                widgets.HTML("<b>Visualize Signal Unmixing Result</b>"),
                widgets.HBox([self._status, self._refresh_btn]),
                self.sample_dropdown,
                self.use_range,
                self.range_box,
                widgets.HBox([self.open_button, self.close_button]),
        ])
        
        save_box = widgets.VBox([
            widgets.HTML("<b>Save Signal Unmixing Result</b>"),
            widgets.HTML("Modifies the tcell channel to the new unmixed channel and saves a new metadata file."),
            self.use_unmix_path,
            self.save_row,
        ])

        self.ui = widgets.VBox([unmix_box, widgets.HTML("<hr>"), visualize_box, widgets.HTML("<hr>"), save_box, widgets.HTML("<hr>"), self.out])

    def display(self):
        widgets.display(self.ui)

    def _toggle_timepoint_inputs(self, change=None):
        show = not self.use_all_timepoints.value
        disp = None if show else 'none'
        self.tp_n.layout.display = disp
        self.tp_n.disabled = not show
        self.tp_n.value = 0 if self.use_all_timepoints.value else self.tp_n.value

    def _persist_params(self):
        pc = self.metadata_loader.behav3d_parameters.setdefault("signal_unmixing", {})
        pc["sink_channel"] = int(self.sink_channel.value)
        pc["source_channels"] = str(self.source_channels.value)
        pc["train_z"] = str(self.train_z.value)
        pc["train_t"] = str(self.train_t.value)
        pc["bg_percentile"] = int(self.bg_percentile.value)
        pc["median_size"] = int(self.median_size.value)
        pc["gaussian_sigma"]   = int(self.gaussian_sigma.value)
        pc["training_log"] = bool(self.training_log.value)
        pc["tp_n"] = int(self.tp_n.value)
        pc["use_all_timepoints"] = bool(self.use_all_timepoints.value)
        pc["timepoints"] = int(self.tp_n.value)
        
        sum_scale = {}
        for row in self.sum_scale_widgets:
            filename = row['label'].value
            input_str = row['input'].value
            try:
                parts = [int(s.strip()) for s in input_str.split(",") if s.strip().isdigit()]
                if len(parts) == 2:
                    sum_scale[filename] = parts
            except Exception:
                pass
        pc["sum_scale"] = sum_scale

        if hasattr(self.metadata_loader, "behav3d_parameters_path"):
            yaml.safe_dump(
                self.metadata_loader.behav3d_parameters,
                self.metadata_loader.behav3d_parameters_path.open("w"),
                sort_keys=False
            )

    def _on_btn_unmix_clicked(self, _):
        with self.out:
            self.out.clear_output()
            self._lock(True)
            previous_pc = deepcopy(_cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing", {}))
            self._persist_params()
            try:
                odir = Path(self.metadata_loader.output_dir).expanduser()
                odir.mkdir(parents=True, exist_ok=True)
                sum_scale_per_file = {}
                for row in self.sum_scale_widgets:
                    filename = row['label'].value
                    val = row['input'].value
                    try:
                        parts = [int(s.strip()) for s in val.split(",") if s.strip().isdigit()]
                        if len(parts) == 2:
                            sum_scale_per_file[filename] = parts
                    except Exception:
                        pass

                self.spinner_unmix.layout.display = None

                new_md = signal_unmixing(
                    metadata=self.metadata_loader.metadata,
                    output_dir=str(odir),
                    sink_channel=int(self.sink_channel.value),
                    source_channels=[int(s.strip()) for s in str(self.source_channels.value).split(",") if s.strip().isdigit()],
                    sum_scale=sum_scale_per_file,
                    train_z=[int(s.strip()) for s in str(self.train_z.value).split(",") if s.strip().isdigit()],
                    train_t=[int(s.strip()) for s in str(self.train_t.value).split(",") if s.strip().isdigit()],
                    bg_percentile=int(self.bg_percentile.value),
                    median_size=int(self.median_size.value),
                    gaussian_sigma=int(self.gaussian_sigma.value),
                    training_log=bool(self.training_log.value),
                    timepoints=int(self.tp_n.value),
                    previous_pc=previous_pc,
                )
                if new_md is not None:
                    self.metadata_loader.metadata = new_md
                    new_md.to_csv(self.metadata_loader.metadata_csv_path, index=False)
                print("✅ Unmixing finished.")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_unmix.layout.display = "none"
                self._lock(False)

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
        if self.use_range.value:
            start_t, end_t = int(self.start_t.value), int(self.end_t.value)
            if end_t < start_t:
                start_t, end_t = end_t, start_t
            time_range = (start_t, end_t)
        else:
            time_range = None

        with self.out:
            self.out.clear_output()
            try:
                print(f"Opening viewer for: {row['sample_name']}")
                self._viewer = visualize_unmix(
                    metadata_row=row,
                    timepoint_range=time_range,
                    channel_colors=self.channel_colors,
                )
            except Exception as e:
                print(f"Error: {e}")
            self.close_button.layout.display = "inline-block" if self._viewer is not None else "none"

    def _ensure_metadata_ready(self):
        df = getattr(self.metadata_loader, "metadata", None)
        if not isinstance(df, pd.DataFrame) or df.empty:
            raise RuntimeError("Metadata not loaded yet.")

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
            self.sample_dropdown.options = []
            self.sample_dropdown.value = None
            self.sample_dropdown.disabled = True
            self.open_button.disabled = True
            return

        desired = _cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing.sample_name", None)
        self.sample_dropdown.options = sample_names
        self.sample_dropdown.value = desired if (desired in sample_names) else sample_names[0]
        self.sample_dropdown.disabled = False
        self.open_button.disabled = False

    def _lock(self, state: bool):
        for w in [
            self.sample_dropdown, self.sink_channel, self.sum_scale_input, self.apply_button,
            self.source_channels, self.gaussian_sigma, self.median_size, self.bg_percentile,
            self.btn_save, self.btn_unmix, self.use_all_timepoints, self.tp_n, self.use_range, self.train_t, self.end_t,
            self.start_t, self.train_z, self.training_log, self.open_button, self.close_button, self.range_box, self.save_row, self.use_unmix_path
        ]:
            if hasattr(w, 'disabled'):
                w.disabled = state

    def _on_open_clicked(self, _):
        self._lock(True)
        self.open_viewer()
        self._lock(False)

    def _on_close_clicked(self, _):
        with self.out:
            try:
                if self._viewer is not None:
                    self._viewer.close()
                    self._viewer = None
            finally:
                self.close_button.layout.display = "none"
                
    def _on_refresh_clicked(self, _):
        with self.out:
            self.out.clear_output()
            try:
                self._maybe_build_from_loader()
            except Exception as e:
                print(f"Refresh failed: {e}")

    def _on_save_clicked(self, _):
        with self.out:
            self.out.clear_output()
            try:
                self._persist_params()
                metadata = self.metadata_loader.metadata
                for idx, sample in metadata.iterrows():
                    unmix_path = Path(sample["signal_unmixing_image_path"]).expanduser()
                    image = load_image(unmix_path)
                    new_tcell_channel = image.shape[1] - 1
                    print(f"Sample '{sample['sample_name']}': setting tcell_channel to {new_tcell_channel}") 
                    metadata.at[idx, 'tcell_channel'] = new_tcell_channel
                    if self.use_unmix_path:
                        if sample["signal_unmixing_image_path"] != metadata.at[idx, 'raw_image_path']:
                            metadata.at[idx, 'original_raw_image_path'] = sample["raw_image_path"]
                            metadata.at[idx, 'raw_image_path'] = sample["signal_unmixing_image_path"]
                
                if metadata is not None:
                    self.metadata_loader.metadata = metadata
                    metadata.to_csv(self.metadata_loader.metadata_csv_path, index=False)
                print(f"✅ Saved updated metadata")
            except Exception as e:
                print(f"Error saving metadata: {e}")
                traceback.print_exc()
