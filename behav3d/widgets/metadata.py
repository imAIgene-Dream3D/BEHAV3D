import ipywidgets as widgets
from IPython.display import display, clear_output
from pathlib import Path
import pandas as pd
import re
from .utils import (
    _load_config, 
    detect_cell_type_category
)
from behav3d.core.metadata import (
    load_behav3d_metadata,
    check_behav3d_metadata,
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata
)

from traitlets import Bool
from copy import deepcopy
import yaml
from behav3d.widgets.utils import (
    _DEFAULT_CONFIG,
    CONFIG,
     _load_config,
     _cfg_get
)

class MetadataBuilder(widgets.VBox):
    """
    Interactive widget for creating/editing BEHAV3D metadata CSV files.
    Supports:
    - Dynamic number of samples
    - Organoid, immune, and other cell type populations
    - Dead channel configuration (global setting)
    - Fill-down from sample 1
    - Load existing CSV for editing
    """
    def __init__(self, output_dir_picker=None, **kwargs):
        super().__init__(**kwargs)
        self.output_dir_picker = output_dir_picker
        
        # Population configuration
        self.n_organoid_types = 0
        self.n_immune_types = 0
        self.n_other_types = 0
        self.organoid_names = []
        self.immune_names = []
        self.other_names = []
        
        # Dead channel configuration
        self.include_dead_channel = False
        
        # Sample data storage
        self.sample_forms = []
        
        # Store loaded dataframe for re-population
        self.loaded_df = None
        
        self._build_ui()
    
    def _build_ui(self):
        """Build the initial configuration UI"""
        header = widgets.HTML('<h3>Metadata Builder</h3>')
        instructions = widgets.HTML(
            '<p>Define the number of samples and configure cell type populations. '
            'Fields marked with * are required.</p>'
        )
        
        self.n_samples_input = widgets.IntText(
            value=1,
            description='Number of samples*:',
            style={'description_width': '150px'}
        )
        
        self.btn_load_csv = widgets.Button(
            description='Load Existing CSV',
            button_style='info',
            icon='upload',
            tooltip='Load an existing metadata CSV to edit'
        )
        self.btn_load_csv.on_click(self._on_load_csv)
        
        self.btn_configure = widgets.Button(
            description='Configure Cell Types',
            button_style='primary',
            icon='cog'
        )
        self.btn_configure.on_click(self._on_configure_populations)
        
        self.main_container = widgets.VBox()
        
        self.children = [
            header,
            instructions,
            widgets.HBox([self.n_samples_input, self.btn_load_csv]),
            self.btn_configure,
            self.main_container
        ]
    
    def _on_load_csv(self, btn):
        """Load existing CSV and populate the form"""
        if self.output_dir_picker is None:
            self.main_container.children = [widgets.HTML('<p style="color:red">⚠️ No output directory picker provided</p>')]
            return
        
        output_dir = Path(self.output_dir_picker.value)
        if not output_dir.exists():
            self.main_container.children = [widgets.HTML('<p style="color:red">⚠️ Output directory does not exist</p>')]
            return
        
        csv_path = output_dir / 'metadata.csv'
        if not csv_path.exists():
            self.main_container.children = [widgets.HTML(f'<p style="color:red">⚠️ No metadata.csv found in {output_dir}</p>')]
            return
        
        try:
            df = pd.read_csv(csv_path)
            self._populate_from_dataframe(df)
            self.main_container.children = [widgets.HTML('<p style="color:green">✅ CSV loaded successfully! Scroll down to edit.</p>')]
        except Exception as e:
            self.main_container.children = [widgets.HTML(f'<p style="color:red">⚠️ Error loading CSV: {e}</p>')]
    
    def _populate_from_dataframe(self, df):
        """Populate the form from an existing DataFrame"""
        self.loaded_df = df.copy()
        
        organoid_types = detect_organoid_types_from_metadata(df)
        immune_types = detect_immune_cell_types_from_metadata(df)
        other_types = detect_other_cell_types_from_metadata(df)
               
        self.n_organoid_types = len(organoid_types)
        self.n_immune_types = len(immune_types)
        self.n_other_types = len(other_types)
        self.organoid_names = organoid_types
        self.immune_names = immune_types
        self.other_names = other_types
        
        # Check for dead channel
        if 'dead_channel' in df.columns and df['dead_channel'].notna().any():
            self.include_dead_channel = True
        else:
            self.include_dead_channel = False
        
        self.n_samples_input.value = len(df)
        self._build_data_entry_form()
    
    def _populate_form_fields(self, df):
        """Populate form fields from DataFrame"""
        for idx, row in df.iterrows():
            if idx >= len(self.sample_forms):
                break
            
            form = self.sample_forms[idx]
            
            for field_name, widget in form['basic'].items():
                if field_name in row and pd.notna(row[field_name]):
                    try:
                        if isinstance(widget, widgets.IntText):
                            widget.value = int(row[field_name])
                        elif isinstance(widget, widgets.FloatText):
                            widget.value = float(row[field_name])
                        else:
                            widget.value = str(row[field_name])
                    except (ValueError, TypeError):
                        widget.value = str(row[field_name]) if row[field_name] != '' else ''
            
            for cell_type, fields_dict in form['cell_types'].items():
                if cell_type in self.organoid_names:
                    prefix = 'or'
                elif cell_type in self.immune_names:
                    prefix = 'im'
                elif cell_type in self.other_names:
                    prefix = 'ot'
                else:
                    continue
                
                line_cond_col = f'{prefix}_{cell_type}_line_condition'
                if line_cond_col in row and pd.notna(row[line_cond_col]) and str(row[line_cond_col]).strip() != '':
                    value = str(row[line_cond_col]).strip()
                    parts = value.split('_', 1)
                    if 'line' in fields_dict:
                        fields_dict['line'].value = parts[0] if parts[0] else ''
                    if 'condition' in fields_dict and len(parts) > 1:
                        fields_dict['condition'].value = parts[1] if parts[1] else ''
                    elif 'condition' in fields_dict:
                        fields_dict['condition'].value = ''
                
                for path_field in ['segments_image_path', 'tracks_image_path', 'tracks_csv_path']:
                    col_name = f'{prefix}_{cell_type}_{path_field}'
                    if col_name in row and pd.notna(row[col_name]) and path_field in fields_dict:
                        fields_dict[path_field].value = str(row[col_name]).strip()
            
            # Populate dead channel fields
            if self.include_dead_channel:
                if 'dead_channel' in row and pd.notna(row['dead_channel']) and str(row['dead_channel']).strip() != '':
                    try:
                        form['dead_channel']['number'].value = int(row['dead_channel'])
                    except (ValueError, TypeError):
                        pass
                if 'dead_mask_path' in row and pd.notna(row['dead_mask_path']) and 'mask_path' in form['dead_channel']:
                    form['dead_channel']['mask_path'].value = str(row['dead_mask_path']).strip()
    
    def _on_configure_populations(self, btn):
        """Show population configuration UI"""
        org_label = widgets.HTML('<h4>Organoid Populations</h4>')
        self.n_organoid_input = widgets.IntText(
            value=self.n_organoid_types,
            description='Number of types:',
            style={'description_width': '120px'}
        )
        
        immune_label = widgets.HTML('<h4>Immune Cell Populations</h4>')
        self.n_immune_input = widgets.IntText(
            value=self.n_immune_types,
            description='Number of types:',
            style={'description_width': '120px'}
        )
        
        other_label = widgets.HTML('<h4>Other Cell Types</h4>')
        self.n_other_input = widgets.IntText(
            value=self.n_other_types,
            description='Number of types:',
            style={'description_width': '120px'}
        )
        
        # Dead channel checkbox (moved here from sample form)
        dead_label = widgets.HTML('<h4>Dead Channel</h4>')
        self.include_dead_checkbox = widgets.Checkbox(
            description='Include dead channel?',
            value=self.include_dead_channel,
            style={'description_width': '150px'}
        )
        
        btn_confirm = widgets.Button(description='Next: Name Cell Types', button_style='success')
        btn_confirm.on_click(self._show_cell_type_naming)
        
        self.main_container.children = [
            org_label,
            self.n_organoid_input,
            immune_label,
            self.n_immune_input,
            other_label,
            self.n_other_input,
            dead_label,
            self.include_dead_checkbox,
            btn_confirm
        ]
    
    def _show_cell_type_naming(self, btn):
        """Show UI to name each cell type"""
        self.n_organoid_types = self.n_organoid_input.value
        self.n_immune_types = self.n_immune_input.value
        self.n_other_types = self.n_other_input.value
        
        self.include_dead_channel = self.include_dead_checkbox.value
        
        naming_widgets = []
        self.organoid_name_inputs = []
        self.immune_name_inputs = []
        self.other_name_inputs = []
        
        if self.n_organoid_types > 0:
            naming_widgets.append(widgets.HTML('<h4>Name Your Organoid Types</h4>'))
            for i in range(self.n_organoid_types):
                default = self.organoid_names[i] if i < len(self.organoid_names) else f'organoid{i+1}'
                w = widgets.Text(
                    value=default,
                    description=f'Organoid {i+1}:',
                    placeholder='e.g., organoid1, organoidWT',
                    style={'description_width': '100px'}
                )
                self.organoid_name_inputs.append(w)
                naming_widgets.append(w)
        
        if self.n_immune_types > 0:
            naming_widgets.append(widgets.HTML('<h4>Name Your Immune Cell Types</h4>'))
            for i in range(self.n_immune_types):
                default = self.immune_names[i] if i < len(self.immune_names) else f'tcell{i+1}' if i == 0 else f'immune{i+1}'
                w = widgets.Text(
                    value=default,
                    description=f'Immune {i+1}:',
                    placeholder='e.g., tcell, nk, macro',
                    style={'description_width': '100px'}
                )
                self.immune_name_inputs.append(w)
                naming_widgets.append(w)
        
        if self.n_other_types > 0:
            naming_widgets.append(widgets.HTML('<h4>Name Your Other Cell Types</h4>'))
            for i in range(self.n_other_types):
                default = self.other_names[i] if i < len(self.other_names) else f'other{i+1}'
                w = widgets.Text(
                    value=default,
                    description=f'Other {i+1}:',
                    placeholder='e.g., tumor, fibroblast',
                    style={'description_width': '100px'}
                )
                self.other_name_inputs.append(w)
                naming_widgets.append(w)
        
        btn_next = widgets.Button(description='Create Data Entry Form', button_style='success')
        btn_next.on_click(self._on_names_confirmed)
        naming_widgets.append(btn_next)
        
        self.main_container.children = naming_widgets
    
    def _on_names_confirmed(self, btn):
        """Collect names and build data entry form"""
        self.organoid_names = [w.value.strip() for w in self.organoid_name_inputs]
        self.immune_names = [w.value.strip() for w in self.immune_name_inputs]
        self.other_names = [w.value.strip() for w in self.other_name_inputs]
        
        self._build_data_entry_form()
    
    def _build_data_entry_form(self):
        """Build the main data entry form for all samples"""
        n_samples = self.n_samples_input.value
        self.sample_forms = []
        
        form_widgets = [widgets.HTML('<h3>Sample Data Entry</h3>')]
        
        btn_fill_down = widgets.Button(
            description='Fill All from Sample 1',
            button_style='warning',
            icon='arrow-down',
            tooltip='Copy all values from Sample 1 to all other samples (except sample names)'
        )
        btn_fill_down.on_click(self._on_fill_down)
        form_widgets.append(btn_fill_down)
        
        for i in range(n_samples):
            sample_form = self._create_sample_form(i)
            self.sample_forms.append(sample_form)
            form_widgets.append(sample_form['widget'])
        
        btn_save = widgets.Button(
            description='Save Metadata to CSV',
            button_style='success',
            icon='save'
        )
        btn_save.on_click(self._on_save_metadata)
        
        self.save_output = widgets.Output()
        
        form_widgets.extend([btn_save, self.save_output])
        self.main_container.children = form_widgets
        
        if self.loaded_df is not None:
            self._populate_form_fields(self.loaded_df)
    
    def _on_fill_down(self, btn):
        """Copy all values from Sample 1 to all other samples (except sample_name)"""
        if len(self.sample_forms) < 2:
            return
        
        source_form = self.sample_forms[0]
        
        for i in range(1, len(self.sample_forms)):
            target_form = self.sample_forms[i]
            
            # Copy basic fields (except sample_name)
            for field_name, src_widget in source_form['basic'].items():
                if field_name != 'sample_name':
                    target_form['basic'][field_name].value = src_widget.value
            
            # Copy cell type fields
            for cell_type in source_form['cell_types']:
                for field_name, src_widget in source_form['cell_types'][cell_type].items():
                    target_form['cell_types'][cell_type][field_name].value = src_widget.value
            
            # Copy dead channel fields
            if self.include_dead_channel:
                if 'mask_path' in source_form['dead_channel'] and 'mask_path' in target_form['dead_channel']:
                    target_form['dead_channel']['mask_path'].value = source_form['dead_channel']['mask_path'].value
                target_form['dead_channel']['number'].value = source_form['dead_channel']['number'].value
        
        with self.save_output:
            self.save_output.clear_output()
            print('✅ Filled all samples from Sample 1!')
    
    def _create_sample_form(self, sample_idx):
        """Create form widgets for a single sample"""
        form_data = {
            'basic': {},
            'cell_types': {},
            'dead_channel': {},
            'widget': None
        }
        
        header = widgets.HTML(f'<h4 style="background:#e1f5fe;padding:8px;">Sample {sample_idx + 1}</h4>')
        
        # Basic information widgets
        form_data['basic']['sample_name'] = widgets.Text(
            description='Sample name*:',
            placeholder='e.g., Sample001',
            style={'description_width': '150px'},
            layout=widgets.Layout(width='400px')
        )
        
        form_data['basic']['exp_nr'] = widgets.IntText(
            description='Exp number*:',
            value=1,
            style={'description_width': '150px'}
        )
        
        form_data['basic']['well'] = widgets.Text(
            description='Well*:',
            placeholder='e.g., A1',
            style={'description_width': '150px'}
        )
        
        form_data['basic']['raw_image_path'] = widgets.Text(
            description='Raw image path*:',
            placeholder='/path/to/image.czi/.ims/.lif/.liff/.tiff/.zarr',
            style={'description_width': '150px'},
            layout=widgets.Layout(width='600px')
        )
        
        form_data['basic']['dimension_order'] = widgets.Text(
            description='Dimension order:',
            placeholder='Optional - e.g., TCZYX',
            style={'description_width': '150px'}
        )
        
        # Imaging parameters
        form_data['basic']['pixel_distance_xy'] = widgets.FloatText(
            description='Pixel xy (μm)*:',
            value=0.5,
            style={'description_width': '150px'}
        )
        
        form_data['basic']['pixel_distance_z'] = widgets.FloatText(
            description='Pixel z (μm)*:',
            value=2.0,
            style={'description_width': '150px'}
        )
        
        form_data['basic']['distance_unit'] = widgets.Text(
            description='Distance unit*:',
            value='μm',
            style={'description_width': '150px'}
        )
        
        form_data['basic']['time_interval'] = widgets.FloatText(
            description='Time interval*:',
            value=1.0,
            style={'description_width': '150px'}
        )
        
        form_data['basic']['time_unit'] = widgets.Text(
            description='Time unit*:',
            value='s',
            style={'description_width': '150px'}
        )
        
        # Dead channel widgets
        dead_channel_widgets = []
        dead_number = None
        if self.include_dead_channel:
            dead_channel_label = widgets.HTML('<h5>Dead Channel Configuration</h5>')
            
            dead_mask_path = widgets.Text(
                description='Dead mask path:',
                placeholder='Optional - path to dead cell mask (filled by pipeline)',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='600px')
            )
            
            # Dead channel number - IntText
            dead_number = widgets.IntText(
                description='Dead channel:',
                value=0,
                style={'description_width': '120px'}
            )
            
            form_data['dead_channel']['mask_path'] = dead_mask_path
            form_data['dead_channel']['number'] = dead_number
            
            dead_channel_widgets = [dead_channel_label, dead_mask_path, dead_number]
        
        # Cell type configuration widgets
        cell_type_widgets = []
        
        for org_name in self.organoid_names:
            cell_type_widgets.append(widgets.HTML(f'<h5>{org_name.capitalize()}</h5>'))
            fields = {}
            fields['line'] = widgets.Text(
                description='Line*:',
                placeholder='e.g., 32',
                style={'description_width': '100px'}
            )
            fields['condition'] = widgets.Text(
                description='Condition*:',
                placeholder='e.g., ko, wt',
                style={'description_width': '100px'}
            )
            fields['segments_image_path'] = widgets.Text(
                description='Segments path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_image_path'] = widgets.Text(
                description='Tracks img path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_csv_path'] = widgets.Text(
                description='Tracks csv path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            form_data['cell_types'][org_name] = fields
            cell_type_widgets.extend([
                fields['line'], 
                fields['condition'],
                fields['segments_image_path'],
                fields['tracks_image_path'],
                fields['tracks_csv_path']
            ])
        
        for immune_name in self.immune_names:
            cell_type_widgets.append(widgets.HTML(f'<h5>{immune_name.capitalize()}</h5>'))
            fields = {}
            fields['line'] = widgets.Text(
                description='Line*:',
                placeholder='e.g., CD8',
                style={'description_width': '100px'}
            )
            fields['condition'] = widgets.Text(
                description='Condition*:',
                placeholder='e.g., activated',
                style={'description_width': '100px'}
            )
            fields['segments_image_path'] = widgets.Text(
                description='Segments path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_image_path'] = widgets.Text(
                description='Tracks img path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_csv_path'] = widgets.Text(
                description='Tracks csv path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            form_data['cell_types'][immune_name] = fields
            cell_type_widgets.extend([
                fields['line'], 
                fields['condition'],
                fields['segments_image_path'],
                fields['tracks_image_path'],
                fields['tracks_csv_path']
            ])
        
        for other_name in self.other_names:
            cell_type_widgets.append(widgets.HTML(f'<h5>{other_name.capitalize()}</h5>'))
            fields = {}
            fields['line'] = widgets.Text(
                description='Line*:',
                placeholder='e.g., WT',
                style={'description_width': '100px'}
            )
            fields['condition'] = widgets.Text(
                description='Condition*:',
                placeholder='e.g., control',
                style={'description_width': '100px'}
            )
            fields['segments_image_path'] = widgets.Text(
                description='Segments path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_image_path'] = widgets.Text(
                description='Tracks img path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_csv_path'] = widgets.Text(
                description='Tracks csv path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            form_data['cell_types'][other_name] = fields
            cell_type_widgets.extend([
                fields['line'], 
                fields['condition'],
                fields['segments_image_path'],
                fields['tracks_image_path'],
                fields['tracks_csv_path']
            ])
        
        # Build form layout
        form_children = [
            header,
            widgets.HTML('<h5>Basic Information</h5>'),
            form_data['basic']['sample_name'],
            form_data['basic']['exp_nr'],
            form_data['basic']['well'],
            form_data['basic']['raw_image_path'],
            form_data['basic']['dimension_order'],
            widgets.HTML('<h5>Imaging Parameters</h5>'),
            form_data['basic']['pixel_distance_xy'],
            form_data['basic']['pixel_distance_z'],
            form_data['basic']['distance_unit'],
            form_data['basic']['time_interval'],
            form_data['basic']['time_unit'],
        ]
                
        # Add dead channel configuration
        form_children.extend(dead_channel_widgets)
        
        # Add cell type configuration
        form_children.append(widgets.HTML('<h5>Cell Type Configuration</h5>'))
        form_children.extend(cell_type_widgets)
        form_children.append(widgets.HTML('<hr>'))
        
        form_widget = widgets.VBox(form_children)
        form_data['widget'] = form_widget
        return form_data
    
    def _on_save_metadata(self, btn):
        """Save all form data to CSV"""
        if self.output_dir_picker is None:
            with self.save_output:
                self.save_output.clear_output()
                print('⚠️ No output directory picker provided')
            return
        
        output_dir = Path(self.output_dir_picker.value)
        if not output_dir.exists():
            with self.save_output:
                self.save_output.clear_output()
                print(f'⚠️ Output directory does not exist: {output_dir}')
            return
        
        rows = []
        for form in self.sample_forms:
            row = {}
                        
            # Save basic fields
            for field_name, widget in form['basic'].items():
                row[field_name] = widget.value
            
            # Save dead channel fields (only if dead channel is included)
            if self.include_dead_channel:
                if 'mask_path' in form['dead_channel']:
                    row['dead_mask_path'] = form['dead_channel']['mask_path'].value.strip()
                row['dead_channel'] = form['dead_channel']['number'].value
            
            # Save cell type fields
            for cell_type, fields in form['cell_types'].items():
                if cell_type in self.organoid_names:
                    prefix = 'or'
                elif cell_type in self.immune_names:
                    prefix = 'im'
                elif cell_type in self.other_names:
                    prefix = 'ot'
                else:
                    continue
                
                line = fields['line'].value.strip()
                condition = fields['condition'].value.strip()
                col_name = f'{prefix}_{cell_type}_line_condition'
                if line and condition:
                    row[col_name] = f'{line}_{condition}'
                elif line:
                    row[col_name] = line
                elif condition:
                    row[col_name] = f'_{condition}'
                else:
                    row[col_name] = ''
                
                row[f'{prefix}_{cell_type}_segments_image_path'] = fields.get('segments_image_path', widgets.Text()).value.strip()
                row[f'{prefix}_{cell_type}_tracks_image_path'] = fields.get('tracks_image_path', widgets.Text()).value.strip()
                row[f'{prefix}_{cell_type}_tracks_csv_path'] = fields.get('tracks_csv_path', widgets.Text()).value.strip()
            
            rows.append(row)
        
        df = pd.DataFrame(rows)
        
        # Reorder columns in the desired order
        ordered_columns = []
        
        # 1. Basic columns first
        basic_order = ['sample_name', 'exp_nr', 'well', 'raw_image_path', 'dimension_order',
                       'pixel_distance_xy', 'pixel_distance_z', 'distance_unit', 
                       'time_interval', 'time_unit']
        for col in basic_order:
            if col in df.columns:
                ordered_columns.append(col)
        
        # 2. Dead channel columns (if included)
        if self.include_dead_channel:
            if 'dead_mask_path' in df.columns:
                ordered_columns.append('dead_mask_path')
            if 'dead_channel' in df.columns:
                ordered_columns.append('dead_channel')
        
        # 4. Cell type columns: or_* then im_* then ot_*
        for org_name in self.organoid_names:
            for suffix in ['_line_condition', '_segments_image_path', '_tracks_image_path', '_tracks_csv_path']:
                col = f'or_{org_name}{suffix}'
                if col in df.columns:
                    ordered_columns.append(col)
        
        for immune_name in self.immune_names:
            for suffix in ['_line_condition', '_segments_image_path', '_tracks_image_path', '_tracks_csv_path']:
                col = f'im_{immune_name}{suffix}'
                if col in df.columns:
                    ordered_columns.append(col)
        
        for other_name in self.other_names:
            for suffix in ['_line_condition', '_segments_image_path', '_tracks_image_path', '_tracks_csv_path']:
                col = f'ot_{other_name}{suffix}'
                if col in df.columns:
                    ordered_columns.append(col)
        
        # Add any remaining columns not in the ordered list
        for col in df.columns:
            if col not in ordered_columns:
                ordered_columns.append(col)
        
        df = df[ordered_columns]
        
        csv_path = output_dir / 'metadata.csv'
        df.to_csv(csv_path, index=False)
        
        with self.save_output:
            self.save_output.clear_output()
            print(f'✅ Metadata saved to: {csv_path}')
            print(f'{len(df)} samples, {len(df.columns)} columns')
            if self.include_dead_channel:
                print(f'Dead channel: column "dead_channel" and "dead_mask_path" included')
            display(df)


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
        func=False,
        button_description: str = "Load metadata",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.metadata = None
        self.output_dir = None
        self.metadata_csv_path = None

        self.func = func
        self.file_picker = metadata_path_picker
        self.output_dir_picker = output_dir_picker

        # initialize pickers with JSON defaults if empty
        try:
            if hasattr(self.file_picker, "value") and not self.file_picker.value:
                self.file_picker.value = _cfg_get(CONFIG, "paths.metadata_csv", "") or ""
        except Exception:
            pass
        try:
            if hasattr(self.output_dir_picker, "value") and not self.output_dir_picker.value:
                self.output_dir_picker.value = _cfg_get(CONFIG, "paths.output_dir", "") or ""
        except Exception:
            pass
        
        self.button = widgets.Button(description=button_description, button_style='success', layout=widgets.Layout(
        width="fit-content",   # size to content
        flex="0 0 auto"        # don't stretch in HBox/VBox
    ))
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
            
            self.behav3d_parameters_path = Path(self.output_dir, "behav3d_parameters.yml")
            
            if self.behav3d_parameters_path.exists():
                self.behav3d_parameters = _load_config(path = self.behav3d_parameters_path)
            else:
                self.behav3d_parameters = deepcopy(_DEFAULT_CONFIG)
                
            self.behav3d_parameters["paths"]["metadata_csv"] = str(self.metadata_csv_path)
            self.behav3d_parameters["paths"]["output_dir"]   = str(self.output_dir)
            yaml.safe_dump(self.behav3d_parameters, self.behav3d_parameters_path.open("w"), sort_keys=False)
            
            check_behav3d_metadata(self.metadata, self.func)
            print("✅ Checks passed!")
            display(self.metadata)

            
    # Button handler
    def _on_click(self, _):
        if self._busy:
            return  # re-entrancy guard prevents double execution
        self._busy = True
        self.load()
        # try:
        #     self.load()
        # except Exception as e:
        #     with self.out:
        #         print(f"❌ Error: {e}")
        # finally:
        self._busy = False
