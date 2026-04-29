"""
GroupBuilder widget for creating cell type groups in BEHAV3D.

Allows users to select existing detected cell types and group them under a single name.
Groups are stored as additional metadata columns with the _merged suffix.
"""

import ipywidgets as widgets
import pandas as pd
from IPython.display import display
from pathlib import Path

from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
    detect_merged_cell_types_from_metadata,
)


def _get_cell_type_category(cell_type, metadata):
    """Determine category (organoid/immune/other) for a cell type."""
    organoid_types = detect_organoid_types_from_metadata(metadata)
    immune_types = detect_immune_cell_types_from_metadata(metadata)
    other_types = detect_other_cell_types_from_metadata(metadata)
    
    if cell_type in organoid_types:
        return 'organoid'
    elif cell_type in immune_types:
        return 'immune'
    elif cell_type in other_types:
        return 'other'
    return None


class GroupBuilder(widgets.VBox):
    """
    Interactive widget for creating cell type groups.
    
    Groups multiple detected cell types under a single name, which is then
    added to metadata as a merged type. Groups are persisted to the metadata CSV.
    
    Parameters
    ----------
    metadata_loader : object
        Object with attributes:
        - metadata : pandas.DataFrame
        - metadata_csv_path : str or None
        - output_dir : str
    """
    
    def __init__(self, metadata_loader, **kwargs):
        super().__init__(**kwargs)
        self.metadata_loader = metadata_loader
        self._build_ui()
    
    def _build_ui(self):
        """Build the widget UI."""
        # Header
        header = widgets.HTML('<h4>Group Builder</h4>')
        instructions = widgets.HTML(
            '<p>Select existing cell types to group, give the group a name, '
            'and click <b>Create Group</b> to persist it to metadata.</p>'
        )
        
        # Get available cell types
        metadata = self.metadata_loader.metadata
        if metadata is None or metadata.empty:
            error_msg = widgets.HTML(
                '<div style="color:red;"><b>Error:</b> No metadata loaded.</div>'
            )
            self.children = [header, error_msg]
            return
        
        organoid_types = detect_organoid_types_from_metadata(metadata)
        immune_types = detect_immune_cell_types_from_metadata(metadata)
        other_types = detect_other_cell_types_from_metadata(metadata)
        existing_groups = detect_merged_cell_types_from_metadata(metadata)
        
        all_cell_types = sorted(organoid_types + immune_types + other_types)
        
        if not all_cell_types:
            warning_msg = widgets.HTML(
                '<div style="color:orange;"><b>Warning:</b> No cell types found in metadata.</div>'
            )
            self.children = [header, warning_msg]
            return
        
        # Cell type selector (multi-select checkboxes)
        selector_label = widgets.HTML('<b>Select cell types to group:</b>')
        self.cell_type_checkboxes = {
            ct: widgets.Checkbox(value=False, description=ct, indent=False)
            for ct in all_cell_types
        }
        checkbox_items = widgets.VBox([
            self.cell_type_checkboxes[ct] for ct in all_cell_types
        ])
        
        # Group name input
        group_name_label = widgets.HTML('<b>Group name*:</b>')
        self.group_name_input = widgets.Text(
            placeholder='e.g., tcells, tumors, etc.',
            description='',
            style={'description_width': '0px'},
            layout=widgets.Layout(width='300px')
        )
        
        # Create button
        self.create_btn = widgets.Button(
            description='Create Group',
            button_style='success',
            tooltip='Create the group and add it to metadata',
            icon='check'
        )
        self.create_btn.on_click(self._on_create_group)
        
        # Status output
        self.status_output = widgets.HTML('')
        
        # Layout
        self.children = [
            header,
            instructions,
            selector_label,
            checkbox_items,
            widgets.VBox([group_name_label, self.group_name_input]),
            self.create_btn,
            self.status_output
        ]
    
    def _validate_group_name(self, group_name):
        """Validate group name format and uniqueness."""
        group_name = str(group_name).strip()
        
        if not group_name:
            return False, "Group name cannot be empty."
        
        # Check for invalid characters
        if not all(c.isalnum() or c in '_-' for c in group_name):
            return False, "Group name can only contain alphanumeric characters, hyphens, and underscores."
        
        # Check if group already exists
        existing_groups = detect_merged_cell_types_from_metadata(self.metadata_loader.metadata)
        if group_name in existing_groups:
            return False, f"Group '{group_name}' already exists in metadata."
        
        # Check if it conflicts with existing cell types
        organoid_types = detect_organoid_types_from_metadata(self.metadata_loader.metadata)
        immune_types = detect_immune_cell_types_from_metadata(self.metadata_loader.metadata)
        other_types = detect_other_cell_types_from_metadata(self.metadata_loader.metadata)
        all_cell_types = organoid_types + immune_types + other_types
        
        if group_name in all_cell_types:
            return False, f"Group name '{group_name}' conflicts with an existing cell type."
        
        return True, ""
    
    def _get_selected_cell_types(self):
        """Get list of selected cell types from checkboxes."""
        return [ct for ct, cb in self.cell_type_checkboxes.items() if cb.value]
    
    def _on_create_group(self, btn):
        """Callback for Create Group button."""
        try:
            # Validate inputs
            group_name = self.group_name_input.value.strip()
            is_valid, error_msg = self._validate_group_name(group_name)
            if not is_valid:
                self.status_output.value = f'<div style="color:red;"><b>Error:</b> {error_msg}</div>'
                return
            
            selected_types = self._get_selected_cell_types()
            if not selected_types:
                self.status_output.value = '<div style="color:red;"><b>Error:</b> Select at least one cell type.</div>'
                return
            
            # Get metadata
            md = self.metadata_loader.metadata
            if md is None or md.empty:
                self.status_output.value = '<div style="color:red;"><b>Error:</b> No metadata loaded.</div>'
                return
            
            # Determine dominant category (organoid > immune > other)
            organoid_types = detect_organoid_types_from_metadata(md)
            immune_types = detect_immune_cell_types_from_metadata(md)
            other_types = detect_other_cell_types_from_metadata(md)
            
            selected_organoid = [t for t in selected_types if t in organoid_types]
            selected_immune = [t for t in selected_types if t in immune_types]
            selected_other = [t for t in selected_types if t in other_types]
            
            if selected_organoid:
                prefix = 'or_'
            elif selected_immune:
                prefix = 'im_'
            else:
                prefix = 'ot_'
            
            # Create merged column names
            merged_name = f"{group_name}_merged"
            col_prefix = f"{prefix}{merged_name}"
            
            # Add columns to metadata for each suffix pattern (line_condition, segments_image_path, etc.)
            suffixes = [
                'line_condition',
                'segments_image_path',
                'tracks_image_path',
                'tracks_csv_path'
            ]
            
            for suffix in suffixes:
                col_name = f"{col_prefix}_{suffix}"
                if col_name not in md.columns:
                    md[col_name] = pd.NA
            
            # Optionally save metadata to CSV
            csv_path = getattr(self.metadata_loader, 'metadata_csv_path', None)
            if csv_path:
                md.to_csv(csv_path, index=False)
                saved_msg = f"<br/>Metadata saved to: <code>{csv_path}</code>"
            else:
                saved_msg = ""
            
            # Display success message
            success_html = (
                f'<div style="color:green;"><b>✓ Group created:</b> <code>{merged_name}</code></div>'
                f'<div style="margin-top:8px; font-size:0.9em;">'
                f'<b>Members:</b> {", ".join(selected_types)}<br/>'
                f'<b>Columns added:</b><br/>'
            )
            
            for suffix in suffixes:
                col_name = f"{col_prefix}_{suffix}"
                success_html += f'&nbsp;&nbsp;- <code>{col_name}</code><br/>'
            
            success_html += f'</div>{saved_msg}'
            self.status_output.value = success_html
            
            # Print to console for reference
            print(f"✓ Group '{merged_name}' created from: {selected_types}")
            print(f"  Columns added: {[f'{col_prefix}_{s}' for s in suffixes]}")
            
            # Clear inputs
            self.group_name_input.value = ''
            for cb in self.cell_type_checkboxes.values():
                cb.value = False
        
        except Exception as e:
            self.status_output.value = (
                f'<div style="color:red;"><b>Error:</b> {str(e)}</div>'
            )
            print(f"Error in GroupBuilder._on_create_group: {e}")
            import traceback
            traceback.print_exc()
