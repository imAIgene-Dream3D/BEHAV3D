"""GroupBuilder widget for creating cell type groups in BEHAV3D.

Groups merge several already-*filtered* cell types' track-features CSVs
into a single pseudo cell-type population (``{name}_merged``) for
downstream analysis (Death Dynamics, Interaction Analysis, ...). Group
membership is recorded in ``behav3d_parameters.yml`` under
``cell_type_groups`` — metadata.csv is never touched, so this widget
should be used *after* Filtering, not before Feature Extraction.
"""

import ipywidgets as widgets
import yaml
from pathlib import Path

from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
)
from behav3d.analysis.grouping import (
    create_cell_type_group,
    create_group_tracked_segments,
    filtered_track_features_csv,
    group_cell_type_id,
    group_tracked_segments_status,
    list_cell_type_groups,
    remove_cell_type_group,
    validate_group_name,
)


class GroupBuilder(widgets.VBox):
    """
    Interactive widget for creating cell type groups.

    Merges each selected (already-filtered) cell type's filtered
    track-features CSV into a single population, selectable as a
    pseudo cell-type in downstream analysis panels.

    Parameters
    ----------
    metadata_loader : object
        Object with attributes:
        - metadata : pandas.DataFrame
        - output_dir : str
        - behav3d_parameters : dict
        - behav3d_parameters_path : Path (where to persist config)
    """

    def __init__(self, metadata_loader, **kwargs):
        super().__init__(**kwargs)
        self.metadata_loader = metadata_loader
        self._build_ui()

    def _params(self) -> dict:
        return getattr(self.metadata_loader, "behav3d_parameters", {}) or {}

    def _persist_params(self):
        params = self._params()
        params_path = getattr(self.metadata_loader, "behav3d_parameters_path", None)
        if params_path is None:
            out_dir = getattr(self.metadata_loader, "output_dir", None)
            if out_dir:
                params_path = Path(out_dir) / "behav3d_parameters.yml"
        if not params_path:
            return None
        try:
            with open(params_path, "w", encoding="utf-8") as f:
                yaml.safe_dump(params, f, sort_keys=False)
            return str(params_path)
        except Exception as e:
            self.status_output.value = f'<div style="color:red;"><b>Error:</b> Could not save config: {e}</div>'
            return None

    def _build_ui(self):
        """Build the widget UI."""
        header = widgets.HTML('<h4>Group Builder</h4>')
        instructions = widgets.HTML(
            '<p>Select cell types that have already been <b>filtered</b> above, '
            'give the group a name, and click <b>Create Group</b>. The group '
            'merges each member\'s filtered track-features CSV into one '
            'population, usable as a cell type in the analysis panels below '
            '(e.g. Death Dynamics, Interaction Analysis). Metadata is not '
            'modified.</p>'
        )

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
        self._all_cell_types = sorted(organoid_types + immune_types + other_types)

        if not self._all_cell_types:
            warning_msg = widgets.HTML(
                '<div style="color:orange;"><b>Warning:</b> No cell types found in metadata.</div>'
            )
            self.children = [header, warning_msg]
            return

        existing_groups = list_cell_type_groups(self._params())
        existing_html = (
            "<i>Existing groups:</i> " + ", ".join(
                f"{gid} ({', '.join(members)})" for gid, members in existing_groups.items()
            )
            if existing_groups else ""
        )
        self.existing_groups_label = widgets.HTML(existing_html)

        output_dir = getattr(self.metadata_loader, "output_dir", None)
        selector_label = widgets.HTML('<b>Select cell types to group (must be filtered already):</b>')
        self.cell_type_checkboxes = {}
        checkbox_rows = []
        for ct in self._all_cell_types:
            has_filtered = bool(
                output_dir and filtered_track_features_csv(output_dir, ct).exists()
            )
            desc = ct if has_filtered else f"{ct} (not yet filtered)"
            cb = widgets.Checkbox(value=False, description=desc, indent=False, disabled=not has_filtered)
            checkbox_rows.append(cb)
            self.cell_type_checkboxes[ct] = cb
        # Existing groups can be inspected/removed via the same checklist.
        self.group_checkboxes = {}
        for gid in existing_groups:
            desc = f"{gid} (group)"
            if output_dir:
                status = group_tracked_segments_status(output_dir, gid, metadata)
                if status["total"] > 0:
                    desc += f"  [tracked: {status['built']}/{status['total']}]"
            cb = widgets.Checkbox(value=False, description=desc, indent=False)
            checkbox_rows.append(cb)
            self.group_checkboxes[gid] = cb
        checkbox_items = widgets.VBox(checkbox_rows)

        group_name_label = widgets.HTML('<b>Group name*:</b>')
        self.group_name_input = widgets.Text(
            placeholder='e.g., tcells, tumors, etc.',
            description='',
            style={'description_width': '0px'},
            layout=widgets.Layout(width='300px')
        )

        self.build_tracked_checkbox = widgets.Checkbox(
            value=False,
            description='Also build tracked segments (slower)',
            indent=False,
            tooltip=(
                'Rebuilds a per-sample tracked label image + tracks CSV for '
                'the group by relabeling each member\'s existing tracked '
                'images. Needed to see the group in Visualization/'
                'backprojection. This rewrites a label image per sample and '
                'can take a while.'
            ),
        )

        self.create_btn = widgets.Button(
            description='Create Group',
            button_style='success',
            tooltip='Merge the selected filtered cell types into a new group',
            icon='check'
        )
        self.create_btn.on_click(self._on_create_group)

        self.remove_btn = widgets.Button(
            description='Remove Group',
            button_style='danger',
            tooltip='Remove the selected group from the config (keeps the merged CSV on disk)',
            icon='trash'
        )
        self.remove_btn.on_click(self._on_remove_group)

        self.preview_btn = widgets.Button(
            description='📊 Track Length Distribution',
            tooltip='Preview the selected group\'s track-length distribution',
        )
        self.preview_btn.on_click(self._on_preview_lengths)

        self.build_tracked_btn = widgets.Button(
            description='🧩 Build Tracked Segments',
            tooltip=(
                'Build (or rebuild) per-sample tracked label images + tracks '
                'CSVs for the selected existing group. Select exactly one '
                '"(group)" entry above.'
            ),
        )
        self.build_tracked_btn.on_click(self._on_build_tracked)

        self.status_output = widgets.HTML('')
        self.preview_output = widgets.Output()

        # Progress bar for the (blocking, in-kernel) tracked-segments build.
        # ipywidgets flushes .value updates to the front end as they happen
        # even mid-call, so this renders live despite running synchronously.
        self.tracked_progress = widgets.IntProgress(
            min=0, max=1, value=0, description='',
            layout=widgets.Layout(width='400px', visibility='hidden'),
            bar_style='info',
        )
        self.tracked_progress_label = widgets.HTML('')

        self.children = [
            header,
            instructions,
            self.existing_groups_label,
            selector_label,
            checkbox_items,
            widgets.VBox([group_name_label, self.group_name_input]),
            self.build_tracked_checkbox,
            widgets.HBox([self.create_btn, self.remove_btn, self.preview_btn, self.build_tracked_btn]),
            widgets.HBox([self.tracked_progress, self.tracked_progress_label]),
            self.status_output,
            self.preview_output,
        ]

    def _tracked_progress_cb(self, current, total, label):
        total = max(1, int(total))
        current = max(0, min(int(current), total))
        self.tracked_progress.layout.visibility = 'visible'
        if self.tracked_progress.max != total:
            self.tracked_progress.max = total
        self.tracked_progress.value = current
        self.tracked_progress_label.value = f'<span style="font-size:0.9em;">{label} ({current}/{total})</span>'

    def _run_tracked_build(self, output_dir, group_id):
        """Run create_group_tracked_segments with a live progress bar.

        Blocking (no background-thread infra in the notebook layer), but
        the IntProgress widget still updates live during the call.
        """
        try:
            written = create_group_tracked_segments(
                output_dir, group_id, self.metadata_loader.metadata,
                log=print,
                progress_cb=self._tracked_progress_cb,
            )
            return written, None
        except Exception as e:
            return None, e
        finally:
            self.tracked_progress.layout.visibility = 'hidden'
            self.tracked_progress_label.value = ''

    def _get_selected_cell_types(self):
        """Get list of selected member cell types from checkboxes."""
        return [ct for ct, cb in self.cell_type_checkboxes.items() if cb.value]

    def _get_selected_groups(self):
        return [gid for gid, cb in getattr(self, "group_checkboxes", {}).items() if cb.value]

    def _on_create_group(self, btn):
        """Callback for Create Group button."""
        try:
            group_name = self.group_name_input.value.strip()
            existing_groups = list_cell_type_groups(self._params())
            is_valid, error_msg = validate_group_name(group_name, self._all_cell_types, existing_groups)
            if not is_valid:
                self.status_output.value = f'<div style="color:red;"><b>Error:</b> {error_msg}</div>'
                return

            selected_types = self._get_selected_cell_types()
            if not selected_types:
                self.status_output.value = '<div style="color:red;"><b>Error:</b> Select at least one cell type.</div>'
                return

            output_dir = getattr(self.metadata_loader, "output_dir", None)
            if not output_dir:
                self.status_output.value = '<div style="color:red;"><b>Error:</b> No output directory set.</div>'
                return

            group_id = group_cell_type_id(group_name)
            try:
                out_csv = create_cell_type_group(output_dir, group_name, selected_types)
            except (FileNotFoundError, ValueError) as e:
                self.status_output.value = f'<div style="color:red;"><b>Error:</b> {e}</div>'
                return

            params = self._params()
            groups = params.setdefault("cell_type_groups", {})
            groups[group_id] = list(selected_types)
            saved_to = self._persist_params()

            success_html = (
                f'<div style="color:green;"><b>✓ Group created:</b> <code>{group_id}</code></div>'
                f'<div style="margin-top:8px; font-size:0.9em;">'
                f'<b>Members:</b> {", ".join(selected_types)}<br/>'
                f'<b>Merged CSV:</b> {out_csv}<br/>'
            )
            if saved_to:
                success_html += f'<b>Config saved to:</b> {saved_to}<br/>'

            print(f"✓ Group '{group_id}' created from: {selected_types}")
            print(f"  Merged CSV: {out_csv}")

            if self.build_tracked_checkbox.value:
                success_html += '<i>Building tracked segments… this may take a while.</i><br/>'
                self.status_output.value = success_html + '</div>'
                print("  Building tracked segments...")
                written, err = self._run_tracked_build(output_dir, group_id)
                if err is None:
                    success_html += (
                        f'<b style="color:green;">✓ Tracked segments built for '
                        f'{len(written)} sample(s).</b><br/>'
                    )
                    print(f"  ✓ Tracked segments built for {len(written)} sample(s).")
                else:
                    success_html += (
                        f'<b style="color:red;">⚠ Could not build tracked segments:</b> {err}<br/>'
                    )
                    print(f"  ⚠ Could not build tracked segments: {err}")

            success_html += '</div>'
            self.status_output.value = success_html

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

    def _on_remove_group(self, btn):
        """Callback for Remove Group button."""
        try:
            selected = self._get_selected_groups()
            if len(selected) != 1:
                self.status_output.value = '<div style="color:red;"><b>Error:</b> Select exactly one existing group to remove.</div>'
                return
            group_id = selected[0]
            params = self._params()
            remove_cell_type_group(params, group_id)
            saved_to = self._persist_params()

            summary = f'<div style="color:green;"><b>✓ Group removed:</b> <code>{group_id}</code></div>'
            if saved_to:
                summary += f'<div style="font-size:0.9em;">Config saved to: {saved_to}</div>'
            self.status_output.value = summary
            print(f"✓ Group '{group_id}' removed from config (merged CSV kept on disk).")

        except Exception as e:
            self.status_output.value = f'<div style="color:red;"><b>Error:</b> {str(e)}</div>'
            print(f"Error in GroupBuilder._on_remove_group: {e}")
            import traceback
            traceback.print_exc()

    def _on_build_tracked(self, btn):
        """Callback for Build Tracked Segments button (existing group)."""
        try:
            selected = self._get_selected_groups()
            if len(selected) != 1:
                self.status_output.value = '<div style="color:red;"><b>Error:</b> Select exactly one existing group to build tracked segments for.</div>'
                return
            group_id = selected[0]
            output_dir = getattr(self.metadata_loader, "output_dir", None)
            if not output_dir:
                self.status_output.value = '<div style="color:red;"><b>Error:</b> No output directory set.</div>'
                return

            self.status_output.value = (
                f'<div><i>Building tracked segments for <code>{group_id}</code>… '
                'this may take a while.</i></div>'
            )
            print(f"Building tracked segments for '{group_id}'...")
            written, err = self._run_tracked_build(output_dir, group_id)
            if err is None:
                status = group_tracked_segments_status(
                    output_dir, group_id, self.metadata_loader.metadata
                )
                self.status_output.value = (
                    f'<div style="color:green;"><b>✓ Tracked segments built for '
                    f'{group_id}:</b> {len(written)} sample(s) written '
                    f'({status["built"]}/{status["total"]} total now have tracked segments).</div>'
                )
                print(f"✓ Tracked segments built for '{group_id}': {len(written)} sample(s) written.")
            elif isinstance(err, (FileNotFoundError, ValueError)):
                self.status_output.value = f'<div style="color:red;"><b>Error:</b> {err}</div>'
            else:
                self.status_output.value = f'<div style="color:red;"><b>⚠ Could not build tracked segments:</b> {err}</div>'
                print(f"⚠ Could not build tracked segments: {err}")

        except Exception as e:
            self.status_output.value = f'<div style="color:red;"><b>Error:</b> {str(e)}</div>'
            print(f"Error in GroupBuilder._on_build_tracked: {e}")
            import traceback
            traceback.print_exc()

    def _on_preview_lengths(self, btn):
        """Callback for the track-length distribution preview button."""
        from behav3d.widgets.utils import build_track_length_distribution_figure
        from IPython.display import display

        self.preview_output.clear_output(wait=True)
        selected = self._get_selected_groups()
        with self.preview_output:
            if len(selected) != 1:
                print("Select exactly one existing group to preview.")
                return
            group_id = selected[0]
            output_dir = getattr(self.metadata_loader, "output_dir", None)
            if not output_dir:
                print("No output directory set.")
                return
            from behav3d.analysis.grouping import merged_track_features_csv
            csv_path = merged_track_features_csv(output_dir, group_id)
            fig = build_track_length_distribution_figure(
                csv_path, title=f"Track length distribution — {group_id}"
            )
            if fig is None:
                print(f"Could not build a track-length distribution for '{group_id}' ({csv_path}).")
                return
            display(fig)
