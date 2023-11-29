import pandas as pd


def transform_csv_to_behav3d_input(
    original_csv,
    pos_x_column_name,
    pos_y_column_name,
    pos_z_column_name,
    pos_t_column_name,
    segment_id_column_name,
    track_id_column_name
    ):
    
    df_tracks = pd.read_csv(original_csv)
    rename_dict={
        pos_x_column_name:"position_x",
        pos_y_column_name:"position_y",
        pos_z_column_name:"position_z",
        pos_t_column_name:"position_t",
        segment_id_column_name:"SegmentID",
        track_id_column_name:"TrackID"
    }
    df_tracks=df_tracks.rename(columns=rename_dict)
    