import pandas as pd
from pathlib import Path

def transform_csv_to_behav3d_input(
    original_csv,
    pos_x_column_name,
    pos_y_column_name,
    pos_z_column_name,
    pos_t_column_name,
    segment_id_column_name,
    track_id_column_name
    ):
    
    original_csv=Path(original_csv)
    df_tracks = pd.read_csv(original_csv)
    df_out_path = Path(original_csv.parent, str(original_csv.stem)+"_BEHAV3D_format.csv")
    # Dictionary that transforms supplied columns to BEHAV3D names
    rename_dict={
        track_id_column_name:"TrackID",
        segment_id_column_name:"SegmentID",
        pos_t_column_name:"position_t", 
        pos_z_column_name:"position_z",
        pos_y_column_name:"position_y",
        pos_x_column_name:"position_x",   
    }
    
    # Rename, select and order the columns based on 'rename_dict'
    df_tracks=df_tracks.rename(columns=rename_dict)
    df_tracks=df_tracks[list(rename_dict.values())]
    df_tracks=df_tracks[rename_dict.values()]
    
    # Set the position_t to always start from index 0
    if df_tracks["position_t"].min()==1:
        df_tracks["position_t"]-=1
        
    df_tracks.to_csv(df_out_path, sep=",", index=False)   
    
    
    