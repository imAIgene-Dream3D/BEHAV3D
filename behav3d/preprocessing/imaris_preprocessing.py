import pandas as pd
import argparse
from pathlib import Path

def run_imaris_preprocessing(
    df_positions_path, 
    df_organoid_distances_path, 
    df_dead_dye_means_path, 
    output_path
    ):
    df_pos_imaris = pd.read_csv(df_positions_path, skiprows=3)
    df_pos_imaris=df_pos_imaris[["TrackID", "ID", "Time", "Position X", "Position Y", "Position Z"]]
    df_orgdist_imaris = pd.read_csv(df_organoid_distances_path, skiprows=3)
    df_orgdist_imaris=df_orgdist_imaris[["TrackID", "ID", "Time", "Intensity Min"]]
    df_dye_imaris = pd.read_csv(df_dead_dye_means_path, skiprows=3)
    df_dye_imaris=df_dye_imaris[["TrackID", "ID", "Time", "Intensity Mean"]]

    df_imaris = pd.merge(df_pos_imaris, df_orgdist_imaris)
    df_imaris = pd.merge(df_imaris, df_dye_imaris)

    df_imaris=df_imaris.rename(
        columns={
            "Position X": "position_x",
            "Position Y": "position_y",
            "Position Z": "position_z",
            "Intensity Min": "organoid_distance",
            "Intensity Mean": "dead_dye_mean",
            "Time": 'position_t',
            "ID": "SegmentID"    
            }
    )
    df_imaris["position_t"]=df_imaris["position_t"]-1

    if output_path==None:
        output_path = df_positions_path.replace("_Position", "_tracks")

    df_imaris.to_csv(output_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
    parser.add_argument('-p', '--positions', type=str, help='path to a imaris stats file containing the 3D position (*_dt_Position)', required=True)
    parser.add_argument('-x', '--organoid_distances', type=str, help='path to a imaris stats file containing the distance to organoids', required=True)
    parser.add_argument('-d', '--dead_dye_means', type=str, help='path to a imaris stats file containing mean dead dye intensity', required=True)
    parser.add_argument('-o', '--output', type=str, help='path to the output csv in BEHAV3D format', required=False, default=None)
    args = parser.parse_args()
    run_imaris_preprocessing(
        args.positions, 
        args.organoid_distances, 
        args.dead_dye_means, 
        args.output
    )