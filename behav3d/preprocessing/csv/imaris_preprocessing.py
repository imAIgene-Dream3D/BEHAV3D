import pandas as pd
import argparse
from pathlib import Path

def run_imaris_preprocessing(
    df_positions_path, 
    df_organoid_distances_path,
    df_dead_dye_means_path, 
    output_path,
    df_tcell_distances_path=None,
    ):
    df_pos_imaris = pd.read_csv(df_positions_path, skiprows=3)
    df_pos_imaris=df_pos_imaris[["TrackID", "ID", "Time", "Position X", "Position Y", "Position Z"]]
    df_pos_imaris=df_pos_imaris.rename(columns={
            "Position X": "position_x",
            "Position Y": "position_y",
            "Position Z": "position_z",
            }
        )
    df_orgdist_imaris = pd.read_csv(df_organoid_distances_path, skiprows=3)
    df_orgdist_imaris=df_orgdist_imaris[["TrackID", "ID", "Time", "Shortest Distance to Surfaces"]]
    df_orgdist_imaris=df_orgdist_imaris.rename(columns={"Shortest Distance to Surfaces": "organoid_distance"})
    df_imaris = pd.merge(df_pos_imaris, df_orgdist_imaris)
    
    if df_tcell_distances_path:
        df_tcelldist_imaris = pd.read_csv(df_tcell_distances_path, skiprows=3)
        df_tcelldist_imaris=df_tcelldist_imaris[["TrackID", "ID", "Time", "Shortest Distance to Surfaces"]]
        df_tcelldist_imaris=df_tcelldist_imaris.rename(columns={"Shortest Distance to Surfaces": "tcell_distance"})
        df_imaris = pd.merge(df_imaris, df_tcelldist_imaris)
        
    df_dye_imaris = pd.read_csv(df_dead_dye_means_path, skiprows=3)
    df_dye_imaris=df_dye_imaris[["TrackID", "ID", "Time", "Intensity Mean"]]
    df_dye_imaris=df_dye_imaris.rename(columns={"Intensity Mean": "mean_dead_dye"})
    df_imaris = pd.merge(df_imaris, df_dye_imaris)

    df_imaris=df_imaris.rename(
        columns={
            "Time": 'position_t',
            "ID": "SegmentID"    
            }
    )
    df_imaris["position_t"]=df_imaris["position_t"]-1

    if output_path==None:
        output_path = df_positions_path.replace("_Position", "_tracks")

    df_imaris.to_csv(output_path, index=False)
    
def batch_imaris_preprocessing(
    folder,
    dead_dye_channel,
    organoid_surfaces_name,
    tcell_surfaces_name=None
    ):
    folder=Path(folder)
    assert folder.is_dir(), f"Given path is not a directory: {folder}"
    ims_samples = [d.resolve() for d in folder.iterdir() if d.is_dir()]
    
    def get_imaris_path(pattern, folder):
        folder = Path(folder)
        file_matches = list(folder.glob(pattern))
        assert len(file_matches)!=0, f"No {pattern} found for {folder.name}"
        assert len(file_matches)==1, f"Multiple {pattern} csv's found for {folder.name}"
        match = file_matches[0]
        return(match)
    
    for ims_sample in ims_samples:
        pos_patt = "*_Position.csv"
        pos_path = get_imaris_path(pos_patt, ims_sample)
        
        pos_patt = f"*_Intensity_Mean_Ch={dead_dye_channel}_Img=1.csv"
        deaddye_path = get_imaris_path(pos_patt, ims_sample)
        
        if not tcell_surfaces_name:
            tcell_dist_path=None
        else:
            pos_patt = f"*_Shortest_Distance_to_Surfaces_Surfaces={tcell_surfaces_name}.csv"
            tcell_dist_path = get_imaris_path(pos_patt, ims_sample)
             
        pos_patt = f"*_Shortest_Distance_to_Surfaces_Surfaces={organoid_surfaces_name}.csv"
        organoid_dist_path = get_imaris_path(pos_patt, ims_sample)
        
        outpath = f"{ims_sample}_behav3d_preprocessed.csv"
        run_imaris_preprocessing(
            df_positions_path=pos_path, 
            df_organoid_distances_path=organoid_dist_path,
            df_tcell_distances_path=tcell_dist_path,
            df_dead_dye_means_path=deaddye_path, 
            output_path=outpath
        )

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