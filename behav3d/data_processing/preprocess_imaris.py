import pandas as pd
import argparse
from pathlib import Path
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
parser.add_argument('-p', '--positions', type=str, help='path to a imaris stats file containing the 3D position (*_dt_Position)', required=True)
parser.add_argument('-x', '--organoid_distances', type=str, help='path to a imaris stats file containing the distance to organoids', required=True)
parser.add_argument('-d', '--dead_dye_means', type=str, help='path to a imaris stats file containing mean dead dye intensity', required=True)
parser.add_argument('-o', '--output', type=str, help='path to the output csv in BEHAV3D format', required=False, default=None)
args = parser.parse_args()

imaris_pos_csv = args.positions
imaris_orgdist_csv = args.organoid_distances
imaris_dye_csv = args.dead_dye_means

output = args.output

imaris_pos_csv="/Users/samdeblank/Documents/1.projects/BHVD_BEHAV3D/BEHAV3D-ilastik/test/imaris_run/old_behav3d/AIM_MB2_Exp58_Img003_donor899_MA quantified_Statistics/AIM_MB2_Exp58_Img003_donor899_MA quantified_Position.csv"
imaris_orgdist_csv="/Users/samdeblank/Documents/1.projects/BHVD_BEHAV3D/BEHAV3D-ilastik/test/imaris_run/old_behav3d/AIM_MB2_Exp58_Img003_donor899_MA quantified_Statistics/AIM_MB2_Exp58_Img003_donor899_MA quantified_Intensity_Min_Ch=7_Img=1.csv"
imaris_dye_csv="/Users/samdeblank/Documents/1.projects/BHVD_BEHAV3D/BEHAV3D-ilastik/test/imaris_run/old_behav3d/AIM_MB2_Exp58_Img003_donor899_MA quantified_Statistics/AIM_MB2_Exp58_Img003_donor899_MA quantified_Intensity_Mean_Ch=4_Img=1.csv"


df_pos_imaris = pd.read_csv(imaris_pos_csv, skiprows=3)
df_pos_imaris=df_pos_imaris[["TrackID", "ID", "Time", "Position X", "Position Y", "Position Z"]]
df_orgdist_imaris = pd.read_csv(imaris_orgdist_csv, skiprows=3)
df_orgdist_imaris=df_orgdist_imaris[["TrackID", "ID", "Time", "Intensity Min"]]
df_dye_imaris = pd.read_csv(imaris_dye_csv, skiprows=3)
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

if output==None:
    output = imaris_pos_csv.replace("_Position", "_tracks")

df_imaris.to_csv(output, index=False)
