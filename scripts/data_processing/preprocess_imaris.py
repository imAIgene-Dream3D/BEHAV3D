import pandas as pd
import argparse
from pathlib import Path
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
parser.add_argument('-i', '--input', type=str, help='path to a imaris stats file containing the 3D position (*_dt_Position)', required=True)
parser.add_argument('-o', '--output', type=str, help='path to the output csv in BEHAV3D format', required=False, default=None)
args = parser.parse_args()

imaris_csv = args.input
output = args.output

imaris_csv="/Users/samdeblank/Documents/1.projects/BHVD_BEHAV3D/BEHAV3D-ilastik/test/imaris_run/2021-06-10_WT1_n1(1)_10T_[ims1_2021-09-09T14-14-13.594]_dt_Position.csv"
df_imaris = pd.read_csv(imaris_csv, skiprows=3)
df_imaris=df_imaris.rename(
    columns={
        "Position X": "position_x",
        "Position Y": "position_y",
        "Position Z": "position_z",
        "Time": 'position_t',
        "ID": "SegmentID"    
        }
)
df_imaris["position_t"]=df_imaris["position_t"]-1
df_imaris = df_imaris[['TrackID', 'SegmentID', 'position_t','position_z', 'position_y', 'position_x']]

if output==None:
    output = imaris_csv.replace("_dt_Position", "_tracks")

df_imaris.to_csv(output, index=False)
