# https://imagej.net/plugins/trackmate/scripting/trackmate-detectors-trackers-keys
import sys
import imagej
import scyjava as sj
import pandas as pd
import numpy as np
import psutil
from pathlib import Path
import yaml
from tifffile import imread, imwrite
from skimage.measure import regionprops_table

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
args = parser.parse_args()

def main(config):
    with open("/Users/samdeblank/Library/CloudStorage/OneDrive-PrinsesMaximaCentrum/github/BEHAV3D-ilastik/scripts/data_processing/config.yml", "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    
    sample_name = config['sample_name']
    output_dir = config['output_dir']

    element_size_x=config['element_size_x']
    element_size_y=config['element_size_y'] 
    element_size_z=config['element_size_z']
    element_size_unit=config['element_size_unit']
    
    ### Track the data using TrackMate
    tcell_segments_path = Path(output_dir, f"{sample_name}_tcell_segments.tiff")   
    df_tracks=run_trackmate(
        image_path=str(tcell_segments_path),
        element_size_x=element_size_x,
        element_size_y=element_size_y,
        element_size_z=element_size_z,
        element_size_unit=element_size_unit
        )
    # Add 1 to every track_id so 0 is not a track in the image (should be background)
    df_tracks["track_id"]=df_tracks["track_id"]+1
    df_tracks = df_tracks.sort_values(by=["track_id","position_t"])

    tracks_out_path = Path(output_dir, f"{sample_name}_tracks.csv")
    df_tracks.to_csv(tracks_out_path, sep=",", index=False)
    ### Assign the tracks to existing segments
    # Loop through spots, link to segments in the image and replace label with track_id
    tcell_segments = imread(tcell_segments_path)
    
    df_centroids = []
    for t, tcell_stack in enumerate(tcell_segments):
        properties=pd.DataFrame(regionprops_table(label_image=tcell_stack, properties=['label', f'centroid']))
        properties["position_t"]=t
        df_centroids.append(properties)
    df_centroids = pd.concat(df_centroids)
    df_centroids["position_z"]=df_centroids["centroid-0"]*element_size_z
    df_centroids["position_y"]=df_centroids["centroid-1"]*element_size_y
    df_centroids["position_x"]=df_centroids["centroid-2"]*element_size_x
    
    tcell_tracked_out_path= Path(output_dir, f"{sample_name}_tcells_tracked.tiff")
    tcells_tracked = np.zeros_like(tcell_segments)
    for _, row in df_tracks.iterrows():
        t,z,y,x = int(row["position_t"]),row["position_z"], row["position_y"], row["position_x"]
        corr_seg = df_centroids[
            (df_centroids["position_x"]==x) &
            (df_centroids["position_y"]==y) &
            (df_centroids["position_z"]==z) &
            (df_centroids["position_t"]==t) 
        ]["label"].iloc[0]
        assert corr_seg!=0, f"Position of center segment corresponds to background (0), which is an error"
        tcells_tracked[t,:,:,:][tcell_segments[t,:,:,:]==corr_seg]=row["track_id"]
        # im_track = im_track[im==corr_seg]=row["track_id"]
    imwrite(
        tcell_tracked_out_path,
        tcells_tracked
    ) 

def run_trackmate(
    image_path,
    element_size_x, 
    element_size_y, 
    element_size_z,
    element_size_unit
    ):
    available_75perc_memory=int(psutil.virtual_memory().available*0.75/(1024**3))
    sj.config.add_options(f'-Xmx{available_75perc_memory}g')
    ij = imagej.init('sc.fiji:fiji', mode='headless')
    IJ = ij.IJ
    WindowManager=ij.WindowManager
    # from ij import IJ
    # from ij import WindowManager
    
    Model = sj.jimport('fiji.plugin.trackmate.Model')
    Settings = sj.jimport('fiji.plugin.trackmate.Settings')
    TrackMate = sj.jimport('fiji.plugin.trackmate.TrackMate')
    SelectionModel = sj.jimport('fiji.plugin.trackmate.SelectionModel')
    SparseLAPTrackerFactory = sj.jimport('fiji.plugin.trackmate.tracking.jaqaman.SparseLAPTrackerFactory')
    LabelImageDetectorFactory = sj.jimport('fiji.plugin.trackmate.detection.LabelImageDetectorFactory')
    FeatureFilter = sj.jimport('fiji.plugin.trackmate.features.FeatureFilter')
    Logger = sj.jimport('fiji.plugin.trackmate.Logger')

    print(image_path)
    imp = IJ.openImage(image_path)
    IJ.run(imp, "Properties...", f"unit={element_size_unit} pixel_width={element_size_x} pixel_height={element_size_y} voxel_depth={element_size_z}")
    # imp.show()
    
    #----------------------------
    # Create the model object now
    #----------------------------
    
    # Some of the parameters we configure below need to have
    # a reference to the model at creation. So we create an
    # empty model now.
    
    model = Model()
    model.setLogger(Logger.IJ_LOGGER)
    
    settings = Settings(imp)
    
    # Configure detector - We use the Strings for the keys
    settings.detectorFactory = LabelImageDetectorFactory()
    settings.detectorSettings = settings.detectorFactory.getDefaultSettings()
    
    # Configure tracker - We want to allow merges and fusions
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = settings.trackerFactory.getDefaultSettings() # almost good enough
    settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = True
    settings.trackerSettings['ALLOW_TRACK_MERGING'] = True
    
    # Add ALL the feature analyzers known to TrackMate. They will 
    # yield numerical features for the results, such as speed, mean intensity etc.
    settings.addAllAnalyzers()
    
    # Configure track filters - We want to get rid of the two immobile spots at
    # the bottom right of the image. Track displacement must be above 10 pixels.
    
    # filter2 = FeatureFilter('TRACK_DISPLACEMENT', 10, True)
    # settings.addTrackFilter(filter2)
    
    trackmate = TrackMate(model, settings)
    ok = trackmate.checkInput()
    if not ok:
        sys.exit(str(trackmate.getErrorMessage()))
    
    ok = trackmate.process()
    if not ok:
        sys.exit(str(trackmate.getErrorMessage()))

    # The feature model, that stores edge and track features.
    fm = model.getFeatureModel()

    keys_df_spots = [
        "spot_id",
        "track_id",
        "position_t",
        "position_z",
        "position_y",
        "position_x",
    ]
    df_spots = pd.DataFrame(columns=keys_df_spots)
    # Iterate over all the tracks that are visible.
    for trackid in model.getTrackModel().trackIDs(True):
        track = model.getTrackModel().trackSpots(trackid)
        for spot in track:
            sid = spot.ID()
            # q=spot.getFeature('QUALITY')
            spot_info = {
                "spot_id":spot.ID(),
                "track_id":trackid,
                "position_t":spot.getFeature("POSITION_T"),
                "position_z":spot.getFeature("POSITION_Z"),
                "position_y":spot.getFeature("POSITION_Y"),
                "position_x":spot.getFeature("POSITION_X"),
            }
            new_row = pd.DataFrame(spot_info, index=[0])
            df_spots = pd.concat([df_spots, new_row], ignore_index=True)
            
    return(df_spots)

if __name__ == "__main__":
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    main(config)