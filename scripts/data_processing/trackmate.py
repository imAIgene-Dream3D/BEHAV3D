# https://imagej.net/plugins/trackmate/scripting/trackmate-detectors-trackers-keys
import sys
import imagej
import scyjava as sj
import pandas as pd
import numpy as np

def run_trackmate(image_path):
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

    # image_path="/Users/samdeblank/surfdrive/Documents/1.projects/BHVD_BEHAV3D/BEHAV3D-ilastik/test/testdata_medium_finalobjects.tiff"
    print(image_path)
    imp = IJ.openImage(image_path)
    imp.show()
    
    #----------------------------
    # Create the model object now
    #----------------------------
    
    # Some of the parameters we configure below need to have
    # a reference to the model at creation. So we create an
    # empty model now.
    
    model = Model()
    # Send all messages to ImageJ log window.
    model.setLogger(Logger.IJ_LOGGER)
    
    #------------------------
    # Prepare settings object
    #------------------------
    
    settings = Settings(imp)
    
    # Configure detector - We use the Strings for the keys
    settings.detectorFactory = LabelImageDetectorFactory()
    settings.detectorSettings = settings.detectorFactory.getDefaultSettings()
    # settings.detectorSettings = {
    #     # 'DO_SUBPIXEL_LOCALIZATION' : True,
    #     # 'RADIUS' : 2.5,
    #     'TARGET_CHANNEL' : 1,
    #     'SIMPLIFY_CONTOURS' : False,
    #     # 'DO_MEDIAN_FILTERING' : False,
    # }  
    
    # Configure spot filters - Classical filter on quality
    filter1 = FeatureFilter('QUALITY', 30, True)
    settings.addSpotFilter(filter1)
    
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
    
    
    #-------------------
    # Instantiate plugin
    #-------------------
    
    trackmate = TrackMate(model, settings)
    
    #--------
    # Process
    #--------
    
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