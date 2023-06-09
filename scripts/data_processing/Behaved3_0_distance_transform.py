# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:30:42 2022
This script takes H5 file with segmented organoids and perfroms distance transformation generating another channel with image transform
@author: imaging
"""


import napari
import numpy as np
import h5py
from scipy.ndimage.morphology import distance_transform_edt
import matplotlib.pyplot as plt
# Import the segmented organoids
filename="E:/malieva/BEHAV3D3.0/20210803_WT1_n3_1_test-dataset_1_org_Object-Identities.h5"
f = h5py.File(filename, 'r')

# Import the original image
filename2="E:/malieva/BEHAV3D3.0/test_data/20210803_WT1_n3_1_test.h5"
f2 = h5py.File(filename2, 'r')


#Get the HDF5 keys
list(f.keys())
data = f['exported_data'][:]
# To view all the channels, need to loop over and add to napari


viewer=napari.Viewer() #optional

mask_org= np.ma.masked_where(data!=0, data)

## Perform distance transform of reverse mask, since distance transform only works on inside the objects
dist_org=distance_transform_edt(mask_org.mask)

plt.imshow(dist_org[0,0,20])
plt.savefig('result.png', bbox_inches='tight')

# Stack processed channels together
merged=np.stack((f2['dataset_1'][:,0], f2['dataset_1'][:,1],f2['dataset_1'][:,2],dist_org[:,0]), axis=1)## Sam can you fix this to have flexibility for different number of channels

## Save file to H5
h5f = h5py.File('E:/malieva/BEHAV3D3.0/test_data/20210803_WT1_n3_1_test_dist_transform.h5', 'w')
h5f.create_dataset('dataset_2', data=merged)
h5f.close()
