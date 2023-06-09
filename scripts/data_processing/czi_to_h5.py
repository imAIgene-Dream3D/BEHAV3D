# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 09:32:35 2022
Script to convert czi to h5
@author: akrashen
"""


from aicsimageio import AICSImage
import h5py
input_file="E:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/AIM_ALLImmune/3.Analysis/behaved3.0setup/test data/20210803_WT1_n3(1).czi"
output_file="E:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/AIM_ALLImmune/3.Analysis/behaved3.0setup/test data/20210803_WT1_n3_1_test.h5"
img = AICSImage(input_file)

## Save file to H5
h5f = h5py.File(output_file, 'w')
h5f.create_dataset('dataset_1', data=img.data)
h5f.close()
