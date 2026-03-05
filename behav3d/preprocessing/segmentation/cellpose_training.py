"""Cellpose training utilities.

This module provides some functions that train a Cellpose model for segmenting
organoids and Immune cells/Other cell types for BEHAV3D.

"""
import os
from cellpose import models, io, train, transforms
import numpy as np
import cv2

def reassign_labels(mask):
    labels_prev=np.unique(mask)
    new_mask=np.zeros_like(mask)
    for i in range(labels_prev.size):
        new_mask[mask==labels_prev[i]]=i
    return new_mask

def generate_2D_planes(image_dir, imf, maskf, dim_order, channels, anisotropy, crop_size=512, n_z_planes=10):
    """Generates the images for training Cellpose, 10 Z-planes per volume with default size 512x512 pixels.

    Parameters
    ----------
    image_dir:
        Root output directory for the images.
    imf:
        Filter for images. Images have the name NAME_OF_IMAGE+_imf.
    maskf:
        Filter for masks. Masks have the name NAME_OF_IMAGE+_maskf.
    dim_order:
        Dimension order for your images.
    channels:
        Channels for Cellpose.
    anisotropy:
        Pixel size in Z/Pixel size in XY. For making the volumes isotropic before training
    """
    assert (imf!='' or maskf!=''), f"Image or mask filter must be non empty"
    dir_2D=os.path.join(image_dir,'2D_cellpose'+maskf)
    os.makedirs(dir_2D, exist_ok=True)
    
    # find images
    image_names = io.get_image_files(image_dir, maskf, imf=imf,look_one_level_down=False)
    np.random.seed(0)
    pm = [(0, 1, 2, 3), (2, 0, 1, 3), (1, 0, 2, 3)]
    npm = ["YX", "ZY", "ZX"]
    for name in image_names:
        name0 = os.path.splitext(os.path.split(name)[-1])[0]
        img0 = io.imread(name)
        mask0 = io.imread(name.replace(imf,maskf))[None,...]
        # Make the images to have z-axis first and channel last
        if (dim_order=='ZYX'):
            img0 = img0.transpose(1,2,3,0)
        elif (dim_order=='YXZ'):
            img0 = img0.transpose(3,1,2,0)
        # Select the useful channels
        if (len(channels)>1):
            img0 = img0[...,channels]
        else:
            img0 = np.concatenate((img0[...,channels],np.zeros_like(img0[...,channels])),axis=-1)
        mask0 = mask0.transpose(1,2,3,0)
        for p in range(len(npm)):
            img = img0.transpose(pm[p]).copy()
            mask = mask0.transpose(pm[p]).copy()
            Ly, Lx = img.shape[1:3]
            id_z=np.where(np.max(mask[...,0],axis=(1,2))>0)[0]
            id_z=id_z[np.random.permutation(id_z.size)]
            imgs = img[id_z[:n_z_planes]]
            masks = mask[id_z[:n_z_planes]]
            if anisotropy > 1.0 and p > 0:
                imgs = transforms.resize_image(imgs, Ly=int(anisotropy * Ly), Lx=Lx)
                masks = transforms.resize_image(masks, Ly=int(anisotropy * Ly), Lx=Lx,interpolation=cv2.INTER_NEAREST)
            for k, img in enumerate(imgs):
                mask=masks[k,...]
                img=img.transpose(2,0,1)
                mask=reassign_labels(mask[:,:,0])
                mask=mask.astype('uint16')
                ly = 0 if Ly - crop_size <= 0 else np.random.randint(0, Ly - crop_size)
                lx = 0 if Lx - crop_size <= 0 else np.random.randint(0, Lx - crop_size)
                io.imsave(os.path.join(dir_2D, f'{name0}_{npm[p]}_{k}_img.tif'),
                        img[:,ly:ly + crop_size, lx:lx + crop_size].squeeze())
                io.imsave(os.path.join(dir_2D, f'{name0}_{npm[p]}_{k}_masks.tif'),
                        mask[ly:ly + crop_size, lx:lx + crop_size].squeeze())

def train_cellpose(train_dir, test_dir, nchan, n_epochs, model_type, model_name):
    """Train Cellpose.

    Parameters
    ----------
    train_dir:
        Root output directory for the train images.
    test_dir:
        Root output directory for the test images.
    channels:
        Channels to use for segmentation
    n_epochs:
        Number of epochs to train
    model_type:
        Model for Cellpose ('cyto3' for organoids, 'nuclei' for T cells).
    model_name:
        Name for saving the model
    """
    io.logger_setup()
    
    # 1. Get dataloaders
    images, labels, image_names, test_images, test_labels, image_names_test = io.load_train_test_data(train_dir, test_dir, image_filter='_img',
                                    mask_filter='_masks', look_one_level_down=False)
    
    # 2. Retrain the Cellpose generalist model
    if (nchan<3):
        model = models.CellposeModel(gpu=True, pretrained_model=False, model_type=model_type, nchan=2)
    else:
        model = models.CellposeModel(gpu=True, pretrained_model=False, nchan=nchan)
    
    model_path, train_losses, test_losses = train.train_seg(model.net,
                                train_data=images, train_labels=labels,
                                channels=None, normalize=True, 
                                test_data=test_images, test_labels=test_labels,
                                weight_decay=1e-4, SGD=True, learning_rate=0.1,
                                min_train_masks=1,
                                n_epochs=n_epochs, model_name=model_name)


