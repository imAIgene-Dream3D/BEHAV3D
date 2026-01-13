import h5py

def load_h5(path, substruct):
    """
    Loading .h5 images
    """
    h5file = h5py.File(name=path, mode="r")
    img = h5file[substruct][:]
    return(img)
