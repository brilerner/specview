import numpy as np

def read_h5(path):
    import h5py
    if path:
        h5 = h5py.File(path, 'r')
    else:
        import warnings
        warnings.warn('Must input filepath.')
        h5 = None
    return h5

def get_img_arr(h5):
    return h5['ImageData']['Image'][()].squeeze()

def get_extra_settings(h5):
    import json
    dct = h5['Acquisition1']['PhysicalData']['ExtraSettings'][()][0] # is string
    dct = json.loads(dct) # convert json dict to python dict
    return dct

def get_sem_resolution(h5):
    """
    Output shape --> (res_x, _resy)
    """
    dct = get_extra_settings(h5)
    return np.array(dct['SEM E-beam']['resolution'][0])

def get_wavelengths(h5):
    return h5['Acquisition2']['ImageData']['DimensionScaleC'][()]