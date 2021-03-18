import holoviews as hv
from holoviews import opts
import numpy as np
from cl_utils.odemis_proc import get_img_arr, get_wavelengths

def get_img_opts(shape):
        y,x = shape
        colorbar_opts = {'height':100}
        frame_width = 250
        frame_height = round(frame_width*(y/x))
        extra_border = round(frame_width*(1 - y/x))
        def hook(plot, element):
            plot.handles['plot'].min_border_top = extra_border//2
            plot.handles['plot'].min_border_bottom = extra_border//2

        return opts.Image(hooks=[hook], colorbar=True, frame_height=frame_height, frame_width=frame_width,\
                   colorbar_opts=colorbar_opts)

def match_aspect_opts(shape):
    y,x = shape
    if y < x:
        frame_width = 250
        frame_height = round(frame_width*(y/x))
    else:
        frame_height = 250
        frame_width = round(frame_height*(x/y))
    return opts.Image(frame_height=frame_height, frame_width=frame_width)

def get_bounds(h5, key):
    arr = get_img_arr(h5[key])
    if key == 'Acquisition2': # aggregate over CL
        arr = np.mean(arr, axis=0)
    yscale = h5[key]['ImageData']['DimensionScaleY'][()]
    xscale = h5[key]['ImageData']['DimensionScaleX'][()]
    ylen = arr.shape[0]*xscale
    xlen = arr.shape[1]*xscale
    bounds = -xlen/2, -ylen/2, xlen/2, ylen/2
    bounds = [i*1e6 for i in bounds]
    return bounds

def normspan(arr):
    arr -= np.min(arr)
    arr = arr/np.max(arr)
    return arr

def plot_scan(h5, key, normed=True):
    arr = get_img_arr(h5[key])
    title = h5[key]['PhysicalData']['Title'][()]
    if key == 'Acquisition2': # aggregate over CL
        arr = np.sum(arr, axis=0)
        title = 'CL (aggregated)'
    if normed:
        arr = normspan(arr)
    bounds = get_bounds(h5, key)
    img_opts = opts.Image(title=title, xlabel='x (um)', ylabel='y (um)', align='center')
    return hv.Image(arr, bounds=bounds).opts(match_aspect_opts(arr.shape)).opts(img_opts)

def plot_cl_full_dist(h5):
    key = 'Acquisition2'
    arr = get_img_arr(h5[key])
    sums = np.array([ np.sum(i) for i in arr])*1e-6
#     title = h5[key]['PhysicalData']['Title'][()]
    title = 'CL Spectrum (aggregated)'
    wavelengths = get_wavelengths(h5)*1e9
    gridstyle = {'minor_xgrid_line_color': 'lightgray'}
    curve_opts = opts.Curve(title=title, xlabel='Wavelength (nm)', ylabel='Counts (/million)', \
                            align='center', frame_height=100, show_grid=True, gridstyle=gridstyle)
    return hv.Curve((wavelengths, sums)).opts(curve_opts)