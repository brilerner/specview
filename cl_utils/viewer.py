# hv.extension('bokeh')
import holoviews as hv
from holoviews import opts, dim, streams
hv.extension('bokeh')
import numpy as np
import panel as pn
pn.extension()
import param
import xarray as xr

from scipy.ndimage import zoom

# local pkgs
from cl_utils.odemis_proc import get_wavelengths, get_img_arr, read_h5, get_sem_resolution

def make_rect_coords(xside=10, yside=10, xoffset=0, yoffset=0):
    x0 = 0 + xoffset
    x1 = x0 + xside
    y0 = 0 + yoffset
    y1 = y0 + yside
    return {('x','y'):[(x0,y0), (x1,y0), (x1,y1), (x0,y1)]}

def bare_hook(plot, element):
    plot.handles['xaxis'].visible = False
    plot.handles['yaxis'].visible = False


def conv_coords(data):
    data = zip(data['x0'], data['x1'], data['y0'], data['y1'])
    return [{('x','y'):[(x0,y0), (x1,y0), (x1,y1), (x0,y1)]} for x0,x1,y0,y1 in data]





class Spectrum(param.Parameterized):
    
    substrate = param.Array()
    xunits = param.String(default='nm')
    path = param.Filename()
    rect_coords = param.List()
    color_cycle = ['red', 'green', 'blue', 'orange', 'purple', 'cyan', \
                   'magenta', 'maroon', 'teal', 'navy', 'brown', 'pink']
    
    def __init__(self, **params):
        super().__init__(**params)
        self.h5 = read_h5(self.path)
        if self.h5:
            self.process_h5()
            self.set_con_dmap()
            self.set_cl_dmap()
            self.set_cl_agg_dmap()
            self.polys = hv.Polygons(self.rect_coords)
            # self.box_stream = streams.BoxEdit(source=self.polys,  styles={'fill_color': self.color_cycle})
        
    def process_h5(self):
        self.cl_arr = get_img_arr(self.h5['Acquisition2'])
        self.con_arr = zoom(get_img_arr(self.h5['Acquisition1']), 1/get_sem_resolution(self.h5))
#         self.con_arr = get_img_arr(self.h5['Acquisition1'])

#     @param.depends('rect_coords', watch=True) # doesn't work to link...
    def view(self):
        
        polys = self.polys
        polys = hv.Polygons(self.rect_coords)
        # box_stream = self.box_stream
        box_stream = streams.BoxEdit(source=polys,  styles={'fill_color': self.color_cycle})
        self.box_stream = box_stream
        hlines = hv.HoloMap({i: hv.VLine(i) for i in self.wl_vals}, self.x_coord)
        self.dmap_roi = hv.DynamicMap(self.roi_curves, streams=[box_stream])
        gridstyle = {'minor_xgrid_line_color': 'lightgray'}
        self.dmap_roi.opts(ylabel='Intensity', frame_height=215, show_grid=True, gridstyle=gridstyle)
#         self.dmap_roi.opts(ylim=(0,1000))
#         self.dmap_roi = self.dmap_roi.redim(CL=hv.Dimension('CL', soft_range=(400, 1000)))
        
        con_slider = pn.widgets.FloatSlider(start=0, end=1, value=1, name='Concurrent electrons alpha')
        con = self.dmap_con.apply.opts(alpha=con_slider.param.value)
        cl_agg_slider = pn.widgets.FloatSlider(start=0, end=1, value=0, name='CL aggregated alpha')
        cl_agg = self.dmap_cl_agg.apply.opts(alpha=cl_agg_slider.param.value)
        
        left = self.dmap_cl*con*cl_agg\
             *polys
        left.opts(xaxis='bare', yaxis='bare')
        
        right = self.dmap_roi*hlines
        right.opts(show_legend=False)
#         right.opts(legend_position='bottom_right')
        
        # set correct dimensions- pretty sure these are correct, but need to confirm
        y,x = self.con_arr.shape
        if y<x:
            diff = round((x-y)/2)
            left.opts(ylim=(0-diff, y+diff))
        else:
            diff = round((y-x)/2)
            left.opts(xlim=(0-diff, y+diff))
            
        dash_opts = [
            opts.Image(frame_width=350, frame_height=350, colorbar=True, colorbar_position='bottom',\
                       cmap='gray', hooks=[bare_hook]),
            opts.Curve(frame_width=350, frame_height=350, framewise=True), #width=400,
            opts.Polygons(fill_alpha=0.2, line_color='white'), 
            opts.Rectangles(fill_alpha=0.2, line_color='white'), 
            opts.VLine(color='black'),
            opts.Overlay(toolbar='right')
            ]
        
        plots = left+right
        plots.opts(dash_opts).opts(toolbar='left')
#         plots.opts(active_tools=['box_edit']) # doesn't work... nor does putting in plot-type specific opts
        sliders = pn.Column(con_slider, cl_agg_slider, width=300)
        dash = pn.Column(sliders, plots, background='white')

        return dash



    def roi_curves(self, data):
        if not data or not any(len(d) for d in data.values()):
            return hv.NdOverlay({0: hv.Curve([], self.x_coord, 'CL').opts(color='red')})
            # return hv.NdOverlay({'substrate': hv.Curve([], self.xlab, 'CL').opts(color='red')})

        
#         self.rect_coords = conv_coords(data)
        curves = {}
        self.data = data
        data = zip(data['x0'], data['x1'], data['y0'], data['y1'])
        for i, (x0, x1, y0, y1) in enumerate(data):
            pass
#             if i == 0:
#                 i = 'substrate'
            selection = self.ds_cl.select(x=(x0, x1), y=(y0, y1))
            curves[i] = hv.Curve(selection.aggregate(self.x_coord, np.mean)).opts(xlabel=self.xlab)
#         curves[0].
        self.curves = curves
        return hv.NdOverlay(curves).opts(opts.Curve(color= hv.Cycle(self.color_cycle)))
        # return hv.NdOverlay({0: hv.Curve([], self.xlab, 'CL').opts(color='red')})

### Beginnings of code to get physical dims (0 would be in middle of first pixel, but scale should be correct)
#     yscale = h5[key]['ImageData']['DimensionScaleY'][()]
#     xscale = h5[key]['ImageData']['DimensionScaleX'][()]
    # set [xscale*i for i in range()]
    # need to change ylim in self.view
    

    def set_con_dmap(self):
        con_y, con_x = self.con_arr.shape
        
        x_rng = range(con_x)
        y_rng = range(con_y)

# doesn't work...
#         xres, yres = get_sem_resolution(self.h5) # correct order?
#         x_rng = np.linspace(0, con_x/xres, con_x)
#         y_rng = np.linspace(0, con_y/yres, con_y)
        
        self.ds_con = hv.Dataset(xr.DataArray(self.con_arr[::-1,:],
                         coords={'x': x_rng, 'y': y_rng},
                         dims=('y','x'), name='Con'))
        self.dmap_con = self.ds_con.to(hv.Image)#.opts(width=500, height=500)
        if np.min(self.con_arr) == 0:
            rng = np.min([i for i in self.con_arr.flat if i>0]), np.max(self.con_arr)
            self.dmap_con.opts(clim = rng, clipping_colors={'min': 'transparent'})

    def set_cl_agg_dmap(self):
        arr = np.sum(self.cl_arr, axis=0)
        y, x = arr.shape
        self.ds_cl_agg = hv.Dataset(xr.DataArray(arr[::-1,:],
                         coords={'x': range(x), 'y': range(y)},
                         dims=('y','x'), name='CL Agg'))
        self.dmap_cl_agg = self.ds_cl_agg.to(hv.Image) # not really a dmap...


    
    def set_cl_dmap(self):    
        # wl_vals = f['Acquisition2']['ImageData']['DimensionScaleC'][()]*1e9 # nm
        self.wl_vals = get_wavelengths(self.h5)*1e9
        if self.xunits == 'nm':
            self.xlab = 'Wavelength (nm)'
            self.x_coord = 'w'
        elif self.xunits == 'eV':
            self.wl_vals = 1240/self.wl_vals
            self.xlab = 'eV'
            self.x_coord = 'eV'
        nwls, y, x = self.cl_arr.shape
        raw_xrds = xr.DataArray(
                        self.cl_arr[:,::-1,:],
                        coords={(self.x_coord):self.wl_vals, 'x': range(x), 'y': range(y)},
                        dims=(self.x_coord,'y','x'), name='CL'
                    )
        if not type(self.substrate) == type(None):
            raw_xrds = raw_xrds/self.substrate
#         raw_xrds[self.x_coord]['long_name'] = self.xlab # doesn't seem to hold
        self.ds_cl = hv.Dataset(raw_xrds)
        self.dmap_cl = self.ds_cl.to(hv.Image, ['x', 'y'], dynamic=True)
        
    def get_avg_substrate(self):
        df = self.dmap_roi.dframe()
        avg_substrate = df[df.Element==0].iloc[:,-1]
#         avg_substrate = df[df.Element=='substrate'].iloc[:,-1]
        self.avg_substrate =  avg_substrate.to_numpy().reshape(len(self.wl_vals),1,1) 
        return self.avg_substrate