import holoviews as hv
from holoviews import opts, dim, streams
import numpy as np
import panel as pn
pn.extension()
import param
import xarray as xr

from cl_utils.odemis_proc import get_wavelengths, get_img_arr, read_h5, get_sem_resolution
from hooks import y_bare_hook, no_border_hook


def hdf2xds(h5_member):
    
    arr = get_img_arr(h5_member)
    title = h5_member['PhysicalData']['Title'][()]
    x_rng, y_rng = get_ranges(h5_member)
    coords= [('y', y_rng), ('x', x_rng)]
    
    # add wl coord for spectrum
    if title == 'Spectrum':
        wl_vals = h5_member['ImageData']['DimensionScaleC'][()]*1e9
        coords.insert(0, ('w', wl_vals))
    
    return xr.DataArray(arr, coords, name=title)

def xds2img(xds):
    return hv.Image(xds, ['x', 'y'])

def get_percentiles(arr, lower_perc=2, upper_perc=98):
    
    lower_bound = np.percentile(arr, lower_perc)    
    upper_bound = np.percentile(arr, upper_perc)
    
    return lower_bound, upper_bound
    
def get_ranges(h5_member, mult_scale=1e6): # GENERALIZE INPUT?
    arr = get_img_arr(h5_member)
    ny, nx = arr.shape[-2:] # need to exlude wl dim for Acq 2
    
    img_data =  h5_member['ImageData']
    yscale = img_data['DimensionScaleY'][()]*mult_scale # units are m
    xscale = img_data['DimensionScaleX'][()]*mult_scale # units are m
    yoff = img_data['YOffset'][()]*mult_scale
    xoff = img_data['XOffset'][()]*mult_scale# units are m
    
    ylen = ny*yscale
    xlen = nx*xscale
    x0 = -xlen/2 + xoff
    x1 = xlen/2 + xoff
    y0 = -ylen/2 + yoff
    y1 = ylen/2 + yoff
    
    x_rng = np.arange(x0+xscale/2, x1, xscale)
    y_rng = np.arange(y0+yscale/2, y1, yscale)[::-1] # for proper xarray coords
    
    return x_rng, y_rng

def get_extents(xds): # MAKE INPUT XRNG, YRNG
    x0, x1 = min(xds.x), max(xds.x)
    y0, y1 = min(xds.y), max(xds.y)
    dx = np.abs(x1-x0)
    dy = np.abs(y1-y0)
    cushion = (dx-dy)/2
    if dx>dy:
        y0 -= cushion
        y1 += cushion
    else:
        x0 -= cushion
        x1 += cushion

    return tuple(map(float, [x0, y0, x1, y1]))

def transp_hook(plot, element):
    plot.handles['plot'].border_fill_alpha = 0

css = '''
.bk.panel-widget-box {
  background: #f0f0f0;
  border-radius: 5px;
  border: 1px black solid;
}
'''

pn.extension(raw_css=[css])

class Spectrum_Base(param.Parameterized):
    
    band = param.Range(default=(500, 501), label='Band (nm)')#, step=20)# can use floats too
    center  = param.Number(default=500, step=1, doc="center of band") # Parameter works too, but can't adjust value in box
    bw  = param.Number(default=100, step=1, doc="center of band", label='BW')
    outlier_percentiles = param.Range(default=(2,98), bounds=(0,100), doc="outlier percentiles", label='Threshold')
    cl_alpha = param.Number(default=0.5, bounds=(0,1), label='Alpha')
    rect_coords = param.List()
    roi_toggle = param.Selector(objects=["CL", "Trans"])
    color_cycle = param.List(['red', 'green', 'blue', 'orange', 'purple', 'cyan', \
                   'magenta', 'maroon', 'teal', 'navy', 'brown', 'pink'])




class Spectrum(Spectrum_Base):
    
    

    def __init__(self, h5, **params):
        # super().__init__(**params)
        
        self.process_h5(h5)

        self.polys =  hv.Polygons(self.rect_coords)
        self.box_stream = streams.BoxEdit(source=self.polys,  styles={'fill_color': self.color_cycle})
        

    ### INITIAL PROCESSING ###

    def process_h5(self, h5): # SHOULD I TURN DATASETS INTO A CLASS
        # set datasets
        self.survey_xds = hdf2xds(h5['Acquisition0'])
        self.concurrent_xds = hdf2xds(h5['Acquisition1'])
        self.cl_xds = hdf2xds(h5['Acquisition2'])

        # set values
        self.agg_spec = self.cl_xds.sum(('x','y'))
        self.raw_wls = self.cl_xds.w.values # MAKE INTO PARAMETER?

        # set parameters
        self.set_wl_bounds()
        # self.set_band_slice()

        # set images
        self.survey_img = xds2img(self.survey_xds)
        x0, y0, x1, y1 = get_extents(self.survey_xds)
        self.survey_img = self.survey_img.redim(x=dict(range=(x0,x1)), y=dict(range=(y0,y1)))
        self.concurrent_img = xds2img(self.concurrent_xds)

        
    def set_wl_bounds(self):
        min_wl, max_wl = min(self.raw_wls), max(self.raw_wls)
        center = (max_wl+min_wl)/2
        min_wl, max_wl, center = round(min_wl), round(max_wl), round(center)
        self.param.band.bounds = (min_wl, max_wl)
        self.center = center

    
    ### DYNAMIC SETTING ###

    @param.depends('band', watch=True)
    def set_text_from_band(self):
        center = np.average(self.band)
        bw = np.diff((self.band))[0]
        center, bw = np.round(center,1), np.round(bw,1)
#         print(center, bw)
        self.center = center
#         print(self.center)
        self.bw = bw
        
    @param.depends('center', 'bw', watch=True)
    def set_band_from_text(self):
#         print(self.center)
#         print('set text')
        # self.ls.append('set text')
        lower = self.center - self.bw/2
        upper = self.center + self.bw/2
        
        lower_bound, upper_bound = self.param.band.bounds
        if (lower < lower_bound) or (upper >= upper_bound):
            self.set_literal_input()
        else:
            self.band = lower, upper
    
    def plot_cl_img(self, band, outlier_percentiles):

        # set band slice
        self.band_slice = self.cl_xds.sel(w=slice(*band)).mean('w')

        # plot img
        img = xds2img(self.band_slice)
        lower_bound, upper_bound = get_percentiles(self.band_slice, *outlier_percentiles)
        return img.redim.range(Spectrum=(lower_bound, upper_bound))

        
    def band_selection(self, band):
        agg_spec_curve = hv.Curve(self.agg_spec).opts(hooks=[y_bare_hook])#, no_border_hook])
        band_vspan = hv.VSpan(*band).opts(color='gray', alpha=0.5)
        band_opts = opts.Overlay(toolbar=None, hooks=[transp_hook], xlabel='',\
                                    frame_width=150, frame_height=50)
        return (agg_spec_curve*band_vspan).opts(band_opts)     


    ### ROI ####

    def roi_curves(self, data): ### CHANGE Spectrum?
        
        if not data or not any(len(d) for d in data.values()):
            default_curve = hv.Curve([], 'Spectrum', 'CL').opts(color='red')
            return hv.NdOverlay({0: default_curve}).opts(show_legend=False) # code breaks without using a curve in overlay
        
        curves = {}
        data = zip(data['x0'], data['x1'], data['y0'], data['y1'])
        for i, (x0, x1, y0, y1) in enumerate(data):
            selection = self.cl_xds.sel(x=slice(x0, x1), y=slice(y1, y0))
            self.selection = selection
            selection_avg = selection.mean(['x','y'])
            if self.roi_toggle == 'Trans':
                if i == 0:
                    substrate = selection_avg.copy()
                selection_avg /= substrate
            curves[i] = hv.Curve(selection_avg)
        
        return hv.NdOverlay(curves).opts(opts.Curve(color= hv.Cycle(self.color_cycle)))

    def get_roi(self):

        box_stream = self.box_stream
        
        roi = hv.DynamicMap(self.roi_curves, streams=[box_stream])
        roi.opts(xlim=(self.raw_wls[0], self.raw_wls[-1]))
#         self.dmap_roi.opts(ylim=(0,1000))
#         self.dmap_roi = self.dmap_roi.redim(CL=hv.Dimension('CL', soft_range=(400, 1000)))
        line = hv.VLine(self.center)
        roi = roi*line
        roi_opts = [opts.Curve(frame_width=275, frame_height=125, framewise=True, \
                    xlabel='Wavelength (nm)', ylabel='Intensity', \
                    show_grid=True, gridstyle = {'minor_xgrid_line_color': 'lightgray'}, toolbar='above'), #width=400,
                    opts.VLine(color='black', line_dash='dashed', line_width=1),]
                    # opts.NdOverlay(toolbar=None, show_legend=False,)]
        roi.opts(roi_opts)
        return roi


    def view(self):
        
 
        

        scan_opts = [opts.Image(frame_width=350, frame_height=350, colorbar=True, colorbar_position='bottom',\
                cmap='gray', xlabel='units in µm', ylabel=''), #hooks=[bare_hook]
            opts.Polygons(fill_alpha=0.2, line_color='white'), 
            opts.Rectangles(fill_alpha=0.2, line_color='white'),
            opts.Overlay(toolbar='above')]

        
        # cl

        streams = dict(band=self.param.band, outlier_percentiles=self.param.outlier_percentiles)
        self.cl_img = hv.DynamicMap(self.plot_cl_img, streams=streams)
        self.cl_img = self.cl_img.apply.opts(alpha=self.param.cl_alpha)
        # img_ov = self.cl_img*self.polys
        img_ov = self.survey_img*self.concurrent_img*(self.cl_img*self.polys)
        img_ov.opts(scan_opts)


        # slider_row = pn.Row('Con', self.alpha_slider, 'CL', width=100)

        # roi
        roi_toggle_widg = pn.Param(self.param.roi_toggle, 
                            widgets={'roi_toggle': {'type':pn.widgets.RadioButtonGroup, 'width':150}})
        roi_col = pn.Column(
                    pn.Column(roi_toggle_widg, align='center', width_policy='max'), # how to align??
                    pn.Row(self.get_roi),
                     css_classes=['panel-widget-box'])


        w = 175


        # band
        self.band_selec = hv.DynamicMap(self.band_selection, streams=dict(band=self.param.band))
        band_indicators = pn.Row(self.param.center, self.param.bw, width=w)
        band_widg =  pn.Param(self.param.band, widgets={'band':{'width':150}} )
        band_lay = pn.Column(band_indicators, pn.Row(band_widg, margin = (0, 0, -25, 0)), \
                        pn.pane.HoloViews(self.band_selec, linked_axes=False), 
                        css_classes=['panel-widget-box'], width=w)


        sidebar = pn.Column(
                    pn.Column('## CL Opts',  background='#44e3a1', width_policy='max'), 
                    pn.Column(self.param.cl_alpha, self.param.outlier_percentiles, css_classes=['panel-widget-box'], width=w),
                    pn.Column(band_lay, css_classes=['panel-widget-box']) ,
                    background='white', width=w)


        
        return pn.Row(
                sidebar,
                img_ov,
                pn.Column(pn.Spacer(width=200), roi_col),
                background='white')


    
#     def pixel_dist(self):
#         hv.Histogram(np.histogram(self.band_slice, bins=100))

#         plots.opts(active_tools=['box_edit']) # doesn't work... nor does putting in plot-type specific opts
        # widg =  pn.Param( self.param.center, widgets={'center':{'type':pn.widgets.NumberInput, 'width':300}} )

  

    # def test(self, band, center):
    #     band_slice = self.cl_xds.sel(w=slice(*band)).mean('w')
    #     img = xds2img(band_slice)
    #     return img
    
    # def tview(self):
    #     streams = dict(band=self.param.band, center=self.param.center)
    #     dm = hv.DynamicMap(self.test, streams=streams)
    #     return pn.Column(self.param.band, dm)
    

        

     
 