import holoviews as hv
from holoviews import opts, dim, streams
import numpy as np
import panel as pn
pn.extension()
import param
import xarray as xr

from cl_utils.odemis_proc import get_wavelengths, get_img_arr, read_h5, get_sem_resolution
from cl_utils.hooks import y_bare_hook, no_border_hook


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
# change css for accordion?
# .widget-accordion .bk-btn-group button {
#   font-size: 45pt;
# }
pn.extension(raw_css=[css])

class BasePanViewer(param.Parameterized):
    threshold = param.Range(default=(2,98), bounds=(0,100), doc="outlier percentiles", label='Threshold')
    alpha = param.Number(default=0.5, bounds=(0,1), label='Alpha')

    def __init__(self, acq, **params):
        """
        acq : acquisition in form of h5 dictionary
        """
        super().__init__(**params)
        self.xds = hdf2xds(acq)
        self.img_data = self.xds
        self.title = acq['PhysicalData']['Title'][()]
        self.acq = acq
        self.streams = dict(threshold=self.param.threshold)
        self.extents  = get_extents(self.xds)

    def proc_img(self, img_data, threshold):
        img = xds2img(img_data).opts(colorbar=True)
        lower_bound, upper_bound = get_percentiles(img_data, *threshold)
        img = img.redim.range(**{self.title : (lower_bound, upper_bound)})

        # set extents
        x0, y0, x1, y1 = self.extents
        img = img.redim(x=dict(range=(x0,x1)), y=dict(range=(y0,y1)))

        # set opts
        img_opts = opts.Image(frame_width=350, frame_height=350, colorbar=True, colorbar_position='bottom',\
                cmap='gray', xlabel='units in µm', ylabel='')
        img.opts(img_opts)
        
        return img

    def plot_img(self, threshold):
        return self.proc_img(self.img_data, threshold)

    def get_dmap(self):
        dmap = hv.DynamicMap(self.plot_img, streams=self.streams)
        return dmap.apply.opts(alpha=self.param.alpha)

    def view(self):
        return pn.Row(self.param, self.get_dmap(), background='white')

class BaseSpecViewer(BasePanViewer):
    
    band = param.Range(default=(500, 501), bounds=(200,800), label='Band (nm)')#, step=20)# can use floats too
    center  = param.Number(default=500, step=1, doc="center of band") # Parameter works too, but can't adjust value in box
    bw  = param.Number(default=100, step=1, doc="center of band", label='BW', precedence=-1)

    def __init__(self, acq, **params):
        """
        acq : acquisition in form of h5 dictionary
        """
        super().__init__(acq, **params)
        self.streams['band'] = self.param.band
        self.set_wl_bounds()

        # self.agg_spec = self.Spectrum.xds.sum(('x','y'))



    def get_dmap(self, band_override=None):
        def plot_img(threshold, band):
            if band_override:
                band = band_override
            img_data = self.xds.sel(w=slice(*band)).mean('w')
            return self.proc_img(img_data, threshold)
        dmap = hv.DynamicMap(plot_img, streams=self.streams)
        return dmap.apply.opts(alpha=self.param.alpha)
    
    def set_wl_bounds(self):
        wls = self.xds.w.values
        min_wl, max_wl = min(wls), max(wls)
        center = (max_wl+min_wl)/2
        min_wl, max_wl, center = round(min_wl), round(max_wl), round(center)
        self.param.band.bounds = (min_wl, max_wl)
        self.band = center - self.bw/2, center + self.bw/2
        self.center = center

#     @param.depends('band', watch=True)
#     def set_text_from_band(self):
#         center = np.average(self.band)
#         bw = np.diff((self.band))[0]
#         center, bw = np.round(center,1), np.round(bw,1)
# #         print(center, bw)
#         self.center = center
# #         print(self.center)
#         self.bw = bw
        
#     @param.depends('center', 'bw', watch=True)
#     def set_band_from_text(self):
# #         print(self.center)
# #         print('set text')
#         # self.ls.append('set text')
#         lower = self.center - self.bw/2
#         upper = self.center + self.bw/2
        
#         lower_bound, upper_bound = self.param.band.bounds
#         if (lower < lower_bound) or (upper >= upper_bound):
#             self.set_literal_input()
#         else:
#             self.band = lower, upper

    # def band_selection(self, band):
    #     agg_spec_curve = hv.Curve(self.agg_spec).opts(hooks=[y_bare_hook])#, no_border_hook])
    #     band_vspan = hv.VSpan(*band).opts(color='gray', alpha=0.5)
    #     band_opts = opts.Overlay(toolbar=None, hooks=[transp_hook], xlabel='',\
    #                                 frame_width=150, frame_height=50)
    #     return (agg_spec_curve*band_vspan).opts(band_opts)  
    # 
    # for view method:
        # self.band_selec = hv.DynamicMap(self.band_selection, streams=dict(band=self.param.band))
        # band_indicators = pn.Row(self.param.center, self.param.bw, width=w)
        # band_widg =  pn.Param(self.param.band, widgets={'band':{'width':150}} )
        # band_lay = pn.Column(band_indicators, pn.Row(band_widg, margin = (0, 0, -25, 0)), \
        #                 pn.pane.HoloViews(self.band_selec, linked_axes=False), 
        #                 css_classes=['panel-widget-box'], width=w)
  


class BaseSpecViewerROI(BaseSpecViewer):
    roi_toggle = param.Selector(objects=["CL", "Trans"])
    color_cycle = param.List(['red', 'green', 'blue', 'orange', 'purple', 'cyan', \
                   'magenta', 'maroon', 'teal', 'navy', 'brown', 'pink'])
    rect_coords = param.List()

    
    def __init__(self, acq, **params):
        """
        acq : acquisition in form of h5 dictionary
        """
        super().__init__(acq, **params)
        self.polys =  hv.Polygons(self.rect_coords)
        self.box_stream = streams.BoxEdit(source=self.polys,  styles={'fill_color': self.color_cycle})

    def roi_curves(self, data): ### CHANGE Spectrum?
            
        if not data or not any(len(d) for d in data.values()):
            default_curve = hv.Curve([], 'Spectrum', 'CL').opts(color='red') # How can I change 'Spectrum', 'CL'
            return hv.NdOverlay({0: default_curve}).opts(show_legend=False) # code breaks without using a curve in overlay
        
        curves = {}
        data = zip(data['x0'], data['x1'], data['y0'], data['y1'])
        for i, (x0, x1, y0, y1) in enumerate(data):
            selection = self.xds.sel(x=slice(x0, x1), y=slice(y1, y0))
            selection_avg = selection.mean(['x','y'])
            if self.roi_toggle == 'Trans': # apparently param knows when this changes without having to make it a 'stream' var
                if i == 0:
                    substrate = selection_avg.copy()
                selection_avg /= substrate
            curves[i] = hv.Curve(selection_avg)
        
        color_cycle_opts = opts.Curve(color= hv.Cycle(self.color_cycle))
        return hv.NdOverlay(curves).opts(color_cycle_opts)

    def get_polys(self): # needs to be separate?
        poly_opts = opts.Polygons(fill_alpha=0.2, line_color='white')
        return self.polys.opts(poly_opts)

    def get_roi(self):

        roi_opts = opts.Curve(frame_width=275, frame_height=125, framewise=True,
                                xlabel='Wavelength (nm)', ylabel='Intensity',
                                show_grid=True, gridstyle = {'minor_xgrid_line_color': 'lightgray'}, 
                                toolbar='above')
        roi = hv.DynamicMap(self.roi_curves, streams=[self.box_stream]).opts(roi_opts)
        roi.opts(xlim=(self.param.band.bounds)) # CHANGE raw_wls

        line_opts = opts.VLine(color='black', line_dash='dashed', line_width=1)
        line = hv.VLine(self.center).opts(line_opts)

        return roi*line


    def get_roi_panel(self):

        roi_toggle = pn.Param(self.param.roi_toggle, 
                            widgets={'roi_toggle': {'type':pn.widgets.RadioButtonGroup, 'width':150}})
        return pn.Column(
                    pn.Column(roi_toggle, align='center', width_policy='max'), # how to align??
                    pn.Row(self.get_roi),
                     css_classes=['panel-widget-box'])

    def view(self):
        
        # img panel
        img_params = pn.Column( self.param.alpha, self.param.threshold, self.param.band)
        img_panel = pn.Row(
                        img_params,
                        self.get_dmap()*self.get_polys())

        # roi panel


        return pn.Row(
                    img_panel,
                    self.get_roi_panel(),
                    background='white')

    

class Spectrum(param.Parameterized):

    def __init__(self, h5, **params):
        super().__init__(**params)
        self.process_h5(h5)

    ### INITIAL PROCESSING ###

    def process_h5(self, h5): # SHOULD I TURN DATASETS INTO A CLASS
        self.Survey = BasePanViewer(h5['Acquisition0'])
        self.Concurrent = BasePanViewer(h5['Acquisition1'])
        self.Spectrum = BaseSpecViewerROI(h5['Acquisition2'])

    def view(self):
        
        img_ov = hv.Overlay([
                            self.Survey.get_dmap(),
                            self.Concurrent.get_dmap(),
                            self.Spectrum.get_dmap()*self.Spectrum.get_polys()]) # how to call polys?
        img_ov.opts(toolbar='above')
        img_ov = img_ov.collate()

        sidebar_tuples = []
        for class_ in [self.Spectrum, self.Concurrent, self.Survey]:

            use_params = [class_.param.alpha, class_.param.threshold] # use super?
            if class_.title == 'Spectrum':
                use_params.append(class_.param.band)
                class_.alpha = 0.5
            else: class_.alpha = 1
            
            sidebar_tuples.append((class_.title, pn.Column(*use_params)))

        sidebar = pn.Accordion(*sidebar_tuples, active=[0])#, style={'font-size': "10"})
        return pn.Row(
                    sidebar,
                    img_ov,
                    self.Spectrum.get_roi_panel(),
                    background='white')


def pick_class(acq):
    arr = get_img_arr(acq)
    if len(arr.shape) == 2:
        return BasePanViewer(acq)
    else: return BaseSpecViewer(acq)

def Splitter(acq, wl_ranges = [(400,500), (500,600), (600,700)]):

    acq_class = BaseSpecViewer(acq)
    thresh = acq_class.threshold
    
    plots = []
    for wl in wl_ranges:
        plot = acq_class.get_dmap(band_override=wl)
        plots.append(plot)
    
    params = [acq_class.param.threshold, acq_class.param.alpha]
    return pn.Column(
                pn.Column(*params),
                pn.Row(*plots),
                background='white'
    )
#     def pixel_dist(self):
#         hv.Histogram(np.histogram(self.band_slice, bins=100))

#         plots.opts(active_tools=['box_edit']) # doesn't work... nor does putting in plot-type specific opts
        # widg =  pn.Param( self.param.center, widgets={'center':{'type':pn.widgets.NumberInput, 'width':300}} )

  

    # def test(self, band, center):
    #     band_slice = self.Spectrum.xds.sel(w=slice(*band)).mean('w')
    #     img = xds2img(band_slice)
    #     return img
    
    # def tview(self):
    #     streams = dict(band=self.param.band, center=self.param.center)
    #     dm = hv.DynamicMap(self.test, streams=streams)
    #     return pn.Column(self.param.band, dm)
    

        

     
 