"""
The classes here provide visualizations of Odemis hdf5 files containing panchromatic and hyperspectral acquisitions.
"""

import holoviews as hv
from holoviews import opts, dim, streams
import numpy as np
import panel as pn
import param

from cl_utils.utils import hdf2xds, get_ranges, stretch_arr, get_percentiles, get_extents
from cl_utils.hooks import y_bare_hook, transp_hook

css = '''
.bk.panel-widget-box {
  background: #f0f0f0;
  border-radius: 5px;
  border: 1px black solid;
}
'''
pn.extension(raw_css=[css])

class BasePanViewer(param.Parameterized):
    """
    Base class for generating dashboards from acquisitons drawn from Odemis hdf5 files. 

    acq : hdf5 group
        Typically selected from Odemis hdf5 file.
    """
    thresholds = param.Range(default=(2,98), bounds=(0,100), doc="outlier percentiles", label='thresholds')
    alpha = param.Number(default=0.5, bounds=(0,1), label='Alpha')

    def __init__(self, acq, **params):
        super().__init__(**params)
        self.xds = hdf2xds(acq)
        self.img_data = self.xds.cts
        self.title = self.xds.name
        self.acq = acq
        self.streams = dict(thresholds=self.param.thresholds)
        self.extents  = get_extents(self.xds.x, self.xds.y)

    def img_proc(self, img_data, thresholds):
        """
        Image processing steps.
        """
        # hacky way to make sure HoloViews displays thresholds correctly in overlay
        img_data = stretch_arr(img_data)
        lower_bound, upper_bound = get_percentiles(img_data, *thresholds)
        img_data = (img_data - lower_bound)/(upper_bound - lower_bound)
        
        img = hv.Image(img_data, ['x', 'y'])
        img = img.redim.range(**{'cts' : (lower_bound, upper_bound)})

        # set extents
        x0, y0, x1, y1 = self.extents
        img = img.redim(x=dict(range=(x0,x1)), y=dict(range=(y0,y1)))

        # set opts
        img_opts = opts.Image(frame_width=350, frame_height=350,
                                cmap='gray', xlabel='units in µm', ylabel='')
        return img.opts(img_opts)        

    def get_dmap(self):
        """
        Convert the current image into a parameterized DynamicMap.
        """
        def plot_img(thresholds):
            return self.img_proc(self.img_data, thresholds)
        dmap = hv.DynamicMap(plot_img, streams=self.streams)
        return dmap.apply.opts(alpha=self.param.alpha)

    def view(self):
        return pn.Row(self.param, self.get_dmap(), background='white')

class BaseSpecViewer(BasePanViewer):
    """ 
    Extends the functionality of BasePanViewer for Odemis acquisitions with a spectral dimension.

    acq : hdf5 group
        Typically selected from Odemis hdf5 file.
    """
    band = param.Range(default=(500, 501), bounds=(200,800), label='Band (nm)')
    center  = param.Number(default=500, step=1, doc="center of band", precedence=-1)
    bw  = param.Number(default=100, step=1, doc="center of band", label='BW', precedence=-1)

    def __init__(self, acq, **params):
        super().__init__(acq, **params)
        self.streams['band'] = self.param.band
        self.set_wl_bounds()

    def get_dmap(self, band_override=None):
        """
        Convert the current image into a parameterized DynamicMap.
        Optionally, set band_override to generate the DynamicMap for a specific band.
        
        Parameters
        ----------
        band_override : tuple
            
        Returns
        -------
        parameterized DynamicMap
        """
        def plot_img(thresholds, band):
            if band_override:
                band = band_override
            img_data = self.xds.cts.sel(wl=slice(*band)).mean('wl')
            return self.img_proc(img_data, thresholds)
        dmap = hv.DynamicMap(plot_img, streams=self.streams).opts(normalize=False)
        return dmap.apply.opts(alpha=self.param.alpha)
    
    def set_wl_bounds(self):
        """ 
        Set the default wavelength bounds.
        """
        wls = self.xds['wl'].values
        min_wl, max_wl = min(wls), max(wls)
        center = (max_wl+min_wl)/2
        min_wl, max_wl, center = round(min_wl), round(max_wl), round(center)
        self.param.band.bounds = (min_wl, max_wl)
        self.band = center - self.bw/2, center + self.bw/2
        self.center = center

    @param.depends('band', watch=True)
    def set_text_from_band(self):
        """ """
        center = np.average(self.band)
        bw = np.diff((self.band))[0]
        center, bw = np.round(center,1), np.round(bw,1)
        self.center = center
        self.bw = bw
        
    @param.depends('center', 'bw', watch=True)
    def set_band_from_text(self):
        """ """
        lower = self.center - self.bw/2
        upper = self.center + self.bw/2

    def get_band_dmap(self):
        """ 
        Get a DynamicMap which is an Overlay of the aggregated spectrum with a shaded box defining the current band selection.
        """
        agg_spec = self.xds.sum(('x','y'))
        def band_selection(band):
            agg_spec_curve = hv.Curve(agg_spec)
            band_vspan = hv.VSpan(*band).opts(color='gray', alpha=0.5)
            band_opts = opts.Overlay(toolbar=None, hooks=[transp_hook, y_bare_hook], xlabel='',\
                                        frame_height=50, responsive=True)
            return (agg_spec_curve*band_vspan).opts(band_opts)  
        return hv.DynamicMap(band_selection, streams=dict(band=self.param.band))

    def get_band_panel(self):
        """ 
        Get an organized panel containing the band controls and visualization.
        """
        input_params = pn.Row(self.param.center, self.param.bw)
        return  pn.Column(
                    input_params, 
                    pn.Column(
                        pn.Row(self.param.band, margin = (0, 0, -25, 0)),
                        pn.pane.HoloViews(self.get_band_dmap(), linked_axes=False)
                        ), 
                    css_classes=['panel-widget-box'],
                    )

    def view(self):
        controls = pn.Column(
                        self.param.alpha,
                        self.param.thresholds, 
                        self.get_band_panel()
                         )
        return pn.Row(controls, self.get_dmap(), background='white')
    
class BaseSpecViewerROI(BaseSpecViewer):
    """ 
    Adds ROI selection and associated spectrum visualization onto BaseSpecViewer.

    acq : hdf5 group
        Typically selected from Odemis hdf5 file.
    """
    roi_toggle = param.Selector(objects=["CL", "Trans"])
    color_cycle = param.List(['red', 'green', 'blue', 'orange', 'purple', 'cyan', \
                   'magenta', 'maroon', 'teal', 'navy', 'brown', 'pink'], precedence=-1)
    rect_coords = param.List(precedence=-1)

    def __init__(self, acq, **params):
        super().__init__(acq, **params)
        self.polys =  hv.Polygons(self.rect_coords)
        self.box_stream = streams.BoxEdit(source=self.polys,  styles={'fill_color': self.color_cycle})

    def roi_curves(self, data): 
        """
        Callable for dynamically streaming curves to ROI chart based on shape selections.
        """
        if not data or not any(len(d) for d in data.values()):
            default_curve = hv.Curve([], 'Spectrum', 'CL').opts(color='red') 
            return hv.NdOverlay({0: default_curve}).opts(show_legend=False) # code breaks without using a curve in ndoverlay
        
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

    def get_polys(self): 
        """ 
        Get the Polygons plot object that enables shape drawing.
        """
        poly_opts = opts.Polygons(fill_alpha=0.2, line_color='white')
        return self.polys.opts(poly_opts)

    def get_roi(self):
        """ 
        Get the ROI chart.
        """
        roi_opts = opts.Curve(frame_width=275, frame_height=125, framewise=True,
                                xlabel='Wavelength (nm)', ylabel='Intensity',
                                show_grid=True, gridstyle = {'minor_xgrid_line_color': 'lightgray'}, 
                                toolbar='above')
        roi = hv.DynamicMap(self.roi_curves, streams=[self.box_stream]).opts(roi_opts)
        roi.opts(xlim=(self.param.band.bounds)) 

        line_opts = opts.VLine(color='black', line_dash='dashed', line_width=1)
        line = hv.VLine(self.center).opts(line_opts)

        return roi*line

    def get_roi_panel(self):
        """ 
        Get a panel composed of the ROI chart and associated toggle between intensity and transmission.
        """
        roi_toggle = pn.Param(self.param.roi_toggle, 
                            widgets={'roi_toggle': {'type':pn.widgets.RadioButtonGroup, 'width':150}})
        return pn.Column(
                    pn.Column(roi_toggle, align='center', width_policy='max'),
                    pn.Row(self.get_roi),
                     css_classes=['panel-widget-box'])

    def view(self):
        img_params = pn.Column(self.param.alpha, self.param.thresholds, self.param.band)
        img_panel = pn.Row(
                        img_params,
                        self.get_dmap()*self.get_polys()
                        )
        return pn.Row(
                    img_panel,
                    self.get_roi_panel(),
                    background='white'
                    )

class Spectrum(param.Parameterized):
    """ 
    Generates a dashboard containing the first three acquisitions of an Odemis hdf5 file with associated controls and ROI interactivity.

    h5 : in-memory hdf5 file
        Typically an Odemis hdf5 file.
    """
    def __init__(self, h5, **params):
        super().__init__(**params)
        self.process_h5(h5)


    def process_h5(self, h5): 
        self.Survey = BasePanViewer(h5['Acquisition0'])
        self.Concurrent = BasePanViewer(h5['Acquisition1'])
        self.Spectrum = BaseSpecViewerROI(h5['Acquisition2'])

    def get_img_overlay(self):
        img_ov = hv.Overlay([
                    self.Survey.get_dmap(),
                    self.Concurrent.get_dmap(),
                    self.Spectrum.get_dmap()*self.Spectrum.get_polys()]) 
        img_ov.opts(toolbar='above', normalize=False)
        return img_ov.collate()

    def view(self):
        sidebar_tuples = []
        for class_ in [self.Spectrum, self.Concurrent, self.Survey]:
            use_params = [class_.param.alpha, class_.param.thresholds]
            if class_.title == 'Spectrum':
                use_params.append(class_.get_band_panel())
                class_.alpha = 0.5
            else: class_.alpha = 1         
            sidebar_tuples.append((class_.title, pn.Column(*use_params)))

        sidebar = pn.Accordion(*sidebar_tuples, active=[0])
        return pn.Row(
                    sidebar,
                    pn.pane.HoloViews( self.get_img_overlay(), linked_axes=False),
                    self.Spectrum.get_roi_panel(),
                    background='white')







     
 