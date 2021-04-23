"""
The classes here combine the viewer classes into complex dashboards.
"""

import holoviews as hv
import panel as pn
pn.extension()
import param

from specview.viewers import BasePanViewer, BaseSpecViewerROI

class MultiView(param.Parameterized):
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
        return img_ov.collate().opts(toolbar='above')


    def view(self):
        sidebar_tuples = []
        for class_ in [self.Spectrum, self.Concurrent, self.Survey]:
            class_.title = {
                            'Spectrum' : 'Spectrum',
                            'Secondary electrons concurrent' : 'Concurrent',
                            'Secondary electrons survey' : 'Survey'
                            }[class_.title]
            if class_.title == 'Spectrum':
                class_.alpha = 0.5
            else: class_.alpha = 1         
            sidebar_tuples.append((class_.title, class_.get_controls()))

        sidebar = pn.Accordion(*sidebar_tuples, active=[0], width=160)

        return pn.Row(
                    sidebar,
                    pn.pane.HoloViews( self.get_img_overlay(), linked_axes=False),
                    self.Spectrum.get_roi_panel(),
                    background='white')










     
 