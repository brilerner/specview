"""
HoloViews hooks for use with bokeh backend.
"""

def y_bare_hook(plot, element):
    plot.handles['yaxis'].visible = False

def transp_hook(plot, element):
    plot.handles['plot'].border_fill_alpha = 0