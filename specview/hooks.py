# for bokeh

def no_border_hook(plot, element):
    plot.handles['plot'].min_border_left = 0
    plot.handles['plot'].min_border_right = 0
#     plot.handles['plot'].min_border_top = 0
#     plot.handles['plot'].min_border_bottom = 0

def y_bare_hook(plot, element):
    plot.handles['yaxis'].visible = False

def transp_hook(plot, element):
    plot.handles['plot'].border_fill_alpha = 0