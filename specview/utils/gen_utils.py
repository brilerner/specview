"""
Utilities for general usage (unreliant on hdf5 file structure).
"""

import numpy as np

def make_rect_coords(xside=10, yside=10, xoffset=0, yoffset=0):
    """
    Parameters
    ----------
    xside : float
         (Default value = 10)
    yside : float
         (Default value = 10)
    xoffset : float
         (Default value = 0)
    yoffset : float
         (Default value = 0)

    Returns
    -------
     dict of tuples of corners : {('x','y'):[(x0,y0), (x1,y0), (x1,y1), (x0,y1)]}
    """
    x0 = 0 + xoffset
    x1 = x0 + xside
    y0 = 0 + yoffset
    y1 = y0 + yside
    return {('x','y'):[(x0,y0), (x1,y0), (x1,y1), (x0,y1)]}

def stretch_arr(img, min_val=0, max_val=1):
    """
    Parameters
    ----------
    img : array-like
        
    min_val : float
         (Default value = 0)
    max_val : float
         (Default value = 1)

    Returns
    -------

    """
    min_ = img.min()
    max_ = img.max()

    a = (max_val - min_val) / (max_ - min_)
    b = -a * min_ + min_val

    return a * img + b

def get_extents(x, y): 
     """
     Calculates square extents for the plot while retaining proper image proportions.

     Parameters
     ----------
     x : array-like

     y : arrray-like


     Returns
     -------
     (x0, y0, x1, y1) : tuple
          The image extents.
     """
     x0, x1 = min(x), max(x)
     y0, y1 = min(y), max(y)
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

def get_percentiles(arr, lthresh=2, uthresh=98):
    """
    Parameters
    ----------
    arr : array-like
        
    lthresh : float
         (Default value = 2)
    uthresh : float
         (Default value = 98)

    Returns
    -------
    lower_percentile, upper_percentile : tuple
    """
    return np.percentile(arr, lthresh), np.percentile(arr, uthresh)