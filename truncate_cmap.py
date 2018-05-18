import pylab as plt
import matplotlib.colors as colors
import numpy as np

def truncate_colormap(cmap, minval=0, maxval=1, n=100):
    return colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f},)'.format(n=cmap.name,a=minval,b=maxval),
        cmap(np.linspace(minval, maxval, n)))
