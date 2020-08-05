"""
===================
mod_eff.py
Erica  Lastufka 9.3.18
===================

Calculate modulation efficiency from the moire patterns from the X-ray lamp test images

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from PIL import Image
from matplotlib import cm
from skimage import feature, exposure
from skimage.transform import  (hough_line, hough_line_peaks,
                               probabilistic_hough_line)
import pickle
from scipy.misc import imrotate
from scipy.ndimage.interpolation import rotate
import time
import glob
from scipy import optimize
from plot_divergent import plot_moire_fancy
import edge_and_hough as eh

#list of dictionaries of the window numbers, pitches (mm), and nominal angles (deg)
global moire

moire=[{'number':11,'period_px':99., 'angle': 31.,'amp':.10},
    {'number':21,'period_px':58., 'angle': -35.,'amp':.10},
    {'number':12,'period_px':43., 'angle': -49.,'amp':.10},
    {'number':22,'period_px':22., 'angle': 42,'amp':.10},
    {'number':31,'period_px':60., 'angle': 49.,'amp':.10},
    {'number':41,'period_px':40., 'angle': -50.,'amp':.10},
    {'number':32,'period_px':18., 'angle': -42.,'amp':.10},
    {'number':42,'period_px':22., 'angle': -42.,'amp':.10},
    {'number':33,'period_px':40., 'angle': 39.,'amp':.10},
    {'number':43,'period_px':35., 'angle': -39.,'amp':.10},
    {'number':34,'period_px':22., 'angle': -45.,'amp':.10},
    {'number':44,'period_px':18., 'angle': 45.,'amp':.10}]

def im2ndarray(filen):
    '''Convert .tif image to numpy.ndarray'''
    imraw = Image.open(filen) #.rotate(45)#.convert('1').convert('L')
    im=np.array(imraw)
    return im

def remove_edges(im,ncrop=50,plot=False):
    '''remove edges from images with edges of the window in them'''
    imshape=np.shape(im)
    im=im[ncrop:imshape[0]-ncrop,ncrop:imshape[1]-ncrop]
    percentage = np.float((np.shape(im)[0]*np.shape(im)[1]))/np.float((imshape[0]*imshape[1]))
    print "Percentage cropped:", (1.-percentage)*100.
    #print np.shape(im),imshape
    if plot:
        fig,ax=plt.subplots()
        ax.imshow(im,cmap=cm.gray)
        fig.show()
    return im

def contrast_stretch(im,plow=2,phigh=98):
    pl = np.percentile(im, plow)
    ph = np.percentile(im, phigh)
    im2 = exposure.rescale_intensity(im, in_range=(pl, ph))
    return im2

def mod_eff(fname,plot=False,cstretch=False):
    im=im2ndarray(fname)
    im=remove_edges(im,ncrop=30,plot=plot)
    if cstretch:
        im=contrast_stretch(im)
    immax=np.max(im)
    immin=np.min(im)
    mean=np.mean(im)
    ratio_max=mean/immax
    ratio_min=mean/immin
    ratio_ratios=ratio_max/ratio_min
    print ratio_ratios
    return ratio_ratios

if __name__ == '__main__':
    fnames=glob.glob('win*.tif')
    lines=[]
    for fn in fnames:
        rr=mod_eff(fn,plot=False,cstretch=True)
        lines.append(fn+', ' +str(rr)+'\n')
    with open('mod_eff_output_cstretch.txt','w') as f:
        f.writelines(lines)
