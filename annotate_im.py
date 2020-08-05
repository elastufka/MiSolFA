
 #######################################
# both_bproj.py
# Erica Lastufka 15/5/19

#Description: make back projection images of MiSolFA and STIX
#######################################

#######################################
# Usage:

# for default output: python imager_widget.py
######################################

from PIL import Image
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from CB_colors import CB_color_cycle as CB

def annotate_moire_im(EM1=True, EM2=False):
    if EM1:
        ims=glob.glob('/Users/wheatley/Documents/Solar/MiSolFA/plots/EM1seg*_tcat_moire.png')
        labels=[['11: 90$\mu$m','21: 90$\mu$m','12: 22.5$\mu$m','22: 22.5$\mu$m'],['31: 45$\mu$m','41: 45$\mu$m','32: 18$\mu$m','42: 18$\mu$m'],['33: 30$\mu$m','43: 30$\mu$m','34: 15$\mu$m','44: 15$\mu$m']]
        coords=[[25,35],[25,550],[700,35],[700,550]]
    elif EM2:
        ims=glob.glob('/Users/wheatley/Documents/Solar/MiSolFA/plots/EM2seg*_tcat_moire.png')
        ims.sort()
        labels=[['11: 90$\mu$m','21: 90$\mu$m','12: 45$\mu$m','22: 45$\mu$m'],['31: 30$\mu$m','41: 30$\mu$m','32: 22.5$\mu$m','42: 22.5$\mu$m'],['33: 18$\mu$m','43: 18$\mu$m','34: 15$\mu$m','44: 15$\mu$m']]
        coords=[[25,35],[25,575],[450,35],[450,575]]

    font = {'family': 'sans-serif',
        'color':  'white',
        'weight': 'normal',
        'size': 16,
        }
    font2 = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

    for im, nums in zip(ims,labels):
        fig,ax=plt.subplots()
        ii=Image.open(im)
        ax.imshow(ii, origin='upper')
        for c,n in enumerate(nums):
            if c !=1:
                ax.text(coords[c][0],coords[c][1],n,fontdict=font)#,fontcolor='white')
            else:
                ax.text(coords[c][0],coords[c][1],n,fontdict=font2)#,fontcolor='white')
        plt.axis('off')
        fig.show()

