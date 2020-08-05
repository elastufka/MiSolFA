
 #######################################
# both_bproj.py
# Erica Lastufka 15/5/19

#Description: make back projection images of MiSolFA and STIX
#######################################

#######################################
# Usage:

# for default output: python imager_widget.py
######################################

from scipy.io import readsav
from scipy.ndimage import rotate
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pickle
import os
import numpy as np
from CB_colors import CB_color_cycle as CB

def stix_plot(fname1='/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/stx_scenario_EL/stx-map.sav'):
    aa=readsav(fname1,python_dict=True)
    data=(aa['map']['data'][0]*.68)/48.
    cmap=plt.get_cmap('jet')
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n='jet', a=0., b=.8),
        cmap(np.linspace(0, 0.8, 100)))
    fig,ax=plt.subplots()
    p=ax.imshow(data,norm=matplotlib.colors.Normalize(vmin=0,vmax=1),cmap=new_cmap)
    cbar=fig.colorbar(p)
    cbar.ax.set_ylabel('Normalized Flux, 8-12 keV')
    ax.set_xlabel('X (arcsec)')
    ax.set_ylabel('Y (arcsec)')
    ax.set_ylim([0,131])
    ax.set_yticks([5,25,45,65,85,105,125])
    ax.set_xticks([5,25,45,65,85,105,125])
    ax.set_yticklabels(['-60','-40','-20','0','20','40','60'])
    ax.set_xticklabels(['-60','-40','-20','0','20','40','60'])
    ax.set_xlim([0,131])
    #plt.title("Effective area of grids "+ data.substrate['Material'] + ' '+ data.substrate['Thickness'] + '$\mu$m, ' + data.slits['Material'] + ' '+data.slits['Thickness'] +'$\mu$m with ' +'$\mu$m attenuator')

    fig.show()

def msfa_plot(fname1='/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/stx_scenario_EL_MISOFA/msfa-map.sav'):
    aa=readsav(fname1,python_dict=True)
    data=(aa['map']['data'][0])/48.
    cmap=plt.get_cmap('jet')
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n='jet', a=0., b=.8),
        cmap(np.linspace(0, 0.8, 100)))

    fig,ax=plt.subplots()
    p=ax.imshow(rotate(data,45.),norm=matplotlib.colors.Normalize(vmin=0,vmax=1),cmap=new_cmap)
    #cbar=fig.colorbar(p)
    #cbar.ax.set_ylabel('Normalized Flux, 4-12 keV')

    ax.set_xlabel('X (arcsec)')
    ax.set_ylabel('Y (arcsec)')
    ax.set_ylim([62,122])
    ax.set_yticks([72,82,92,102,112])
    ax.set_xticks([62,72,82,92,102,112,122])
    ax.set_yticklabels(['-40','-20','0','20','40'])
    ax.set_xticklabels(['-60','-40','-20','0','20','40','60'])
    ax.set_xlim([62,122])
    #plt.title("Effective area of grids "+ data.substrate['Material'] + ' '+ data.substrate['Thickness'] + '$\mu$m, ' + data.slits['Material'] + ' '+data.slits['Thickness'] +'$\mu$m with ' +'$\mu$m attenuator')

    fig.show()


