
 #######################################
# gridwinder_plots.py
# Erica Lastufka 04-05-18

#Description: make transmission plot for Stefan's grids
#######################################

#######################################
# Usage:

# for default output: python imager_widget.py
######################################

from numpy import arange, sin, pi
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pickle
import os
import numpy as np
import glob        

def make_plot(xlog=True):
    '''Assumes a series of Imager objects are stored in the given location. Should be easily adaptable to values read from a spreadsheet however'''
    objects=glob.glob('/Users/wheatley/Documents/Solar/gridwinder/gridwinder*mm.p')
    print objects
    fig,ax=plt.subplots()
    labels=['1.8mm','1.9mm','2.2mm','2.8mm','3mm','4.3mm','4.6mm','4.9mm']
    colors=['k','c','m','g','r','y','b','0.75']
    for o,label,c in zip(objects,labels,colors):
        obj=pickle.load(open(o,'rb'))
        energy=list(obj.slits['Energy'])
        eff_area=list(obj.eff_area)
        sub=list(obj.substrate['eff_area'])
        slats=list(obj.slits['eff_area'])
        #print sub[-1],slats[-1]
        energy.extend(list(np.arange(433,1000,2)))
        eff_area.extend(list(eff_area[-1]+np.zeros(284)))
        sub.extend(list(sub[-1]+np.zeros(284)))
        slats.extend(list(slats[-1]+np.zeros(284)))
        slatt=np.array(slats)*(.5*np.array(sub))
        #print len(energy),len(sub)
   
        ax.plot(energy, .5*np.array(sub)+slatt, color=c,label=label,linewidth='2')
        #ax.plot(energy, .5*np.array(sub), color="g",linewidth='1')
        ax.plot(energy, slatt, color=c,linestyle='--',linewidth='2')

    ax.grid(which='both')
    ax.set_xlabel('Energy (keV)')
    ax.set_ylabel('Transmission')
    ax.set_ylim([0.001,2])
    ax.set_yscale('log')
    if xlog:
        ax.set_xscale('log')
        ax.set_xlim([0,1000])
    else:
        ax.set_xlim([0,120])
    plt.title("Transmission probability for a single grid "+ data.substrate['Material'] + ' '+ data.substrate['Thickness'] + '$\mu$m, ' + data.slits['Material'])
    #l1=ax.legend(["total","slits (C)","slats (C + Au)"],title='Minimum Configuration',bbox_to_anchor=(.99,.47),fontsize='medium', frameon=False)
    ax.legend(loc="lower right")
    #l2=ax.legend([p1,p2,p3],["total","slits (C)","slats (C+Au)"],title='Engineering Model',bbox_to_anchor=(.955,.25),fontsize='medium', frameon=False)
    #plt.gca().add_artist(l1)
    
    fig.show()


