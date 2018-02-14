
 #######################################
# thermal_test_profile.py
# Erica Lastufka 1/2/18  

#Description: Plot the profile used for the thermal test
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
        

def plot_fig():
    timevec=list(np.linspace(0,1,13))[0:7]+range(1,33,1)+[32.5,33,34]
    tempvec=np.zeros(len(timevec))
    print len(timevec),len(tempvec)
    tempvec[0:5]=np.nan
    tempvec[5]=20
    tempvec[6]=20 #now start the cycle
    cycle=[t for t in [70,70,-30,-30]]
    for n in [7,11,15,19,23,27,31,35]:
        tempvec[n:n+4]=cycle

    tempvec[-3:]=20
        
    fig,ax=plt.subplots(figsize=(14,5))
    ax.plot(timevec, tempvec, color="r",label="qualif.",linewidth='2')
    ax.plot(timevec, 70+np.zeros(len(timevec)),"k--")
    ax.plot(timevec, 20+np.zeros(len(timevec)),"k--")
    ax.plot(timevec, -30+np.zeros(len(timevec)),"k--")
    for t in range(5,len(timevec)):
        #print tempvec[t]/120.+.333
        ax.axvline(x=timevec[t],ymax=tempvec[t]/120.+.33,color='k', linestyle='--')
    #ax.grid(which='both')
    ax.set_xlabel('$t$ (h)')
    plt.xticks(range(1,35,1))
    plt.yticks([-30,0,20,70])
    ax.set_ylabel('T ($\circ$C)')
    ax.set_ylim([-40,80])
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_xlim([0,34])
    #plt.title("Effective area of grids "+ data.substrate['Material'] + ' '+ data.substrate['Thickness'] + '$\mu$m, ' + data.slits['Material'] + ' '+data.slits['Thickness'] +'$\mu$m with ' +'$\mu$m attenuator')
    #ax.legend()
    #l2=ax.legend([p1,p2,p3],["total","slits (C)","slats (C+Au)"],title='Engineering Model',bbox_to_anchor=(.955,.25),fontsize='medium', frameon=False)
    #plt.gca().add_artist(l1)
    
    fig.show()
    return timevec,tempvec


