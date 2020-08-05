
 #######################################
# simulator_gui.py
# Erica Lastufka 10/2/17

#Description: Make the widget for getting effective area plots for the grids
#######################################

#######################################
# Usage:

# for default output: python imager_widget.py
######################################

from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pickle
import os
import numpy as np
import pandas as pd

from CB_colors import CB_color_cycle as CB

def transm_plot(fname1='/Users/wheatley/Documents/Solar/papers/thesis_fig/stix_grids_4_rounded.csv',fname2='/Users/wheatley/Documents/Solar/MiSolFA/grids/transmEM1.csv',fname3='/Users/wheatley/Documents/Solar/MiSolFA/grids/transmEM2.csv',EM='1'):
    stix=pd.read_csv(fname1)
    if EM == '1':
        msfa=pd.read_csv(fname2)
    else:
        msfa=pd.read_csv(fname3)

    #STIX
    sftransm=stix['Front Slit width [mm]']/stix['Front Pitch [mm]']
    srtransm=stix['Rear Slit width']/stix['Rear Pitch']
    bftransm=1.-stix['Front Bridge width [mm]']/stix['Front Bridge pitch [mm]']
    brtransm=1.-stix['Rear Bridge width']/stix['Rear Bridge pitch']
    spitch=np.mean([stix['Front Pitch [mm]'],stix['Rear Pitch']],axis=0)*1000.
    swin=stix['\xef\xbb\xbfWindow']
    for i,w,sf,sr in zip(stix.index,swin,sftransm,srtransm):
        if w in [11,12,13,17,18,19]:
            n=float(np.shape(np.where(swin==w))[1])
            if n == 3:
                m=6
            else:
                m=3
            sftransm[i]=(stix['Front Slit width [mm]'][i]/m)/(stix['Front Pitch [mm]'][i]/(n))
            srtransm[i]=(stix['Rear Slit width'][i]/m)/(stix['Rear Pitch'][i]/(n))
            spitch[i]=spitch[i]/n#np.mean([stix['Front Pitch [mm]'][i]/n,stix['Rear Pitch'][i]/n],axis=0)*1000.
            bftransm[i]=1.-(stix['Front Bridge width [mm]'][i]/m)/(stix['Front Bridge pitch [mm]'][i]/n)
            brtransm[i]=1.-(stix['Rear Bridge width'][i]/m)/(stix['Rear Bridge pitch'][i]/n)
            print n,m,sftransm[i],spitch[i]
    stot=sftransm*srtransm


    #MiSolFA
    mftransm=msfa['Front Width']/msfa['Front Pitch']
    mrtransm=msfa['Rear Width']/msfa['Rear Pitch']
    mpitch=np.mean([msfa['Front Pitch'],msfa['Rear Pitch']],axis=0)
    mtot=mftransm*mrtransm
    print mtot

    fig,ax=plt.subplots()
    ax.scatter(spitch, sftransm, color=CB[0],marker='+',s=60,label="STIX Front")
    ax.scatter(spitch, srtransm, color=CB[1],marker='+',s=60,label="STIX Rear")
    ax.scatter(spitch, bftransm, color=CB[2],marker='+',s=60,label="STIX Bridges")
    ax.scatter(spitch, brtransm, color=CB[2],marker='+',s=60)
    ax.scatter(spitch, stot, color='k',marker='+',s=60,label="STIX Total")

    ax.scatter(mpitch, mftransm, color=CB[0],label="MiSolFA EM"+EM+" Front")
    ax.scatter(mpitch, mrtransm, color=CB[1],label="MiSolFA EM"+EM+" Rear")
    ax.scatter(mpitch, mtot, color='k',label="MiSolFA Total")

    ax.set_xlabel('Pitch ($\mu$m)')
    ax.set_ylabel('Transmission')
    ax.set_ylim([0,1])
    #ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim([10,1000])
    #plt.title("Effective area of grids "+ data.substrate['Material'] + ' '+ data.substrate['Thickness'] + '$\mu$m, ' + data.slits['Material'] + ' '+data.slits['Thickness'] +'$\mu$m with ' +'$\mu$m attenuator')
    #l1=ax.legend(["total","slits (C)","slats (C + Au)"],title='Nominal Configuration',bbox_to_anchor=(.99,.47),fontsize='medium', frameon=False)
    #l2=ax.legend([p1,p2,p3],["total","slits (C)","slats (C+Au)"],title='Acheived Thickness',bbox_to_anchor=(.955,.25),fontsize='medium', frameon=False)
    #plt.gca().add_artist(l1)
    ax.legend(loc='upper left')
    fig.show()





