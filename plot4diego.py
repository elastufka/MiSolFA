
 #######################################
# simulator_gui.py
# Erica Lastufka 10/2/17

#Description: Make the widget for getting effective area plots for the grids
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
from transmission import read_data
from CB_colors import CB_color_cycle as CB

def eff_area(fname1='test1.p', fname2='EM.p'):
    obj=pickle.load(open(fname1,'rb'))
    em=pickle.load(open(fname2,'rb'))

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

    emenergy=list(em.slits['Energy'])
    emeff_area=list(em.eff_area)
    emsub=list(em.substrate['eff_area'])
    emslats=list(em.slits['eff_area'])
    emenergy.extend(list(np.arange(433,1000,2)))
    emeff_area.extend(list(emeff_area[-1]+np.zeros(284)))
    emsub.extend(list(emsub[-1]+np.zeros(284)))
    emslats.extend(list(emslats[-1]+np.zeros(284)))
    emslatt=np.array(emslats)*(.5*np.array(emsub))

    fig,ax=plt.subplots()
    ax.plot(energy, .5*np.array(sub)+slatt, color=CB[0],label="total",linewidth='2')
    ax.plot(energy, .5*np.array(sub), color=CB[1],label="slits (C)",linewidth='2')
    ax.plot(energy, slatt, color=CB[2],label="slats (C + Au)",linewidth='2')
    p1,=ax.plot(emenergy, .5*np.array(emsub)+emslatt, color=CB[0],linestyle="--",label="total",linewidth='2')
    p2,=ax.plot(emenergy, .5*np.array(emsub),color=CB[1],linestyle="--",label="slits (C)",linewidth='2')
    p3,=ax.plot(emenergy,emslatt,color=CB[2],linestyle="--",label="slats (C + Au)",linewidth='2')
    #ax.plot(emenergy, .1+np.zeros(784), 'k--')
    #ax.plot(emenergy, .2+np.zeros(784), 'k--')
    ax.grid(which='both')
    ax.set_xlabel('Energy (keV)')
    ax.set_ylabel('Transmission')
    ax.set_ylim([0.001,1])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim([0,1000])
    #plt.title("Effective area of grids "+ data.substrate['Material'] + ' '+ data.substrate['Thickness'] + '$\mu$m, ' + data.slits['Material'] + ' '+data.slits['Thickness'] +'$\mu$m with ' +'$\mu$m attenuator')
    l1=ax.legend(["total","slits (C)","slats (C + Au)"],title='Nominal Configuration',bbox_to_anchor=(.99,.47),fontsize='medium', frameon=False)
    l2=ax.legend([p1,p2,p3],["total","slits (C)","slats (C+Au)"],title='Acheived Thickness',bbox_to_anchor=(.955,.25),fontsize='medium', frameon=False)
    plt.gca().add_artist(l1)

    fig.show()


def transm_all(materials=['Au','W','Au80Sn20'], thickness='P(xi) (d=250 um)'):
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/calculations/data')
    fig,ax=plt.subplots()
    cc=CB[:3]#['c','m','k']
    for m,c in zip(materials,cc):
        pd = read_data(m+'.csv')
        ax.plot(pd['E (keV)'], np.array(pd[thickness]), color=c,label=m,linewidth='2')
    ax.grid(which='both')
    ax.set_xlabel('Energy (keV)')
    ax.set_ylabel('Transmission')
    ax.set_ylim([0.001,2])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim([1,1000])
    ax.legend(loc='upper left')
    #plt.title("Effective area of grids "+ data.substrate['Material'] + ' '+ data.substrate['Thickness'] + '$\mu$m, ' + data.slits['Material'] + ' '+data.slits['Thickness'] +'$\mu$m with ' +'$\mu$m attenuator')
   # l1=ax.legend(["total","slits (C)","slats (C + Au)"],title='Minimum Configuration',bbox_to_anchor=(.99,.47),fontsize='medium', frameon=False)
    #l2=ax.legend([p1,p2,p3],["total","slits (C)","slats (C+Au)"],title='Engineering Model',bbox_to_anchor=(.955,.25),fontsize='medium', frameon=False)
    #plt.gca().add_artist(l1)

    fig.show()



