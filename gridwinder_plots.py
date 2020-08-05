
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
import imager

def get_objects(slit_material,thicknesses,sub='C'):
    '''Assume substrate is C'''
    objects=[]
    for th in thicknesses:
        int_th=int(float(th[:-2])*1000.)
        o=imager.Imager(slits=[slit_material,str(int_th)],substrate=[sub,str(int_th)],filled=True)
        objects.append(o)
    return objects

def make_plot(objects=False, xlog=True):
    '''Assumes a series of Imager objects are stored in the given location. Should be easily adaptable to values read from a spreadsheet however'''
    if not objects:
        objects=glob.glob('/Users/wheatley/Documents/Solar/gridwinder/gridwinder*mm.p')
        print objects
    fig,ax=plt.subplots()
    labels=['1.8mm','1.9mm','2.2mm','2.8mm','3mm','4.3mm','4.6mm','4.9mm']
    colors=['k','c','m','g','r','y','b','0.75']
    for o,label,c in zip(objects,labels,colors):
        if not objects:
            obj=pickle.load(open(o,'rb'))
        else:
            obj=o
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
    #plt.title("Transmission probability for a single grid "+ data.substrate['Material'] + ' '+ data.substrate['Thickness'] + '$\mu$m, ' + data.slits['Material'])
    #l1=ax.legend(["total","slits (C)","slats (C + Au)"],title='Minimum Configuration',bbox_to_anchor=(.99,.47),fontsize='medium', frameon=False)
    ax.legend(loc="lower right")
    #l2=ax.legend([p1,p2,p3],["total","slits (C)","slats (C+Au)"],title='Engineering Model',bbox_to_anchor=(.955,.25),fontsize='medium', frameon=False)
    #plt.gca().add_artist(l1)

    fig.show()

def make_plot_future(material='Steel',xlog=True):
    '''Assumes a series of Imager objects are stored in the given location. Should be easily adaptable to values read from a spreadsheet however'''
    labels=['1.8mm','1.9mm','2.2mm','2.8mm','3mm','4.3mm','4.6mm','4.9mm']
    objects=get_objects(material,labels)
    #objects=glob.glob('/Users/wheatley/Documents/Solar/gridwinder/gridwinder*mm.p')
    #print objects
    fig,ax=plt.subplots()
    #labels=['1.8mm','1.9mm','2.2mm','2.8mm','3mm','4.3mm','4.6mm','4.9mm']
    colors=['k','c','m','g','r','y','b','0.75']
    for obj,label,c in zip(objects,labels,colors):
        #obj=pickle.load(open(o,'rb'))
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
    #plt.title("Transmission probability for a single grid "+ obj.slits['Material'] + ' '+ obj.substrate['Thickness'] + '$\mu$m, ' + obj.slits['Material'])
    plt.title('Tungsten slats, Carbon substrate')
    #l1=ax.legend(["total","slits (C)","slats (C + Au)"],title='Minimum Configuration',bbox_to_anchor=(.99,.47),fontsize='medium', frameon=False)
    #ax.legend(loc="lower right")
    #l2=ax.legend([p1,p2,p3],["total","slits (C)","slats (C+Au)"],title='Engineering Model',bbox_to_anchor=(.955,.25),fontsize='medium', frameon=False)
    #plt.gca().add_artist(l1)

    fig.show()

def make_plot_proposal(material=['Steel','W'],xlog=True):
    '''Assumes a series of Imager objects are stored in the given location. Should be easily adaptable to values read from a spreadsheet however'''
    #titles=['Steel slats, Carbon substrate','Tungsten slats, Carbon substrate']
    allslats=[]
    energies=[]
    for i,m in enumerate(material):
        objects=get_objects(m,['2.0mm'])
        fig,ax=plt.subplots()
        obj=objects[0]
        energy=list(obj.slits['Energy'])
        #eff_area=list(obj.eff_area)
        sub=list(obj.substrate['eff_area'])
        slats=list(obj.slits['eff_area'])
        #print sub[-1],slats[-1]
        energy.extend(list(np.arange(433,1000,2)))
        energies.append(energy)
        #eff_area.extend(list(eff_area[-1]+np.zeros(284)))
        sub.extend(list(sub[-1]+np.zeros(284)))
        slats.extend(list(slats[-1]+np.zeros(284)))
        allslats.append(slats)
        slatt=np.array(slats)*(.5*np.array(sub))
        #print len(energy),len(sub)

    ax.plot(energy, np.array(sub),label='C',linewidth='2')
    ax.plot(energy, np.array(allslats[0]), label='Steel', color="g",linewidth='2')
    ax.plot(energy, np.array(allslats[1]),color='c',linewidth='2',label='W')
    ax.legend(loc='upper left')
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
    #plt.title("Transmission probability for a single grid "+ obj.slits['Material'] + ' '+ obj.substrate['Thickness'] + '$\mu$m, ' + obj.slits['Material'])
        #plt.title(titles[i])
    #l1=ax.legend(["total","slits (C)","slats (C + Au)"],title='Minimum Configuration',bbox_to_anchor=(.99,.47),fontsize='medium', frameon=False)
    #ax.legend(loc="lower right")
    #l2=ax.legend([p1,p2,p3],["total","slits (C)","slats (C+Au)"],title='Engineering Model',bbox_to_anchor=(.955,.25),fontsize='medium', frameon=False)
    #plt.gca().add_artist(l1)

    fig.show()

def make_plot_HXI(material=['Au','W'],thick='1.0mm',xlog=True,title=False):
    '''Assumes a series of Imager objects are stored in the given location. Should be easily adaptable to values read from a spreadsheet however'''
    #titles=['Steel slats, Carbon substrate','Tungsten slats, Carbon substrate']
    allslats,allsubs=[],[]
    energies=[]
    subs=['C','Si']
    for i,m in enumerate(material):
        objects=get_objects(m,[thick],sub=subs[i])
        fig,ax=plt.subplots()
        obj=objects[0]
        energy=list(obj.slits['Energy'])
        #eff_area=list(obj.eff_area)
        sub=list(obj.substrate['eff_area'])
        slats=list(obj.slits['eff_area'])
        #print sub[-1],slats[-1]
        energy.extend(list(np.arange(433,1000,2)))
        energies.append(energy)
        #eff_area.extend(list(eff_area[-1]+np.zeros(284)))
        sub.extend(list(sub[-1]+np.zeros(284)))
        slats.extend(list(slats[-1]+np.zeros(284)))
        allslats.append(slats)
        allsubs.append(sub)
        slatt=np.array(slats)*(.5*np.array(sub))
        #print len(energy),len(sub)

    ax.plot(energy, np.array(allsubs[0]),label='C',linewidth='2')
    ax.plot(energy, np.array(allsubs[1]),label='Si',linewidth='2')
    ax.plot(energy, 2*np.array(allslats[0]), label='Au', color="m",linewidth='2')
    ax.plot(energy, 2*np.array(allslats[1]),color='c',linewidth='2',label='W')
    ax.plot(energy, np.array(allsubs[0])-2*np.array(allslats[0]), label='C-Au', color="k",linewidth='2')
    ax.plot(energy, np.array(allsubs[0])-2*np.array(allslats[1]),color='r',linewidth='2',label='Si-W')
    ax.legend(loc='lower right')
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
    #plt.title("Transmission probability for a single grid "+ obj.slits['Material'] + ' '+ obj.substrate['Thickness'] + '$\mu$m, ' + obj.slits['Material'])
        #plt.title(titles[i])
    #l1=ax.legend(["total","slits (C)","slats (C + Au)"],title='Minimum Configuration',bbox_to_anchor=(.99,.47),fontsize='medium', frameon=False)
    #ax.legend(loc="lower right")
    #l2=ax.legend([p1,p2,p3],["total","slits (C)","slats (C+Au)"],title='Engineering Model',bbox_to_anchor=(.955,.25),fontsize='medium', frameon=False)
    #plt.gca().add_artist(l1)
    if title:
        ax.set_title(title)

    fig.show()

def make_plot_HXI_linear(material=['Au','W'],thick=['.25mm','1.0mm'],xlog=True,title=False):
    '''Assumes a series of Imager objects are stored in the given location. Should be easily adaptable to values read from a spreadsheet however'''
    #titles=['Steel slats, Carbon substrate','Tungsten slats, Carbon substrate']
    allslats,allsubs=[],[]
    energies=[]
    subs=['C','Si']
    for i,m in enumerate(material):
        for t in thick:
            objects=get_objects(m,[t],sub=subs[i])
            obj=objects[0]
            energy=list(obj.slits['Energy'])
            #eff_area=list(obj.eff_area)
            sub=list(obj.substrate['eff_area'])
            slats=list(obj.slits['eff_area'])
            #print sub[-1],slats[-1]
            energy.extend(list(np.arange(433,1000,2)))
            energies.append(energy)
            #eff_area.extend(list(eff_area[-1]+np.zeros(284)))
            sub.extend(list(sub[-1]+np.zeros(284)))
            slats.extend(list(slats[-1]+np.zeros(284)))
            allslats.append(slats)
            allsubs.append(sub)
            slatt=np.array(slats)*(.5*np.array(sub))
        #print len(energy),len(sub)
    print np.shape(allslats),np.shape(allsubs)
    fig,ax=plt.subplots()

    #ax.plot(energy, np.array(allsubs[0]),label='C',linewidth='2')
    #ax.plot(energy, np.array(allsubs[1]),label='Si',linewidth='2')
    #ax.plot(energy, 2*np.array(allslats[0]), label='Au', color="m",linewidth='2')
    ax.plot(energy, 1.0-2*np.array(allslats[2]),color='g',linewidth='2',label='.25 mm W')
    ax.plot(energy, np.array(allsubs[0])-2*np.array(allslats[0]), label='.25 mm C-Au', color="k",linewidth='2')
    ax.plot(energy, 1.0-2*np.array(allslats[3]),'m--',linewidth='2',label='1mm W')
    ax.plot(energy, np.array(allsubs[1])-2*np.array(allslats[1]), "b--", label='1mm C-Au',linewidth='2')
    #ax.plot(energy, np.array(allsubs[0])-2*np.array(allslats[1]),color='r',linewidth='2',label='Si-W')
    ax.legend(loc='lower center')
    ax.grid(which='both')
    ax.set_xlabel('Energy (keV)')
    ax.set_ylabel('Transmission')
    ax.set_ylim([0,1.1])
    ax.set_xlim([3,400])
    #ax.set_yscale('log')
    if xlog:
        ax.set_xscale('log')
        #ax.set_xlim([0,1000])
    else:
        ax.set_xlim([0,120])
    #plt.title("Transmission probability for a single grid "+ obj.slits['Material'] + ' '+ obj.substrate['Thickness'] + '$\mu$m, ' + obj.slits['Material'])
        #plt.title(titles[i])
    #l1=ax.legend(["total","slits (C)","slats (C + Au)"],title='Minimum Configuration',bbox_to_anchor=(.99,.47),fontsize='medium', frameon=False)
    #ax.legend(loc="lower right")
    #l2=ax.legend([p1,p2,p3],["total","slits (C)","slats (C+Au)"],title='Engineering Model',bbox_to_anchor=(.955,.25),fontsize='medium', frameon=False)
    #plt.gca().add_artist(l1)
    if title:
        ax.set_title(title)

    fig.show()
