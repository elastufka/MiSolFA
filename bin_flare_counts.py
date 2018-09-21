#######################################
# bin_flare_counts.py
# Erica Lastufka 11/11/2016  

#Description: Bin flare counts in different ways, plot error, see how it looks. don't judge me on this piece of code and its messiness
#######################################

#######################################
# Usage:

# for default output: python transmission.py
######################################

#import glob
#import os
import numpy as np
#from numpy import genfromtxt
import matplotlib.pyplot as plt
#import Image
#import csv
import scipy.constants as sc
from transmission import read_data
#import grid_eff_area
import pickle

def bin_flare_counts():
    plotdata =pickle.load(open('EMplotdata.p','rb'))
    #plotdata = evector,eff_area,sub,slit,substrate_properties,slit_properties,ntnt,thth,fac, attenuator
    energy = plotdata[0]
    eff_area = plotdata[1]
    sub = plotdata[2]
    slit = plotdata[3]
    substrate_properties = plotdata[4]
    slit_properties = plotdata[5]
    ntnt = plotdata[10] #these are the ones after applying the grids
    thth = plotdata[9]
    fac = plotdata[8]

    title= substrate_properties['material'] + ' '+ substrate_properties['thickness'] + '$\mu$m substrate with '+ slit_properties['material'] + ' '+ slit_properties['thickness'] + '$\mu$m grids'

    test, binEdges=np.histogram(energy, bins=np.logspace(0, 2.4, 32)) #these are the energy bins. Now we need to fill them with the data ....
    binEdges[np.where(binEdges[1:] - binEdges[:-1] < 1.0)] =0 #if a bin is less than 1 keV wide, get rid of it

    maybe = np.digitize(energy, binEdges) # now each data point has an index of which bin it belongs in
    #put the data in bins
    total_nt=np.linspace(1,32,32)
    total_th=np.linspace(1,32,32)
    
    for i in range(0,31):
        foo = np.where(maybe == i)
        if foo:
            #print np.shape(eff_area),foo, np.sum(eff_area[foo])
            total_nt[i] = np.sum(ntnt[foo])#sum the values of test at those points and assign it to output vector
            total_th[i] = np.sum(thth[foo])
            #print total[i]
    factor = 12#*30 #12 detectors, 2 grids, 30s integration time
    total_nt = total_nt*factor
    total_th = total_th*factor
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    menStd_nt     = 1/np.sqrt(total_nt[0:31])
    menStd_th     = 1/np.sqrt(total_th[0:31])
    infs = np.where(menStd_th == float("inf")) and np.where(menStd_th > 10**2)
    menStd_th[infs] = 0
    infs = np.where(menStd_nt == float("inf")) and np.where(menStd_nt > 10**2)
    menStd_nt[infs] = 0
    
    width      = (binEdges[1:]-binEdges[:-1])

    evector, nbinEdges=np.histogram(energy, bins=np.linspace(0, 250, 32)) #these are the energy bins. Now we need to fill them with the data ....
    index = np.digitize(energy, nbinEdges) # now each data point has an index of which bin it belongs in
    nbinEdges[np.where(nbinEdges[1:] - nbinEdges[:-1] < 1.0)] =0 #if a bin is less than 1 keV wide, get rid of it
    #put the data in bins
    normal_nt=np.linspace(1,32,32)
    normal_th=np.linspace(1,32,32)
    for i in range(0,31):
        foo = np.where(index == i)
        #print foo
        if foo != -1:
            #print np.shape(eff_area),foo, np.sum(eff_area[foo])
            normal_nt[i] = np.sum(ntnt[foo])#sum the values of test at those points and assign it to output vector
            normal_th[i] = np.sum(thth[foo])
            #print total[i]
    normal_nt = normal_nt*factor
    normal_th = normal_th*factor
    nbincenters = 0.5*(nbinEdges[1:]+nbinEdges[:-1])
    nmenStd_nt     = 1/np.sqrt(normal_nt[0:31])
    nmenStd_th     = 1/np.sqrt(normal_th[0:31])
    nwidth      = (nbinEdges[1:]-nbinEdges[:-1])
    infs = np.where(nmenStd_th == float("inf")) and np.where(nmenStd_th > 10**2)
    nmenStd_th[infs] = 0
    infs = np.where(nmenStd_nt == float("inf")) and np.where(nmenStd_nt > 10**2)
    nmenStd_nt[infs] = 0
    #print nmenStd_th

    xlim=[0,80]
    ylim=[10,10**6]
    return xlim,ylim,binEdges,total_nt,total_th,width,menStd_nt,menStd_th,nbinEdges,normal_nt,normal_th,nwidth,nmenStd_nt,nmenStd_th,ntnt,thth,factor, title

def plot_binned_counts(inputs,nfigures=2):
    xlim=inputs[0]
    ylim=inputs[1]
    binEdges=inputs[2]
    total_nt=inputs[3]
    total_th=inputs[4]
    width=inputs[5]
    menStd_nt=inputs[6]
    menStd_th=inputs[7]
    nbinEdges=inputs[8]
    normal_nt=inputs[9]
    normal_th=inputs[10]
    nwidth=inputs[11]
    nmenStd_nt=inputs[12]
    nmenStd_th=inputs[13]
    ntnt=inputs[14]
    thth=inputs[15]
    factor=inputs[16]
    title=inputs[17]
    
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    if nfigures !=2:
        ax1 = fig.add_subplot(131)
        
    ax1.step(binEdges[0:31], total_nt[0:31])#, width=width, color='r',yerr=menStd_nt)#, alpha=0.6)
    ax1.step(binEdges[0:31], total_th[0:31])#, bottom=total_nt[0:31], width=width, color='g',yerr=menStd_th)#, alpha=0.6)
    plt.suptitle(title)
    plt.title('Log binning')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts $s^{-1} keV^{-1}$')
    #ax1.xscale('log')
    ax1.set_xlim(xlim)
    plt.yscale('log', nonposy='clip')
    ax1.set_ylim(ylim)

    ax2 = fig.add_subplot(122)
    if nfigures !=2:
        ax2 = fig.add_subplot(132)

    #ax2.hist() #histogram with normal binning
    ax2.bar(nbinEdges[0:31], normal_nt[0:31], width=nwidth, color='r',yerr=nmenStd_nt,alpha=0.6)
    ax2.bar(nbinEdges[0:31], normal_th[0:31], bottom=normal_nt[0:31],width=nwidth, color='g',yerr=nmenStd_th,alpha=0.6)
    ax2.set_xlim(xlim)
    plt.yscale('log', nonposy='clip')
    ax2.set_ylim(ylim)
    plt.title('Linear binning')
    plt.xlabel('Energy (keV)')
    #plt.ylabel('Counts')

    if nfigures !=2:
        ax3 = fig.add_subplot(133)
        ax3.plot(energy, ntnt*factor, color='r')
        ax3.plot(energy, thth*factor, color='g')
        plt.title('Input spectrum')
        plt.xlabel('Energy (keV)')
        plt.ylabel('Counts')
        #ax1.xscale('log')
        ax3.set_xlim(xlim)
        plt.yscale('log', nonposy='clip')
        ax3.set_ylim([1,100000])
   
    fig.show()
    print 'non-thermal to thermal ratio: ', (np.sum(ntnt)/np.sum(thth))#/np.sum(ntnt+thth))*100
    
def plot_binned_counts1(inputs):
    xlim=inputs[0]
    ylim=inputs[1]
    binEdges=inputs[2]
    total_nt=inputs[3]
    total_th=inputs[4]
    print np.max(total_nt)
    print np.max(total_th)
    width=inputs[5]
    menStd_nt=inputs[6]
    menStd_th=inputs[7]
    nbinEdges=inputs[8]
    normal_nt=inputs[9]
    normal_th=inputs[10]
    nwidth=inputs[11]
    nmenStd_nt=inputs[12]
    nmenStd_th=inputs[13]
    ntnt=inputs[14]
    thth=inputs[15]
    factor=inputs[16]
    title=inputs[17]
    bincenters=[(binEdges[i+1]-binEdges[i])/2.+binEdges[i] for i in range(0,len(binEdges)-1)]
    fig,ax1 = plt.subplots()
        
    ax1.step(binEdges[0:31], total_nt[0:31],label='nonthermal',linewidth=3,color='b')#, width=width, color='r',yerr=menStd_nt)#, alpha=0.6)
    #ax1.errorbar(bincenters, total_nt[1:32],yerr=menStd_nt,fmt='.',color='b',linewidth=2)
    ax1.step(binEdges[0:31], total_th[0:31],label='thermal',linewidth=3,color='k')#, bottom=total_nt[0:31], width=width, color='g',yerr=menStd_th)#, alpha=0.6)
    ax1.step(binEdges[0:31], total_th[0:31]+total_nt[0:31],label='total',linewidth=2,color='r')#, bottom=total_nt[0:31], width=width, color='g',yerr=menStd_th)#, alpha=0.6)
    ax1.errorbar(bincenters[-10:], total_nt[22:32],yerr=menStd_nt[-10:],fmt='.',color='r',linewidth=2)
    #plt.suptitle(title)
    #plt.title('Log binning',fontsize='large')
    plt.xlabel('Energy (keV)',fontsize='large')
    plt.ylabel('Counts $s^{-1} keV^{-1}$',fontsize='large')
    plt.legend(fontsize='large')
    #ax1.xscale('log')
    ax1.set_xlim(xlim)
    plt.yscale('log', nonposy='clip')
    ax1.set_ylim([1,100000])

    fig.show()

    return inputs
#if __name__ == "__main__":
#    inputs=bin_flare_counts()
#    plot_binned_counts1(inputs)
