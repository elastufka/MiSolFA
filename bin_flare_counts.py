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

def plot_binned_counts():
    plotdata =pickle.load(open('plotdata.p','rb'))
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
    factor = 12*30 #12 detectors, 2 grids, 30s integration time
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

    evector, nbinEdges=np.histogram(energy, bins=np.linspace(0, 150, 32)) #these are the energy bins. Now we need to fill them with the data ....
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
    ylim=[1,10**6]
    
    fig = plt.figure()
    ax1 = fig.add_subplot(121)

    ax1.bar(binEdges[0:31], total_nt[0:31], width=width, color='r',yerr=menStd_nt, alpha=0.6)
    ax1.bar(binEdges[0:31], total_th[0:31], width=width, color='g',yerr=menStd_th, alpha=0.6)
    plt.suptitle(title)
    plt.title('Log binning')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts')
    #ax1.xscale('log')
    ax1.set_xlim(xlim)
    plt.yscale('log', nonposy='clip')
    ax1.set_ylim(ylim)

    ax2 = fig.add_subplot(122)
    #ax2.hist() #histogram with normal binning
    ax2.bar(nbinEdges[0:31], normal_nt[0:31], width=nwidth, color='r',yerr=nmenStd_nt,alpha=0.6)
    ax2.bar(nbinEdges[0:31], normal_th[0:31], width=nwidth, color='g',yerr=nmenStd_th,alpha=0.6)
    ax2.set_xlim(xlim)
    plt.yscale('log', nonposy='clip')
    ax2.set_ylim(ylim)
    plt.title('Linear binning')
    plt.xlabel('Energy (keV)')
    #plt.ylabel('Counts')

    #ax3 = fig.add_subplot(133)
    #ax3.plot(energy, ntnt, color='r')
    #ax3.plot(energy, thth, color='g')
    #plt.title('Flare counts, log binning')
    #plt.xlabel('Energy (keV)')
    #plt.ylabel('Counts')
    ##ax1.xscale('log')
    #ax3.set_xlim([0,25])
    #plt.yscale('log', nonposy='clip')
    #ax3.set_ylim([1,10**3])
   
    fig.show()

if __name__ == "__main__":
    plot_binned_counts()
