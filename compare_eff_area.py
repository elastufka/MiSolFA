 #######################################
# compare_eff_area.py
# Erica Lastufka 30/11/2016  

#Description: compare effective area of different grid configurations
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
import grid_eff_area as ga
import pickle

def compare_attenuators():
    a1 = ['Al','50']
    a2 = ['Al','100']
    a3 = ['Be','300']
    
    data1=ga.grid_eff_area(['Au','200'],['C','500'],distribution='uniform',attenuator=a1, filled=True)
    data2=ga.grid_eff_area(['Au','200'],['C','500'],distribution='uniform',attenuator=a2, filled=True)
    data3=ga.grid_eff_area(['Au','200'],['C','500'],distribution='uniform',attenuator=a3, filled=True)

    energy1 = data1[0]
    eff_area1 = data1[1]
    attenuator1 = data1[6]

    energy2 = data2[0]
    eff_area2 = data2[1]
    attenuator2 = data2[6]

    energy3 = data3[0]
    eff_area3 = data3[1]
    attenuator3 = data3[6]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(energy1, eff_area1, color="r",label="Al 50 $\mu$m",linewidth='2')
    ax1.plot(energy2, eff_area2, color="g",label="Al 100 $\mu$m",linewidth='2')
    ax1.plot(energy3, eff_area3, color="b",label="Be 300 $\mu$m",linewidth='2')
    ax1.plot(energy1, attenuator1, color="r",linewidth='2', linestyle='dashed')
    ax1.plot(energy2, attenuator2, color="g",linewidth='2', linestyle='dashed')
    ax1.plot(energy3, attenuator3, color="b",linewidth='2', linestyle='dashed')
    
    plt.xlabel('Energy (keV)')
    plt.ylabel('Effective area (cm$^2$)')
    ax1.set_ylim([0,1])
    ax1.set_xlim([0,50])
    plt.title("Effective area of grids C 500 $\mu$m, Au 200$\mu$m" )
    ax1.legend(loc='lower right',fontsize='medium')

    fig.show()

def compare_attenuators_flare():
    a1 = ['Al','50']
    a2 = ['Al','100']
    a3 = ['Be','300']
    
    data1=ga.grid_eff_area(['Au','200'],['C','500'],distribution='flare',attenuator=a1, filled=True)
    data2=ga.grid_eff_area(['Au','200'],['C','500'],distribution='flare',attenuator=a2, filled=True)
    data3=ga.grid_eff_area(['Au','200'],['C','500'],distribution='flare',attenuator=a3, filled=True)

    energy1 = data1[0]
    eff_area1 = data1[1]
    attenuator1 = data1[6]

    energy2 = data2[0]
    eff_area2 = data2[1]
    attenuator2 = data2[6]

    energy3 = data3[0]
    eff_area3 = data3[1]
    attenuator3 = data3[6]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.semilogy(energy1, eff_area1, color="r",label="Al 50 $\mu$m",linewidth='2')
    ax1.semilogy(energy2, eff_area2, color="g",label="Al 100 $\mu$m",linewidth='2')
    ax1.semilogy(energy3, eff_area3, color="b",label="Be 300 $\mu$m",linewidth='2')
    #ax1.semilogy(energy1, attenuator1, color="r",linewidth='2', linestyle='dashed')
    #ax1.semilogy(energy2, attenuator2, color="g",linewidth='2', linestyle='dashed')
    #ax1.semilogy(energy3, attenuator3, color="b",linewidth='2', linestyle='dashed')
    
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts $s^{-1} keV^{-1}$')
    ax1.set_ylim([1,10000])
    ax1.set_xlim([0,50])
    plt.title("Counts with grids C 500 $\mu$m, Au 200$\mu$m" )
    ax1.legend(loc='lower right',fontsize='medium')

    fig.show()

def compare_detectors():
    d1 = ['1000']
    d2 = ['2000']    
    d3 = ['3000']
    
    data1=ga.grid_eff_area(['Au','200'],['C','500'],distribution='uniform',attenuator=['Be','300'], filled=True, detector=d1)
    data2=ga.grid_eff_area(['Au','200'],['C','500'],distribution='uniform',attenuator=['Be','300'], filled=True, detector=d2)
    data3=ga.grid_eff_area(['Au','200'],['C','500'],distribution='uniform',attenuator=['Be','300'], filled=True, detector=d3)

    energy1 = data1[0]
    eff_area1 = data1[1]
    detector1 = data1[6]

    energy2 = data2[0]
    eff_area2 = data2[1]
    detector2 = data2[6]

    energy3 = data3[0]
    eff_area3 = data3[1]
    detector3 = data3[6]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(energy1, eff_area1, color="r",label="CdTe 1 mm",linewidth='2')
    ax1.plot(energy2, eff_area2, color="g",label="CdTe 2 mm",linewidth='2')
    ax1.plot(energy3, eff_area3, color="b",label="CdTe 3 mm",linewidth='2')
    ax1.plot(energy1, detector1, color="r",linewidth='2', linestyle='dashed')
    ax1.plot(energy2, detector2, color="g",linewidth='2', linestyle='dashed')
    ax1.plot(energy3, detector3, color="b",linewidth='2', linestyle='dashed')
    
    plt.xlabel('Energy (keV)')
    plt.ylabel('Effective area (cm$^2$)')
    ax1.set_ylim([0,1])
    ax1.set_xlim([50,150])
    plt.title("Effective area of grids C 500 $\mu$m, Au 200$\mu$m" )
    ax1.legend(loc='lower right',fontsize='medium')

    fig.show()

