 #######################################
# grid_eff_area.py
# Erica Lastufka 24/10/2016  

#Description: Calculate effective area as a function of energy for MiSolFA (and other) grids
#######################################

#######################################
# Usage:

# for default output: python transmission.py
######################################

import glob
import os
import numpy as np
#from numpy import genfromtxt
import matplotlib.pyplot as plt
#import Image
#import csv
import scipy.constants as sc
#from transmission import read_data #so I can read the .csv files into python dictionaries
#am I really limited to the energy range given by the data?

def input_distribution(energy): 
    '''Returns uniform (or other) distribution of photons as a function of energy'''
    Nphotons = 10e20
    counts = Nphotons*data_dict[key]
    return counts

def slits(thickness, material=material):
    '''Define geometry, thickness, material of slits, return them as a Python dictionary'''
    if not material:
        material='vacuum'
    slit_properties={"material":material, "thickness":thickness, "percent_area":0.5, "transmission":[]}
    return slit_properties

def substrate(thickness, material):
    '''Define geometry, thickness, material, return them as a Python dictionary'''
    substrate_area= 12 #cm^2
    substrate_properties={"material":material, "thickness":thickness,"substrate_area":substrate_area,"transmission":[]}
    return substrate_properties

def calc_transmission(component, energy):
    '''Calculate transmission for either the slits or the substrate'''

    return transmission

def grid_eff_area(slits, substrate, distribution,instrument=instrument):
    '''Calculate effective area for the entire grid. Return that plus whatever properties are important for labeling the plot.'''
    if not instrument:
        instrument='MiSolFA'
    slit_properties=slits(slit_thickness, material=slit_material)
    substrate_properties=substrate(substrate_thickness, substrate_material)
        
    return eff_area

def plot_eff_area(eff_area,energy):
    '''Make the plot'''
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.ylog(energy, eff_area, color="b", label="Planck")
    #ax1.loglog(energy, 10E10*MB_Np, color="r",label="Maxwell-Boltzmann")
    plt.xlabel('Energy (keV)')
    plt.ylabel('Effective area (cm^2)')
    #ax1.set_ylim([1,10e20])
    ax1.set_xlim([1,1000])
    plt.title("Effective area of " + grid parameters + "grids")
    ax1.legend(loc='lower right',fontsize='medium')

    fig.show()


    #if __name__ == "__main__":
