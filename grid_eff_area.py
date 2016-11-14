 #######################################
# grid_eff_area.py
# Erica Lastufka 24/10/2016  

#Description: Calculate effective area as a function of energy for MiSolFA (and other) grids
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
import q
import pickle

@q

def input_distribution(energy,distribution): 
    '''Returns uniform (or other) distribution of photons as a function of energy'''
    if distribution == 'Planck':
        from transmission import Planck
        Nphotons = Planck(energy)
        prob = Nphotons/max(Nphotons)
    if distribution == 'MB':
        from transmission import MaxwellBoltzmann
        prob = MaxwellBoltzmann(energy)/np.max(MaxwellBoltzmann(energy))
    if distribution == 'flare':
        import flare_XR_dist as fd
        dist = fd.flare_XR_dist() #[energy, thermal part, non-thermal part]
        prob = np.interp(energy,dist[0],dist[1])
        #print np.shape(dist[0]),np.shape(dist[2]),np.shape(dist[3])
        ntnt = np.interp(energy,dist[0],dist[2])
        thth = np.interp(energy,dist[0],dist[3])

        return prob, ntnt, thth

def slits(thickness, material='vaccuum'):
    '''Define geometry, thickness, material of slits, return them as a Python dictionary'''
    transmission = 1.
    elementf = '/Users/wheatley/Documents/Solar/MiSolFA/calculations/data/' + material + '.csv'
    data_dict = read_data(elementf)
    slit_properties={"material":material, "thickness":thickness, "percent_area":0.5, "transmission":data_dict['P(xi) (d='+thickness+' um)'], "energy":data_dict['E (keV)']}
    return slit_properties

def substrate(thickness, material):
    '''Define geometry, thickness, material, return them as a Python dictionary'''
    substrate_area= 1. #cm^2
    elementf = '/Users/wheatley/Documents/Solar/MiSolFA/calculations/data/' + material + '.csv'
    data_dict = read_data(elementf)    
    substrate_properties={"material":material, "thickness":thickness,"substrate_area":substrate_area,"transmission":data_dict['P(xi) (d='+thickness+' um)'],"energy":data_dict['E (keV)']}
    return substrate_properties

def grid_eff_area(slit_input, substrate_input, distribution='uniform',attenuator=False):
    '''Calculate effective area for the entire grid. Return that plus whatever properties are important for labeling the plot. Inputs are arrays [material, thickness] for slits and substrate'''
    slit_thickness = slit_input[1]
    slit_material = slit_input[0]
    substrate_thickness = substrate_input[1]
    substrate_material = substrate_input[0]
    #print slit_thickness, slit_material,substrate_thickness,substrate_material
    
    slit_properties = slits(slit_thickness, slit_material)
    substrate_properties = substrate(substrate_thickness, substrate_material)

    #need to interpolate everything to the same grid. Let's say 100 points for fun

    #first find the limits of the energies
    emin = np.min([np.min(slit_properties['energy']),np.min(substrate_properties['energy'])])
    emax = np.max([np.max(slit_properties['energy']),np.max(substrate_properties['energy'])])
    print emin,emax
    evector = np.linspace(emin, 1000, num=500, endpoint=True)
    #print np.min(evector),np.max(evector)
    slit_transmission_interp = np.interp(evector,slit_properties['energy'],slit_properties['transmission'])
    substrate_transmission_interp = np.interp(evector,substrate_properties['energy'],substrate_properties['transmission'])

    detector_dict = read_data('/Users/wheatley/Documents/Solar/MiSolFA/calculations/data/CdTe.csv')
    detector = 1-np.interp(evector,detector_dict['E (keV)'],detector_dict['P(xi) (d=1000 um)']) #1mm CdTe - detector response is # of photons STOPPED by the detector, ie, 1-T_dectector
    if attenuator:
        data_dict = read_data('/Users/wheatley/Documents/Solar/MiSolFA/calculations/data/'+attenuator[0]+'.csv')
        ff = data_dict['P(xi) (d='+attenuator[1] +' um)']
        fac = np.interp(evector,data_dict['E (keV)'],ff) * detector
    else:
        fac = 1.*detector

    eff_area = (substrate_properties['substrate_area']*substrate_transmission_interp - slit_properties['percent_area']*slit_transmission_interp)*fac #area of substrate * transmission % + area of slits * transmission %? can this be greater than 1?
    sub=substrate_properties['substrate_area']*substrate_transmission_interp*fac
    slit=slit_properties['percent_area']*slit_transmission_interp*fac
    
    if distribution != 'uniform' and distribution != 'flare':
         factor =  input_distribution(evector,distribution)
         eff_area = (substrate_properties['substrate_area']*substrate_transmission_interp - slit_properties['percent_area']*slit_transmission_interp) * factor*fac
         sub=substrate_properties['substrate_area']*substrate_transmission_interp*factor*fac
         slit=slit_properties['percent_area']*slit_transmission_interp*factor*fac

    plotdata = evector,eff_area,sub,slit,substrate_properties,slit_properties,fac

    if distribution == 'flare':
         factor,ntnt,thth =  input_distribution(evector,distribution)
         eff_area = (substrate_properties['substrate_area']*substrate_transmission_interp - slit_properties['percent_area']*slit_transmission_interp) * factor*fac
         sub=substrate_properties['substrate_area']*substrate_transmission_interp*factor*fac
         slit=slit_properties['percent_area']*slit_transmission_interp*factor*fac
         plotdata = evector,eff_area,sub,slit,substrate_properties,slit_properties,ntnt,thth,fac
         
    return plotdata

def plot_eff_area(plotdata):
    '''Make the plot'''
    energy = plotdata[0]
    eff_area = plotdata[1]
    sub = plotdata[2]
    slit = plotdata[3]
    substrate_properties = plotdata[4]
    slit_properties = plotdata[5]
    attenuator = plotdata[6]
     
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax1.ylog(energy, eff_area, color="b", label="Planck")
    ax1.plot(energy, eff_area, color="r",label="total",linewidth='2')
    ax1.plot(energy, sub, color="g",label="substrate",linewidth='2')
    ax1.plot(energy, slit, color="b",label="slits",linewidth='2')
    ax1.plot(energy, attenuator, color="m",label="attenuator+detector (1mm CdTe)",linewidth='2')
    #ax1.plot(energy, plotdata[7], color="c",label="detector",linewidth='2')
    
    
    plt.xlabel('Energy (keV)')
    plt.ylabel('Effective area (cm$^2$)')
    ax1.set_ylim([0,1.5])
    ax1.set_xlim([0,150])
    plt.title("Effective area of grids "+ substrate_properties['material'] + ' '+ substrate_properties['thickness'] + '$\mu$m, ' + slit_properties['material'] + ' '+ slit_properties['thickness']+'$\mu$m' )
    ax1.legend(loc='upper right',fontsize='medium')

    ax1.plot()

    fig.show()
    return fig

def plot_flare_counts(plotdata):
    '''Make the plot for if you have a small flare'''
    #print np.shape(plotdata)
    energy = plotdata[0]
    eff_area = plotdata[1]
    sub = plotdata[2]
    slit = plotdata[3]
    substrate_properties = plotdata[4]
    slit_properties = plotdata[5]
    ntnt = plotdata[6]
    thth = plotdata[7]
    fac = plotdata[8]

    slit_transmission_interp = np.interp(energy,slit_properties['energy'],slit_properties['transmission'])
    substrate_transmission_interp = np.interp(energy,substrate_properties['energy'],substrate_properties['transmission'])
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax1.ylog(energy, eff_area, color="b", label="Planck")
    ax1.loglog(energy, eff_area, color="r",label="total",linewidth='2')
    ax1.loglog(energy, (substrate_properties['substrate_area']*substrate_transmission_interp - slit_properties['percent_area']*slit_transmission_interp)*thth*fac, color="g",label="thermal")
    ax1.loglog(energy, (substrate_properties['substrate_area']*substrate_transmission_interp - slit_properties['percent_area']*slit_transmission_interp)*ntnt*fac, color="b",label="non-thermal")
   # print np.mean(eff_area),np.min(eff_area),np.max(eff_area)
    
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts $s^{-1} keV^{-1}$')
    ax1.set_ylim([10e-4,10e4])
    ax1.set_xlim([1,1000])
    plt.title("Flare counts, "+ substrate_properties['material'] + ' '+ substrate_properties['thickness'] + '$\mu$m, ' + slit_properties['material'] + ' '+ slit_properties['thickness']+'$\mu$m' )

    #plt.title("Effective area of grids")
    ax1.legend(loc='upper right',fontsize='medium')

    ax1.plot()

    fig.show()

    #save current plot data in pickle file
    pickle.dump(plotdata, open('plotdata.p','wb'))

    return fig


if __name__ == "__main__":
    data=grid_eff_area(['Au','200'],['C','200'],distribution='uniform',attenuator=['Al','100'])
    #pprint.pprint(data)
    plot_eff_area(data)
    data=grid_eff_area(['Au','200'],['C','200'],distribution='flare',attenuator=['Al','100'])
    plot_flare_counts(data)
