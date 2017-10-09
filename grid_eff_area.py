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
#import q
import pickle

#@q

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

def interpolate_transmission(slit_input, substrate_input, attenuator=False, filled=False, detector=False):
    '''Interpolate everything on to one energy range'''
    slit_thickness = slit_input[1]
    slit_material = slit_input[0]
    substrate_thickness = substrate_input[1]
    substrate_material = substrate_input[0]
    #print slit_thickness, slit_material,substrate_thickness,substrate_material
    
    slit_properties = slits(slit_thickness, slit_material)
    substrate_properties = substrate(substrate_thickness, substrate_material)

    #first find the limits of the energies
    emin = np.min([np.min(slit_properties['energy']),np.min(substrate_properties['energy'])])
    emax = np.max([np.max(slit_properties['energy']),np.max(substrate_properties['energy'])])
    print emin,emax
    evector = np.linspace(emin, 1000, num=500, endpoint=True)
    #print np.min(evector),np.max(evector)
    slit_transmission_interp = np.interp(evector,slit_properties['energy'],slit_properties['transmission'])
    substrate_transmission_interp = np.interp(evector,substrate_properties['energy'],substrate_properties['transmission'])

    #are the slits filled with polymer? If so, calculate the transmission through the material. This material covers 50% of the area of the front part of grid (layer on top of the substrate). Think of this layer as two separate layers, each receiving 50% of the incoming light. Then, the transmission through the total layer is: .5*T_1 * .5*T_2. So the effective area will actually increase fractionally if the slits are filled with a material. Does that make sense though? I would expect a small decrease. Think about this on Monady.
    if filled:
        filling_dict = read_data('/Users/wheatley/Documents/Solar/MiSolFA/calculations/data/Polymer.csv')
        filling = np.interp(evector,filling_dict['E (keV)'],filling_dict['P(xi) (d='+slit_thickness +' um)'])
        filling = filling*0.5 #because it has half the effective area
    else:
        filling = 1.0

    #now for the detector
    if detector:
        detector_dict = read_data('/Users/wheatley/Documents/Solar/MiSolFA/calculations/data/CdTe.csv')
        key='P(xi) (d='+detector[0]+' um)'
        detector = 1-np.interp(evector,detector_dict['E (keV)'],detector_dict[key]) #1mm CdTe - detector response is # of photons STOPPED by the detector, ie, 1-T_dectector
    else:
        detector_dict = read_data('/Users/wheatley/Documents/Solar/MiSolFA/calculations/data/CdTe.csv')
        detector = 1-np.interp(evector,detector_dict['E (keV)'],detector_dict['P(xi) (d=1000 um)']) #1mm CdTe - detector response is # of photons STOPPED by the detector, ie, 1-T_dectector

    #and also the attenuator
    if attenuator:
        data_dict = read_data('/Users/wheatley/Documents/Solar/MiSolFA/calculations/data/'+attenuator[0]+'.csv')
        ff = data_dict['P(xi) (d='+attenuator[1] +' um)']
        fac = np.interp(evector,data_dict['E (keV)'],ff) * detector**2
    else:
        fac = 1.*detector**2 #since there is a front and rear grid, we have to take into account those effects ^2
    
    return substrate_properties, slit_properties, fac, substrate_transmission_interp, slit_transmission_interp, evector, filling
        

def grid_eff_area(slit_input, substrate_input, distribution='uniform',attenuator=False, filled=False, detector=False):
    '''Calculate effective area for the entire grid. Return that plus whatever properties are important for labeling the plot. Inputs are arrays [material, thickness] for slits and substrate'''
    griddata=interpolate_transmission(slit_input, substrate_input,attenuator=attenuator, filled=filled,detector=detector)
    substrate_properties = griddata[0]
    slit_properties = griddata[1]
    fac = griddata[2]
    substrate_transmission_interp = griddata[3]
    slit_transmission_interp = griddata[4]
    evector = griddata[5]
    filling = griddata[6]
    
    eff_area = (substrate_properties['substrate_area']*substrate_transmission_interp - (slit_properties['percent_area']*slit_transmission_interp * filling))*fac #area of substrate * transmission % + area of slits * transmission %? can this be greater than 1?
    sub=substrate_properties['substrate_area']*substrate_transmission_interp*fac
    slit=(slit_properties['percent_area']*slit_transmission_interp * filling)*fac
    
    if distribution != 'uniform' and distribution != 'flare':
         factor =  input_distribution(evector,distribution)
         eff_area =(substrate_properties['substrate_area']*substrate_transmission_interp - (slit_properties['percent_area']*slit_transmission_interp*filling)) * factor*fac
         sub=substrate_properties['substrate_area']*substrate_transmission_interp*factor*fac
         slit=(slit_properties['percent_area']*slit_transmission_interp * filling)*factor*fac

    plotdata = evector,eff_area,sub,slit,substrate_properties,slit_properties,fac,attenuator

    if distribution == 'flare':
         factor,ntnt,thth =  input_distribution(evector,distribution)
         eff_area = (substrate_properties['substrate_area']*substrate_transmission_interp - (slit_properties['percent_area']*slit_transmission_interp*filling)) * factor*fac
         sub=substrate_properties['substrate_area']*substrate_transmission_interp*factor*fac
         slit=(slit_properties['percent_area']*slit_transmission_interp*filling)*factor*fac
         plotdata = evector,eff_area,sub,slit,substrate_properties,slit_properties,ntnt,thth,fac, attenuator
         
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
    atype=plotdata[7]
     
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax1.ylog(energy, eff_area, color="b", label="Planck")
    ax1.plot(energy, eff_area, color="r",label="total",linewidth='2')
    ax1.plot(energy, sub, color="g",label="substrate",linewidth='2')
    ax1.plot(energy, slit, color="b",label="slits",linewidth='2')
    if atype:
        ax1.plot(energy, attenuator, color="m",label="attenuator+\ndetector (1mm CdTe)",linewidth='2')
    else:
        atype= ['0','0']
    #ax1.plot(energy, plotdata[7], color="c",label="detector",linewidth='2')
    
    
    plt.xlabel('Energy (keV)')
    plt.ylabel('Effective area (cm$^2$)')
    ax1.set_ylim([0,1])
    ax1.set_xlim([0,150])
    plt.title("Effective area of grids "+ substrate_properties['material'] + ' '+ substrate_properties['thickness'] + '$\mu$m, ' + slit_properties['material'] + ' '+ slit_properties['thickness']+'$\mu$m with ' + atype[0] + ' '+ atype[1]+'$\mu$m attenuator')
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
    atype=plotdata[9]

    slit_transmission_interp = np.interp(energy,slit_properties['energy'],slit_properties['transmission'])
    substrate_transmission_interp = np.interp(energy,substrate_properties['energy'],substrate_properties['transmission'])

    thermal_counts=((substrate_properties['substrate_area']*substrate_transmission_interp - slit_properties['percent_area']*slit_transmission_interp)*thth*fac) #thermal counts affected by the grid
    nonthermal_counts=((substrate_properties['substrate_area']*substrate_transmission_interp - slit_properties['percent_area']*slit_transmission_interp)*ntnt*fac) #nonthermal counts affected by the grid
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)    #ax1.ylog(energy, eff_area, color="b", label="Planck")
    ax1.semilogy(energy, eff_area, color="r",label="total",linewidth='2')
    ax1.semilogy(energy, thermal_counts, color="g",label="thermal")
    ax1.semilogy(energy, nonthermal_counts, color="b",label="non-thermal")
   # print np.mean(eff_area),np.min(eff_area),np.max(eff_area)
    
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts $s^{-1} keV^{-1}$')
    ax1.set_ylim([1,10000])
    ax1.set_xlim([0,150])
    plt.title("Expected flare counts for grids "+ substrate_properties['material'] + ' '+ substrate_properties['thickness'] + '$\mu$m, ' + slit_properties['material'] + ' '+ slit_properties['thickness']+'$\mu$m'# with ' + atype[0] + ' '+ atype[1]+'$\mu$m attenuator')

    #plt.title("Effective area of grids")
    ax1.legend(loc='upper right',fontsize='medium')
    #ax1.plot()

    fig.show()

    plotdata=[energy,eff_area,sub,slit,substrate_properties,slit_properties,ntnt,thth,fac,thermal_counts,nonthermal_counts]
    #save current plot data in pickle file
    pickle.dump(plotdata, open('plotdata.p','wb'))

    return fig


#if __name__ == "__main__":
    #data=grid_eff_area(['Au','200'],['C','500'],distribution='uniform',attenuator=['Al','50'], filled=True)
    #pprint.pprint(data)
    #plot_eff_area(data)
    #data=grid_eff_area(['Au','200'],['C','500'],distribution='flare',attenuator=['Al','50'], filled=True)
    #plot_flare_counts(data)
    #data=grid_eff_area(['Au','200'],['Si','500'],distribution='uniform',attenuator=['Al','50'], filled=False)
    #pprint.pprint(data)
    #plot_eff_area(data)

