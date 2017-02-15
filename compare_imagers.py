 #######################################
# compare_imagers.py
# Erica Lastufka 7/2/2107

#Description: Compare behaviour of various possible configurations for the imager
#######################################

#######################################
# Usage:

######################################

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import pickle
import itertools
import imager

def PSI_Imager(picklename):
    if not picklename:
        im = imager(attenuator=[],slits=[],substrate=[],filled=[])
    else:
        im = pickle.load(open(picklename),'rb')
    return im

def uworks_Imager(picklename):
    if not picklename:
        im = imager(attenuator=[],slits=[],substrate=[],filled=[])
    else:
        im = pickle.load(open(picklename),'rb')
    return im

def compare_attenuators():
    one=Imager(attenuator=['W','300'])
    two=Imager(attenuator=['Steel','300'])
    three=Imager(attenuator=['Ti','300'])

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    
    ax1.plot(one.slits['Energy'],one.eff_area, color='r',linewidth='2')
    ax1.plot(two.slits['Energy'],two.eff_area, color='g',linewidth='2')
    ax1.plot(three.slits['Energy'],three.eff_area, color='b',linewidth='2')
    ax1.set_xlabel('Energy (keV)')
    ax1.set_ylabel('Effective area (cm$^2$)')
    ax1.set_title('Effective area of grids with ' +'300$\mu$m attenuator')
    
    ax2.plot(one.slits['Energy'],one.attenuator['Transmission'], color='r',label='W',linewidth='2')
    ax2.plot(two.slits['Energy'],two.attenuator['Transmission'], color='g',label='Steel',linewidth='2') #these 2 look weird! why do they go above .5? Is that right?
    ax2.plot(three.slits['Energy'],three.attenuator['Transmission'], color='b',label='Ti',linewidth='2')
    ax2.set_xlabel('Energy (keV)')
    ax2.set_ylabel('Transmission')
    ax2.set_title('Tranmission through 300$\mu$m attenuators')
    ax2.legend(loc='lower right',fontsize='medium')
    ax2.plot()
    fig.show()

def find_closest_match(first_element='W',elements=['Ti','Steel'],first_thickness=['100','200','500','1000'],match_thickness=['600','700','800','900','1000','1500','2000','2500','3000','3500','4000','4500','5000','6000','7000','8000','9000','10000']):
    import pandas as pd
    os.chdir('../../calculations/data')
    data0=pd.read_csv(first_element+'.csv',sep=',',header=0)
    filenames,arr=[],[]
    for element in elements: filenames.append(element + '.csv')
    for file in filenames: arr.append([])
    for i,file in enumerate(filenames):
        arr[i] = pd.read_csv(file,sep=',',header=0)

    #interoplate
    evector = np.linspace(2.08, 433, num=500, endpoint=True)
    data1= [[] for i in range(len(first_thickness))]
    #arr1= [[] for i in range(len(match_thickness))]
    #arr2= [[] for i in range(len(match_thickness))] #because I'm lazy
    closest_match1,closest_match2,match1,match2 = [],[],[],[]
    for i,thick in enumerate(first_thickness):    
        data1[i]=np.interp(evector,data0['E (keV)'],data0['P(xi) (d=' + thick + ' um)'])
        print np.min(data1[i]),np.max(data1[i])
        arr1= [[] for k in range(len(match_thickness))]
        arr2= [[] for k in range(len(match_thickness))]
        diff1 = [[] for k in range(len(match_thickness))]
        diff2 = [[] for k in range(len(match_thickness))]        
        for j,t in enumerate(match_thickness):
            arr1[j]=np.interp(evector,arr[0]['E (keV)'],arr[0]['P(xi) (d=' + t + ' um)'])   
            arr2[j]=np.interp(evector,arr[1]['E (keV)'],arr[1]['P(xi) (d=' + t + ' um)'])
            diff1[j] = arr1[j][0:56] - data1[i][0:56] # E=50 keV at index 56
            diff2[j] = arr2[j][0:56] - data1[i][0:56]
            #print np.min(arr1[j]),np.max(arr1[j]),np.min(diff1[j]),np.max(diff1[j])
        closest_match1.append(np.where(np.mean(np.abs(diff1), axis=1) == np.min(np.mean(np.abs(diff1), axis=1)))) #need to get the index somehow - don't think this is doing it right however
        closest_match2.append(np.where(np.mean(np.abs(diff2), axis=1) == np.min(np.mean(np.abs(diff2), axis=1))))
        match1.append(match_thickness[closest_match1[i][0][0]])
        match2.append(match_thickness[closest_match2[i][0][0]])
    #print np.shape(match1),np.shape(match2)
    fig, axes = plt.subplots(2, 2, figsize=(10, 12))
    #plt.subplots_adjust(wspace=0,hspace=0)
    
    for ax,thick,i in zip(axes.flat,first_thickness,range(0,len(first_thickness))):
        ax.semilogy(data0['E (keV)'],data0['P(xi) (d=' + thick + ' um)'],color='k',label='W ' +thick + 'um', linewidth='2')
        ax.semilogy(arr[0]['E (keV)'],arr[0]['P(xi) (d=' + match1[i] + ' um)'],color='r',label='Ti '+ match1[i] + 'um',linewidth='2') #plot the closest matches
        ax.semilogy(arr[1]['E (keV)'],arr[1]['P(xi) (d=' + match2[i] + ' um)'],color='g',label='Steel '+ match2[i] + 'um',linewidth='2')
        ax.set_xlabel('E (keV)')
        ax.set_ylabel('Transmission')
        ax.set_xlim([0,100])
        ax.set_ylim([10**-10,10**0])
        ax.legend(loc='lower right')
    plt.suptitle('Materials and thicknesses with closest transmission various thicknesses of W',size=16)
    plt.show()
