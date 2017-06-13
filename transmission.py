#######################################
# transmission.py
# Erica Lastufka 12/10/2016  

#Description: Calculate transmission for different materials at different energies
#######################################

#######################################
# Usage:

# for default output: python transmission.py
######################################

import glob
import os
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
#import Image
import csv
import scipy.constants as sc

def read_data(elementf): 
    '''Reads data in from csv files and puts it in a Python dictionary'''
    #dir='/home/erica/Documents/Solar/MiSolFA/calculations/data'
    #os.chdir(dir)
    filename=elementf

    with open(filename, 'rb') as f:
        reader = csv.reader(f)
        header=next(reader) #get headers
    #print header
    data=genfromtxt(filename,delimiter=',') #column 0 will be NaN because it's text
    datadim = np.shape(data) #get dimensions

    data_dict={header[i]:data[1:datadim[0],i] for i in range(0,len(header))}
    #data_dict = {header[0]:data[1:datadim[0],0], header[1]:data[1:datadim[0],1],header[2]:data[1:datadim[0],2],header[3]:data[1:datadim[0],3],header[4]:data[1:datadim[0],4],header[5]:data[1:datadim[0],5],header[6]:data[1:datadim[0],6]}
    return data_dict

def Wien_approx(energy):
    '''Use Wien approximation to calculate solar spectrum. See notes in wiki.'''
    kev_to_joules = 10E-3*sc.e
    T_sun = 5700.
    energy_joules= energy*kev_to_joules
    Nphotons = (2*(energy_joules/(sc.h*sc.c))**2)*np.exp(-1.*energy_joules/(sc.k*T_sun))
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.loglog(energy, Nphotons, marker="s")
    #fig.show()

    return Nphotons

def Planck(energy):
    '''Use Planck function to calculate solar spectrum. See notes in wiki.'''
    kev_to_joules = 10E-3*sc.e
    #T_sun = 5700.
    T_sun=20.*10E6 #for a flare
    energy_joules= energy*kev_to_joules
    Nphotons = (2*(energy_joules/(sc.h*sc.c))**2)*(1/(np.exp(energy_joules/(sc.k*T_sun)) -1))

    return Nphotons

def MaxwellBoltzmann(energy):
    '''Use Maxwell Boltzmann distribution to calculate solar spectrum. See notes in wiki.'''
    kev_to_joules = 10E-3*sc.e
    #T_sun = 5700.
    T_sun=20.*10E6 #for a flare
    energy_joules= energy*kev_to_joules
    prob = (2*(energy_joules/sc.pi)**(.5))*(1/(sc.k*T_sun))**(3/2)*(np.exp(-1*energy_joules/(sc.k*T_sun)))
    #Nphotons = (1/(sc.k*T_sun))*(np.exp(-1*energy_joules/(sc.k*T_sun)))
    return prob

def Planckfig(Nphotons,energy):
    '''Make a figure of the Planck function '''
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.loglog(energy, Nphotons*energy, marker="o")
    plt.xlabel('Energy (keV)')
    plt.ylabel('Intensity (W/m^2 per photon energy)')
    ax1.set_ylim([0,10e15])
    ax1.set_xlim([1,100])
    plt.title("Planck distribution ")
    fig.show()

    return fig

def DistributionFig(P_Np,MB_prob,energy):
    '''Make a figure of the Planck function compared to the Maxwell Boltzmann distrubution'''
    #need to normalize both distributions
    P_dist = (P_Np)/np.trapz(P_Np)
    MB_dist = MB_prob/np.trapz(MB_prob)
    print np.trapz(P_dist), np.trapz(MB_dist)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.semilogx(energy, P_dist, color="b", label="Planck")
    ax1.semilogx(energy, MB_dist, color="r",label="Maxwell-Boltzmann")
    plt.xlabel('Energy (keV)')
    plt.ylabel('Normalized Distribution')
    ax1.set_ylim([0, 0.015])
    #ax1.set_xlim([1,100000000])
    plt.title("Planck distribution vs. Maxwell Boltzmann distribution, T=20 million K")
    ax1.legend(loc='upper left',fontsize='medium')

    fig.show()

    #return fig

def calc_counts(data_dict,key): 
    '''calculate number of counts for given thickness'''
    Nphotons = Planck(data_dict['E (keV)'])
    counts = Nphotons*data_dict[key]
    return counts

def plot_data(element,data_dict,thickness='thickness'):
    '''Plot E vs T and E vs counts for given element and thickness'''
    colors=['b','g','r','c','m']
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.xlabel('Energy (keV)')
    plt.ylabel('Transmission')
    ax1.set_ylim([10e-50,10e0])
    ax1.set_xlim([0,50])

    #ax2 = plt.subplot(211)
    #plt.xlabel('Energy (keV)')
    #plt.ylabel('Counts')
    #ax2.set_ylim([1,10e20])
    #ax2.set_xlim([1,100])


    if thickness != 'thickness': #to compare different elements - now it is a list of dictionaries,not just a dictionary. So get the keys
        keys=data_dict[0].keys()
        plt.title(thickness+ ' $\mu$m thickness')
        for item in range(0,np.size(data_dict)): #find the data with the right thickness
            for key in keys:
                if thickness in key:
                    x1=data_dict[item]['E (keV)'] 
                    ax1.semilogy(x1, data_dict[item][key], marker="s", c=colors[item],label=data_dict[item]['element'])
                    #ax2.loglog(x1,calc_counts(data_dict[item], key), marker="s", c=colors[item])
                    #fig.show()
      
      
    else: 
        keys=data_dict.keys()
        x1 = data_dict['E (keV)']
        n=0
        plt.title(element)
        for key in keys:
            if key.startswith('P'):
                print key
                x1=data_dict['E (keV)']
                ax1.semilogy(x1, data_dict[key], marker="s", c=colors[n],label=key)
                #ax2.loglog(x1,calc_counts(data_dict, key),marker="s", c=colors[n])
                n=n+1
                
    ax1.legend(loc='lower right',fontsize='medium')

    return fig

def plot_cover(element,data_dict,thickness=[.1,.15,.2], log=True): #thickness in cm
    '''Plot E vs T and E vs counts for given element and thickness'''
    colors=['b','g','r','c','m']
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.xlabel('Energy (keV)')
    plt.ylabel('Transmission')
    #ax1.set_ylim([10e-10,10e0])
    ax1.set_xlim([0,80])

    plt.title('Transmission of ' + element + ' at various (total) thicknesses')
    for i,thick in enumerate(thickness):
        pxi = np.exp(data_dict['a (cm^2/g)']*data_dict['rho (g/cm^3)'][0]*thick*-1) #calculate P(xi)
        print np.min(pxi),np.max(pxi)
        x1=data_dict['E (keV)']
        if log:
            ax1.semilogy(x1, pxi, marker="s", c=colors[i],label=str(thick)+' cm')
            ax1.set_ylim([10e-10,10e0])
            loc= 'lower right'
        else:
             ax1.plot(x1, pxi, marker="s", c=colors[i],label=str(thick)+' cm')
             ax1.set_ylim([0,1])
             loc='upper left'
                 
    ax1.legend(loc=loc,fontsize='medium')
    fig.show()
    
    return fig

def plot_cover_compare(elements,densities,thickness,log=True):
    '''Plot E vs T and E vs counts for given elements and thicknesses'''
    colors=['b','g','r','c','m']
    markers=['s','o','v','*']
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.xlabel('Energy (keV)')
    plt.ylabel('Transmission')
    #ax1.set_ylim([10e-10,10e0])
    ax1.set_xlim([0,80])

    plt.title('Transmission of ' +', '.join(el for el in elements)  + ' at various (total) thicknesses')
    for j,el in enumerate(elements):
        #get the data dictionary
        data_dict=read_data(el+'.csv')
        if not'rho (g/cm^3)' in data_dict.keys():
            data_dict['rho (g/cm^3)']=[densities[j]]
        color=colors[j]
        for i,thick in enumerate(thickness):
            pxi = np.exp(data_dict['a (cm^2/g)']*data_dict['rho (g/cm^3)'][0]*thick*-1) #calculate P(xi)
            print np.min(pxi),np.max(pxi)
            x1=data_dict['E (keV)']
            if log:
                ax1.semilogy(x1, pxi, marker=markers[i], c=color,label=str(thick)+' cm '+el)
                ax1.set_ylim([10e-10,10e0])
                loc= 'lower right'
            else:
                ax1.plot(x1, pxi, marker=markers[i], c=color,label=str(thick)+' cm '+el)
                ax1.set_ylim([0,1])
                loc='upper left'
                 
    ax1.legend(loc=loc,fontsize='medium')
    fig.show()
    
    return fig


#convert figure to data
def fig2data ( fig ):
    """
    from http://www.icare.univ-lille1.fr/node/1141
    @brief Convert a Matplotlib figure to a 4D numpy array with RGBA channels and return it
    @param fig a matplotlib figure
    @return a numpy 3D array of RGBA values
    """
    # draw the renderer
    fig.canvas.draw ( )
 
    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    buf = np.fromstring( fig.canvas.tostring_argb(), dtype=np.uint8 ).reshape(w, h, 4)
    #buf.shape = ( w*5, h*5,4 )
    # canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
    buf = np.roll ( buf, 3, axis = 2 )
    return buf

#convert figure to image
def fig2image(fig):
    """
    @brief Convert a Matplotlib figure to a PIL Image in RGBA format and return it
    @param fig a matplotlib figure
    @return a Python Imaging Library ( PIL ) image
    """

    buf = fig2data ( fig )
    w, h, d = buf.shape
    return Image.frombytes( "RGBA", ( w ,h ), buf.tostring( ) )

def compare_thickness(elements=0, show=0,ymin=10**-15):
    '''Make plots for various thicknesses and output to pdf'''

    if elements is 0:
         elements=['W','Au','Au80Sn20','Si', 'C', 'Polymer']

     #assume 6 elements
    newfig, axes = plt.subplots(3, 2, figsize=(10, 12),
                         subplot_kw={'xticks': [], 'yticks': []})
    plt.subplots_adjust(wspace=0,hspace=0)

    if show:
        for element,ax in zip(elements,axes.flat):
            data_dict=read_data(element+'.csv')
            fig=plot_data(element,data_dict)
            ax.set_ylim([.000001,.1])
            fig.show()

    if not show:
         for element, ax in zip(elements,axes.flat):
             data_dict=read_data(element+'.csv')
             fig=plot_data(element,data_dict)
             image= fig2image(fig)
             ax.imshow(image,interpolation='none')

     #newfig.show()
    newfig.savefig('compare_thickness.pdf',dpi=300)

def compare_elements(elements=0, thickness=0, show=0,ymin=10**-15): #thickness is also an array
    '''Make plots for various elements and output to pdf'''
    if elements is 0:
        elements=['W','Au','Au80Sn20']
    if thickness is 0:
        thickness=['100','150','200','250','300']

    #make a new data_dict
    list_dict = []

    for element in elements:
         data_dict=read_data(element+'.csv')
         data_dict['element']=element
         list_dict.append(data_dict)
    
    newfig, axes = plt.subplots(3, 2, figsize=(10, 12),
                         subplot_kw={'xticks': [], 'yticks': []})
    plt.subplots_adjust(wspace=0,hspace=0)

    if show:
        for thick,ax in zip(thickness,axes.flat):
            fig = plot_data('foo',list_dict,thickness=thick)
            ax.set_ylim([.0000001,.1])
            fig.show()

    if not show:
        
        for thick,ax in zip(thickness,axes.flat):
            fig = plot_data('foo',list_dict,thickness=thick)
            image= fig2image(fig)
            ax.imshow(image,interpolation='none')

        planckfig=Planckfig(Planck(list_dict[0]['E (keV)']),list_dict[0]['E (keV)'])#let's plot the Planck disttribution in the final spot   
        planckimage=fig2image(planckfig)
        axes.flat[5].imshow(planckimage,interpolation='none')

        newfig.savefig('compare_elements.pdf',dpi=300)

    
#if __name__ == "__main__":

os.chdir('../../calculations/data')
#compare_thickness(['W','Ti','Steel'], show=1)
#compare_elements(['W','Ti','Steel'], ['100','150','200','250','300'],show=1)
#compare_elements()
#compare_thickness()
#data_dict=read_data('Au.csv')
#en=np.linspace(1,100000,200)

#DistributionFig(Planck(en), MaxwellBoltzmann(en), en)
