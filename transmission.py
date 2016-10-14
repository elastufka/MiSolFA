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
import Image
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


    data_dict = {header[0]:data[1:datadim[0],0], header[1]:data[1:datadim[0],1],header[2]:data[1:datadim[0],2],header[3]:data[1:datadim[0],3],header[4]:data[1:datadim[0],4],header[5]:data[1:datadim[0],5],header[6]:data[1:datadim[0],6]}
    return data_dict

def Wien_approx(energy):
    '''Use Wien approximation to calculate solar spectrum. See notes in wiki.'''
    kev_to_joules = 10E-3*sc.e
    T_sun = 5280.
    energy_joules= energy*kev_to_joules
    Nphotons = (2*(energy_joules/(sc.h*sc.c))**2)*np.exp(-1.*energy_joules/(sc.k*T_sun))
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.loglog(energy, Nphotons, marker="s")
    #fig.show()

    return Nphotons

def Planck(energy):
    '''Use Wien approximation to calculate solar spectrum. See notes in wiki.'''
    kev_to_joules = 10E-3*sc.e
    T_sun = 5280.
    energy_joules= energy*kev_to_joules
    Nphotons = (2*(energy_joules/(sc.h*sc.c))**2)*(1/(np.exp(energy_joules/(sc.k*T_sun)) -1))
    
    #fig = plt.figure()
    #ax1 = fig.add_subplot(111)
    #ax1.loglog(energy, Nphotons, marker="o")
    #fig.show()

    return Nphotons

def Planckfig(Nphotons,energy):
    '''Make a figure of the Planck function '''
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.loglog(energy, Nphotons*energy, marker="o")
    plt.xlabel('Energy (keV)')
    plt.ylabel('Intensity (W/m^2 per photon energy)')
    #ax1.set_ylim([0,10e30])
    ax1.set_xlim([1,100])
    plt.title("Planck distribution")
    fig.show()

    return fig

def calc_counts(data_dict,key): 
    '''calculate number of counts for given thickness'''
    Nphotons = Planck(data_dict['E (keV)'])
    counts = Nphotons*data_dict[key]
    return counts

def plot_data(element,data_dict,thickness='thickness'):
    '''Plot E vs T and E vs counts for given element and thickness'''
    colors=['b','g','r','c','m']
    fig = plt.figure()
    ax1 = fig.add_subplot(212)
    plt.xlabel('Energy (keV)')
    plt.ylabel('Transmission')
    ax1.set_ylim([10e-50,10e0])
    ax1.set_xlim([1,100])
    ax1.legend(loc='lower right',fontsize='medium')

    ax2 = plt.subplot(211)
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts')
    ax2.set_ylim([1,10e20])
    ax2.set_xlim([1,100])


    if thickness != 'thickness': #to compare different elements - now it is a list of dictionaries,not just a dictionary. So get the keys
        keys=data_dict[0].keys()
        plt.title(thickness+ ' um thickness')
        for item in range(0,np.size(data_dict)): #find the data with the right thickness
            for key in keys:
                if thickness in key:
                    x1=data_dict[item]['E (keV)'] 
                    ax1.loglog(x1, data_dict[item][key], marker="s", c=colors[item],label=data_dict[item]['element'])
                    ax2.loglog(x1,calc_counts(data_dict[item], key), marker="s", c=colors[item])
                    #fig.show()
      
      
    else: 
        keys=data_dict.keys()
        x1 = data_dict['E (keV)']
        n=0
        plt.title(element)
        for key in keys:
            if key.startswith('P'):
                print key
                ax1.loglog(x1, data_dict[key], marker="s", c=colors[n],label=key)
                ax2.loglog(x1,calc_counts(data_dict, key),marker="s", c=colors[n])
                n=n+1     

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

def compare_thickness(elements=0):
    '''Make plots for various thicknesses and output to pdf'''

    if elements is 0:
         elements=['W','Au','Au80Sn20','Si','Be','Al']

     #assume 6 elements
    newfig, axes = plt.subplots(3, 2, figsize=(10, 12),
                         subplot_kw={'xticks': [], 'yticks': []})
    plt.subplots_adjust(wspace=0,hspace=0)

    for element, ax in zip(elements,axes.flat):
         data_dict=read_data(element+'.csv')
         fig=plot_data(element,data_dict)
         image= fig2image(fig)
         ax.imshow(image,interpolation='none')

     #newfig.show()
    newfig.savefig('compare_thickness.pdf',dpi=300)

def compare_elements(elements=0, thickness=0): #thickness is also an array
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

    for thick,ax in zip(thickness,axes.flat):
        fig = plot_data('foo',list_dict,thickness=thick)
        image= fig2image(fig)
        ax.imshow(image,interpolation='none')

    planckfig=Planckfig(Planck(list_dict[0]['E (keV)']),list_dict[0]['E (keV)'])#let's plot the Planck disttribution in the final spot   
    planckimage=fig2image(planckfig)
    axes.flat[5].imshow(planckimage,interpolation='none')

    newfig.savefig('compare_elements.pdf',dpi=300)

#if __name__ == "__main__":

#compare_thickness(['W','Au','Au80Sn20','Si','Be','Al'])
#compare_elements(['W','Au','Au80Sn20'], ['100','150','200','250','300'])
compare_elements()
compare_thickness()
