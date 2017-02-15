 #######################################
# imager.py
# Erica Lastufka 7/2/2107

#Description: Class for the MiSolFA imager object
#######################################

#######################################
# Usage:
#im=Imager() #for the defaults
#im=Imager(widget=True) #to input via GUI
######################################

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import pickle
import itertools

class Imager(object):
    '''Create the imager object. Read in the material transmission data and interpolate it to vectors of the same length. Calculate the transmission of each individual component and assign it to attributes.'''
    def __init__(self, widget=False, attenuator = False, slits= ['Au','200'], substrate = ['C','500'], filled = False, diffraction=False,detector=False,distribution='uniform',psi=False, uworks=False):
        if psi==True:
            try:
                self = pickle.load(open('psi_imager.p'),'rb') #can I do this? I'd need to skip to the end then
                widget=None
            except IOError:
                attenuator =['Al','100']
                slits=['Au','200']
                substrate=['Si','500']
                filled='Si'
                detector=['CdTe','1000']            

        if uworks ==True:
            try:
                self = pickle.load(open('uworks.p'),'rb') #can I do this?
                widget=None
            except IOError:
                attenuator =['Al','100']
                slits=['Au','200']
                substrate=['C','500']
                filled='Polymer'
                detector=['CdTe','1000']                

        if widget == False:
            self.attenuator = attenuator
            if attenuator:
                self.attenuator = {'Material':attenuator[0],'Thickness':attenuator[1],'Transmission':1, 'Energy':[], 'eff_area':[]}      
            self.slits = {'Material':slits[0], 'Thickness':slits[1],'Percent_area':0.5,'Transmission':[], 'Energy':[], 'eff_area':[]}
            self.substrate = {'Material':substrate[0], 'Thickness':substrate[1],'Substrate_area':1,'Transmission':[], 'Energy':[], 'eff_area':[]}
            self.filled = filled
            if filled:
                self.filled = {'Material':filled,'Thickness':slits[1],'Transmission':1, 'Energy':[], 'eff_area':[]}
            self.detector = detector
            if detector:
                self.detector = {'Material':detector[0],'Thickness':detector[1],'Transmission':0, 'Energy':[], 'percent_area':.5}
            self.distribution = distribution
            self.diffraction = diffraction
        elif widget == True: #run the widget, get the inputs from there
            import imager_widget as iw
            done=False
            while not done:
                w = iw.Imager_Widget()
                w.identify_choices()
                for att in ['Substrate_Material','Substrate_Thickness','Slits_Material','Slits_Thickness']:
                    try:
                        foo=getattr(w,att)
                        done=True
                    except AttributeError:
                        print 'You forgot to select the '+att+'! Please choose again.'
                        done=False
                        break
     
            self.attenuator = attenuator
            if w.Attenuator_Material!='None': # attenuator chosen
                self.attenuator = {'Material':w.Attenuator_Material,'Thickness':w.Attenuator_Thickness,'Transmission':1, 'Energy':[], 'eff_area':[]} 
            self.slits = {'Material':w.Slits_Material, 'Thickness':w.Slits_Thickness,'Percent_area':0.5,'Transmission':[], 'Energy':[]}
            self.substrate = {'Material':w.Substrate_Material, 'Thickness':w.Substrate_Thickness,'Substrate_area':1,'Transmission':[], 'Energy':[]}
            self.filled = filled
            try:
                if w.Filled_choices=='Y':
                    self.filled = {'Material':w.Substrate_Material,'Thickness':w.Slits_Thickness,'Transmission':1, 'Energy':[], 'percent_area':.5}
            except AttributeError: #assume no if nothingselected
                pass
            self.diffraction = diffraction
            self.detector = detector
            if w.Detector_Material!='None': # detector chosen
                self.detector = {'Material':w.Detector_Material,'Thickness':w.Detector_Thickness,'Transmission':0, 'Energy':[], 'eff_area':[]}

        #now read the csv's to fill in the transmission info
        import pandas as pd
        filenames =[self.substrate['Material']+'.csv',self.slits['Material']+'.csv']
        thicknesses = [self.substrate['Thickness'],self.slits['Thickness']]
        array=[0,0,0]
        if self.attenuator:
            filenames.append(self.attenuator['Material']+'.csv')
            thicknesses.append(self.attenuator['Thickness'])
            array[0] = 1
        if self.detector:
            filenames.append(self.detector['Material']+'.csv')
            thicknesses.append(self.detector['Thickness'])
            array[1] = 1
        if self.filled:
            filenames.append(self.filled['Material']+'.csv') 
            thicknesses.append(self.slits['Thickness']) #thickness is the same as slit thickness
            array[2] = 1
            
        energy = [[] for i in range(len(filenames))]
        transmission=[[] for i in range(len(filenames))]
        tinterp = [[] for i in range(len(filenames))]
        os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/calculations/data')
        for filename,thickness,index in zip(filenames,thicknesses,range(0,len(filenames))):
            data=pd.read_csv(filename,sep=',', header=0)
            energy[index] = data['E (keV)']
            try:
                transmission[index] = data['P(xi) (d='+thickness+' um)']
                #print transmission[index][0:10]
            except KeyError:
                print 'the key '+ 'P(xi) (d='+thickness+' um)'+ ' does not exist'
                #start over?

        os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/code/MiSolFA')

        #interpolate everything to the same energy range

        evector = np.linspace(2.08, 433, num=500, endpoint=True) #check that this spans the correct range
        #print evector[0:10]
        for i in range(0,len(filenames)):
            tinterp[i]=np.interp(evector,map(float,energy[i]),map(float,transmission[i]))

        #make any adjustments and assign back to self
        self.substrate['Energy']=evector
        self.slits['Energy']=evector   
        self.substrate['Transmission']=tinterp[0]
        self.slits['Transmission']=tinterp[1]
        i=2
        for a,p in zip(array,['attenuator','detector','filled']): #should check that this works correctly
            if a == 1:
                adict = getattr(self,p)
                setattr(self,p,{'Transmission':tinterp[i],'Energy': evector, 'Material':adict['Material'],'Thickness':adict['Thickness']})
                i=i+1
                
        self.substrate['eff_area']=self.substrate['Substrate_area']*self.substrate['Transmission']
        self.slits['eff_area']=self.slits['Percent_area']*self.slits['Transmission']
        if attenuator:
            self.substrate['eff_area']=self.substrate['eff_area']*self.attenuator['Transmission']
            self.slits['eff_area']=self.slits['eff_area']*self.attenuator['Transmission']
        if filled:
            self.slits['eff_area']=self.slits['eff_area']*self.filled['Transmission']
        if detector:
            self.substrate['eff_area']=self.substrate['eff_area']*(1-self.detector['Transmission'])**2
            self.slits['eff_area']=self.slits['eff_area']*(1-self.detector['Transmission'])**2
        self.eff_area = self.substrate['eff_area']- self.slits['eff_area']
                 
    def plot_eff_area(self): #should I make a plotting widget? possibly....
        '''Plot the effective area of the selected configuration.'''
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        #ax1.ylog(energy, eff_area, color="b", label="Planck")
        ax1.plot(self.slits['Energy'], self.eff_area, color="r",label="total",linewidth='2')
        ax1.plot(self.slits['Energy'], self.substrate['eff_area'], color="g",label="substrate",linewidth='2')
        ax1.plot(self.slits['Energy'], self.slits['eff_area'], color="b",label="slits",linewidth='2')
        try:
            ax1.plot(self.slits['Energy'], self.attenuator['Transmission'], color="m",label="attenuator+\ndetector (1mm CdTe)",linewidth='2') #make the label better
        except TypeError:
            pass
        plt.xlabel('Energy (keV)')
        plt.ylabel('Effective area (cm$^2$)')
        ax1.set_ylim([0,1])
        ax1.set_xlim([0,150])
        plt.title("Effective area of grids "+ self.substrate['Material'] + ' '+ self.substrate['Thickness'] + '$\mu$m, ' + self.slits['Material'] + ' '+self.slits['Thickness'] +'$\mu$m with ' +'$\mu$m attenuator')
        ax1.legend(loc='upper right',fontsize='medium')
        ax1.plot()
        fig.show()
        
    def plot_flare_counts(self,xlim=150,ylim=1000):
        '''Make the plot for if you have a small flare'''
        import flare_XR_dist as fd
        dist = fd.flare_XR_dist() #[energy, thermal part, non-thermal part]
        prob = np.interp(self.slits['Energy'],dist[0],dist[1])
        ntnt = np.interp(self.slits['Energy'],dist[0],dist[2])
        thth = np.interp(self.slits['Energy'],dist[0],dist[3])

        total_counts= prob*self.eff_area
        thermal_counts= thth*self.eff_area 
        nonthermal_counts= ntnt*self.eff_area
    
        fig = plt.figure()
        ax1 = fig.add_subplot(111)    #ax1.ylog(energy, eff_area, color="b", label="Planck")
        ax1.semilogy(self.slits['Energy'], total_counts, color="r",label="total",linewidth='2')
        ax1.semilogy(self.slits['Energy'], thermal_counts, color="g",label="thermal")
        ax1.semilogy(self.slits['Energy'], nonthermal_counts, color="b",label="non-thermal")
    
        plt.xlabel('Energy (keV)')
        plt.ylabel('Counts $s^{-1} keV^{-1}$')
        ax1.set_ylim([1,ylim])
        ax1.set_xlim([0,xlim])
        plt.title("Expected flare counts for grids "+ self.substrate['Material'] + ' '+ self.substrate['Thickness'] + '$\mu$m, ' + self.slits['Material'] + ' '+self.slits['Thickness'] +'$\mu$m with ' +'$\mu$m attenuator')
        ax1.legend(loc='upper right',fontsize='medium')
        fig.show()

    def print_config(self):
        '''print the relevant information about the imager'''
        print 'Substrate:  ' + self.substrate['Thickness'] + ' um of ' + self.substrate['Material'] 
        print 'Slits:      ' + self.slits['Thickness'] + ' um of ' + self.slits['Material'] 
        print 'Attenuator: ' + self.attenuator['Thickness'] + ' um of ' + self.attenuator['Material'] 
        print 'Detector:   ' + self.detector['Thickness'] + ' um of ' + self.detector['Material'] 
  
    def export2pickle(self, picklename):
        ''' Saves the data object in a .p file'''
        import pickle
        pickle.dump(self, open(picklename, 'wb'))

    def export2idl(self, idlname):
        ''' Saves the data object in a .p file'''
        import pidly
        print 'foo'
        
    def export2csv(self, csvname):
        ''' Saves the data object in a .p file'''
        d=self.__dict__
        columns,headers=[],[]
        for key in d.keys():
            try:
                for k in d[key].keys():
                    headers.append(key+'.'+k)
                    if type(d[key][k]) == int or type(d[key][k]) == float:
                        columns.append(str(d[key][k]))
                    else:
                        columns.append(d[key][k])
            except AttributeError:
                headers.append(key)
                if type(d[key]) == bool: #or int? Or will this not happen... prob not
                    columns.append(str(d[key]))
                else:
                    columns.append(d[key])

        import csv
        from itertools import izip_longest
        with open(csvname,'wb') as f:
            writer=csv.writer(f)
            writer.writerow(headers)
            for col in izip_longest(*columns): # works but it doesn't keep strings together... splits up by chars
                writer.writerow(col)



