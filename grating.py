"""
===================
grating.py
Erica  Lastufka 15.5.11
===================
Given the ideal grid parameters, simulate the transmission images/profiles I should get out of them. For a SINGLE grating.

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from PIL import Image
from matplotlib import cm
import pickle
from scipy.misc import imrotate
from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit
from scipy import interpolate
import time
import glob
import random
import itertools
import analyze_general as ag


#list of dictionaries of the window numbers, pitches (mm), and nominal angles (deg)
global EMwindows #should these be inheirited from an Imager object if available? or just a file somewhere?
global EMwindowsr
global EMparams
global QMwindows
global QMwindowsr
global QMparams
global cbp

cbp = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

EMwindows=[{'number':11,'pitch':89.6752,  'nominal angle':-44.79357,'phase':0.,'height':250,'slit_width':.5},
                {'number':21,'pitch':90.326,  'nominal angle': 45.20793,'phase':0.,'height':250,'slit_width':.5},
                {'number':12,'pitch':22.4797,'nominal angle': 44.94825,'phase':0.,'height':250,'slit_width':.5},
                {'number':22,'pitch':22.5203,'nominal angle':-45.05184,'phase':0.,'height':250,'slit_width':.5},
                {'number':31,'pitch':45.0814, 'nominal angle':-45.10378,'phase':0.,'height':250,'slit_width':.5},
                {'number':41,'pitch':44.9187, 'nominal angle': 44.8966,'phase':0.,'height':250,'slit_width':.5},
                {'number':32,'pitch':18.013, 'nominal angle': 45.04146,'phase':0.,'height':250,'slit_width':.5},
                {'number':42,'pitch':17.987, 'nominal angle':-44.95859,'phase':0.,'height':250,'slit_width':.5},
                {'number':33,'pitch':29.9639,  'nominal angle':-44.93102,'phase':0.,'height':250,'slit_width':.5},
                {'number':43,'pitch':30.0362,  'nominal angle': 45.06914,'phase':0.,'height':250,'slit_width':.5},
                {'number':34,'pitch':14.991, 'nominal angle': 44.96549,'phase':0.,'height':250,'slit_width':.5},
                {'number':44,'pitch':15.009, 'nominal angle':-45.03455,'phase':0.,'height':250,'slit_width':.5}]

EMwindowsr=[{'number':11,'pitch':90.326,  'nominal angle':45.20793,'phase':0.,'height':250,'slit_width':.5},
                {'number':21,'pitch':89.6752,  'nominal angle': -44.79357,'phase':0.,'height':250,'slit_width':.5},
                {'number':12,'pitch':22.5203,'nominal angle': -45.05184,'phase':0.,'height':250,'slit_width':.5},
                {'number':22,'pitch':22.4797,'nominal angle':44.94825,'phase':0.,'height':250,'slit_width':.5},
                {'number':31,'pitch':44.9187, 'nominal angle':44.8966,'phase':0.,'height':250,'slit_width':.5},
                {'number':41,'pitch':45.0814, 'nominal angle': -45.10378,'phase':0.,'height':250,'slit_width':.5},
                {'number':32,'pitch':17.987, 'nominal angle': -44.95859,'phase':0.,'height':250,'slit_width':.5},
                {'number':42,'pitch':18.013, 'nominal angle':45.04146,'phase':0.,'height':250,'slit_width':.5},
                {'number':33,'pitch':30.0362,  'nominal angle':45.06914,'phase':0.,'height':250,'slit_width':.5},
                {'number':43,'pitch':29.9639,  'nominal angle': -44.93102,'phase':0.,'height':250,'slit_width':.5},
                {'number':34,'pitch':15.009, 'nominal angle': -45.03455,'phase':0.,'height':250,'slit_width':.5},
                {'number':44,'pitch':14.991, 'nominal angle':44.96549,'phase':0.,'height':250,'slit_width':.5}] #dectector side angle for ease, windows swapped...

EMparams={'11':{'tol':9, 'htol':6,'sigma':2,'r_or_f_tol':9,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':235},
          '21':{'tol':9, 'htol':6,'sigma':2,'r_or_f_tol':9,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':235},
          '12':{'tol':4, 'htol':5,'sigma':1.5,'r_or_f_tol':4,'csigma':3,'ll':[50,100,150,200,250],'spread':1.5,'n_hough':800},
          '22':{'tol':4, 'htol':5,'sigma':1.5,'r_or_f_tol':4,'csigma':3,'ll':[50,100,150,200,250],'spread':1.5,'n_hough':800},
          '31':{'tol':6, 'htol':6,'sigma':2,'r_or_f_tol':6,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':571},
          '41':{'tol':6, 'htol':6,'sigma':2,'r_or_f_tol':6,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':571},
          '32':{'tol':3, 'htol':3,'sigma':1.5,'r_or_f_tol':3,'csigma':3,'ll':[50,100,150,200,250],'spread':1.4,'n_hough':800},
          '42':{'tol':3, 'htol':3,'sigma':1.5,'r_or_f_tol':3,'csigma':3,'ll':[50,100,150,200,250],'spread':1.4,'n_hough':800},
          '33':{'tol':5, 'htol':5,'sigma':1.5,'r_or_f_tol':5,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':800},
          '43':{'tol':5, 'htol':5,'sigma':1.5,'r_or_f_tol':5,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':800},
          '34':{'tol':2, 'htol':2,'sigma':1.,'r_or_f_tol':2,'csigma':2,'ll':[20,50,75,100,150],'spread':1.,'n_hough':800},
          '44':{'tol':2, 'htol':2,'sigma':1.,'r_or_f_tol':2,'csigma':2,'ll':[20,50,75,100,150],'spread':1.,'n_hough':800}}

QMwindows=[{'number':11,'pitch':89.6,  'nominal angle':-44.77,'phase':0.,'height':250,'slit_width':.5},
                {'number':21,'pitch':90.3,  'nominal angle': 45.23,'phase':0.,'height':250,'slit_width':.5},
                {'number':12,'pitch':44.9,'nominal angle': 44.90,'phase':0.,'height':250,'slit_width':.5},
                {'number':22,'pitch':45.1,'nominal angle':-45.10,'phase':0.,'height':250,'slit_width':.5},
                {'number':31,'pitch':30.0, 'nominal angle':-45.08,'phase':0.,'height':250,'slit_width':.5},
                {'number':41,'pitch':30.1, 'nominal angle': 44.92,'phase':0.,'height':250,'slit_width':.5},
                {'number':32,'pitch':22.5, 'nominal angle': 45.08,'phase':0.,'height':250,'slit_width':.5},
                {'number':42,'pitch':22.5, 'nominal angle':-44.92,'phase':0.,'height':250,'slit_width':.5},
                {'number':33,'pitch':18.0,  'nominal angle':-44.93,'phase':0.,'height':250,'slit_width':.5},
                {'number':43,'pitch':18.0,  'nominal angle': 45.07,'phase':0.,'height':250,'slit_width':.5},
                {'number':34,'pitch':15.0, 'nominal angle': 44.94,'phase':0.,'height':250,'slit_width':.5},
                {'number':44,'pitch':15.0, 'nominal angle':-45.06,'phase':0.,'height':250,'slit_width':.5}]

QMwindowsr=[{'number':11,'pitch':90.3,  'nominal angle':45.23,'phase':0.,'height':250,'slit_width':.5},
                {'number':21,'pitch':89.6,  'nominal angle': -44.77,'phase':0.,'height':250,'slit_width':.5},
                {'number':12,'pitch':45.1,'nominal angle': -45.10,'phase':0.,'height':250,'slit_width':.5},
                {'number':22,'pitch':44.9,'nominal angle':44.90,'phase':0.,'height':250,'slit_width':.5},
                {'number':31,'pitch':30.1, 'nominal angle':44.92,'phase':0.,'height':250,'slit_width':.5},
                {'number':41,'pitch':30.0, 'nominal angle': -45.08,'phase':0.,'height':250,'slit_width':.5},
                {'number':32,'pitch':22.5, 'nominal angle': -44.92,'phase':0.,'height':250,'slit_width':.5},
                {'number':42,'pitch':22.5, 'nominal angle':45.08,'phase':0.,'height':250,'slit_width':.5},
                {'number':33,'pitch':18.0,  'nominal angle':45.09,'phase':0.,'height':250,'slit_width':.5},
                {'number':43,'pitch':18.0,  'nominal angle':-44.93,'phase':0.,'height':250,'slit_width':.5},
                {'number':34,'pitch':15.0, 'nominal angle':-45.06,'phase':0.,'height':250,'slit_width':.5},
                {'number':44,'pitch':15.0, 'nominal angle':44.94,'phase':0.,'height':250,'slit_width':.5}]

# QMparams={'11':{'tol':9, 'htol':6,'sigma':2,'r_or_f_tol':9,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':235},
#           '21':{'tol':9, 'htol':6,'sigma':2,'r_or_f_tol':9,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':235},
#           '32':{'tol':4, 'htol':5,'sigma':1.5,'r_or_f_tol':4,'csigma':3,'ll':[50,100,150,200,250],'spread':1.5,'n_hough':800},
#           '42':{'tol':4, 'htol':5,'sigma':1.5,'r_or_f_tol':4,'csigma':3,'ll':[50,100,150,200,250],'spread':1.5,'n_hough':800},
#           '12':{'tol':6, 'htol':6,'sigma':2,'r_or_f_tol':6,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':571},
#           '22':{'tol':6, 'htol':6,'sigma':2,'r_or_f_tol':6,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':571},
#           '33':{'tol':3, 'htol':3,'sigma':1.5,'r_or_f_tol':3,'csigma':3,'ll':[50,100,150,200,250],'spread':1.4,'n_hough':800},
#           '43':{'tol':3, 'htol':3,'sigma':1.5,'r_or_f_tol':3,'csigma':3,'ll':[50,100,150,200,250],'spread':1.4,'n_hough':800},
#           '31':{'tol':5, 'htol':5,'sigma':1.5,'r_or_f_tol':5,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':800},
#           '41':{'tol':5, 'htol':5,'sigma':1.5,'r_or_f_tol':5,'csigma':3,'ll':[50,100,150,200,250],'spread':2.,'n_hough':800},
#           '34':{'tol':2, 'htol':2,'sigma':1.5,'r_or_f_tol':2,'csigma':2,'ll':[20,50,75,100,150],'spread':1.,'n_hough':800},
#           '44':{'tol':2, 'htol':2,'sigma':1.5,'r_or_f_tol':2,'csigma':2,'ll':[20,50,75,100,150],'spread':1.,'n_hough':800}}

QMparams={'11':{'tol':9, 'htol':6,'sigma':2,'r_or_f_tol':9,'csigma':3,'ll':[15,20,30],'spread':2.,'n_hough':235},
          '21':{'tol':9, 'htol':6,'sigma':2,'r_or_f_tol':9,'csigma':3,'ll':[15,20,30],'spread':2.,'n_hough':235},
          '32':{'tol':4, 'htol':5,'sigma':1.5,'r_or_f_tol':4,'csigma':3,'ll':[15,20,30],'spread':1.5,'n_hough':800},
          '42':{'tol':4, 'htol':5,'sigma':1.5,'r_or_f_tol':4,'csigma':3,'ll':[15,20,30],'spread':1.5,'n_hough':800},
          '12':{'tol':6, 'htol':6,'sigma':2,'r_or_f_tol':6,'csigma':3,'ll':[15,20,30],'spread':2.,'n_hough':571},
          '22':{'tol':6, 'htol':6,'sigma':2,'r_or_f_tol':6,'csigma':3,'ll':[15,20,30],'spread':2.,'n_hough':571},
          '33':{'tol':3, 'htol':3,'sigma':1.5,'r_or_f_tol':3,'csigma':3,'ll':[15,20,30],'spread':1.4,'n_hough':800},
          '43':{'tol':3, 'htol':3,'sigma':1.5,'r_or_f_tol':3,'csigma':3,'ll':[15,20,30],'spread':1.4,'n_hough':800},
          '31':{'tol':5, 'htol':5,'sigma':1.5,'r_or_f_tol':5,'csigma':3,'ll':[15,20,30],'spread':2.,'n_hough':800},
          '41':{'tol':5, 'htol':5,'sigma':1.5,'r_or_f_tol':5,'csigma':3,'ll':[15,20,30],'spread':2.,'n_hough':800},
          '34':{'tol':2, 'htol':2,'sigma':1.,'r_or_f_tol':2,'csigma':2,'ll':[15,20,30],'spread':1.,'n_hough':800},
          '44':{'tol':2, 'htol':2,'sigma':1.,'r_or_f_tol':2,'csigma':2,'ll':[15,20,30],'spread':1.,'n_hough':800}}

#EM
# Xparams={'11':{'tolvec':[9,9,7,7,7,5,7,7,7,9,9], 'sigma':2,'psigma':2},
#           '21':{'tolvec':[9,9,7,7,7,5,7,7,7,9,9], 'sigma':2,'psigma':2},
#           '32':{'tolvec':[5,5,5,7,7,9,7,7,7,5,5,5], 'sigma':1.5,'psigma':1.5},
#           '42':{'tolvec':[5,5,5,7,7,9,7,7,7,5,5,5], 'sigma':1.5,'psigma':1.5},
#           '12':{'tolvec':[5,5,5,6,6,6,5,5,5], 'sigma':1.25,'psigma':2},
#           '22':{'tolvec':[5,5,5,6,6,6,5,5,5], 'sigma':1.25,'psigma':2},
#           '33':{'tolvec':[5,5,5,6,6,6,5,5,5], 'sigma':1.25,'psigma':1.5},
#           '43':{'tolvec':[5,5,5,6,6,6,5,5,5], 'sigma':1.25,'psigma':1.5},
#           '31':{'tolvec':[5,5,5,5,6,6,6,5,5,5,5], 'sigma':.75,'psigma':1.5},
#           '41':{'tolvec':[5,5,5,5,6,6,6,5,5,5,5], 'sigma':.75,'psigma':1.5},
#           '34':{'tolvec':[4,4,4,5,5,5,4,4,4], 'sigma':.75,'psigma':1.},
#           '44':{'tolvec':[4,4,4,5,5,5,4,4,4], 'sigma':.75,'psigma':1.}}

#QM
Xparams={'11':{'tolvec':[9,9,7,7,7,5,7,7,7,9,9], 'sigma':2,'psigma':2},
          '21':{'tolvec':[9,9,7,7,7,5,7,7,7,9,9], 'sigma':2,'psigma':2},
          '12':{'tolvec':[5,5,5,7,7,9,7,7,7,5,5,5], 'sigma':1.5,'psigma':1.5},
          '22':{'tolvec':[5,5,5,7,7,9,7,7,7,5,5,5], 'sigma':1.5,'psigma':1.5},
          '31':{'tolvec':[5,5,5,5,6,6,6,5,5,5,5], 'sigma':1.25,'psigma':2},
          '41':{'tolvec':[5,5,5,5,6,6,6,5,5,5,5], 'sigma':1.25,'psigma':2},
          '32':{'tolvec':[5,5,5,6,6,6,5,5,5], 'sigma':1.25,'psigma':1.5},
          '42':{'tolvec':[5,5,5,6,6,6,5,5,5], 'sigma':1.25,'psigma':1.5},
          '33':{'tolvec':[5,5,5,5,6,6,6,5,5,5,5], 'sigma':.75,'psigma':1.5},
          '43':{'tolvec':[5,5,5,5,6,6,6,5,5,5,5], 'sigma':.75,'psigma':1.5},
          '34':{'tolvec':[4,4,4,5,5,5,4,4,4], 'sigma':.75,'psigma':1.},
          '44':{'tolvec':[4,4,4,5,5,5,4,4,4], 'sigma':.75,'psigma':1.}}

class Grating():

    def __init__(self,win,side=1.0,nominal_thickness=False,EM=True,source=False,analyze=False,data_path=False,calc_measured=False,rename=False,Apr=False,May=True, Sept=False,gparams=False):
        '''Construct the simulated grating'''
        self.size=[1.1,1.1] #cm
        if Apr:
            EM=True
            analyze='Xray'
            data_path='/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/SLS_Apr2018'
            opath='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/EMassembly_2017_11_15_FrontGrid'
            if side==-1.0:
                #opath='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/EMassembly_2017_11_15_RearGrid'
                opath='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw469sub2765_2018_01_31'
            #this is reversed! May is when we did the EM front grid. April was the rear grid.
            #sort out later, the results are what's important for now
        if May:
            analyze='Xray'
            data_path='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018'
            if EM:
                opath='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw469sub2765_2018_01_31'
                side=-1.0
            else:
                opath='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw501sub3437_2018_05_04'
        if Sept:
            EM=False
            analyze='Xray'
            data_path='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018'
            #if side == 1.0: #front grid
            opath='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw561sub3501_2018_06_26'
            #else: #rear grid
            #opath='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw501sub3501_2018_06_26'

        if EM:
            self.model='EM'
            fgparams=EMwindows
            rgparams=EMwindowsr
        elif not gparams:
            self.model='QM'
            fgparams=QMwindows
            rgparams=QMwindowsr

        if side==1.0:
            if not gparams:
                gparams=fgparams
            sidestr='front'
        else:
            if not gparams:
                gparams=rgparams
            sidestr='rear'

        windex=[gp['number'] for gp in gparams].index(win)
        self.nominal={"pitch":gparams[windex]['pitch'],"orientation":gparams[windex]['nominal angle'],"phase":gparams[windex]['phase'],"height":gparams[windex]['height'],"slit_width":gparams[windex]['slit_width']}
        self.win=win
        self.side=side
        self.source=source #contains the pix2um conversion...
        self.analyze=analyze #can call methods using self.analyze.ExampleMethod() ? perhaps...how can I make the attributes of self visible to the Analyze object?
        #self.method=method
        self.data=Data(empty=True)
        if data_path:
            xpath=data_path
            #opath=data_path
        else:
            if EM:
                xpath='EMmodel/'+sidestr
            else:
                xpath='QMmodel/'+sidestr
        try:
            self.data.Xdata=XData(path=xpath,logfile=glob.glob(xpath+'/*.log')[0],imlist=glob.glob(xpath+'/'+str(win)+'*.tif')) #think more about how this should work later
        except IndexError: #it's transmission so there's one more layer of folders
            self.data.Xdata=XData(path=xpath) #think more about how this should work later
        if May or Apr or Sept:
            try:
                self.data.Odata=OData(path=opath,logfile=glob.glob(opath+'/*.log')[0],imlist=glob.glob(opath+'/'+str(win)+'*.tif')) #think more about how this should work later
            except IndexError: #it's transmission so there's one more layer of folders
                self.data.Odata=OData(path=opath) #think more about how this should work later

        if analyze:
            self.data.type=analyze #optical or transmission?
            self.results=''
            import analyze_general as analyze
            if self.data.type == 'optical':
                import analyze_optical as ao
            elif self.data.type == 'Xray':
                import analyze_xray as aXray
            self.__get_and_store_data__(rename=rename)

        #now construct the simulated nominal grating. Assume it's NOT rotated for starters
        #self.nominal_grating=self.calc_3D_nominal()

        ###probably the wrong place for this###
        #if requested construct the simulated measured grating. First have to set the measured attributes, by use of the Analyze methods?
        #if calc_measured:
        #    self.measured_grating=self.calc_3D_measured()

        #if there is data to analyze, get and store the filenames and coordinates of the images
        self.__get_and_store_data__(rename=rename)

        if not self.results:
            self.results={}


    def __get_and_store_data__(self,plot=False,rename=False):
        '''sort different files with different tags into the right places in dictionary. Store metadata (coordinates)'''

        def _sort_by_p(imlist):
            pdict={'p0':'','p1':'','p2':'','p3':'','p4':'','p5':'','p6':''}
            for p in range(0,7):
                key='p'+str(p)
                plist=[im for im in imlist if key in im]
                #print plist
                pdict[key] = plist
                #if there's no p, just return the original list
                if plist ==[]:
                    return {'imlist':imlist}

            #print pdict
            return pdict

        if self.data.type == 'optical' or self.data.Odata != '':
          win=self.win
          opath=self.data.Odata.path
          os.chdir(opath)
          plist=glob.glob(opath+'/win'+str(win)+'*.p')
          import make_mosaic as mm
          if self.model=='EM':
            self.data.Odata.rawim=ag.EM_list(str(win),'5.0',ending='.tif')
            self.data.Odata.corrim=ag.EM_list(str(win),'5.0',ending='_corrected.tif')
            self.data.Odata.edges=ag.EM_list(str(win),'5.0',ending='_corrected_edges.p')
            self.data.Odata.cat_hough=ag.EM_list(str(win),'5.0',ending='corrected_hough_all.p')
            self.data.Odata.coords=mm.get_coordinates(str(win),self.nominal,opath,self.data.Odata.rawim)
          else:
            self.data.Odata.rawim=glob.glob('win'+str(win)+'*_5.0X.tif')
            if self.side !=-1.0:
                self.data.Odata.corrim=ag.EM_list(str(win),'5.0',ending='_corrected.tif')
                self.data.Odata.edges=[im for im in plist if im.endswith('_corrected_edges.p')]
                self.data.Odata.hough=[im for im in plist if im.endswith('_hough_ll*.p')]
            else: #it's may
                self.data.Odata.corrim=ag.EM_list(str(win),'5.0',ending='_May_corrected.tif')
                self.data.Odata.edges=[im for im in plist if im.endswith('_May_corrected_edges.p')]
                self.data.Odata.hough=[im for im in plist if im.endswith('_May_corrected_hough_ll*.p')]

            self.data.Odata.cat_hough=[im for im in plist if im.endswith('_hough_all.p')]
            self.data.Odata.coords=mm.get_coordinates(str(win),self.nominal, opath,self.data.Odata.rawim)

        #elif self.data.type == 'Xray': #lamp or transmission?
        #first get flagfile:
        fdir='/Users/wheatley/Documents/Solar/MiSolFA/code/MiSolFA/'
        if self.side == 1.0:
            tside='Front'
        else:
            tside='Rear'
        ff=fdir+self.model+tside+'flags.txt'
        self.data.Xdata.flagfile=ff
        self.flags_from_file(ff)

        xpath=self.data.Xdata.path
        root,dirs,files=os.walk(xpath).next()
        if dirs !=[]: #transmission. need to re-write imlist
            for d in dirs:
                if d.endswith('_'):
                    dirs.remove(d)
            for d in dirs:
                if d.endswith('_'):
                    dirs.remove(d)
            moire_dirs=[d for d in dirs if d.startswith('moire')]
            try:
                transm_dirs=[d for d in dirs if d.startswith('transm'+self.model) and not d.endswith('_')]
            except AttributeError: #no model specified
                transm_dirs=[d for d in dirs if d.startswith('transm') and not d.endswith('_')]
            self.data.Xdata.logfile=transm_dirs[0]+'.log'
            log_dir='log'
            for mdir in moire_dirs:
                os.chdir(xpath+'/'+mdir)
                mimlist=glob.glob('win'+str(self.win)+'*.tif')
                os.chdir(xpath+'/../')
                self.data.Xdata.moire={'mdir':mdir,'imlist':mimlist}
            for i,tdir in enumerate(transm_dirs): #group by level, type, position
                os.chdir(xpath+'/'+tdir)
                timlist=glob.glob('win'+str(self.win)+'*.tif')
                #print timlist
                oimlist=glob.glob('*window0'+str(self.win)+'*.tif')

                if oimlist != [] and rename: #rename images
                    import analyze_xray as aXray
                    onames,newnames=aXray.rename_files(self.win,oimlist,logfile=False)
                else:
                    onames,newnames=oimlist,timlist
                if onames == [] and newnames ==[]: #there's nothing there. Don't make data dict
                    continue

                try:
                    self.data.Xdata.transm[i]['folder']=tdir
                except IndexError:
                    self.data.Xdata.transm.append({})
                    self.data.Xdata.transm[i]['folder']=tdir
                self.data.Xdata.transm[i]['level00']=_sort_by_p([o for o in onames if o[-5].isdigit()]) #level 0
                if newnames !=[]:
                    self.data.Xdata.transm[i]['level01']=_sort_by_p([o for o in newnames if o[-5].isdigit()]) #level 1 (renamed)
                test2=[o for o in newnames if o.endswith('_flatdark.tif')]
                if test2 !=[]:
                    self.data.Xdata.transm[i]['level02']= _sort_by_p(test2) #level 2
                test3=[o for o in newnames if o.endswith('_corrected.tif')] #level 3
                if test3 !=[]:
                   self.data.Xdata.transm[i]['level03']=_sort_by_p(test3)
                test4=[o for o in newnames if o.endswith('_groupstretch.tif')] #level 4
                if test4 !=[]:
                    self.data.Xdata.transm[i]['level04']=_sort_by_p(test4)
                #has any higher-level processing been done? Check for edges, lines etc
                plist = glob.glob(xpath + '/'+tdir+'/win'+str(self.win)+'*.p')
                os.chdir('/../')
                test5=[o for o in plist if o.endswith('_edges.p')]
                if test5 !=[]:
                    self.data.Xdata.transm[i]['level05']= _sort_by_p(test5)#level 5 - Canny edges
                test6=[o for o in plist if 'hough' in o]
                #print test6
                if test6 !=[]:
                    self.data.Xdata.transm[i]['level06']= _sort_by_p(test6) #level 6

    def flags_from_file(self,flagfile):
        flaglist=[] #[[],[],[],[],[],[],[],[],[],[],[],[]]
        with open(flagfile) as f:
           lines=f.readlines()
        for l in lines:
           if l != '\n':
               ll=l[:-1].split(',')
               #ll[2]=ll[2][:-2]
               if int(ll[0])==self.win:
                   #windex=wins.index(int(ll[0]))
                   flaglist.append(ll[1:]) #need int if just p
        if lines[0] !='\n':
            for fl in flaglist:
                if fl[0] =='p':
                    self.flag_data(p=fl[1])
                else:
                    self.flag_data(p0ang=fl[1])
        if flaglist ==[]:
            for i,t in enumerate(self.data.Xdata.transm):
                self.data.Xdata.transm[i]['flags']=[]



    def flag_data(self,p=True,p0ang=False):
        '''In case of a bad sequence or frame, flag data so that it will not be processed but replaced with a constant value or nan instead depending on method...For now, always flag indivdual values before sequences or else it might get appended to the data lists for some reason'''
        from copy import deepcopy
        if type(p) == int: #flag whole sequence of p's, at all levels. Ex. Win 21 needs p0 and p1 flagged.
            for i,t in enumerate(self.data.Xdata.transm):
                tcopy=deepcopy(t)
                for k in tcopy.keys():
                    if k.startswith('level') and k !='level00':
                        try:
                            self.data.Xdata.transm[i][k]['flagged'].extend(tcopy[k]['p'+str(p)])
                        except KeyError:
                            self.data.Xdata.transm[i][k]['flagged']=tcopy[k]['p'+str(p)]
                        #print self.data.Xdata.transm[i][k]['flagged'],t[k]['p'+str(p)]
        elif type(p0ang) == str: #flag individual image described by p0 and angle (as string)
            for i,t in enumerate(self.data.Xdata.transm):
                try:
                    self.data.Xdata.transm[i]['flags'].append(p0ang)
                except KeyError:
                    self.data.Xdata.transm[i]['flags']=[p0ang]
#             for i,t in enumerate(self.data.Xdata.transm):
#                 tcopy=deepcopy(t)
#                 for k in tcopy.keys():
#                     if k.startswith('level') and k !='level00':
#                         try:
#                             example=tcopy[k]['p0'][0]
#                             suffix=example[example.rfind('_'):]
#                             self.data.Xdata.transm[i][k]['flagged'].extend(['win'+str(self.win)+'_'+p0ang+suffix])
#                         except KeyError:
#                             self.data.Xdata.transm[i][k]['flagged']=['win'+str(self.win)+'_'+p0ang+suffix]

    def trim_optical_edges(self, maxx='10',maxy='13',lpix=50,rpix=50,tpix=50,bpix=50):
        ''' trim edges off of the microscope images '''
        os.chdir(self.data.Odata.path)
        imdata=self.data.Odata.rawim #list of images
        for filen in imdata:
            ag.remove_edges(filen,ag.im2ndarray(filen),maxx=maxx,maxy=maxy,lpix=lpix,rpix=rpix,tpix=tpix,bpix=bpix)

    def parameterize_optical(self,ll=[50,100,150,200,250],coverwrite=False,eoverwrite=False, hoverwrite=False,tol=9,plot=True,stats=True,htol=5,ftags='a',pix2um=1.955,mask45=False,csigma=3,sigma=2.,r_or_f_tol=3,spread=2.0,n_hough=235,filter_nominal=False,example=False,auto=True,angle_only=False):
        '''Reduce data from optical measurements to set of parameters and error: period, orientation,phase,errors'''
        import analyze_optical as ao
        os.chdir(self.data.Odata.path)
        if auto: #use the parameters from the global dictionary
            if self.model=='EM':
                pdict=EMparams[str(self.win)]
            else:
                pdict=QMparams[str(self.win)]
            tol=pdict['tol']
            htol=pdict['htol']
            sigma=pdict['sigma']
            r_or_f_tol=pdict['r_or_f_tol']
            csigma=pdict['csigma']
            spread=pdict['spread']
            n_hough=pdict['n_hough']
            ll=pdict['ll']

        self.results['optical']=ao.analyze_optical(self, ll=ll,coverwrite=coverwrite,hoverwrite=hoverwrite,tol=tol,plot=plot,sigma=sigma,mask45=mask45,r_or_f_tol=r_or_f_tol,stats=plot,htol=htol,spread=spread,n_hough=n_hough,ftags=ftags,filter_nominal=filter_nominal,example=example, angle_only=angle_only)
        #self.results['optical'].nominal=self.nominal

    def parameterize_transmission(self, method='sum',mkwargs=False,plot=True,sigma=False,tolvec=False,angle_only=False,autoflag=True,auto=True,write=True):
        '''Reduce data from transmission measurements to set of parameters and error: period, orientation,height, transmission, errors'''
        import analyze_xray as aXray
        #mkwargs are the kwargs associated with the method
        if autoflag:
            self.flags_from_file(self.data.Xdata.flagfile)
        if auto: #use the parameters from the global dictionary
            pdict=Xparams[str(self.win)]
            psigma=pdict['psigma']
            tolvec=pdict['tolvec']
            sigma=pdict['sigma']

        for j in range(0,len(self.data.Xdata.transm)):
            #first test if window is in this folder or not
            try:
                os.chdir(self.data.Xdata.path + '/'+self.data.Xdata.transm[j]['folder'])
            except KeyError:
                 continue
            if 'level04' not in self.data.Xdata.transm[j].keys():
                continue
            if not mkwargs:
                updates,results = aXray.analyze_window(self,j,method,sigma=sigma,tolvec=tolvec,angle_only=angle_only,psigma=psigma)
                #print 'here'
            else:
                updates,results = aXray.analyze_window(self,j,method,**mkwargs) #updates need to be written to transm
            if len(updates) !=0:
                for i in range(0,len(updates)):
                    self.data.Xdata.transm[j][updates[i]['key']]=updates[i]['data'] #not right yet
            print j
        if write:
            self.results[method]=results
            self.results[method].nominal=self.nominal
        if write and plot:
            if method == 'widths':
                yran=[.1,.6]
            else:
                yran=[.1,1.1]
            self.results[method].scat_allP(self.win,yran=yran)


    def compare_transm_methods(self,yran=False,figname=False,xvmid=False,ex_mid='L',excl=False,overwrite=False):
      ''' Compare the methods'''
      import analyze_xray as aXray
      if self.model=='EM' and self.side == 1.:
        title='EM1 Front'
      elif self.model=='EM' and self.side == -1.:
        title='EM1 Rear'
      elif self.model=='QM' and self.side == 1.:
        title='EM2 Front'
      elif self.model=='QM' and self.side == -1.:
        title='EM2 Rear'

      if 'fdict' not in self.results.keys() or overwrite:
          ddict,fdict=aXray.compare_methods(self.win,self.results['sum'],self.results['widths'],self.nominal['pitch']/2.,yran=yran,figname=figname,title=title,xvmid=xvmid,ex_mid=ex_mid,excl=excl)
          self.results['fdict']=fdict
          self.results['ddict']=ddict
      else:
          aXray.compare_methods_plot(self.win,self.results['ddict'],self.results['fdict'],yran=yran,figname=figname,title=title)


    def fit_profile_to_data_optimize(self,pindex,plot=True,ret=False,height_bounds=False, width_bounds=False,dc_bounds=False,yoff_bounds=False,fit_bounds=False,fix_dc=True,show=True,xvmid=False, oplot_nominal=True,soverwrite=False,from_sum=False):
        period=self.nominal['pitch']
        height=self.nominal['height']

        if from_sum:
            dc=self.results['fdict']['int_y']/self.results['widths'].period_stats['mean']
            self.results['sum'].fit_profile_to_data_optimize(pindex,period,height,plot=plot,ret=ret,height_bounds=height_bounds, width_bounds=width_bounds,dc_bounds=dc_bounds,yoff_bounds=yoff_bounds,fit_bounds=fit_bounds,fix_dc=fix_dc,show=show,xvmid=xvmid, oplot_nominal=oplot_nominal,soverwrite=soverwrite,sum_fac=dc)
        else: #use widths
            self.results['widths'].fit_profile_to_data_optimize(pindex,period,height,plot=plot,ret=ret,height_bounds=height_bounds, width_bounds=width_bounds,dc_bounds=dc_bounds,yoff_bounds=yoff_bounds,fit_bounds=fit_bounds,fix_dc=fix_dc,show=show,xvmid=xvmid, oplot_nominal=oplot_nominal,soverwrite=soverwrite,sum_fac=False)

#     def parameterize_moire():
#         '''Reduce data from moire measurements to set of parameters and error: period, orientation,height, transmission, errors'''
#         self.measured='foo'

#     def calc_moire_nominal(self,imsize=False):
#         '''Calculate nominal moire period and orientation. If requested, also an image (this goes in Subcollimator, not Grating)'''
#         import plot_divergent as pd
#         fp=
#         rp=
#         fa=
#         ra=
#         self.nominal.moire['period'],self.nominal.moire['orientation'] = pd.calc_moire_period(fp,rp,fa,ra,eff=False,quiet=True):

#     def calc_transm_2D_nominal(self, ret=False):
#         '''Calculate and return dictionary of binary 2D array describing the transmission images, given the nominal parameters'''
#         #for each image used to calculate the 'widths' method, find the first edge


#         #for each image, find the number of horizontal gaps (and their size?)

#         #
#         if ret:
#             return dict2d


    def transm_func(xnew,w0,h,dc,yoff): #can I fit the duty cycle too?
        return yoff+(w0-np.abs(h*np.tan(np.deg2rad((xnew-pkth)))))/(w0/dc)

    def calc_nominal_slat_centers(self,optical=True):
        '''1D vector of ideal slat centers, n+ at least 1, to compare with measured values in order to calculate pitch RMS, starting point for phase measurement'''
        if optical: #pixel size = 1.995 um, img size = [640,480]
            pix2um=1.955
            imdim= [5650,5650] #let's round up quite a bit from 5626 (11mm to pix)
        else: #pixel size= ?, img size= 2kx2k?, assume transmission images
            pix2um=1.955
            imdim=[640,480]

        pitch=self.nominal['pitch']/pix2um
        #let's just assume for now that orientation is more or less 45 degrees
        ncenters=imdim[0]*np.sqrt(2.)/pitch
        icenters=int(ncenters)+5
        slat_centers=np.linspace(0,pitch*icenters,num=icenters+1) #this is now a vector of float slat center positions that is slightly larger than necessary...
        self.nominal['slat_centers']=slat_centers


    def calc_nominal_transm(self, th_vec=False,th_offsets=False):
        '''calculate the transmission profile given ideal conditions'''
        if "transm_prof" not in self.nominal.__dict__.keys():
            self.nominal.transm_prof = {}
        #if not th_vec: #get from a file
        for p in range(0,7):
            yvals=transm_func(th_vec,self.nominal['pitch']/self.nominal['slit_width'],self.nominal['height'],self.nominal['slit_width'],0.) #no y-offset
            self.nominal.transm_prof['p'+str(pindex)]=yvals

    def calc_3D_nominal(self):
        '''Calculate and return binary 3D array describing the grating, given the nominal parameters'''
        return self.calc_3D_grating(nominal=True)

    def calc_3D_nominal_error(self, error_params):
        '''Calculate and return binary 3D array describing the grating, given the nominal parameters and an error parameter dictionary specifiying amount and type of error'''
        return calc_3D_grating_error(self.nominal,error_params)

    def calc_3D_measured(self):
        '''Calculate and return binary 3D array describing the grating, given the measured parameters'''
        return calc_3D_grating(self.measured)

    def calc_3D_measured_error():
        '''Calculate and return binary 3D array describing the grating, given the measured parameters and error parameter dictionary'''
        return calc_3D_grating_error(self.measured,error_params)


    def calc_3D_grating(self,nominal=False):
        '''Do the actual calculation, the other two are just wrapper methods'''
        #self.nominal={"pitch":,"orientation":,"phase":,"height":,"slit_width":]}
        #probably need the size too
        if nominal:
            param_dict=self.nominal

        #or else?

        #get the pixel to um conversion
        if not self.source:
            pix2um=1.
        else:
            pix2um=self.source.pix2um

        #intialize array of zeros in the correct size
        arrx=int(self.size[0]*10000.*pix2um) #cm to um to pix
        arry=int(self.size[1]*10000.*pix2um) #cm to um to pix
        arrz=int(param_dict['height']*pix2um) #um to pix

        garr=np.zeros([arrx,arry,arrz])

        #convert quantities to pixels
        period_pix=int(pix2um*param_dict['pitch'])
        angle_pix=param_dict['orientation']
        phase_pix=int(pix2um*param_dict['phase'])
        width_pix=int(pix2um*param_dict['slit_width']) #height doesn't matter...
        slat_pix=period_pix-width_pix

        #calculate the number of rows you can go down before you need to shift to get the orientation
        angle_factor=int(np.round(np.tan(np.deg2rad(angle_pix)),0))

        #put 1's where the grid slats are. Probably easiest to do this row-by row
        #define phase to be the offset of the rising edge of the first slat w.r.t. corner (0,0,*) in the array
        for i in range(0,arry):
            nperiods=arrx/period_pix #hopefully this is integer
            #print nperiods
            if phase_pix >=0:
                garr[phase_pix:phase_pix+slat_pix,i,:]=1
            else:
                garr[0:phase_pix+slat_pix,i,:]=1
            n=1
            while n<nperiods+1:
                start_pix=phase_pix+n*period_pix
                end_pix=start_pix+slat_pix
                #print n, start_pix,end_pix
                garr[start_pix:end_pix,i,:]=1
                #print np.sum(garr[:,i,:])
                n+=1
            #get the orientation by changing the phase pixel location
            if i % angle_factor ==0:
                phase_pix+=1 #of course if this gets bigger than a slit width you have to reset to 0..
                if phase_pix >width_pix:
                    phase_pix=-width_pix
        return garr

    def calc_3D_grating_error(param_dict,error_params):
        '''Introduce error in the grating, given an error parameter dictionary'''
        #start off with an ideal grating array
        errarr=calc_3D_grating(param_dict)

        #what kind of error to introduce? parse the error_params dict
        ekeys=error_params.keys()
        if 'slit_width' in ekeys:
            swdev=error_params['slit_width'][0]
            swlen=error_params['slit_width'][1]        #length scale over which error affects a single slit. -1 means entire slit
            swpercent=error_params['slit_width'][2]
        if 'angle' in ekeys:
            angdev=error_params['angle'][0]
            anglen=error_params['angle'][1]        #length scale over which error affects a single slat's orientation
            angpercent=error_params['angle'][2]
        if 'slat_missing' in ekeys:
            smdev=error_params['slat_missing'][0]
            smlen=error_params['slat_missing'][1]        #length scale
            smpercent=error_params['slat_missing'][2]
        if 'thickness' in ekeys:
            thdev=error_params['thickness'][0]
            thlen=error_params['thickness'][1]        #length scale
            thpercent=error_params['thickness'][2]

        for key in ekeys:
            if key == 'slit_width':
                errarr=_add_slit_width_err(errarr,swdev,swlen,swpercent)
            elif key == 'angle':
                errarr=_add_angle_err(errarr,angdev,anglen,angpercent)
            elif key == 'slat_missing':
                errarr=_add_slat_missing_err(errarr,smdev,smlen,smpercent)
            elif key == 'thickness':
                errarr=_add_thickness_err(errarr,thdev,thlen,thpercent,slat_width=param_dict['slit_width']) #slit width is close enough for now
        return errarr

    ############################ Error Addition Methods (internal) ##################################

    def _random_gen(low, high):
        while True:
            yield random.randrange(low, high)

    def _add_slit_width_error(arr,dev,elen,percent):
        return errarr

    def _add_angle_error(arr,dev,elen,percent):
        '''should affect BOTH sides of slat'''
        return errarr

    def _add_slat_missing_error(arr,dev,elen,percent):
        return errarr

    def _add_thickness_error(arr,dev,elen,percent,slat_width=False):
        #first flatten array
        arr2d=np.sum(arr,axis=2)
        arrnonzero=np.where(arr2d !=0.0)
        nslatpix=np.sum(arr2d,axis=2)
        nomth=np.max(arr2d)
        tmin=nomth-dev
        tmax=nomth+dev
        npixvar=int(tmax)-int(tmin) #number of pixels in variation that is possible
        #varran=np.linspace(int(tmin),int(tmax),npixvar)
        if len == -1: #entire array is affected
            #generate random vector of thicknesses and gaussian smooth it locally? use gaussian with same width as slat if length scale is not defined
            rvec=set()
            gen=_random_gen(tmin,tmax)
            for x in itertools.takewhile(lambda x: len(items) < nslatpix, gen):
                rvec.add(x)
            #now assign values in rvec back to slats in arr2D
            for item,loc in zip(rvec,arrnonzero):
                arr2d[loc]=item
        #3D gaussian smooth to get rid of discontinuities
        #first make it a masked array so that the edges are preserved


        #re-inflate arr2D

        return errarr


    ############################ Display Methods ##################################

    def plot3D(self, nominal=True):
        '''plot the binary nominal or measured array in 3D. Doesn't quite work yet.'''
        if nominal:
            garr=self.nominal_grating
        else:
            garr=self.measured_grating

        arrx,arry,arrz=np.shape(garr)
        X,Y=np.meshgrid(range(0,arrx),range(0,arry))
        Z=garr[:,:,0]

        from mpl_toolkits.mplot3d import Axes3D
        fig=plt.figure()
        ax=fig.add_subplot(111,projection='3d')
        ax.plot_surface(X,Y,Z)
        fig.show()

    def plot2D(self,nominal=True,rotate=False):
        '''plot the binary nominal or measured array in 2D, rotated at a given angle'''
        if nominal:
            garr=self.nominal_grating
        else:
            garr=self.measured_grating

        plotim=np.transpose(garr[:,:,0])
        #if rotate:
        #do something else

        fig,ax=plt.subplots()
        ax.imshow(plotim,cmap=cm.binary,origin='bottom left')
        fig.show()

    #def plot1D(self,x):

class XData():
    def __init__(self, path='.', logfile=False, imlist=False):
        '''Empty object for storing X-ray data in Grating()'''
        #self.type='Xray'
        self.path=path
        self.logfile=logfile
        self.imlist=imlist
        self.moire={}
        self.transm=[{}]

class OData():
    def __init__(self, path='.',logfile=False, imlist=False):
        '''Empty object for storing Optical data in Grating()'''
        self.path=path
        self.logfile=logfile
        self.imlist=imlist
        self.edges=[]
        self.hough=[]
        self.cat_hough=[]
        self.rawim=[]
        #self.type='Optical'

class Data():
    def __init__(self,empty=True):
        '''Empty object for storing data in Grating()'''
        self.Xdata=''
        self.Odata=''




