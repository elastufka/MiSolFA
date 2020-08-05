"""
===================
grid.py
Erica  Lastufka 22.11.18
===================
Combine grating objects into a grid object

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
# global EMwindows #should these be inheirited from an Imager object if available? or just a file somewhere?
# global EMwindowsr
# global QMwindows
# global QMwindowsr
global winnums
global winfloats
winnums=['11','21','12','22','31','41','32','42','33','43','34','44']
winfloats=[11,21,12,22,31,41,32,42,33,43,34,44]


# EMwindows=[{'number':11,'pitch':89.6752,  'nominal angle':-44.79357,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':21,'pitch':90.326,  'nominal angle': 45.20793,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':12,'pitch':22.4797,'nominal angle': 44.94825,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':22,'pitch':22.5203,'nominal angle':-45.05184,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':31,'pitch':45.0814, 'nominal angle':-45.10378,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':41,'pitch':44.9187, 'nominal angle': 44.8966,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':32,'pitch':18.013, 'nominal angle': 45.04146,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':42,'pitch':17.987, 'nominal angle':-44.95859,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':33,'pitch':29.9639,  'nominal angle':-44.93102,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':43,'pitch':30.0362,  'nominal angle': 45.06914,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':34,'pitch':14.991, 'nominal angle': 44.96549,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':44,'pitch':15.009, 'nominal angle':-45.03455,'phase':0.,'height':250,'slit_width':.5}]

# EMwindowsr=[{'number':11,'pitch':90.326,  'nominal angle':45.20793,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':21,'pitch':89.6752,  'nominal angle': -44.79357,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':12,'pitch':22.5203,'nominal angle': -45.05184,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':22,'pitch':22.4797,'nominal angle':44.94825,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':31,'pitch':44.9187, 'nominal angle':44.8966,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':41,'pitch':45.0814, 'nominal angle': -45.10378,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':32,'pitch':17.987, 'nominal angle': -44.95859,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':42,'pitch':18.013, 'nominal angle':45.04146,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':33,'pitch':30.0362,  'nominal angle':45.06914,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':43,'pitch':29.9639,  'nominal angle': -44.93102,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':34,'pitch':15.009, 'nominal angle': -45.03455,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':44,'pitch':14.991, 'nominal angle':44.96549,'phase':0.,'height':250,'slit_width':.5}] #dectector side angle for ease, windows swapped...

# QMwindows=[{'number':11,'pitch':89.6,  'nominal angle':-44.77,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':21,'pitch':90.3,  'nominal angle': 45.23,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':12,'pitch':44.9,'nominal angle': 44.90,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':22,'pitch':45.1,'nominal angle':-45.10,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':31,'pitch':30.0, 'nominal angle':-45.08,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':41,'pitch':30.1, 'nominal angle': 44.92,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':32,'pitch':22.5, 'nominal angle': 45.08,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':42,'pitch':22.5, 'nominal angle':-44.92,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':33,'pitch':18.0,  'nominal angle':-44.93,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':43,'pitch':18.0,  'nominal angle': 45.07,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':34,'pitch':15.0, 'nominal angle': 44.94,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':44,'pitch':15.0, 'nominal angle':-45.06,'phase':0.,'height':250,'slit_width':.5}]

# QMwindowsr=[{'number':11,'pitch':90.3,  'nominal angle':45.23,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':21,'pitch':89.6,  'nominal angle': -44.77,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':12,'pitch':45.1,'nominal angle': -45.10,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':22,'pitch':44.9,'nominal angle':44.90,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':31,'pitch':30.1, 'nominal angle':44.92,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':41,'pitch':30.0, 'nominal angle': -45.08,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':32,'pitch':22.5, 'nominal angle': -44.92,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':42,'pitch':22.5, 'nominal angle':45.08,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':33,'pitch':18.0,  'nominal angle':45.09,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':43,'pitch':18.0,  'nominal angle':-44.93,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':34,'pitch':15.0, 'nominal angle':-45.06,'phase':0.,'height':250,'slit_width':.5},
#                 {'number':44,'pitch':15.0, 'nominal angle':44.94,'phase':0.,'height':250,'slit_width':.5}]
class Grid():

    def __init__(self,side=1.0,EM=True,obj_dict=False,obj_tag=False,combine_results_obj=False):
        '''Construct the simulated grid'''
        #self.size=[1.1,1.1] #cm

        if obj_dict: #this makes things easy
           self.__set_atts_from_dict__(obj_dict)
           self.sub = raw_input('What is the substrate number? eg. subXXXX ')

        elif EM:
            model='EM'
            if side==1.0: #get the pickled objects from correct directory
                self.sub='2737'
                opath='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/EMassembly_2017_11_15_FrontGrid/'
                xpath='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/'
            else: #rear grid
                self.sub='sub2765'
                opath='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw469sub2765_2018_01_31/'
                xpath='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Apr2018/'

        elif not EM: #QM
            model='QM'
            if side ==1.0:
                self.sub='3501'
                opath='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw561sub3501_2018_06_26/'
                xpath='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018'
            else:
                self.sub='3437'
                opath='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw501sub3437_2018_05_04/'
                xpath='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/'

        #now find the objects (currently no optical ones... fix this eventually....)
        oobj=glob.glob(opath+'sub'+self.sub+'_win*.p') #what about tags? deal with that later #[opath+'sub'+self.sub+'_win'+w+'.p' for w in winnums]
        xobj=glob.glob(xpath+'transm'+model+'/sub'+self.sub+'_win*.p')
        xobj.extend(glob.glob(xpath+'transm'+model+'1/sub'+self.sub+'_win*.p'))
        #filter out the ones without the desired tags
        xfobj=[xo for xo in xobj if len(xo) <= 87]
        oobj.sort()
        xfobj.sort()

        #for each grating object...
        if combine_results_obj:
            newobjs=[self.__combine_results_obj___(oob,xob) for oob,xob in zip(oobj,xfobj)] #assuming they are the same length
        else:
            if len(xfobj) > len(oobj):
                pick='xray'
            elif len(xfobj) == len(oobj):
                pick=raw_input('Use optical or Xray object? ')
            else:
                pick='o'
            if pick == 'optical' or pick =='o':
                newobjs=oobj
            else:
                newobjs=xfobj

        new_obj_dict={}
        for i,w in enumerate(winnums):
            new_obj_dict[w]=newobjs[i]

        self.__set_atts_from_dict__(new_obj_dict)

    def __set_atts_from_dict__(self,obj_dict):
           self.win11=pickle.load(open(obj_dict['11'],'rb'))
           self.win12=pickle.load(open(obj_dict['21'],'rb'))
           self.win21=pickle.load(open(obj_dict['12'],'rb'))
           self.win22=pickle.load(open(obj_dict['22'],'rb'))
           self.win31=pickle.load(open(obj_dict['31'],'rb'))
           self.win41=pickle.load(open(obj_dict['41'],'rb'))
           self.win32=pickle.load(open(obj_dict['32'],'rb'))
           self.win42=pickle.load(open(obj_dict['42'],'rb'))
           self.win33=pickle.load(open(obj_dict['33'],'rb'))
           self.win43=pickle.load(open(obj_dict['43'],'rb'))
           self.win34=pickle.load(open(obj_dict['34'],'rb'))
           self.win44=pickle.load(open(obj_dict['44'],'rb'))


    def __combine_results_obj___(self,optical,xray):
        '''put the Optical and X-ray results objects in the same Grating object'''
        xkeys=xray.results.keys()
        for k in keys:
            optical.results[k]=xray.results[k]
        return optical

    def __get_and_store_data__(self,plot=False,rename=False):
        '''sort different files with different tags into the right places in dictionary. Store metadata (coordinates)'''


    def flag_data(self,p=True,p0ang=False):
        '''In case of a bad sequence or frame, flag data so that it will not be processed but replaced with a constant value or nan instead depending on method...For now, always flag indivdual values before sequences or else it might get appended to the data lists for some reason'''

    def trim_optical_edges(self, maxx='10',maxy='13',lpix=50,rpix=50,tpix=50,bpix=50):
        ''' trim edges off of the microscope images '''
        os.chdir(self.data.Odata.path)
        imdata=self.data.Odata.rawim #list of images
        for filen in imdata:
            ag.remove_edges(filen,ag.im2ndarray(filen),maxx=maxx,maxy=maxy,lpix=lpix,rpix=rpix,tpix=tpix,bpix=bpix)

    def parameterize_optical(self,ll=[50,100,150,200,250],coverwrite=False,hoverwrite=False,tol=9,plot=True,stats=True,htol=5,ftags='a',pix2um=1.955,mask45=False,sigma=2.,r_or_f_tol=3,filter_nominal=False):
        '''Reduce data from optical measurements to set of parameters and error: period, orientation,phase,errors'''
        import analyze_optical as ao
        os.chdir(self.data.Odata.path)
        self.results['optical']=ao.analyze_optical(self, ll=ll,coverwrite=coverwrite,hoverwrite=hoverwrite,tol=tol,plot=plot,sigma=sigma,mask45=mask45,r_or_f_tol=r_or_f_tol,stats=plot,htol=htol,ftags=ftags,filter_nominal=filter_nominal)
        self.results['optical'].nominal=self.nominal

    def parameterize_transmission(self, method='sum',mkwargs=False,plot=True,sigma=False,tolvec=False):
        '''Reduce data from transmission measurements to set of parameters and error: period, orientation,height, transmission, errors'''
        import analyze_xray as aXray
        #mkwargs are the kwargs associated with the method
        for j in range(0,len(self.data.Xdata.transm)):
            #first test if window is in this folder or not
            try:
                os.chdir(self.data.Xdata.path + '/'+self.data.Xdata.transm[j]['folder'])
            except KeyError:
                 continue
            if 'level01' not in self.data.Xdata.transm[j].keys():
                continue
            if not mkwargs:
                updates,results = aXray.analyze_window(self.win,self.model,self.data.Xdata.transm[j],method,sigma=sigma,tolvec=tolvec)
                #print 'here'
            else:
                updates,results = aXray.analyze_window(self.win,self.model,self.data.Xdata.transm[j],method,**mkwargs) #updates need to be written to transm
            if len(updates) !=0:
                for i in range(0,len(updates)):
                    self.data.Xdata.transm[j][updates[i]['key']]=updates[i]['data'] #not right yet
            print j
        self.results[method]=results
        self.results[method].nominal=self.nominal
        if plot:
            if method == 'widths':
                yran=[.1,.6]
            else:
                yran=[.1,1.1]
            self.results[method].scat_allP(self.win,yran=yran)

    def compare_transm_methods(self,yran=False,figname=False):
        ''' Compare the methods'''
        import analyze_xray as aXray
        aXray.compare_methods(self.win,self.results['binary'],self.results['sum'],self.results['widths'],yran=yran,figname=figname)

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
        return calc_3D_grating(self.nominal)

    def calc_3D_nominal_error(self, error_params):
        '''Calculate and return binary 3D array describing the grating, given the nominal parameters and an error parameter dictionary specifiying amount and type of error'''
        return calc_3D_grating_error(self.nominal,error_params)

    def calc_3D_measured(self):
        '''Calculate and return binary 3D array describing the grating, given the measured parameters'''
        return calc_3D_grating(self.measured)

    def calc_3D_measured_error(self):
        '''Calculate and return binary 3D array describing the grating, given the measured parameters and error parameter dictionary'''
        return calc_3D_grating_error(self.measured,error_params)


    def calc_3D_grating(param_dict):
        '''Do the actual calculation, the other two are just wrapper methods'''
        #self.nominal={"pitch":,"orientation":,"phase":,"height":,"slit_width":]}
        #probably need the size too

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
        angle_factor=int(np.round(np.tan(np.deg2rad(angle)),0))

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


    ############################ Display and Output Methods ##################################

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

    def plot_optical_diffs(self): #need to re-work this to work with data stored within the object
        means,noms,yerr=[],[],[]
        x=np.linspace(1,12,12)
        for w in winnums:
            os.chdir(self.win11.data.Odata.path)
            stats=pickle.load(open('win'+w+'_width_stats5.0Xa.p','rb'))
            #stats=pickle.load(open('win'+str(w)+'_angle_stats_5.0.p','rb'))
            means.append(stats['mean'])
            noms.append(dummy.nominal['pitch'])
            #noms.append(dummy.nominal['orientation'])
            yerr.append(stats['stddev'])
        fig,ax=plt.subplots()
        ax.errorbar(x, np.array(means)-np.array(noms), yerr=np.transpose(np.array(yerr)),
                fmt='o', ecolor='g', capthick=2)
        fig.suptitle('Sub '+self.sub)
        ax.set_xticks(x)
        ax.set_xticklabels(winfloats)
        ax.set_xlabel('Window Number')
        ax.set_ylabel('Measured Period Difference ($\mu$m)')
        ax.set_xlim([0,13])
        #ax.set_ylim([-.5,.5])
        fig.show()

        fig,ax=plt.subplots()
        ax.errorbar(x, np.array(means)-np.array(noms), yerr=np.transpose(np.array(yerr)),
                fmt='o', ecolor='g', capthick=2)
        fig.suptitle('Sub '+self.sub)
        ax.set_xticks(x)
        ax.set_xticklabels(winfloats)
        ax.set_xlabel('Window Number')
        ax.set_ylabel('Measured Orientation Difference (degrees)')
        ax.set_xlim([0,13])
        #ax.set_ylim([-.5,.5])
        fig.show()

    #def optical_to_csv():
    #    '''print out optical analysis properties to csv file'''

    #def xray_to_csv():
    #    '''print out Xray analysis properties to csv file'''

