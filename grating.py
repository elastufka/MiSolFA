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


#list of dictionaries of the window numbers, pitches (mm), and nominal angles (deg)
global EMwindows #should these be inheirited from an Imager object if available? or just a file somewhere?
global EMwindowsr
global QMwindows
global QMwindowsr

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

EMwindowsr=[{'number':11,'pitch':90.326,  'nominal angle':-45.20793,'phase':0.,'height':250,'slit_width':.5},
                {'number':21,'pitch':89.6752,  'nominal angle': 44.79357,'phase':0.,'height':250,'slit_width':.5},
                {'number':12,'pitch':22.5203,'nominal angle': 45.05184,'phase':0.,'height':250,'slit_width':.5},
                {'number':22,'pitch':22.4797,'nominal angle':-44.94825,'phase':0.,'height':250,'slit_width':.5},
                {'number':31,'pitch':44.9187, 'nominal angle':-44.8966,'phase':0.,'height':250,'slit_width':.5},
                {'number':41,'pitch':45.0814, 'nominal angle': 45.10378,'phase':0.,'height':250,'slit_width':.5},
                {'number':32,'pitch':17.987, 'nominal angle': 44.95859,'phase':0.,'height':250,'slit_width':.5},
                {'number':42,'pitch':18.013, 'nominal angle':-45.04146,'phase':0.,'height':250,'slit_width':.5},
                {'number':33,'pitch':30.0362,  'nominal angle':-45.06914,'phase':0.,'height':250,'slit_width':.5},
                {'number':43,'pitch':29.9639,  'nominal angle': 44.93102,'phase':0.,'height':250,'slit_width':.5},
                {'number':34,'pitch':15.009, 'nominal angle': 45.03455,'phase':0.,'height':250,'slit_width':.5},
                {'number':44,'pitch':14.991, 'nominal angle':-44.96549,'phase':0.,'height':250,'slit_width':.5}] #dectector side angle for ease, windows swapped...

QMwindows=[{'number':11,'pitch':89.6,  'nominal angle':-44.77,'phase':0.,'height':250,'slit_width':.5},
                {'number':21,'pitch':90.3,  'nominal angle': 45.23,'phase':0.,'height':250,'slit_width':.5},
                {'number':12,'pitch':44.9,'nominal angle': 44.90,'phase':0.,'height':250,'slit_width':.5},
                {'number':22,'pitch':45.1,'nominal angle':-45.10,'phase':0.,'height':250,'slit_width':.5},
                {'number':31,'pitch':30.0, 'nominal angle':-45.08,'phase':0.,'height':250,'slit_width':.5},
                {'number':41,'pitch':30.01, 'nominal angle': 44.92,'phase':0.,'height':250,'slit_width':.5},
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
class Grating():

    def __init__(self,win,side=1,nominal_thickness=False,EM=True,source=False,analyze=False,data_path=False,calc_measured=False,rename=False,Apr=False,May=True, Sept=False):
        '''Construct the simulated grating'''
        self.size=[1.1,1.1] #cm
        if Apr:
            EM=True
            analyze='Xray'
            data_path='/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/SLS_Apr2018'
        if May:
            analyze='Xray'
            data_path='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018'
        if Sept:
            EM=False
            analyze='Xray'
            data_path='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018'
        if EM:
            self.model='EM'
            fgparams=EMwindows
            rgparams=EMwindowsr
        else:
            self.model='QM'
            fgparams=QMwindows
            rgparams=QMwindowsr
        if side==1:
            gparams=fgparams
            sidestr='front'
        else:
            gparams=rgparams
            sidestr='rear'

        windex=[gp['number'] for gp in gparams].index(win)
        self.nominal={"pitch":gparams[windex]['pitch'],"orientation":gparams[windex]['nominal angle'],"phase":gparams[windex]['phase'],"height":gparams[windex]['height'],"slit_width":gparams[windex]['slit_width']}
        self.win=win
        self.source=source #contains the pix2um conversion...
        self.analyze=analyze #can call methods using self.analyze.ExampleMethod() ? perhaps...how can I make the attributes of self visible to the Analyze object?
        #self.method=method
        self.data=Data(empty=True)
        if data_path:
            xpath=data_path
            opath=data_path
        else:
            if EM:
                xpath='EMmodel/'+sidestr
                opath='EMmodel/'
            else:
                xpath='QMmodel/'+sidestr
                opath='QMmodel/'
        try:
            self.data.Xdata=XData(path=xpath,logfile=glob.glob(xpath+'/*.log')[0],imlist=glob.glob(xpath+'/'+str(win)+'*.tif')) #think more about how this should work later
        except IndexError: #it's transmission so there's one more layer of folders
            self.data.Xdata=XData(path=xpath) #think more about how this should work later
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
            #print pdict
            return pdict

        if self.data.type == 'optical':
            opath=self.data.Odata.path
            plist=glob.glob(opath+'/*.p')
            self.data.Odata.edges==[im for im in plist if im.endswith('_edges.p')]
            self.data.Odata.hough==[im for im in plist if im.endswith('_hough_ll*.p')]
            self.data.Odata.cat_hough==[im for im in plist if im.endswith('_hough_all.p')]
        elif self.data.type == 'Xray': #lamp or transmission?
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
                transm_dirs=[d for d in dirs if d.startswith('transm'+self.model) and not d.endswith('_')]
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
                    test6=[o for o in plist if o.endswith('_hough.p')]
                    print test6
                    if test6 !=[]:
                        self.data.Xdata.transm[i]['level06']= _sort_by_p(test6) #level 6

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
                tcopy=deepcopy(t)
                for k in tcopy.keys():
                    if k.startswith('level') and k !='level00':
                        try:
                            example=tcopy[k]['p0'][0]
                            suffix=example[example.rfind('_'):]
                            self.data.Xdata.transm[i][k]['flagged'].extend(['win'+str(self.win)+'_'+p0ang+suffix])
                        except KeyError:
                            self.data.Xdata.transm[i][k]['flagged']=['win'+str(self.win)+'_'+p0ang+suffix]


    def parameterize_optical(self):
        '''Reduce data from optical measurements to set of parameters and error: period, orientation,phase,errors'''
        self.measured='foo'

    def parameterize_transmission(self, method='sum',mkwargs=False,plot=True):
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
                updates,results = aXray.analyze_window(self.win,self.model,self.data.Xdata.transm[j],method)
            else:
                updates,results = aXray.analyze_window(self.win,self.model,self.data.Xdata.transm[j],method,**mkwargs) #updates need to be written to transm
            if len(updates) !=0:
                for i in range(0,len(updates)):
                    self.data.Xdata.transm[j][updates[i]['key']]=updates[i]['data'] #not right yet
        self.results[method]=results
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
#         '''Calcualte nominal moire period and orientation. If requested, also an image (this goes in Subcollimator, not Grating)'''
#         import plot_divergent as pd
#         fp=
#         rp=
#         fa=
#         ra=
#         self.nominal.moire['period'],self.nominal.moire['orientation'] = pd.calc_moire_period(fp,rp,fa,ra,eff=False,quiet=True):

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

class Data():
    def __init__(self,empty=True):
        '''Empty object for storing data in Grating()'''
        self.Xdata=''
        self.Odata=''




