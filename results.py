"""
===================
results.py
Erica  Lastufka 13.9.18
===================
Results object to store results, errors, steps performed
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
from analyze_xray import read_logfile
import os
from grating import EMparams, QMparams

global EMreq
global QMreq

EMreq={'11':{'pitch':0.085, 'angle':0.085},
          '21':{'pitch':0.085, 'angle':0.085},
          '12':{'pitch':0.011, 'angle':0.011},
          '22':{'pitch':0.011, 'angle':0.011},
          '31':{'pitch':0.035, 'angle':0.035},
          '41':{'pitch':0.035, 'angle':0.035},
          '32':{'pitch':0.007, 'angle':0.007},
          '42':{'pitch':0.007, 'angle':0.007},
          '33':{'pitch':0.018, 'angle':0.018},
          '43':{'pitch':0.018, 'angle':0.018},
          '34':{'pitch':0.005, 'angle':0.005},
          '44':{'pitch':0.005, 'angle':0.005}}

QMreq={'11':{'pitch':0.095, 'angle':0.095},
          '21':{'pitch':0.095, 'angle':0.095},
          '12':{'pitch':0.038, 'angle':0.038},
          '22':{'pitch':0.038, 'angle':0.038},
          '31':{'pitch':0.023, 'angle':0.023},
          '41':{'pitch':0.023, 'angle':0.023},
          '32':{'pitch':0.012, 'angle':0.012},
          '42':{'pitch':0.012, 'angle':0.012},
          '33':{'pitch':0.0087, 'angle':0.0087},
          '43':{'pitch':0.0087, 'angle':0.0087},
          '34':{'pitch':0.0063, 'angle':0.0063},
          '44':{'pitch':0.0063, 'angle':0.0063}}

class Results():

    def __init__(self,win,kind,flags,method=False,errors=False,data=False,nominal=False,filenames=False,stats=False,params=False):
        '''Make an object to hold results. if optical: kind ={'type':optical,'model':model,'folder':folder,'side':side}'''
        self.win=win
        self.kind=kind
        self.flags=flags
        if self.kind['type'] == 'optical' and self.kind['model']=='EM':
            self.requirements=EMreq[str(win)]
        elif self.kind['type'] == 'optical' and self.kind['model']=='QM':
            self.requirements=QMreq[str(win)]
        if method:
            self.method=method #else it's probably not relevant if not Xray
        if errors:
            self.errors=errors
        if data:
            self.data=data
        else:
            self._load_data()
        if nominal:
            self.nominal=nominal
        #else:
        #    if self.kind['type'] == 'optical' and self.kind['model']=='EM':
        #        self.nominal=
        #    elif self.kind['type'] == 'optical' and self.kind['model']=='QM':
        #        self.nominal=
        if stats:
            self.stats=stats
        else:
            self.stats=self._get_stats() # get_stats  - do stats on data
        if 'theta' in self.data.keys():
            self.angle_stats={'theta0':self.data['theta0'],'theta':self.data['theta'],'theta_std':self.data['theta_std']}
            self.fix_angle_xray()
        if filenames:
            self.filenames=filenames
        else:
            self.filenames={} # gen_filenames from data <-- what filenames are these supossed to be?
        if params:
            self.calc_params=params
        else: #use default params from gratings globals. if optical
            if self.kind['type']=='optical' and self.kind['model']=='EM':
                self.calc_params=EMparams[str(self.win)]
            elif self.kind['type']=='optical' and self.kind['model']=='QM':
                self.calc_params=QMparams[str(self.win)]
        self.tex={}

    def _load_data(self):
        ''' load data, optical or x-ray'''
        folder=self.kind['folder']
        win=self.win
        os.chdir(folder)
        #default ftags='a'
        datafile='win'+str(win)+'_width_data5.0Xa.p'
        self.data=pickle.load(open(datafile,'rb'))
        os.chdir('../..')

    def _get_stats(self,overwrite=False):
        '''get stats from statfiles. if no statfiles or overwrite, get from data?'''
        folder=self.kind['folder']
        win=self.win
        #default ftags='a'
        os.chdir(folder)
        pstatfile='win'+str(win)+'_width_stats5.0Xa.p'
        try:
            self.period_stats=pickle.load(open(pstatfile,'rb'))
        except IOError:
            pass
        astatfile='win'+str(win)+'_angle_stats_5.0.p'
        try:
            self.angle_stats=pickle.load(open(astatfile,'rb'))
        except IOError:
            os.chdir
            hdata=glob.glob('win'+str(win)+'*hough_ll200.p')
            if hdata !=[]:
                theta=ag.hough_hist(hdata,win,0.0,figname=figname,tol=htol,spread=2.,mask45=False,ret=True) #make histogram
                theta0=aXray.read_logfile(data_dict['folder']+'.log', win)
                self.angle_stats={'theta0':theta0,'theta':theta,'theta_std':np.nanstd(theta)}
        os.chdir('../..')

    def fix_angle_xray(self,lim=88.):
        th=np.abs(self.angle_stats['theta'])
        th.sort()
        th=np.array(th)
        try:
            aa=np.where(th < lim)[0][-1]
            newth=th[aa+1:]
            self.angle_stats['theta_std']=np.std(newth)
            self.angle_stats['theta']=newth
        except IndexError:
            newth=self.angle_stats['theta']
        #adjust w.r.t. nominal...
        thmean=np.nanmean(newth)
        th0=self.angle_stats['theta0']
        fac=int(th0/45.)
        #want th_adj to have the same sign as the original...
        if fac % 2 !=0:
            th0_adj=45.*(fac+1) -th0
        else:
            th0_adj=45.*(fac) -th0
        thm_adj=90.-thmean
        nang=self.nominal['orientation']
        if np.sign(thm_adj) != np.sign(nang):
            thm_adj=-1.*thm_adj
        if np.sign(nang) ==1:
            delta_nom = min(np.abs(th0_adj+thm_adj)-nang,np.abs(th0_adj-thm_adj)-nang)
        else:
            delta_nom=max(np.abs(th0_adj+thm_adj)+nang,np.abs(th0_adj-thm_adj)+nang)

        print th0,fac,th0_adj,thm_adj,nang,delta_nom
        self.angle_stats['delta_nom']=delta_nom

    def reload_widthdata_from_files(self,flagged=False):
        #{'xdata':xvec,'ydata':ydata,'mean_ang':ydataa, \
        #                                             'raw_widths_P':sumvecp,'raw_widths_A':sumveca, \
        #                                            'raw_periods_P':sumpp,'raw_periods_A':sumpa, \
        #                                            'mean_widths_P':meanvecp,'mean_widths_A':meanveca, \
        #                                            'mean_periods_P':meanpp,'mean_periods_A':meanpa}

        #DON'T LOAD FLAGGED DATA UNLESS SPECIFIED
        #basic check
        if self.method != 'widths':
            return

        #build filenames
        xvec=self.data['xdata']
        xstr=['{:g}'.format(xv) for xv in xvec]
        print xstr
        p0ang=['p'+str(pp)+'_'+xv for pp in range(0,7) for xv in xstr]
        if not flagged: #if flagged, get rid of it...
            flags=self.flags
            p0new=[p for p in p0ang if p not in flags]
        else:
            p0new=p0ang
        #dfiles=['win'+str(self.win)+'_width_data_'+p+'.p' for p0new]
        #sfiles=['win'+str(self.win)+'_width_stats_'+p+'.p' for p in p0new]
        #print dfiles[:10]
        #print sfiles[:10]
        newparams={}
        allwidths,allperiods=[],[]
        for p in p0ang:
            if p in p0new:
                data=pickle.load(open('win'+str(self.win)+'_width_data_'+p+'.p','rb'))
                aa=pickle.load(open('win'+str(self.win)+'_width_stats_'+p+'.p','rb'))
                try:
                    sig=aa[1]
                    n=aa[2]
                except KeyError:
                    try:
                        sig=self.calc_params['sigma']
                        n=self.calc_params['tolvec']
                        psig=self.calc_params['psigma']
                    except KeyError:
                        #sig,n=self.calc_params[s] #what was s??
                        continue

                #newparams[s]=[sig,n]
                pp=data['period']
                ww=data['widths']
            else:
                pp=np.empty(len(xstr))
                pp[:]=np.NaN
                ww=np.empty(len(xstr))
                ww[:]=np.NaN

            allwidths.append(np.array(ww)) #need to put row of null where flagged
            allperiods.append(np.array(pp))
        allwidthsp=[allwidths[i::len(xvec)] for i in range (0,len(xvec))]#[allwidths[i::7] for i in range(0,7)]
        allwidthsa=[allwidths[i::7] for i in range (0,7)]
        allperiodsp=[allperiods[i::len(xvec)] for i in range (0,len(xvec))]
        allperiodsa=[allperiods[i::7] for i in range (0,7)]

        sumvecp,sumveca,sumpp,sumpa=np.array(allwidthsp),np.array(allwidthsa),np.array(allperiodsp),np.array(allperiodsa)

        meanvecp,meanpp=[],[]
        for p in range(0,7):
            meanvecp.append([np.mean(sumvecp[i,p]) for i in range(0,len(xvec))])
            meanpp.append([np.mean(sumpp[i,p]) for i in range(0,len(xvec))])

        meanvecp=np.array(meanvecp)
        meanveca=np.transpose(meanvecp)#np.mean(sumveca,axis=1)

        meanpp=np.array(meanpp)
        meanpa=np.transpose(meanpp)#np.mean(sumveca,axis=1)
        ydata=meanvecp/meanpp #duty cycle, grouped by p
        ydataa=meanveca/meanpa #duty cycle, grouped by angle
        self.data={'xdata':xvec,'ydata':ydata,'mean_ang':ydataa, \
                                                     'raw_widths_P':sumvecp,'raw_widths_A':sumveca, \
                                                    'raw_periods_P':sumpp,'raw_periods_A':sumpa, \
                                                    'mean_widths_P':meanvecp,'mean_widths_A':meanveca, \
                                                    'mean_periods_P':meanpp,'mean_periods_A':meanpa}
        #self.calc_params=newparams


    def _mask_raw(self,var,p=True,percent=0.5):
        '''mask outliers in raw data'''
        maskedarr=[]
        count=0
        if p: #sort by p
            dx,dy=np.shape(var)
        else:
            dy,dx=np.shape(var)

        for i in range(0,dx):
            maskedrow=[]
            for j in range(0,dy):
                med=np.median(var[i][j])
                marr=np.ma.masked_outside(var[i][j], percent*med,(1.0+(1.0-percent))*med)
                count+=np.ma.count_masked(marr)
                maskedrow.append(marr)
            maskedarr.append(maskedrow)
        return np.array(maskedarr),count

    def inspect_edges(self,pval,ang,ret=False):
        '''tool to look at the edges and cleaned centres'''
        win=self.win
        dfile=glob.glob('win'+str(win)+'_width_data_p'+str(pval)+'_'+ang+'.p')[0]
        efile=glob.glob('win'+str(win)+'_p'+str(pval)+'_'+ang+'_corrected_edges.p')[0]
        edges=pickle.load(open(efile,'rb'))
        ce=ag.clean_centers(edges, tolerance=False,plot=False)
        ag.plot_centers_and_edges(edges,ce,dfile)
        if ret:
            return ce

    def recalc_sc(self,pix2um=1.955):
        '''recalculate slit and slat centers'''
        rising=self.data['rising']
        falling=self.data['falling']
        #from rising_or_falling_final
        slitc,slatc=[],[]
        for rise,fall in zip(rising,falling):
            flist=[0 for f in fall]
            rlist=[1 for r in rise]
            rlist.extend(flist)
            masterlist=rise
            masterlist.extend(fall)
            locs, rorf = (list(t) for t in zip(*sorted(zip(masterlist, rlist))))
            cslit,cslat=[],[]

            for i,loc,code in zip(range(0,len(locs[:-1])),locs[:-1],rorf[:-1]):
                testsum=code+rorf[i+1]
                if testsum == 1: #r and f (but in what order?)
                    if code ==1: #r then f
                        cslat.append((locs[i+1]-loc)/2. + loc) #slat centers, not slit centers
                    else:
                        cslit.append((locs[i+1]-loc)/2. + loc) #slit centers
            slatc.extend(cslat)
            slitc.extend(cslit)
        self.data['slat_centers']=slatc
        self.data['slit_centers']=slitc

        #now calculate periods
        def calc_p(centers):
            period=[]
            for n,c in enumerate(centers[:-1]):
                period.append(centers[n+1]-c)
            return period

        pslat=pix2um*np.array(calc_p(slatc))
        pslit=pix2um*np.array(calc_p(slitc))
        self.data['periods_from_slat_centers']=pslat
        self.data['periods_from_slit_centers']=pslit
        #update stats
        self.period_stats['anum']=len(pslat)
        self.period_stats['amean']=np.mean(pslat)
        self.period_stats['amed']=np.median(pslat)
        self.period_stats['adev']=np.std(pslat)
        self.period_stats['inum']=len(pslit)
        self.period_stats['imean']=np.mean(pslit)
        self.period_stats['imed']=np.median(pslit)
        self.period_stats['idev']=np.std(pslit)


    def filter_xray_widths(self,key,fac=2.,percent=5.,update_stats=True,update_data=True,pix2um=0.65):
        '''wrapper for do_filter_nominal and do-filter_multiples'''
        period=pix2um*self.data[key]
        wkey='mean_widths'+key[-2:]
        width=pix2um*self.data[wkey]
        npitch=self.nominal['pitch']
        dim=np.shape(period)
        fparr=[]
        mparr=np.zeros(dim)
        scount=0
        for i in range(0,dim[0]):
            frow=[]
            for j in range (0,dim[1]):
                fp,mult=ag.do_filter_nominal(period[i][j],npitch,fac=fac)
                frow.append(list(fp))
                mparr[i][j]=np.mean(fp)
                scount=scount+len(list(fp))
            fparr.append(frow)

        newp=np.array(fparr)
        #mean, etc
        newy=np.transpose(width)/mparr

        if update_stats:
            snum=key+'num'
            smean=key+'mean'
            smed=key+'med'
            sdev=key+'dev'
            self.period_stats[snum]=scount
            self.period_stats[smean]=np.mean(mparr)
            self.period_stats[smed]=np.median(mparr)
            self.period_stats[sdev]=np.std(mparr)

        if update_data:
            self.data[key]=newp
            self.data['mean_periods'+key[-2:]]=mparr #have to update mean whatevers too
            self.data['ydata']=newy

    def filter_optical_widths(self,key,fac=2.,percent=5.,update_stats=True,update_data=True):
        '''wrapper for do_filter_nominal and do-filter_multiples'''
        period=self.data[key]
        npitch=self.nominal['pitch']
        fp,mult=ag.do_filter_nominal(period,npitch,fac=fac)
        #fm=ag.do_filter_multiples(mult,npitch,frac=fac)
        #fp.extend(fm)
        newp=fp
        if len(key)>10:
            newkey='filtered_periods_'+key[13:17]
            iora=key[15]
        else:
            newkey='filtered_periods'
        self.data[newkey]=newp
        if update_stats:
            snum=iora+'num'
            smean=iora+'mean'
            smed=iora+'med'
            sdev=iora+'dev'
            self.period_stats[snum]=len(newp)
            self.period_stats[smean]=np.mean(newp)
            self.period_stats[smed]=np.median(newp)
            self.period_stats[sdev]=np.std(newp)
        if update_data:
            self.data[key]=newp

    def filter_outliers_widths(self,percent=0.5):
        '''mask out elements that are multiple times the median value in a raw data set of widths and periods'''
        raww=self.data['raw_widths_P']
        rawp=self.data['raw_periods_P']
        rawwA=self.data['raw_widths_A']
        rawpA=self.data['raw_periods_A']
        #create the masks
        mraww,nmasked_wp=self._mask_raw(raww,percent=percent)
        mrawp,nmasked_pp=self._mask_raw(rawp,percent=percent)
        mrawwA,nmasked_wa=self._mask_raw(rawwA,percent=percent)#,p=False)
        mrawpA,nmasked_pa=self._mask_raw(rawpA,percent=percent)#,p=False)

        #write masked arrays to self, along with a record of what happened
        #take out x-values corresponding to NaN in y-values? or do this later for the fitting?
        mwp=self.data['mean_widths_P']
        mpp=self.data['mean_periods_P']
        mwa=self.data['mean_widths_A']
        mpa=self.data['mean_periods_A']
        yd=self.data['ydata']
        yda=self.data['mean_ang']
        xvec=self.data['xdata']

        changelog={'mean_widths_P_old':mwp,'mean_periods_P_old':mpp,'mean_widths_A_old':mwa,'mean_periods_A_old':mpa,'ydata_old':yd,'mean_ang_old':yda} #keep the old means so I don't have to recalc

        #recalculate means with masked arrays
        meanvecp,meanpp=[],[]
        for p in range(0,7):
            meanvecp.append([np.ma.mean(mraww[i,p]) for i in range(0,len(xvec))])
            meanpp.append([np.ma.mean(mrawp[i,p]) for i in range(0,len(xvec))])
        meanvecp=np.array(meanvecp)
        meanveca=np.transpose(meanvecp)#np.mean(sumveca,axis=1)
        #meanpp=np.array([[np.mean(svp) for svp in sumvp] for sumvp in sumpp])
        #for a in range(0,len(xvec)):
        #    meanpp.append([np.mean(sumpp[i,a]) for i in range(0,7)])

        meanpp=np.array(meanpp)
        meanpa=np.transpose(meanpp)#np.mean(sumveca,axis=1)
        ydata=meanvecp/meanpp #duty cycle, grouped by p
        ydataa=meanveca/meanpa #duty cycle, grouped by angle

        self.data['ydata']=ydata
        self.data['ydataa']=ydataa
        self.data['mean_widths_P']=meanvecp
        self.data['mean_periods_P']=meanpp
        self.data['mean_widths_A']=meanveca
        self.data['mean_periods_A']=meanpa

        nmasked={'widths':nmasked_wp,'periods':nmasked_pp}
        self.meta={'procedure':'mask outliers','nmasked':nmasked,'changelog':changelog}

    def unfilter(self):
        '''restore old unfiltered data'''
        cl=self.meta['changelog']
        self.data['mean_widths_P']=cl['mean_widths_P_old']
        self.data['mean_periods_P']=cl['mean_periods_P_old']
        self.data['mean_widths_A']=cl['mean_widths_A_old']
        self.data['mean_periods_A']=cl['mean_periods_A_old']
        self.data['ydata']=cl['ydata_old']
        self.data['mean_ang']=cl['mean_ang_old']
        self.meta=None

    def interpolate_transm(self,xdata,ydata,yints=False,period=False,height=False,plot=True):
        '''interpolate transmission results to find the y-intercepts'''
        from scipy.interpolate import interp1d
        #ydata=self.data['ydata'][0]
        #xdata=self.data['xdata']
        if period and height:
            thetac=np.rad2deg(np.arctan2(period/2.,height)) #theoretical thetac, with some padding
        f=interp1d(xdata,ydata,fill_value="extrapolate")
        if yints != False:
            xnew=np.linspace(yints[0],yints[1],100)
        else:
            xnew=np.linspace(-5.*thetac,5.*thetac)
        ynew=f(xnew)
        if plot:
            fig,ax=plt.subplots()
            ax.scatter(xdata,ydata)
            ax.plot(xnew,ynew)
            ax.set_ylim([0,1])
            fig.show()
        return xnew,ynew

    def fit_profile_to_data(self,pindex,period,height,plot=True,fixed='height',ret=False,fit_bounds=[1,-1]):
        '''fit a triangular transmission profile to the given one. Return average width and tilt angle'''
        #get data
        ydata=list(self.data['ydata'][pindex])[fit_bounds[0]:fit_bounds[1]]
        xv=self.data['xdata'][fit_bounds[0]:fit_bounds[1]]
        #check for NaN
        ff=np.where(~np.isfinite(ydata))
        if ff[0].size !=0:
            for loc in list(ff[0]):
              ydata.pop(loc)
              xv.pop(loc)
        xvmid=xv.index(0.0)

        #find peak theta of an ideal triangular distribution fit to the data
        fplus=np.polyfit(xv[:xvmid],ydata[:xvmid],1)
        fminus=np.polyfit(xv[xvmid:],ydata[xvmid:],1)
        fpslope=fplus[0]
        fpintercept=fplus[1]
        fpline=fpslope*np.array(xv[:xvmid+1]) + fpintercept
        fmslope=fminus[0]
        fmintercept=fminus[1]
        fmline=fmslope*np.array(xv[xvmid:]) + fmintercept
        lineintx,lineinty=ag.get_intersect(fpslope,fmslope,fpintercept,fmintercept)

        pkth=lineintx #peak theta

        #now find y-intercepts of triangular profile. Use the line fits for this!
        xnew,ypnew=self.interpolate_transm(xv[:xvmid+1],fpline,period,height,plot=False)
        xnew,ymnew=self.interpolate_transm(xv[xvmid:],fmline,period,height,plot=False)
        yintminus=xnew[(np.abs(ypnew - 0.0)).argmin()]
        yintplus=xnew[(np.abs(ymnew - 0.0)).argmin()]
        b=yintplus-yintminus #degrees

        #calculate new ideal profile using peak theta and b
        #ww.append((p/2.) - np.abs(h*np.tan(np.deg2rad(theta))))
        #w0=h*tan(b/2)
        xnew,ynew=self.interpolate_transm(xv,ydata,period,height,plot=False)
        dc=np.max(ynew)
        if fixed == 'height':
            w0=height*np.tan(np.deg2rad(b)/2.)
        elif fixed == 'pitch':
            #if w0 is fixed, then: height = w0/np.tan(np.deg2rad(b)/2.)
            height = (period*dc)/np.tan(np.deg2rad(b)/2.)
            w0=period*dc
        neww=(w0 - np.abs(height*np.tan(np.deg2rad((xnew-pkth)))))
        #print w0,b,pkth
        yvals=neww/(w0/dc)
        #plot everything
        if plot:
            fig,ax=plt.subplots()
            ax.scatter(xv,ydata)
            ax.plot(xnew[:50],ypnew[:50],'--k')
            ax.plot(xnew[50:],ymnew[50:], '--k')
            ax.plot(xnew, yvals, 'r')
            ax.set_ylim([0,1])
            yran=ax.get_ylim()
            ymax=float(yran[1])
            xran=ax.get_xlim()
            xmin=float(xran[0])
            xmax=float(xran[1])
            ax.text(.95*xmin,.9*ymax,'fit width: ' + str(np.round(w0,3))+ ' $\mu$m')
            ax.text(.95*xmin,.8*ymax,'fit pitch: ' + str(np.round(w0/dc,3))+ ' $\mu$m')
            ax.text(.95*xmin,.7*ymax,'fit height: ' + str(np.round(height,3))+ ' $\mu$m')
            ax.text(.95*xmin,.6*ymax,'fit dtheta: ' + str(np.round(pkth,3)) + ' degrees')
            fig.show()

        if ret:
            return xnew,yvals
        else: #save in results object - how do I want to visualize this though?
            valdict={'xnew':xnew, 'fit_width':w0,'fit_height':height,'peak_theta':pkth,'yvals':yvals}
            if fixed == 'height':
                try:
                    self.fit_width['p'+str(pindex)]=valdict
                except NameError:
                    self.fit_width={'p'+str(pindex):valdict}
            else:
                try:
                    self.fit_height['p'+str(pindex)]=valdict
                except NameError:
                    self.fit_height={'p'+str(pindex):valdict}

    def fit_profile_to_data_optimize(self,pindex,period,height,plot=True,ret=False,height_bounds=False, width_bounds=False,dc_bounds=False,yoff_bounds=False,fit_bounds=False,fix_dc=True,show=True,xvmid=False, oplot_nominal=True,soverwrite=False,sum_fac=False):
        '''fit a triangular transmission profile to the given one, leaving both height and width free. Return average width and tilt angle'''
        if not width_bounds:
            try:
                mpitch=self.period_stats['mean']/2.
                mdev=self.period_stats['stddev']/2.
                width_bounds=[mpitch-mdev,mpitch+mdev]
            except AttributeError:
                pass

        #get data
        if fit_bounds:
            ydata=list(self.data['ydata'][pindex])[fit_bounds[0]:fit_bounds[1]]
            xv=self.data['xdata'][fit_bounds[0]:fit_bounds[1]][:]
        else:
            ydata=list(self.data['ydata'][pindex])
            xv=self.data['xdata'][:]#[fit_bounds[0]:fit_bounds[1]]
        #print ydata
        if type(sum_fac)!=bool:
            ydata=list(sum_fac*np.array(ydata))
        #print ydata,sum_fac
        #check for NaN
        ff=np.where(~np.isfinite(ydata))
        if ff[0].size !=0:
            ydata = [j for i, j in enumerate(ydata) if i not in list(ff[0])]
            xv = [j for i, j in enumerate(xv) if i not in list(ff[0])]

        if not xvmid: #if xvmid is pre-defined, define it as 1+real Xmid index in order to include it
            xvmid=xv.index(0.0)
        #print xv[:xvmid],ydata[:xvmid]
        #print xv[xvmid:],ydata[xvmid:]
        #find peak theta of an ideal triangular distribution fit to the data. Don't fit the peak point
        fplus=np.polyfit(xv[:xvmid+1],ydata[:xvmid+1],1)
        fminus=np.polyfit(xv[xvmid+1:],ydata[xvmid+1:],1)
        fpslope=fplus[0]
        fpintercept=fplus[1]
        fpline=fpslope*np.array(xv[:xvmid+1]) + fpintercept
        fmslope=fminus[0]
        fmintercept=fminus[1]
        fmline=fmslope*np.array(xv[xvmid+1:]) + fmintercept
        lineintx,lineinty=ag.get_intersect(fpslope,fmslope,fpintercept,fmintercept)

        pkth=lineintx #peak theta

        #now find y-intercepts of triangular profile. Use the line fits for this!
        xpnew,ypnew=self.interpolate_transm(xv[:xvmid+1],fpline,period=period,height=height,plot=False) #this uses the ideal params...
        xmnew,ymnew=self.interpolate_transm(xv[xvmid+1:],fmline,period=period,height=height,plot=False) #this uses the ideal params...i think
        yintminusidx=(np.abs(ypnew - 0.0)).argmin()
        yintplusidx=(np.abs(ymnew - 0.0)).argmin()
        yintminus=xmnew[yintminusidx]
        yintplus=xpnew[yintplusidx]
        #print yintminus,ymnew[yintminus]
        #print yintplus,ypnew[yintplus]
        b=yintplus-yintminus #degrees

        #calculate new ideal profile using peak theta and b
        #ww.append((p/2.) - np.abs(h*np.tan(np.deg2rad(theta))))
        #w0=h*tan(b/2)
        xnew,ynew=self.interpolate_transm(xv,ydata,yints=[yintminus,yintplus],plot=False)
        #print xnew[0],ynew[0]
        #print xnew[-1],ynew[-1]
        #trim to wher y>=0
        #aa=np.where(ynew[:50] >=0.0)
        #bb=np.where(ynew[50:] >=0.0)
        #print bb[0]
        #xnew=xnew[yintminusidx:yintplusidx]
        #ynew=ynew[yintminusidx:yintplusidx]
        #print len(xnew),len(ynew)
        dc=np.max(ynew)

        def func(xnew,w0,h,dc,yoff): #can I fit the duty cycle too?
            return yoff+(w0-np.abs(h*np.tan(np.deg2rad((xnew-pkth)))))/(w0/dc)

        if fix_dc: #fix dc to equal 0.5 and adjust profiles to be consistent with this
            ynew=ynew-(dc-0.5)
            dc=0.5

        if not width_bounds:
            width_bounds=[0,np.inf]
        if not height_bounds:
            height_bounds=[0,np.inf]
        if not dc_bounds:
            dc_bounds=[.48,.52]
        if not yoff_bounds:
            yoff_bounds=[-.5,.2]
        #print [width_bounds[0],height_bounds[0]],[width_bounds[1],height_bounds[1]]
        popt, pcov=curve_fit(func,xnew,ynew,bounds=([width_bounds[0],height_bounds[0],dc_bounds[0],yoff_bounds[0]],[width_bounds[1],height_bounds[1],dc_bounds[1],yoff_bounds[1]]),absolute_sigma=True)


        neww=func(xnew,*popt)#(w0 - np.abs(height*np.tan(np.deg2rad((xnew-pkth)))))
        #print w0,b,pkth
        yvals=neww#/(w0/dc)
        w0=popt[0]
        height=popt[1]
        dc=popt[2]
        yoff=popt[3]

        #calculate error from difference of data to fit values...sqrt (data-fit)^2
        #error = []
        #perr = np.sqrt(np.diag(pcov)) yvals...
        newpoints=func(xv,*popt)
        error=np.mean(np.abs(ydata-newpoints))
        #print newpoints
        #print ydata
        #print error

        if oplot_nominal: #overplot nominal profile
            wnom=self.nominal['pitch']/2.
            hnom=self.nominal['height']
            ynomvals=func(xnew, wnom,hnom,0.5,yoff)


        #plot everything
        if plot:
            fig,ax=plt.subplots()
            ax.scatter(xv,ydata, label='data')
            ax.plot(xnew,ynew,'--k')
            ax.plot(xmnew,ymnew, '--g', label='- fit')
            ax.plot(xpnew[:50],ypnew[:50], '--b', label='+ fit')
            ax.plot(xnew, yvals, 'r', label= 'fit')
            if oplot_nominal:
                ax.plot(xnew, ynomvals, '--m',label='ideal')
            ax.set_ylim([0,1])
            try:
                ax.set_xlim([xmnew[yintminusidx-5],xpnew[yintplusidx+5]])
            except IndexError:
                pass
            yran=ax.get_ylim()
            ymax=float(yran[1])
            xran=ax.get_xlim()
            xmin=float(xran[0])
            xmax=float(xran[1])
            ax.text(.95*xmin,.9*ymax,'fit width: ' + str(np.round(w0,3))+ ' $\mu$m')
            ax.text(.95*xmin,.8*ymax,'fit pitch: ' + str(np.round(w0/dc,3))+ ' $\mu$m')
            ax.text(.95*xmin,.7*ymax,'fit duty cycle: ' + str(np.round(dc,3)))
            ax.text(.95*xmin,.6*ymax,'fit height: ' + str(np.round(height,3))+ ' $\mu$m')
            ax.text(.95*xmin,.5*ymax,'fit dtheta: ' + str(np.round(pkth,3)) + ' degrees')
            ax.text(.95*xmin,.4*ymax,'fit yoffset: ' + str(np.round(yoff,3)) )
            ax.legend(loc='upper right')
            if show==True:
                fig.show()
            else:
                fig.savefig(show+'.png')
        if ret:
            return xnew,yvals
        else: #save in results object - how do I want to visualize this though?
            valdict={'xnew':xnew, 'fit_width':w0,'fit_height':height,'fit_period':w0/dc,'peak_theta':pkth,'yvals':yvals,'dc':dc,'yoff':yoff,'pitch_nom':period,'error':error}
            bdict={'height_bounds':height_bounds,'width_bounds':width_bounds,'dc_bounds':dc_bounds,'fit_bounds':fit_bounds}
            try:
                self.fit_both['p'+str(pindex)]=valdict
                self.fit_both['bounds']=bdict
                #print 'written'
            except AttributeError:
                self.fit_both={'p'+str(pindex):valdict,'bounds':bdict}
            #    print 'writenew'

    def print_height_stats(self):
        fw,fh,fp,dc,dt,yoff=[],[],[],[],[],[]
        for p in range(0,7):
            pkey='p'+str(p)
            fw.append(self.fit_both[pkey]['fit_width'])
            fh.append(self.fit_both[pkey]['fit_height'])
            fp.append(self.fit_both[pkey]['fit_period'])
            dc.append(self.fit_both[pkey]['dc'])
            dt.append(self.fit_both[pkey]['peak_theta'])
            yoff.append(self.fit_both[pkey]['yoff'])
        bounds=self.fit_both['bounds']
        print "########### Height Stats for win",str(self.win),"##############"
        print "Fit width: ",np.mean(fw)
        print "Fit height: ",np.mean(fh)
        print "Fit period: ",np.mean(fp)
        print "Fit duty cycle: ",np.mean(dc)
        print "Fit dtheta: ",np.mean(dt)
        print "Fit yoffset: ",np.mean(yoff)
        print "Bounds: ",bounds

    def transm_func(self,xnew,w0,h,dc,yoff,pkth=0.): #can I fit the duty cycle too?
        return yoff+(w0-np.abs(h*np.tan(np.deg2rad((xnew-pkth)))))/(w0/dc)

    def calc_nominal_transm(self, th_vec=False,th_offsets=False):
        '''calculate the transmission profile given ideal conditions'''
        if "nominal_transm_prof" not in self.__dict__.keys():
            self.nominal_transm_prof = {}
        if not th_vec:
            th_vec=self.data['xdata']
        for pindex in range(0,7):
            if not th_offsets:
                th_offset=self.fit_both['p'+str(pindex)]['peak_theta']
            yvals=self.transm_func(th_vec,self.nominal['pitch']*self.nominal['slit_width'],self.nominal['height'],self.nominal['slit_width'],0.,pkth=th_offset) #no y-offset
            self.nominal_transm_prof['p'+str(pindex)]=yvals

    def zip_fit_vals(self,indict,ykey):
        '''stack individual vectors into arrays'''
        yvec=[list(indict['p0'][ykey])]
        for p in range(1,6):
            yvec.append(list(indict['p'+str(p)][ykey]))
        print np.shape(yvec)
        return yvec

    def period_vs_angle(self,pix2um=.65,yran=False, fit=True,plot=True):
        rpp=self.data['raw_periods_A'] #plot errorbars...
        mpp=self.data['mean_periods_P']
        xv=np.array(self.data['xdata']) #angles
        labels=['p0','p1','p2','p3','p4','p5','p6']


            #ax.(xv,np.ma.array(mp)*pix2um,label=labels[i])
        if fit: #do polyfit of a quadratic to the data
            z=np.polyfit(xv,np.ma.mean(mpp,axis=0),2)
            fityvals=np.poly1d(z)#a*xv**2. + b*xv + c
            #print the vertex x=-b/2a

        if plot:
            fig,ax=plt.subplots()
            for i,mp in enumerate(mpp):
                yerr=[np.ma.std(rpp[i][j])*pix2um for j in range(0,len(mp))]
                ax.errorbar(xv,np.ma.array(mp)*pix2um,fmt='o',yerr=np.array(yerr),label=labels[i])
            ax.plot(xv,fityvals(xv)*pix2um, '--k', label= 'best fit')
            print 'vertex : (', np.rad2deg(-1.*z[1]/(2.*z[0])),' degrees,', np.max(fityvals(xv))*pix2um, ' um)'
            ax.legend(loc='upper left')
            ax.set_xlim([xv[0]-.5,xv[-1]+.5])
            ax.set_xlabel('tan $\\theta$')
            ax.set_ylabel('Period ($\mu$m)')

            if yran:
                ax.set_ylim(yran)
            fig.show()

        return fityvals(xv) #in pixels

    def period_hist(self,adjust_p0=False,pix2um=.65,stats=True):
        '''if adjust_p0, use calc_p0 to fix'''
        rpp=self.data['raw_periods_A'] #plot errorbars...
        #mpp=self.data['mean_periods_P']
        #xv=np.array(self.data['xdata']) #angles
        #labels=['p0','p1','p2','p3','p4','p5','p6']
        periods=rpp.flatten() #this is now 49x?
        flatperiods=[]

        for p in periods:
            pl=list(p)
            if flatperiods !=[]:
                flatperiods.extend(pl)
            else:
                flatperiods=pl

        flatperiods=pix2um*np.array(flatperiods)
        flatperiods.sort()
        aa=np.where(np.isnan(flatperiods))
        if aa != np.array([]):
            flatperiods=flatperiods[:aa[0][0]]
        bins=np.arange(np.nanmin(flatperiods),np.nanmax(flatperiods),np.nanstd(flatperiods)/5.)
        fig,ax=plt.subplots()
        ax.hist(flatperiods,bins)
        ax.set_yscale('log')
        ax.set_xlim([.9*self.nominal['pitch'],1.1*self.nominal['pitch']])
        ax.set_xlabel('tan $\\theta$')
        fig.show()

        stats={"num":np.shape(flatperiods),"mean":np.nanmean(flatperiods),"med":np.nanmedian(flatperiods),"stddev":np.nanstd(flatperiods)}
        if stats: #print
            print "-------------------STATISTICS FOR WINDOW "+str(self.win)+"---------------------"
            print '              Mean P: ' + str(stats['mean'])
            print '            Median P: ' + str(stats['med'])
            print 'Standard Deviation P: ' + str(stats['stddev'])
            print '               Num P: ' + str(stats['num'])
        self.period_stats=stats


    def calc_p0(self,recalc_ydata=False):
        ''' calulate period at real theta=0 using parabolic fit'''
        fyvals=self.period_vs_angle(fit=True,plot=False)
        #th0=
        p0=np.max(fyvals)

        #put it in the object.... later use it to re-plot and fit the profile
        self.data['period_0']=p0
        if recalc_ydata: #recalculate w_m/p0 using new p0
            ynew=self.data['mean_widths_P']/p0
            self.data['ydata']=np.array(ynew)

######################### make TeX table entry ##############################

    def opt_stats_tex(self, do_print=True):
        ''' for the results table. win | npitch | pmean | pmean - npitch | pdev | preq | nang | amean | amean -nang | adev | areq'''
        def rv(value):
            return str(np.round([value],decimals=4)[0])
        win=str(self.win)
        side=self.kind['side']
        npitch=self.nominal['pitch']
        try:
            nang=self.nominal['nominal angle']
        except KeyError:
            nang=self.nominal['orientation']
        req = self.requirements
        pstats=self.period_stats
        astats=self.angle_stats
        if type(astats['mean']) == np.ma.core.MaskedArray:
            amean=float(astats['mean'].data)
            #adev=float(astats['stddev'].data)
        else:
            amean=astats['mean']
        adev=astats['stddev']
        line=win + ' & ' + str(npitch) +' & '+rv(pstats['amean'])+' & '+rv(pstats['amean']-npitch)+' & '+rv(pstats['adev'])+' & '+ rv(req['pitch'])+ ' & '+ str(nang) +' & '+rv(amean)+' & '+rv(np.abs(amean)-np.abs(nang))+' & '+rv(adev)+' & '+ rv(req['angle'])+' \\\\'
        self.tex['stats']=line
        if do_print:
            print line

    def opt_params_tex(self,do_print=True):
        ''' for the params table'''
        def rv(value):
            return str(np.round([value],decimals=4)[0])

        win=str(self.win)
        params=self.calc_params
        for p in params.keys():
            params[p]=str(params[p])

        line=win + ' & ' + params['csigma'] +' & '+ params['ll'] +' & '+params['spread']+' & '+params['n_hough']+' & '+ params['htol']+ ' & '+ params['sigma'] +' & '+params['tol']+' & '+ params['r_or_f_tol'] +' \\\\'
        self.tex['params']=line
        if do_print:
            print line

    #do the same for X-ray results
    def xray_stats_tex(self, do_print=True):
        ''' for the results table. win | npitch | pmean | pmean - npitch | pdev | preq | nang | amean | amean -nang | adev | areq'''
        def rv(value):
            return str(np.round([value],decimals=4)[0])
        win=str(self.win)
        side=self.kind['side']
        npitch=self.nominal['pitch']
        try:
            nang=self.nominal['nominal angle']
        except KeyError:
            nang=self.nominal['orientation']
        #req = self.requirements
        pstats=self.period_stats
        astats=self.angle_stats
        th0=astats['theta0']
        amean=np.nanmean(astats['theta'])
        adev=astats['theta_std']

        line=win + ' & ' + str(npitch) +' & '+rv(pstats['mean'])+' & '+rv(pstats['mean']-npitch)+' & '+rv(pstats['stddev'])+ ' & '+ str(nang) +' & '+rv(amean) +' & '+str(th0)+' & '+rv(adev)+' \\\\'
        self.tex['stats']=line
        if do_print:
            print line

    def xray_fit_stats_tex(self, do_print=True):
        ''' for the results table. win | npitch | fit pitch | fit width | fit dc | height '''
        def rv(value):
            return str(np.round([value],decimals=4)[0])
        win=str(self.win)
        side=self.kind['side']
        npitch=self.nominal['pitch']
        try:
            nang=self.nominal['nominal angle']
        except KeyError:
            nang=self.nominal['orientation']
        #req = self.requirements
        pstats=self.period_stats
        astats=self.angle_stats
        amean=astats['theta']
        adev=astats['theta_std']
        line=win + ' & ' + str(npitch) +' & '+rv(pstats['mean'])+' & '+rv(pstats['mean']-npitch)+' & '+rv(pstats['stddev'])+ ' & '+ str(nang) +' & '+rv(amean)+' & '+rv(np.abs(amean)-np.abs(nang))+' & '+rv(adev)+' \\\\'
        self.tex['stats']=line
        if do_print:
            print line

######################PLOT FUNCTIONS##########################

    def scat_allP(self,meanfit=True,figname=False,theta_ran=False,print_height=False,yran=False,yran2=False,fit=False,compare_nominal=False,widths_only=False):
        win=self.win
        if fit == False:
            xv=self.data['xdata']
            if widths_only:
                sumvecs=self.data['mean_widths_A']
            else:
                sumvecs=self.data['ydata']
            xvmid=xv.index(0.0)
            lnsty='o-'
        elif fit == 'height': #plot the fits
            xv=list(self.fit_height['p0']['xnew'])
            sumvecs=self.zip_fit_vals(self.fit_height,'yvals')
            xvmid= 50
            lnsty='-'
        elif fit == 'width':
            xv=list(self.fit_width['p0']['xnew'])
            sumvecs=self.zip_fit_vals(self.fit_width,'yvals')
            xvmid = 50
            lnsty='-'
        elif fit == 'both':
            xv=self.data['xdata']
            if widths_only:
                sumvecs=self.data['mean_widths_A']
            else:
                sumvecs=self.data['ydata']
            xvmid=xv.index(0.0)
            lnsty='o-'
            compare_nominal=True

        if compare_nominal:
            xv2=list(self.fit_both['p0']['xnew'])
            sumvecs2=self.zip_fit_vals(self.fit_both,'yvals')
            xvmid2=50

        fitsplus=[]
        fitsminus=[]
        errvecs=[]

        meanvec=np.nanmean(sumvecs,axis=0)
        fplus=np.polyfit(xv[:xvmid],meanvec[:xvmid],1) #formerly xvmid -1
        fminus=np.polyfit(xv[xvmid+1:],meanvec[xvmid+1:],1)
        #print np.min(svsort), np.max(svsort),fplus,fminus
        fitsplus.append(fplus)
        fitsminus.append(fminus)

        labels=['p0','p1','p2','p3','p4','p5','p6']
        cols=['b','c','m','k','r','g','y']
        fig,ax=plt.subplots(2,sharex=True,figsize=[10,8])
        #ax=fig.add_subplot(112)
        for sv,c,lab in zip(sumvecs,cols,labels):
            #if len(xv) <12:
            ax[0].plot(xv,sv,lnsty,color=c,label=lab)
            #else:
            #    ax[0].plot(xv,sv,lnsty,color=c,label=lab)

        ax[0].axvline(color='k',linestyle='dashed')

        try:
            for sv,c in zip(sumvecs2,cols):
                ax[0].plot(xv2,sv,'--',color=c)
            #ax[0].axvline(color='k',linestyle='dashed')
        except NameError:
            pass

        if not yran:
            ax[0].set_ylim([0,1.1])
        else:
            ax[0].set_ylim(yran)
        if meanfit:
            meanplusslope=np.nanmean([fp[0] for fp in fitsplus])
            meanplusintercept=np.nanmean([fp[1] for fp in fitsplus])
            meanplusline=meanplusslope*np.array(xv[:xvmid+1]) + meanplusintercept
            #print meanplusline
            meanminusslope=np.nanmean([fm[0] for fm in fitsminus])
            meanminusintercept=np.nanmean([fm[1] for fm in fitsminus])
            meanminusline=meanminusslope*np.array(xv[xvmid:]) + meanminusintercept
            lineintx,lineinty=ag.get_intersect(meanplusslope,meanminusslope,meanplusintercept,meanminusintercept)
            #print meanminusline
            ax[0].plot(xv[:xvmid+1],meanplusline,'k--')
            ax[0].plot(xv[xvmid:],meanminusline,'k--')
            yran=ax[0].get_ylim()
            ymax=float(yran[1])
            xran=ax[0].get_xlim()
            xmin=float(xran[0])
            xmax=float(xran[1])
            ax[0].text(.75*xmin,.2*ymax,'avg. slope: ' + str(np.round(meanplusslope,3)))
            ax[0].text(.75*xmin,.1*ymax,'avg. intercept: ' + str(np.round(meanplusintercept,3)))
            ax[0].text(.5*xmax,.2*ymax,'avg. slope: ' + str(np.round(meanminusslope,3)))
            ax[0].text(.5*xmax,.1*ymax,'avg. intercept: ' + str(np.round(meanminusintercept,3)))
            ax[0].text(.15*xmin,.3*ymax,'adj. peak: (' + str(np.round(lineintx,3))+ ',' + str(np.round(lineinty,3))+')')

        #ax.set_ylim([0,1.1])
        #ax.set_xlim([xv[0],xv[-1]])
        ax[1].set_xlabel('Angle (degrees)') #should I set the lables to the actual angle values? xaxlabels... but they would be in the wrong place
    #     if widths:
    #        ax[0].set_title('Window '+ win + ' slit width as a function of angle')
    #        ax[0].set_ylabel('Slit width ($\mu$m)')
    #    else:
        ax[0].set_title('Window '+ str(win) + ' analyzed with '+ self.method + ' method')
        if widths_only:
            ax[0].set_ylabel('Slit Width ($\mu$m)')
        else:
            ax[0].set_ylabel('Percent of maximum transmission')
        ax[0].legend(loc='upper right')

        #make the lower subplot showing the error bars
        #ax2=fig.add_subplot(121)
        if meanfit:
            #generate the points for the fit lines
            meanpluspoints=(meanplusslope*np.array(xv[:xvmid+1]) + meanplusintercept).tolist()
            meanminuspoints=(meanminusslope*np.array(xv[xvmid+1:]) + meanminusintercept).tolist()
            meanpluspoints.extend(meanminuspoints)
            #print meanpluspoints
            #should plot the error bars here
            for sv,c,lab in zip(sumvecs,cols,labels):
                yv=sv-np.array(meanpluspoints)
                errvecs.append(yv)
                ax[1].plot(xv,yv,lnsty,color=c,label=lab)
                #ax[1].errorbar(xv,yv,yerr=err)
            ax[1].axhline(color='k',linestyle='dashed')
            if yran2:
                ax[1].set_ylim(yran2)
            else:
                ax[1].set_ylim([np.ma.min(errvecs),np.ma.max(errvecs)])

        if compare_nominal: #compare with the nominal profiles that were calculated
            sumvecs_nom=[self.nominal_transm_prof['p'+str(p)] for p in range(0,7)]
            ax[1].set_ylabel('Nominal Profile Difference')
            for sv,sv_nom,c,lab in zip(sumvecs,sumvecs_nom,cols,labels):
                yv=sv-np.array(sv_nom)
                errvecs.append(yv)
                ax[1].plot(xv,yv,'o-',color=c,label=lab)

        #re-adjust plot sizes
        box1=ax[0].get_position()
        box2=ax[1].get_position()
        ax[0].set_position([box1.x0,box1.y0*.7,box1.width,box1.height*1.5])
        ax[1].set_position([box2.x0,box2.y0,box2.width,box2.height*.7])
        if figname:
            plt.savefig(figname)
            #data=[xv,sumvecs,fitsplus,fitsminus,meanplusslope,meanplusintercept,meanminusslope,meanminusintercept,lineintx,lineinty]
            #pickle.dump(data,open(figname[:-4]+'.p','wb'))
        else:
            fig.show()
        if print_height:
            #print allfitsplus,allfitsminus
            for fp,fm,i in zip(allfitsplus,allfitsminus,range(0,len(allfitsplus))):
                heights,meanheight=calc_height_from_fit(xv,fp[0],fp[1],fm[0],fm[1])
                print 'p',i,' mean height: ', meanheight
        #return xv,sumvecs,errvecs
        try:
            if not self.errors:
                self.errors=errvecs
        except AttributeError:
            self.errors=errvecs

    def fit_hists(width=True, height=False, dtheta=False):
        '''plot histograms of fit parameters'''
        if width:
            yvals='foo'

    def fit_plots(self,figname=False,bounds=False,saveplot=False,ret=False):
        '''plot fit parameters, with errors'''
        xvals = range(0,7)
        pitch=self.fit_both['p0']['pitch_nom']
        fig,ax=plt.subplots(2,3, figsize=[12,12],sharex=True)
        width_yvals=[self.fit_both['p'+str(k)]['fit_width'] for k in xvals]
        ax[0][0].plot(xvals,np.array(width_yvals)-pitch/2.,'o-')
        ax[0][0].set_title('fit_width -nominal (um)')
        yran=ax[0][0].get_ylim()
        ymax=float(yran[-1])
        ax[0][0].text(1,.8*ymax,'mean width: ' + str(np.round(np.mean(width_yvals),3))+ ' $\mu$m')


        period_yvals=[self.fit_both['p'+str(k)]['fit_period'] for k in xvals]
        ax[0][1].plot(xvals,np.array(period_yvals)-pitch,'o-')
        ax[0][1].set_title('fit period -nominal (um)')
        yran=ax[0][1].get_ylim()
        ymax=float(yran[-1])
        ax[0][1].text(1,.8*ymax,'mean period: ' + str(np.round(np.mean(period_yvals),3))+ ' $\mu$m')

        height_yvals=[self.fit_both['p'+str(k)]['fit_height'] for k in xvals]
        #hcov=[self.fit_both['p'+str(k)]['error'][1] for k in xvals]
        #print hcov
        #ax[0][2].errorbar(xvals,height_yvals,yerr=np.array(hcov))
        ax[0][2].plot(xvals,height_yvals,'o-')
        yran=ax[0][2].get_ylim()
        ymin=float(yran[0])
        ax[0][2].text(1,ymin,'mean height: ' + str(np.round(np.mean(height_yvals),3))+ ' $\mu$m')
        #ax[0][2].errorbar(xvals,hcov)
        ax[0][2].set_title('fit_height (um)')

        th_yvals=[self.fit_both['p'+str(k)]['peak_theta'] for k in xvals]
        #th_yvals2=[self.fit_height['p'+str(k)]['peak_theta'] for k in xvals]
        ax[1][0].plot(xvals,th_yvals,'o-')
        #ax[2].plot(xvals,th_yvals2,'o-',label='fixed width')
        ax[1][0].set_title('dtheta (degrees)')

        dc=[self.fit_both['p'+str(k)]['dc'] for k in xvals]
        dccov=[self.fit_both['p'+str(k)]['error'][2] for k in xvals]
        #print dccov
        #ax[1][1].errorbar(xvals,dc,yerr=np.array(dccov))
        #ax[2].plot(xvals,th_yvals2,'o-',label='fixed width')
        ax[1][1].plot(xvals,dc,'o-')
        yran=ax[1][1].get_ylim()
        ymin=float(yran[0])
        ax[1][1].text(1,ymin,'mean duty cycle: ' + str(np.round(np.mean(dc),3)))
        ax[1][1].set_title('Duty Cycle')

        ycov=[self.fit_both['p'+str(k)]['error'][3] for k in xvals]
        yoff=[self.fit_both['p'+str(k)]['yoff'] for k in xvals]
        #print ycov
        #ax[1][2].errorbar(xvals,yoff,yerr=np.array(ycov))
        ax[1][2].plot(xvals,yoff,'o-')
        yran=ax[1][2].get_ylim()
        ymax=float(yran[-1])
        ax[1][2].text(1,.8*ymax,'mean y-offset: ' + str(np.round(np.mean(yoff),3)))

        #ax[2].plot(xvals,th_yvals2,'o-',label='fixed width')
        ax[1][2].set_title('Y-offset')

        if bounds: #set all the boundaries
            ax[0][0].set_ylim(np.array(bounds[0])-pitch/2.)
            #ax[0][1].set_ylim(np.array(bounds[1])-pitch)
            ax[0][2].set_ylim(bounds[2])
            ax[1][0].set_ylim(bounds[3])
            ax[1][1].set_ylim(bounds[4])
            ax[1][2].set_ylim(bounds[5])

        #ax[2].legend()
        for a in ax[1]:
            a.set_xlabel('p')
        for a in ax[0]:
            a.set_xlabel('p')
        fig.suptitle('Window ' +str(self.win))
        if figname:
            plt.savefig(figname)
        #if saveplot !=False: #save plot data
            pickle.dump({'xvals':xvals,'pitch_nom':pitch,'width_yvals':width_yvals,'period_yvals':period_yvals,'height_yvals':height_yvals,'peak_theta':th_yvals,'dc':dc,'yoff':yoff,'bounds':bounds},open(saveplot+'.p','wb'))
        else:
            fig.show()

        vdict={'xvals':xvals,'pitch_nom':pitch,'width_yvals':width_yvals,'period_yvals':period_yvals,'height_yvals':height_yvals,'peak_theta':th_yvals,'dc':dc,'yoff':yoff,'bounds':bounds}
        if ret:
            return vdict


    def plot_angle_hist(self,mask45=False,imtag=False):
        '''what it sounds like. for optical. wrapper for ag.hough_hist'''
        win=self.win
        try:
            nang=self.nominal['nominal angle']
        except KeyError:
            nang=self.nominal['orientation']
        htol=int(self.calc_params['htol'])
        spread=float(self.calc_params['spread'])
        folder=self.kind['folder']
        if not imtag:
            figname='win'+str(win) + 'ang_hist.png'
        else:
            figname='win'+str(win) + 'ang_hist'+imtag+'.png'
        os.chdir(folder)
        hough_all=glob.glob('win'+str(win)+'*_5.0X_corrected_hough_all.p')
        ag.hough_hist(hough_all,win,nang,figname=figname,tol=htol,spread=spread,mask45=mask45)
        os.chdir('../..')


    def plot_period_hist(self,key='slat',filtered=True, rfp=False,nbins=False,xran=False):
        '''what it sounds like. for optical'''
        #from analyze_optical
        #make the histogram
        npitch=self.nominal['pitch']
        if not nbins:
            bins=np.arange(npitch-npitch/10., npitch+npitch/10., npitch/40.)
        else:
            bins=np.linspace(.9*npitch,1.1*npitch,nbins)
        if key == 'slit' or key == 'slat':
            if filtered:
                dkey='filtered_periods_'+key
            else:
                dkey='periods_from_'+key+'_centers'
        else:
            dkey=key

        periods=self.data[dkey]

        if rfp:
            periodr=self.data['rising']
            periodf=self.data['falling']
            from_c=self.data['periods_from_slat_centers']
            from_c2=self.data['periods_from_slit_centers']
            fig,ax=plt.subplots(1,3,sharex=True,sharey=True)
            ax[0].hist(periodr,bins,facecolor='g')
            ax[1].hist(periodf,bins,facecolor='r')
            ax[2].hist(from_c,bins,facecolor='b',alpha=0.6)
            ax[2].hist(from_c2,bins,facecolor='m',alpha=0.6)

            #ax[2].hist(filtered_c,bins,facecolor='b',alpha=0.6)
            #ax[2].hist(filtered_c2,bins,facecolor='m',alpha=0.6)
            ax[0].set_xlim([bins[0],bins[-1]])
            ax[0].set_title('Rising')
            ax[1].set_title('Falling')
            ax[2].set_title('Win'+str(self.win)+' Total')
            ax[0].set_xlabel('Period $\mu$m')
            ax[0].set_ylabel('Counts')
            ax[0].set_yscale('log')
            ax[0].set_ylim([1,10000])
            #figfilename='win'+str(win)+'_group_periods5.0X'+ftags+'.png'
            #fig.savefig(figfilename)
            fig.show()

        else: #just one histogram
            fig,ax=plt.subplots()
            ax.hist(periods,bins,facecolor='g')
            if not xran:
                ax.set_xlim([bins[0],bins[-1]])
            else:
                ax.set_xlim(xran)
            ax.set_title(str(self.win)+ ' ' + dkey)
            ax.set_xlabel('Period $\mu$m')
            ax.set_ylabel('Counts')
            ax.set_yscale('log')
            ax.set_ylim([1,10000])
            #figfilename='win'+str(win)+'_group_periods5.0X'+ftags+'.png'
            #fig.savefig(figfilename)
            fig.show()
