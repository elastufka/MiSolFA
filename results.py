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

class Results():

    def __init__(self,win,kind,method,errors,data,nominal=False,filenames=False,stats=False,params=False):
        '''Make an object to hold results'''
        self.win=win
        self.kind=kind
        self.method=method
        self.errors=errors
        self.data=data
        self.nominal=nominal
        if stats:
            self.stats=stats
        else:
            self.stats={} # get_stats  - do stats on data
        if filenames:
            self.filenames=filenames
        else:
            self.filenames={} # gen_filenames from data
        if params:
            self.calc_params=params

    def reload_widthdata_from_files(self):
        #{'xdata':xvec,'ydata':ydata,'mean_ang':ydataa, \
        #                                             'raw_widths_P':sumvecp,'raw_widths_A':sumveca, \
        #                                            'raw_periods_P':sumpp,'raw_periods_A':sumpa, \
        #                                            'mean_widths_P':meanvecp,'mean_widths_A':meanveca, \
        #                                            'mean_periods_P':meanpp,'mean_periods_A':meanpa}

        #basic check
        if self.method != 'widths':
            return

        #build filenames
        xvec=self.data['xdata']
        xstr=['{:g}'.format(xv) for xv in xvec]
        print xstr
        dfiles=['win'+str(self.win)+'_width_data_p'+str(pp)+'_'+xv+'.p' for pp in range(0,7) for xv in xstr]
        sfiles=['win'+str(self.win)+'_width_stats_p'+str(pp)+'_'+xv+'.p' for pp in range(0,7) for xv in xstr]
        print dfiles[:10]
        print sfiles[:10]
        newparams={}
        allwidths,allperiods=[],[]
        for d,s in zip(dfiles,sfiles):
            data=pickle.load(open(d,'rb'))
            aa=pickle.load(open(d,'rb'))
            try:
                sig=aa[1]
                n=aa[2]
            except KeyError:
                try:
                    sig=self.calc_params['sigma']
                    n=self.calc_params['tolvec']
                except KeyError:
                    sig,n=self.calc_params[s]

            newparams[s]=[sig,n]
            pp=data['period']
            ww=data['widths']

            allwidths.append(np.array(ww))
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
        self.calc_params=newparams


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

    def fit_profile_to_data_optimize(self,pindex,period,height,plot=True,ret=False,height_bounds=False, width_bounds=False,dc_bounds=False,yoff_bounds=False,fit_bounds=[1,-1],fix_dc=True,show=True,xvmid=False):
        '''fit a triangular transmission profile to the given one, leaving both height and width free. Return average width and tilt angle'''
        #get data
        ydata=list(self.data['ydata'][pindex])[fit_bounds[0]:fit_bounds[1]]
        xv=self.data['xdata'][fit_bounds[0]:fit_bounds[1]]
        #check for NaN
        ff=np.where(~np.isfinite(ydata))
        if ff[0].size !=0:
            for loc in list(ff[0]):
              ydata.pop(loc)
              xv.pop(loc)
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
        b=yintplus-yintminus #degre,\es

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

        if fix_dc: #fix dc to equal 0.5 and adjust profiles to be consistent with this
            ynew=ynew-(dc-0.5)
            dc=0.5


        #if fixed == 'height': #get rid of this and use scipy.optimize.curve_fit to fit both height and w0
        #    w0=height*np.tan(np.deg2rad(b)/2.)
        #elif fixed == 'pitch':
        #    #if w0 is fixed, then: height = w0/np.tan(np.deg2rad(b)/2.)
        #    height = (period*dc)/np.tan(np.deg2rad(b)/2.)
        #    w0=period*dc

        def func(xnew,w0,h,dc,yoff): #can I fit the duty cycle too?
            return yoff+(w0-np.abs(h*np.tan(np.deg2rad((xnew-pkth)))))/(w0/dc)

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

        #calculate error from pcov
        error = []
        for i in range(len(popt)):
            try:
                error.append(np.absolute(pcov[i][i])**0.5)
            except:
                error.append( 0.00 )

        neww=func(xnew,*popt)#(w0 - np.abs(height*np.tan(np.deg2rad((xnew-pkth)))))
        #print w0,b,pkth
        yvals=neww#/(w0/dc)
        w0=popt[0]
        height=popt[1]
        dc=popt[2]
        yoff=popt[3]
        #plot everything
        if plot:
            fig,ax=plt.subplots()
            ax.scatter(xv,ydata)
            ax.plot(xnew,ynew,'--k')
            ax.plot(xmnew,ymnew, '--g')
            ax.plot(xpnew[:50],ypnew[:50], '--b')
            ax.plot(xnew, yvals, 'r')
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
            if show==True:
                fig.show()
            else:
                fig.savefig(show+'.png')
        if ret:
            return xnew,yvals
        else: #save in results object - how do I want to visualize this though?
            valdict={'xnew':xnew, 'fit_width':w0,'fit_height':height,'fit_period':w0/dc,'peak_theta':pkth,'yvals':yvals,'dc':dc,'yoff':yoff,'pitch_nom':period,'error':error}
            try:
                self.fit_both['p'+str(pindex)]=valdict
            except AttributeError:
                self.fit_both={'p'+str(pindex):valdict}

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

    def scat_allP(self,meanfit=True,figname=False,theta_ran=False,print_height=False,yran=False,fit=False,compare_nominal=False):
        win=self.win
        if fit == False:
            xv=self.data['xdata']
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

        meanvec=np.ma.mean(sumvecs,axis=0)
        fplus=np.polyfit(xv[:xvmid-1],meanvec[:xvmid-1],1)
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
            meanplusslope=np.mean([fp[0] for fp in fitsplus])
            meanplusintercept=np.mean([fp[1] for fp in fitsplus])
            meanplusline=meanplusslope*np.array(xv[:xvmid+1]) + meanplusintercept
            #print meanplusline
            meanminusslope=np.mean([fm[0] for fm in fitsminus])
            meanminusintercept=np.mean([fm[1] for fm in fitsminus])
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
        ax[1].set_xlabel('tan(Angle)') #should I set the lables to the actual angle values? xaxlabels... but they would be in the wrong place
    #     if widths:
    #        ax[0].set_title('Window '+ win + ' slit width as a function of angle')
    #        ax[0].set_ylabel('Slit width ($\mu$m)')
    #    else:
        ax[0].set_title('Window '+ str(win) + ' analyzed with '+ self.method + ' method')
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
        if not self.errors:
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


