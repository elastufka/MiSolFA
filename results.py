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

    def __init__(self,win,kind,method, errors, data,filenames=False,stats=False):
        '''Make an object to hold results'''
        self.win=win
        self.kind=kind
        self.method=method
        self.errors=errors
        self.data=data
        if stats:
            self.stats=stats
        else:
            self.stats={} # get_stats  - do stats on data
        if filenames:
            self.filenames=filenames
        else:
            self.filenames={} # gen_filenames from data

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

    def scat_allP(self,meanfit=True,figname=False,theta_ran=False,print_height=False,yran=False):
        win=self.win
        xv=self.data['xdata']
        sumvecs=self.data['ydata']
        fitsplus=[]
        fitsminus=[]
        errvecs=[]
        xvmid=xv.index(0.0)
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
            ax[0].plot(xv,sv,'-o',color=c,label=lab)
        ax[0].axvline(color='k',linestyle='dashed')
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
                ax[1].plot(xv,yv,'o-',color=c,label=lab)
                #ax[1].errorbar(xv,yv,yerr=err)
            ax[1].axhline(color='k',linestyle='dashed')
            ax[1].set_ylim([np.ma.min(errvecs),np.ma.max(errvecs)])

            #re-adjust plot sizes
            box1=ax[0].get_position()
            box2=ax[1].get_position()
            ax[0].set_position([box1.x0,box1.y0*.7,box1.width,box1.height*1.5])
            ax[1].set_position([box2.x0,box2.y0,box2.width,box2.height*.7])
        if figname:
            plt.savefig(figname)
            data=[xv,sumvecs,fitsplus,fitsminus,meanplusslope,meanplusintercept,meanminusslope,meanminusintercept,lineintx,lineinty]
            pickle.dump(data,open(figname[:-4]+'.p','wb'))
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
