"""
===================
basic_sum.py
Erica  Lastufka 14.5.18
===================

Calculate the transmission via the sum of all (dark and flat corrected) pixels in a given image

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from PIL import Image
from matplotlib import cm
import pickle

import glob
import process_transmission as pt

def basic_sum(imlist,corrected=False,groupstretch=False,show=True):
    '''make a plot of the sum of the image as a function of angle'''
    sumvec=[]
    xvec=[]
    for im in imlist:
        imarr=pt.im2ndarray(im)
        sumvec.append(float(np.sum(imarr)))
        if corrected:
            xvec.append(float(im[9:-14]))
        elif groupstretch:
            xvec.append(float(im[im.find("_")+4:-17])) 
        else:
            xvec.append(float(im[im.rfind("_")+1:-4]))
    if show:
        fig,ax=plt.subplots()
        ax.scatter(xvec,sumvec/np.max(sumvec),s=25)
        #ax.set_ylim([0,1])
        ax.set_xlabel('Angle')
        fig.show()
    return xvec,sumvec

def binary_sum(imlist,corrected=True,show=True):
    '''classify first by light (1) and dark(0) then do the sum. Use the corrected images'''
    sumvec=[]
    xvec=[]
    for im in imlist:
        imarr=pt.im2ndarray(im)
        sumvec.append(float(np.sum(np.round(imarr,0))))
        if corrected:
            xvec.append(float(im[9:im.rfind("_")]))
        else:
            xvec.append(float(im[im.rfind("_")+1:-4]))
    if show:
        fig,ax=plt.subplots()
        ax.scatter(xvec,sumvec/np.max(sumvec),s=25)
        #ax.set_ylim([0,1])
        ax.set_xlabel('Angle')
        fig.show()
    return xvec,sumvec

def get_intersect(m1,m2,b1,b2):
    '''get intersection point (x,y) between two lines'''
    x=(b2-b1)/(m1-m2)
    y=m1*x+b1
    return x,y

def scat_allP(win,basicsum=False,binarysum=True,meanfit=True,figname=False,theta_ran=False,print_height=True,save=True):
    sumvecs=[]
    errvecs=[]
    fitsplus=[]
    fitsminus=[]
    allfitsplus,allfitsminus=[],[]
    for p in range(0,7):
        files=glob.glob('win'+win+'_p'+str(p)+'*'+'_groupstretch.tif')
        if basicsum:
            xv,sv=basic_sum(files,groupstretch=True,show=False)
        elif binarysum:
            xv,sv=binary_sum(files,corrected=True,show=False)

        xvmid=xv.index(0.0)
        #svsort=(sv[:xvmid][::-1]+sv[xvmid:])/np.max(sv)
        aa=zip(xv,sv/np.float64(sv[xvmid]))
        aa.sort()
        xv.sort()
        svsort=[b[1] for b in aa]
        sumvecs.append(svsort)

    meanvec=np.mean(sumvecs,axis=0)
    fplus=np.polyfit(xv[:xvmid-1],meanvec[:xvmid-1],1)
    fminus=np.polyfit(xv[xvmid+1:],meanvec[xvmid+1:],1)
    #print np.min(svsort), np.max(svsort),fplus,fminus
    fitsplus.append(fplus)
    fitsminus.append(fminus)

    labels=['p0','p1','p2','p3','p4','p5','p6']
    cols=['b','c','m','k','r','g','y']
    fig,ax=plt.subplots(2,sharex=True)
    #ax=fig.add_subplot(112)
    for sv,c,lab in zip(sumvecs,cols,labels):
        ax[0].plot(xv,sv,'-o',color=c,label=lab)
    ax[0].axvline(color='k',linestyle='dashed')
    if meanfit:
        meanplusslope=np.mean([fp[0] for fp in fitsplus])
        meanplusintercept=np.mean([fp[1] for fp in fitsplus])
        meanplusline=meanplusslope*np.array(xv[:xvmid+1]) + meanplusintercept
        #print meanplusline
        meanminusslope=np.mean([fm[0] for fm in fitsminus])
        meanminusintercept=np.mean([fm[1] for fm in fitsminus])
        meanminusline=meanminusslope*np.array(xv[xvmid:]) + meanminusintercept
        lineintx,lineinty=get_intersect(meanplusslope,meanminusslope,meanplusintercept,meanminusintercept)
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
    ax[0].set_title('Window '+ win + ' pixel sum as a function of angle')
    ax[0].set_ylim([0,1.1])
    ax[0].set_ylabel('Percent of maximum transmission')
    ax[0].legend(loc='upper right')

    #make the lower subplot showing the error bars
    #ax2=fig.add_subplot(121)
    #generate the points for the fit lines
    meanpluspoints=(meanplusslope*np.array(xv[:xvmid+1]) + meanplusintercept).tolist()
    meanminuspoints=(meanminusslope*np.array(xv[xvmid+1:]) + meanminusintercept).tolist()
    meanpluspoints.extend(meanminuspoints)
    #should plot the error bars here
    for sv,c,lab,err in zip(sumvecs,cols,labels,errvecs):
        yv=np.array(sv)-np.array(meanpluspoints)
        print np.min(yv),np.max(yv)
        ax[1].plot(xv,yv,'o',color=c,label=lab)
        ax[1].errorbar(xv,yv,yerr=err)
    ax[1].axhline(color='k',linestyle='dashed')

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
        print allfitsplus,allfitsminus
        for fp,fm,i in zip(allfitsplus,allfitsminus,range(0,len(allfitsplus))):
            heights,meanheight=pt.calc_height_from_fit(xv,fp[0],fp[1],fm[0],fm[1])
            print 'p',i,' mean height: ', meanheight
    if save:
        pickle.dump([xv,sumvecs,errvecs],open('win'+win+'basic_sum_data.p','wb'))
    return xv,sumvecs,errvecs

