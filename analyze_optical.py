"""
===================
analyze_optical.py
Erica  Lastufka 20.9.17
===================

Implementations of general Analyze methods specific to optical data. Based on get_all_angles_and_widths.py

Input: Grating object
Output: Result object

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from PIL import Image
from matplotlib import cm
from skimage import feature, exposure
from skimage.transform import  (hough_line, hough_line_peaks,
                               probabilistic_hough_line)
import itk
import pickle
from scipy.misc import imrotate
from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit
import time
import glob
import analyze_general as ag
import results

# def prep_all_images(grating,folder=False):
#     '''Edge trimming, contrast optimizing, etc. of images '''
#     if folder:
#         os.chdir(folder)
#     else:
#         os.chdir(grating.Folder.optical)
#     win=grating.win
#     if grating.Model == "EM":
#         filen=eh.EM_list(str(win),'5.0',ending='.tif')
#         filen=glob.glob('win'+str(win)+'*_5.0X.tif')
#         filep=ag.EM_list(str(win),'5.0',ending='_edges.p')
#         fileh=ag.EM_list(str(win),'5.0',ending='_hough.p')
#     else:
#         filen=glob.glob('win'+str(win)+'*_5.0X.tif')
#         filep=glob.glob('win'+str(win)+'*_5.0X_edges.p')
#         fileh=glob.glob('win'+str(win)+'*_5.0X_hough_all.p')
#     return filen,filep,fileh

def get_Canny_edges(fn,mag=5.0,sigma=3,tag='_edges'):
    '''Wrapper for running Canny_edge for the optical case'''
    edges=ag.Canny_edge(fn, mag=mag, outfilen=fn[:-4]+tag+".p",sigma=sigma)
    return fn[:-4]+tag+".p"

def analyze_optical(grating, folder=False, ll=[50,100,150,200,250],coverwrite=False,hoverwrite=False,tol=9,plot=True,stats=True,htol=5,ftags='a',pix2um=1.955,mask45=False,sigma=2.,r_or_f_tol=3,filter_nominal=False):
    '''main method - run everything'''
    win=grating.win
    nang=grating.nominal['orientation']
    filen=grating.data.Odata.rawim
    filep=grating.data.Odata.edges
    fileh=grating.data.Odata.cat_hough

    if filep==[] or len(filep) != len(filen) or coverwrite: #need to generate the pickles by doing the edges
        filep=[get_Canny_edges(fn) for fn in filen]#.append(eh.Canny_edge(f,mag=5.0,outfilen=f[:-4]+"_ani_edges.p")) #new list of picklefiles
        #update grating object
    if fileh==[]or len(fileh) != len(filen) or hoverwrite: #need to generate the pickles by doing the edges
        fileh=[]
        #for f in filep:
        #    fileh.append(ag.prob_hough(f,nang,line_length=ll[0],spread=2.,overwrite=hoverwrite)) #new list of picklefiles
        for l in ll:
            tag='ll'+str(l)
            for f in filep:
                fileh.append(ag.prob_hough(f,nang,line_length=l,spread=2.,tag=tag,overwrite=hoverwrite)) #new list of picklefiles
    tags=['ll'+str(t) for t in ll]
    if hoverwrite:
        ag.cat_hough(win,weights=[1,1,1,1,1],tags=tags,mag=5.0,EM=False)
    hough_all=glob.glob('win'+str(win)+'*_5.0X_hough_all.p')
    #ag.hough_peek(hough_all[10],filep[10])
         #NEED TO RE-WRITE THIS IN THE RIGHT ORDER!!!
    figname='win'+str(win) + 'ang_hist.png'
    ag.hough_hist(hough_all,win,nang,figname=figname,tol=htol,spread=2.,mask45=mask45) #what is this tol now?

    #a,b,c=get_rotangle(win,filep[0],filen[0],nang,testplot=True)
    periodr,periodf,periods,widths=[],[],[],[]
    for im,ef in zip(filen,filep):
        rang,rotim,rotarr=get_rotangle(win,ef,im,nang) #should probably check if nang and the nominal angle are the same sign...
        if im=='win44_07_08_5.0X.tif':
            tplot=True
            #print sigma
        else:
            tplot=False
        xpeak,xpeak_int=get_xpeaks(rang,rotim,rotarr,im,ef,sigma=sigma,tol=tol,plot=tplot)  #why is sigma staying equal to 2???
        if len(xpeak) == 0: #really small frame... skip it!
            continue
        #         #NOW do slit_or_slat
        #need to mask values of rotim where rotim == 0
        rotim_masked=np.ma.masked_equal(rotim,0.0)
        imsum=np.ma.sum(rotim_masked,axis=0)
        imgrad=np.gradient(np.array(imsum,dtype=float))

        #if the number of xpeaks is way less than expected (say 5-10X less, meaning it's probably an edge frame), activate filter_nominal
        #image size is 480/640
        #number of edges in an image should be about:
        #sqrt(2)* (2*640/number of pixels in a period)
        #= sqrt(2)*(2*640/grating.nominal['pitch']/pix2um)
        nedges=np.sqrt(2.)*(1280./(grating.nominal['pitch']/pix2um))
        if len(xpeak) < nedges/5.:
            filter_nominal = True
        if filter_nominal:
            filter_nominal=grating.nominal['pitch']/pix2um
        rising,falling,rperiod1,fperiod1,tperiod1,width1=ag.rising_or_falling_final(xpeak,xpeak_int,imgrad,n=r_or_f_tol,filter_multiples=True,filter_nominal=filter_nominal)

        #if they are masked arrays, flatten them...
        if type(tperiod1) == np.ma.core.MaskedArray:
            tperiod1=tperiod1.compressed()
        if type(width1) == np.ma.core.MaskedArray:
            width1=width1.compressed()

        #need to multiply by pix2um
        period1r=pix2um*np.array(rperiod1)
        period1f=pix2um*np.array(fperiod1)
        period1s=pix2um*np.array(tperiod1)
        width1s=pix2um*np.array(width1)

        #add to master lists
        periodr.extend(period1r)
        periodf.extend(period1f)
        periods.extend(period1s)
        widths.extend(width1s)

        #pickle rising and falling so it's easier to inspect the data later.
        picklename=im[:-4]+'inspect.p'

        single_data_dict={'rising':rising,'falling':falling,'earr':rotarr,'rotim':rotim_masked,'periodr':period1r,'periodf':period1f,'widths':width1s,'periods':period1s,'rang':rang,'imshape':np.shape(ag.im2ndarray(im))}
        pickle.dump(single_data_dict,open(picklename,'wb'))


    rdict={'periods':periods,'periodr':periodr,'periodf':periodf,'widths':widths}
    if stats: #print out and pickle stats
        avg=np.ma.mean(periods)
        med=np.ma.median(periods)
        stdv=np.ma.std(periods)
        statfile='win'+str(win)+'_width_stats5.0X'+ftags+'.p'
        datafile='win'+str(win)+'_width_data5.0X'+ftags+'.p'
        print "-------------------STATISTICS FOR WINDOW "+str(win)+"---------------------"
        print '             Total Period: ' + str(len(periods))
        print '              Mean Period: ' + str(avg)
        print '            Median Period: ' + str(med)
        print 'Standard Deviation Period: ' + str(stdv)
        print '             Total Widths: ' + str(len(widths))
        print '              Mean Widths: ' + str(np.ma.mean(widths))
        print '            Median Widths: ' + str(np.ma.median(widths))
        print 'Standard Deviation Widths: ' + str(np.ma.std(widths))
        print '          Mean Duty Cycle: ' + str(np.ma.mean(widths)/avg)
        print 'Results saved in ' + statfile
        data={'period':periods,'rising':periodr,'falling':periodf,'widths':widths}
        stdict={'mean':avg,'median':med,'stddev':stdv}
        pickle.dump(stdict,open(statfile,'wb'))
        pickle.dump(data,open(datafile,'wb'))

    if plot:
        #make the histogram
        if type(periods) == np.ma.core.MaskedArray:
            print 'masked arrays'
            periods=periods.compressed()
            periodr=periodr.compressed()
            periodf=periodf.compressed()
        fig,ax=plt.subplots(1,3,sharex=True,sharey=True)
        bins=np.arange(np.min(periods),np.max(periods),np.std(periods)/2.)
        ax[0].hist(periodr,bins,facecolor='g')
        ax[1].hist(periodf,bins,facecolor='r')
        ax[2].hist(periods,bins)
        ax[0].set_xlim([grating.nominal['pitch']-5,grating.nominal['pitch']+5])
        ax[0].set_title('Rising')
        ax[1].set_title('Falling')
        ax[2].set_title('Total')
        ax[0].set_xlabel('Period $\mu$m')
        ax[0].set_ylabel('Counts')
        ax[0].set_yscale('log')
        ax[0].set_ylim([1,10000])
        figfilename='win'+str(win)+'_group_periods5.0X'+ftags+'.png'
        fig.savefig(figfilename)
        fig.show()

    res=results.Results(win,'optical','usual',False,rdict,params={'tol':tol, 'htol':htol,'sigma':sigma,'r_or_f_tol':r_or_f_tol}) #add angle data later


    return res

def get_xpeaks(rang,rotim,rotarr,im,ef,sigma=2.,tol=9,plot=False):
    #now clean edges of the rotated arrays?
    carray=rotate(np.transpose(np.ones(np.shape(ag.im2ndarray(im)))),(rang), reshape=True)
    carr=np.sum(carray,axis=0)
    cleaned_edges=ag.clean_centers(rotarr,cfac=carr, sigma=sigma,plot=plot) #make the cutoff line follow a curve determined by the local maxima? Or smooth the array ... idk
    #then group edges together using group_edges_by_idx (I think)
    xpeak=ag.sum_peaks(cleaned_edges, tol=tol,irange=False)
    xpeak_int=[int(x) for x in xpeak]
    xpeak_int.sort()
    if plot:
        fig,ax=plt.subplots()
        ax.imshow(rotarr,origin='lower',alpha=1,cmap=cm.Reds)
        ax.imshow(rotim,origin='lower',alpha=.5,cmap=cm.Greys)
        fig.show()
        pickle.dump([xpeak,xpeak_int],open('testpeaks.p','wb'))
        #plot_centers_and_edges(fname=im,datafile=,tol=2.0,cfac=False)
    return xpeak,xpeak_int


def get_rotangle(window_num, ef,im,nang,offset=False,testplot=False,recalc=False):
    if offset:
            nang=nang+offset
    edges=pickle.load(open(ef,'rb'))
    ima=ag.im2ndarray(im)
    if not recalc:
        try:
            testang=pickle.load(open(ef[:-8]+'_hough_all_mean.p','rb'))
        except IOError:
            testang=mean_hough_angle(ef[:-8]+'_hough_all.p')
    else:
        testang=mean_hough_angle(ef[:-8]+'_hough_all.p',plot=testplot)
    if not np.isnan(testang):
        rang=testang #otherwise keep as nominal angle
    else:
        rang=nang
        #print nang
    rotim=rotate(np.transpose(ima),(rang), reshape=True)
    rotarr=rotate(np.transpose(edges),(rang),reshape=True)
    if testplot:
        fig,ax=plt.subplots()
        ax.imshow(rotim,origin='lower',alpha=.6,cmap=cm.gray)
        ax.imshow(rotarr,origin='lower',cmap=cm.Reds)
        fig.show()
    return rang,rotim,rotarr

def mean_hough_angle(f_all,save=True, gaussfit=True,plot=False):
    lines=pickle.load(open(f_all,'rb'))
    angles=[]
    for l in lines:
        angles.append(ag.get_angle(l))
    meana=np.mean(angles) #maybe should fit a Gaussian?
    if np.isnan(meana):
        return meana
    if gaussfit:
        thetaax=np.linspace(np.min(angles),np.max(angles),15)
        yhist,xhist=np.histogram(angles,thetaax)
        #xh = np.where(yhist > 0)[0]
        #yh = yhist[xh]
        #print np.max(yhist), np.mean(yhist), np.std(yhist)
        def gaussian(x, a, mean, sigma):
            return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))
        try:
            y_at_xmean=yhist[np.where(xhist>np.mean(angles))[0][0]]
        except IndexError:
            return meana
        try:
            popt, pcov = curve_fit(gaussian, xhist[:-1], yhist, [y_at_xmean, np.mean(angles), np.std(angles)])
            if plot:
                fig,ax=plt.subplots()
                ax.hist(angles,thetaax)
                ax.plot(thetaax, gaussian(thetaax, *popt))
                fig.show()
            #print np.max(gaussian(thetaax, *popt)),np.min(gaussian(thetaax, *popt))
            mean=np.max(gaussian(thetaax, *popt))
        except RuntimeError: #not enough line segments! ex. frames on the edges
            mean=np.mean(angles)


    if save:
        pickle.dump(meana,open(f_all[:-2]+'_mean.p','wb'))
    return meana

# def get_angle(line):
#     '''Get angle of line via tangent. Note that because the top left corner is 0,0 in the background image we multiply x's by -1'''
#     deltax=(line[1][0]-line[0][0])
#     deltay=line[1][1]-line[0][1]
#     theta=np.arctan(float(deltay)/float(deltax))  #np.arctan2(float(deltay)/float(deltax))
#     thetadeg=np.rad2deg(theta) #theta*180./np.pi
#     return thetadeg
