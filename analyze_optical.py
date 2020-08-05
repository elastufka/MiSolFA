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
import scipy
from scipy.misc import imrotate
from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit
import time
import glob
import os
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

def analyze_optical(grating, folder=False, ll=[50,100,150,200,250],coverwrite=False,eoverwrite=False,hoverwrite=False,tol=9,plot=True,stats=True,htol=5,ftags='a',pix2um=1.955,mask45=False,csigma=3,sigma=2.,r_or_f_tol=3,spread=2.0,n_hough=235,filter_nominal=False,from_centers=True, example=True, angle_only=False):
    '''main method - run everything. example means... select one image to show all plots for'''
    win=grating.win
    nang=grating.nominal['orientation']
    filen=grating.data.Odata.rawim
    filec=grating.data.Odata.corrim
    filep=grating.data.Odata.edges
    fileh=grating.data.Odata.cat_hough

    if filec==[] or len(filec) != len(filen) or coverwrite: #if not contrast stretched, do so
        filec=[ag.contrast_stretch(fn) for fn in filen]
        if example:
            foo=ag.contrast_stretch('win'+str(win)+'_05_05_5.0X.tif',show=True,hist=True)

    if grating.side == -1.0 and grating.model=='QM': #it's the May dataset
        filec=fix_May_im(filec)
        if len(filec) == 0: #assume all is well?
            filec=glob.glob('win'+str(win)+'*May_corrected.tif')
        else:
            print filec[10]
        #eoverwrite=True
    if filep==[] or len(filep) != len(filec) or eoverwrite: #need to generate the pickles by doing the edges
        filep=[get_Canny_edges(fn,sigma=csigma) for fn in filec]
        #.append(eh.Canny_edge(f,mag=5.0,outfilen=f[:-4]+"_ani_edges.p")) #new list of picklefiles
        if example:
            foo=ag.Canny_edge('win'+str(win)+'_05_05_5.0X_corrected.tif',plot=True)

        #update grating object
    if fileh==[]or len(fileh) != len(filen) or hoverwrite: #need to generate the pickles by doing the edges
        fileh=[]
        #for f in filep:
        #    fileh.append(ag.prob_hough(f,nang,line_length=ll[0],spread=2.,overwrite=hoverwrite)) #new list of picklefiles
        for l in ll:
            tag='ll'+str(l)
            for f in filep:
                fileh.append(ag.prob_hough(f,nang,line_length=l,spread=spread,n=n_hough,tag=tag,overwrite=hoverwrite)) #new list of picklefiles
    tags=['ll'+str(t) for t in ll]
    if hoverwrite:
        imtags='X_corrected_'
        if grating.model=='QM' and grating.side == -1.0:
            imtags='X_May_corrected_'
        ag.cat_hough(win,weights=[1,1,1,1,1],tags=tags,mag=5.0,EM=False,imtags=imtags)

    if angle_only:
        hough_all=glob.glob('win'+str(win)+'*_5.0X_*corrected_hough_all.p')
        if example:
            hf=glob.glob('win'+str(win)+'_05_05_5.0X_*corrected_hough_all.p')[0]
            ef=glob.glob('win'+str(win)+'_05_05_5.0X_*corrected_edges.p')[0]
            ag.hough_peek(hf,ef,)
        #NEED TO RE-WRITE THIS IN THE RIGHT ORDER!!!
        figname='win'+str(win) + 'ang_hist.png'
        ag.hough_hist(hough_all,win,nang,figname=figname,tol=htol,spread=spread,mask45=mask45) #what is this tol now?
        return #maybe make this return a results object eventually... but for now ok

    #sort files so the right ones are matching!!
    filen.sort()
    filec.sort()
    filep.sort()

    #a,b,c=get_rotangle(win,filep[0],filen[0],nang,testplot=True)
    rising,falling,periodr,periodf,periods,widths, from_c,filtered_c,from_c2,filtered_c2=[],[],[],[],[],[],[],[],[],[]
    for im,ef in zip(filec,filep):
        rang,rotim,rotarr=get_rotangle(win,ef,im,nang) #should probably check if nang and the nominal angle are the same sign...
        if ('win'+str(win)+'_05_05_5.0X') in im: #_corrected.tif':
            tplot=True
            #print sigma
        else:
            tplot=False
        xpeak,xpeak_int=get_xpeaks(rang,rotim,rotarr,im,ef,sigma=sigma,tol=tol,plot=tplot)  #why is sigma staying equal to 2???

        #get xpeaks of image rotated by nominal angle, for calculation of rms
        #_,rotim_nom,rotarr_nom=get_rotangle(win,ef,im,nang,nominal=True)
        #xpeak_rms,xpeak_int_rms=get_xpeaks(nang,rotim_nom,rotarr_nom,im,ef,sigma=sigma,tol=tol,plot=tplot)
        #need to get slat centers from this

        if len(xpeak) == 0: #really small frame... skip it!
            continue
        #         #NOW do slit_or_slat
        #need to mask values of rotim where rotim == 0
#         rotim_masked=np.ma.masked_equal(rotim,0.0)
#         imsum=np.ma.sum(rotim_masked,axis=0)
#         imgrad=np.gradient(np.array(imsum,dtype=float))

#         #if the number of xpeaks is way less than expected (say 5-10X less, meaning it's probably an edge frame), activate filter_nominal
#         #image size is 480/640
#         #number of edges in an image should be about:
#         #sqrt(2)* (2*640/number of pixels in a period)
#         #= sqrt(2)*(2*640/grating.nominal['pitch']/pix2um)
#
        nedges=np.sqrt(2.)*(1280./(grating.nominal['pitch']/pix2um))
#         if len(xpeak) < nedges/5.:
#             filter_nominal = True
#         if filter_nominal:
#             filter_nominal=grating.nominal['pitch']/pix2um
#         rising,falling,rperiod1,fperiod1,tperiod1,width1, centers,centers2=ag.rising_or_falling_final(xpeak,xpeak_int,imgrad,n=r_or_f_tol,filter_multiples=True,filter_nominal=filter_nominal,centers=from_centers)
#         if im=='win'+str(win)+'_05_05_5.0X_corrected.tif' and example:
#             plot_peaks_and_points(xpeak,rising, falling, centers,centers2)
#         if from_centers==True:
#             period,filtered_periods,period2,filtered_periods2=ag.slit_widths_from_centers(centers, grating.nominal['pitch'],centers2=centers2, pix2um=pix2um,filter_nominal=filter_nominal)
#             fperiodum=pix2um*np.array(filtered_periods)
#             periodum=pix2um*np.array(period)
#             fperiodum2=pix2um*np.array(filtered_periods2)
#             periodum2=pix2um*np.array(period2)
#         #if they are masked arrays, flatten them...
#         if type(tperiod1) == np.ma.core.MaskedArray:
#             tperiod1=tperiod1.compressed()
#         if type(width1) == np.ma.core.MaskedArray:
#             width1=width1.compressed()

#         #need to multiply by pix2um
#         period1r=pix2um*np.array(rperiod1)
#         period1f=pix2um*np.array(fperiod1)
#         period1s=pix2um*np.array(tperiod1)
#         width1s=pix2um*np.array(width1)


        #do rms difference, add it to stats ...
        #rotim_masked_nom=np.ma.masked_equal(rotim_nom,0.0)
        #imsum_nom=np.ma.sum(rotim_masked_nom,axis=0)
        #imgrad_nom=np.gradient(np.array(imsum_nom,dtype=float))
        #_,_,_,_,_,_,centers_rms,_=ag.rising_or_falling_final(xpeak_rms,xpeak_int_rms,imgrad_nom,n=r_or_f_tol,filter_multiples=True,filter_nominal=filter_nominal,centers=from_centers)
#        prms,phase=get_pitch_rms(grating,centers_rms,plot=tplot)

        rise, fall, period1r,period1f,period1s,width1s,fperiodum,periodum,fperiodum2,periodum2,centers,centers2=peaks2periods(grating,im,rotim,rotarr,rang,pix2um,xpeak,xpeak_int,r_or_f_tol=r_or_f_tol,from_centers=from_centers)

        #add to master lists
        rising.append(rise)
        falling.append(fall)
        periodr.extend(period1r)
        periodf.extend(period1f)
        periods.extend(period1s)
        widths.extend(width1s)
        if from_centers:
            from_c.extend(periodum)
            filtered_c.extend(fperiodum)
            from_c2.extend(periodum2)
            filtered_c2.extend(fperiodum2)

    rdict={'periods':periods,'periodr':periodr,'periodf':periodf,'widths':widths,'periods_from_slat_centers':from_c,'filtered_periods_slat':filtered_c,'periods_from_slit_centers':from_c2,'filtered_periods_slit':filtered_c2}
    #what's the difference between single_data_dict and rdict??

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
        print '             Total Period (slat centers): ' + str(len(from_c))#str(len(filtered_c))
        print '              Mean Period (slat centers): ' + str(np.mean(from_c))#str(np.mean(filtered_c))
        print '            Median Period (slat centers): ' + str(np.median(from_c))#str(np.median(filtered_c))
        print 'Standard Deviation Period (slat centers): ' + str(np.std(from_c)) #str(np.std(filtered_c))
        print '             Total Period (slit centers): ' + str(len(from_c2))
        print '              Mean Period (slit centers): ' + str(np.mean(from_c2))
        print '            Median Period (slit centers): ' + str(np.median(from_c2))
        print 'Standard Deviation Period (slit centers): ' + str(np.std(from_c2))

        print '          Mean Duty Cycle: ' + str(np.ma.mean(widths)/avg)
        print 'Results saved in ' + statfile
        data={'rising':rising,'falling':falling,'period':periods,'prising':periodr,'pfalling':periodf,'widths':widths,'periods_from_slat_centers':from_c,'filtered_periods_slat':filtered_c,'periods_from_slit_centers':from_c2,'filtered_periods_slit':filtered_c2}
        stdict={'mean':avg,'median':med,'stddev':stdv,'wmean':np.ma.mean(widths),'wmed':np.ma.median(widths),'wdev':np.ma.std(widths),'amean':np.mean(from_c),'amed':np.median(from_c),'adev':np.std(from_c),'imean':np.mean(from_c2),'imed':np.median(from_c2),'idev':np.std(from_c2)}
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
        npitch=grating.nominal['pitch']
        bins=np.arange(npitch-npitch/10., npitch+npitch/10., npitch/40.)#np.arange(np.min(from_c),np.max(from_c),np.std(from_c)/2.)
        ax[0].hist(periodr,bins,facecolor='g')
        ax[1].hist(periodf,bins,facecolor='r')
        ax[2].hist(from_c,bins,facecolor='b',alpha=0.6)
        ax[2].hist(from_c2,bins,facecolor='m',alpha=0.6)

        #ax[2].hist(filtered_c,bins,facecolor='b',alpha=0.6)
        #ax[2].hist(filtered_c2,bins,facecolor='m',alpha=0.6)
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

    kind={'type':'optical','model':grating.model,'folder':grating.data.Odata.path,'side':grating.side}
    #results init: (self,win,kind,method=False,errors=False,data=False,nominal=False,filenames=False,stats=False,params=False)
    res=results.Results(win,kind,[],method='optical',errors='usual',data=data,nominal=grating.nominal,stats=stdict,params={'tol':tol, 'htol':htol,'sigma':sigma,'r_or_f_tol':r_or_f_tol, 'll':ll,'csigma':csigma}) #add angle data later


    return res

def peaks2periods(grating,im,rotim,rotarr, rang,pix2um,xpeak,xpeak_int,r_or_f_tol=2,from_centers=True,filter_nominal=False,example=True):
    '''consolidate everything that was in analyze_optical '''
    win=grating.win
    rotim_masked=np.ma.masked_equal(rotim,0.0)
    imsum=np.ma.sum(rotim_masked,axis=0)
    imgrad=np.gradient(np.array(imsum,dtype=float))

    #if the number of xpeaks is way less than expected (say 5-10X less, meaning it's probably an edge frame), activate filter_nominal
    #image size is 480/640
    #number of edges in an image should be about:
    #sqrt(2)* (2*640/number of pixels in a period)
    #= sqrt(2)*(2*640/grating.nominal['pitch']/pix2um)
    nedges=np.sqrt(2.)*(1280./(grating.nominal['pitch']/pix2um))
    #if len(xpeak) < nedges/5.:
    #    filter_nominal = True
    #if filter_nominal:
    #    filter_nominal=grating.nominal['pitch']/pix2um

    rising,falling,rperiod1,fperiod1,tperiod1,width1, centers,centers2=ag.rising_or_falling_final(xpeak,xpeak_int,imgrad,n=r_or_f_tol,filter_multiples=True,filter_nominal=filter_nominal,centers=from_centers)
    if im=='win'+str(win)+'_05_05_5.0X_corrected.tif' and example:
        plot_peaks_and_points(xpeak,rising, falling, centers,centers2)
    if from_centers==True:
        period,filtered_periods,period2,filtered_periods2=ag.slit_widths_from_centers(centers, grating.nominal['pitch'],centers2=centers2, pix2um=pix2um,filter_nominal=filter_nominal)
        fperiodum=pix2um*np.array(filtered_periods)
        periodum=pix2um*np.array(period)
        fperiodum2=pix2um*np.array(filtered_periods2)
        periodum2=pix2um*np.array(period2)
    #if they are masked arrays, flatten them...
    if type(tperiod1) == np.ma.core.MaskedArray:
        tperiod1=tperiod1.compressed()
    if type(width1) == np.ma.core.MaskedArray:
        width1=width1.compressed()
    period1r=pix2um*np.array(rperiod1)
    period1f=pix2um*np.array(fperiod1)
    period1s=pix2um*np.array(tperiod1)
    width1s=pix2um*np.array(width1)

    #pickle rising and falling so it's easier to inspect the data later.
    picklename=im[:-4]+'inspect.p'

#                 from_c.extend(periodum)
#             filtered_c.extend(fperiodum)
#             from_c2.extend(periodum2)
#             filtered_c2.extend(fperiodum2)


    single_data_dict={'rising':rising,'falling':falling,'rotim':rotim_masked,'rotarr':rotarr,'periodr':period1r,'periodf':period1f,'widths':width1s,'periods':period1s,'rang':rang,'imshape':np.shape(ag.im2ndarray(im)), 'periods_from_slat_centers':periodum,'filtered_periods_slat':fperiodum,'periods_from_slit_centers':periodum2,'filtered_periods_slit':fperiodum2}
    pickle.dump(single_data_dict,open(picklename,'wb'))


    return rising, falling,period1r,period1f,period1s,width1s,fperiodum,periodum,fperiodum2,periodum2,centers,centers2

def get_pitch_rms(grating, centers_rms, plot=False): #need centers first
    if 'slat_centers' not in grating.nominal.keys():
        grating.calc_nominal_slat_centers()
    scnom=grating.nominal['slat_centers'] #float vector of slat centers
    #vector function to find avg. diff between position in nominal slat centers and calculated. (how to treat missing slat centers?)
    #from analyze_general: def find_nearest(array, value,ret_idx=True):
    xidx_list=[]
    phases=[]
    for x in centers_rms:
        xidx,nval=ag.find_nearest(scnom,x,ret_idx=True)
        #if index not there, it's a missing slat. get rid of these indexes...
        xidx_list.append(xidx)
        diff=nval-x
        phases.append(diff)

    #use this to get the phase
    #then shift by phase
    scnoml=list(scnom)
    phase=np.mean(phases)
    scnom_shift=[s-phase for s in scnoml]
    scnom_trimmed=[s for s in scnom_shift if scnom_shift.index(s) in xidx_list]
    #should also return missed detections
    missed=len(scnom_shift)-len(centers_rms) #not quite right since scnom is padded, have to trim those

    #finally calculate rms. after trimming vector to fit...
    rms_list=[(x-s)**2. for x,s in zip(centers_rms,scnom_trimmed)]
    #RMSD is the square root of the average of squared errors.
    prms=np.sqrt(np.mean(rms_list))

    if plot:
        print 'phase mean and std:', np.mean(phases),np.std(phases)
        print 'pitch rms (pix):',prms
        print len(scnom_trimmed,centers_rms)
        fig,ax=plt.subplots() #make plot showing measured & theoretical slit centers, print rms and phase
        for tp in scnom_shift:
            ax.axvline(tp,color='g')
        for mp in centers_rms:
            ax.axvline(mp,color='m')
        ax.axvline(tp,color='g',label='nominal slat center')
        ax.axvline(mp,color='m',label='measured slat center')

        at0='phase: '+str(np.round(phase,3))
        at1='pitch rms: '+str(np.round(prms,3))
        ax.text(centers_rms[1],1.75,at0)
        ax.text(centers_rms[1],1.5,at1)
        ax.legend()
        ax.set_xlim([centers_rms[0]-1,centers_rms[-1]+1])
        ax.set_ylim([0,2])
        ax.legend()
        fig.show()

    return prms,phase

def fix_May_im(imlist,midrange=[.4,.6],plot=False):
    '''take images in bad data set and turn them binary based off of if they are grey (slits, 'midrange') or greatly varying black/white (minmax, slats) '''
    xvec=[]
    for im in imlist:
        imarr=ag.im2ndarray(im)
        pname=im[:-13]+'May'+im[-14:]
        if os.path.isfile(pname): #assume this has already been done
            continue
        #print pname
        xvec.append(pname)
        #convert each image to binary representation before summing. Use the mean as the cutoff
        bmax=np.max(imarr)
        bmin=np.min(imarr)
        bmean=np.mean(imarr)
        imgt=np.ma.masked_greater(imarr,(1.+midrange[1])*bmean) #what happens if the value equals the mean?
        #imgtfilled=np.ma.filled(imgt, 1.0)
        imlt=np.ma.masked_less(imarr,(1.-midrange[1])*bmean)
        samask = imgt+imlt
        #imltfilled=np.ma.filled(imlt, 1.0)
        #print midrange[0]*bmean,(1.+(1.-midrange[1]))*bmean
        simask=np.ma.masked_inside(imarr,midrange[0]*bmean,(1.+(1.-midrange[1]))*bmean)
        #sllt=np.ma.masked_less(imarr,midrange[0]*bmean) #have to do this on the original image, not where it's already been masked...
        #simask=slgt+sllt
        #as for the rest?
        #imarr=imltfilled
        imslits=np.ma.filled(simask,0.0)
        imslats=np.ma.filled(samask,255.0)
        imfixed=np.ma.filled(imslits,255.0)

        #pickle
        #pickle.dump(imfixed,open(pname,'wb'))
        savec=scipy.misc.toimage(imfixed,mode='F')
        savec.save(pname)
    if plot:
        fig,ax=plt.subplots(2,2)
        ax[0][0].imshow(imarr,cmap='gray')
        ax[0][1].imshow(imfixed,cmap='gray')
        ax[1][0].imshow(imslats,cmap='gray')
        ax[1][1].imshow(imslits,cmap='gray')
        fig.show()
    return xvec

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

def plot_peaks_and_points(xpeak,rising, falling, centers,centers2):
    '''plot peaks, rising, falling, centers'''
    fig,ax=plt.subplots()
    for xp in xpeak:
        ax.axvline(xp)
    orising=rising[:-len(falling)]
    ax.scatter(orising, np.ones(len(orising)), marker='<',c='b', label='rising')
    ax.scatter(falling, np.ones(len(falling)), marker='>',c='g', label='falling')
    ax.scatter(centers, np.ones(len(centers)), marker='v',c='r', label='slat center')
    ax.scatter(centers2, np.ones(len(centers2)), marker='^',c='m', label='slit center')
    ax.legend()
    ax.set_xlim([xpeak[0]-1,xpeak[-1]+1])
    ax.set_ylim([0,2])
    ax.legend()
    fig.show()

def get_rotangle(window_num, ef,im,nang,offset=False,testplot=False,recalc=False,nominal=False):
    if '/' in ef:
        ef=ef[ef.rfind('/')+1:]
    if offset:
            nang=nang+offset
    edges=pickle.load(open(ef,'rb'))
    ima=ag.im2ndarray(im)
    if not recalc:
#         try:
#             testang=pickle.load(open(ef[:26]+'_hough_all_mean.p','rb'))
#         except IOError:
#             testang=mean_hough_angle(ef[:26]+'_hough_all.p')
#     else:
        if 'May' in ef:
            hall=ef[:30]+'_hough_all.p'
        elif 'corrected' in ef:
            hall=ef[:17]+'corrected_hough_all.p'
        else:
            hall=ef[:26]+'corrected_hough_all.p'
        #try:
        testang=mean_hough_angle(hall,plot=testplot)
        #except IOError:
            #testang=mean_hough_angle(ef[:17]+'corrected_hough_all.p',plot=testplot)

    if not np.isnan(testang):
        rang=testang #otherwise keep as nominal angle
    else:
        rang=nang
        #print nang
    if nominal:
        rang=nang
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
