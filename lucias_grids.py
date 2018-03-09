"""
===================
lucias_grids.py
Erica  Lastufka 2.3.18
===================

Some special adaptations to edge-and_hough.py to deal with the PSI grid images

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from PIL import Image
from matplotlib import cm
from skimage import feature, exposure
from skimage.transform import  (hough_line, hough_line_peaks,
                               probabilistic_hough_line)
import pickle
from scipy.misc import imrotate
from scipy.ndimage.interpolation import rotate
import time
import glob
import edge_and_hough as eh
    

def crop(filen,cropx,cropy,save=False):
    '''remove edges from images with edges of the window in them'''
    im=Image.open(filen)
    imcrop=im[cropx[0]:cropx[1],cropy[0]:cropy[1]]
    if save:
        #
    return imcrop

def group_edges_by_idx(rotarr,mask,nomp,mod=3,plot=False,tolerance=.6,from_line = False):
    '''Group edges based on mask index'''
    rising=list(np.where(mask[0,:] == 1.0)[0])
    falling=list(np.where(mask[0,:] == -1.0)[0])
    rmean,fmean=[],[]
    #for r in rising:
    #    subarr=rotarr[:,r-mod:r+mod]
    #    tarr = thin_to_binary(subarr,r-mod)
    #    #print r-mod
    #    if np.isfinite(np.mean(tarr)):
    #        total_edges=np.sum([1 for i in np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))])
    #        rmean.append(np.mean(tarr[np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))]))

    #for f in falling:
    #    subarr=rotarr[:,f-mod:f+mod]
    #    tarr = thin_to_binary(subarr,f-mod)
    #    if np.isfinite(np.mean(tarr)):
    #        total_edges=np.sum([1 for i in np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))])
    #        fmean.append(np.mean(tarr[np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))]))

    #rising.extend(falling)

    #if plot:
    #    fig,ax=plt.subplots()
    #    ax.scatter(rising,rmean+fmean)
    #    fig.show()

    periodr = [rmean[j+1]-rmean[j] for j in range(0,len(rmean)-1)]
    #periodr = [p for p in periodr if p > tolerance*nomp and p < (2.-tolerance)*nomp]
    #print rmean
    #print periodr
    periodf= [fmean[j+1]-fmean[j] for j in range(0,len(fmean)-1)]
    #periodf = [p for p in periodf if p > tolerance*nomp and p < (2.-tolerance)*nomp]
    #print np.mean(periodr)*1.955,np.mean(periodf)*1.955#, np.mean(periodr.extend(periodf))*1.955
            
    return periodr,periodf


def get_period_by_grouping(window_num,mag=5.0,plot=True,side=1.0,pix2um=1.955,stats=True,EM=True,tolerance=0.6,mosaic=True,offset=False, fitangle=True):
    '''Put everything together for a single window'''
    #first get the files
    if EM:
        edgef=EM_list(str(window_num),str(mag),ending='_edges.p')
    else:    
        edgef=glob.glob('win'+str(window_num)+'_*' + str(mag)+'X*_edges.p')
    if EM:
        imf=EM_list(str(window_num),str(mag),ending='.tif')
    else:
        imf=glob.glob('win'+str(window_num)+'_*' + str(mag)+'X*.tif')
    if mosaic:
        edgef=['win11_mosaic_5.0_edges.p']
        imf=['window11mosaic_5.0.tiff']
    
    rperiods,fperiods=[],[]
    badmasks=0
    if side == 1.0:
        nang=[aa['nominal angle'] for aa in windows if aa['number']==window_num][0]
        if offset:
            nang=nang+offset
        nperiod=[aa['pitch'] for aa in windows if aa['number']==window_num][0]
    else:
        nang=[aa['nominal angle'] for aa in windowsr if aa['number']==window_num][0]
        if offset:
            nang=nang+offset
        nperiod=[aa['pitch'] for aa in windowsr if aa['number']==window_num][0]
    for im,ef in zip(imf,edgef):
        #print im,ef
        edges=pickle.load(open(ef,'rb'))
        ima=im2ndarray(im)
        #ima=remove_edges(im,ima)
        #contrast stretching
        #p2 = np.percentile(ima, 2)
        #p98 = np.percentile(ima, 98)
        #ima = exposure.rescale_intensity(ima, in_range=(p2, p98))
        #immean=np.mean(ima)
        #now make mask
        if fitangle:
            try:
                testang=pickle.load(open(ef[:-8]+'_hough_all_mean.p','rb'))
            except IOError:
                testang=mean_hough_angle(ef[:-8]+'_hough_all.p',side=side)
            if not np.isnan(testang):
                nang=testang #otherwise keep as nominal angle
        print nang
        rotim=rotate(np.transpose(ima),side*nang, reshape=True)
        rotarr=rotate(np.transpose(edges),side*nang,reshape=True)
        #print im,ef
        mask=rising_or_falling(rotim[300:450,:],rotarr[300:450,:],np.mean(rotim[300:450,100:-100]), shape=np.shape(rotarr),imfile=im,plot=plot)
        #group
        if np.shape(mask)[0] > 0:
            periodr,periodf=group_edges_by_idx(rotarr,mask,nperiod/pix2um,tolerance=tolerance,mod=3,plot=plot)
            return periodr,periodf
            rperiods.extend(periodr)
            fperiods.extend(periodf)
        else:
            badmasks+=1
            print im

    print 'bad masks ', badmasks
            
    periods=rperiods
    periods.extend(fperiods)

    rperiods=np.array(rperiods)*pix2um
    fperiods=np.array(fperiods)*pix2um
    periods=np.array(periods)*pix2um
        
    if plot:
        #make the histogram
        fig,ax=plt.subplots(1,3,sharex=True,sharey=True)
        bins=np.arange(np.min(periods),np.max(periods),np.std(periods)/2.)
        ax[0].hist(rperiods,bins,facecolor='g')
        ax[1].hist(fperiods,bins,facecolor='r')
        ax[2].hist(periods,bins)
        ax[0].set_xlim([nperiod-5,nperiod+5])
        ax[0].set_title('Rising')
        ax[1].set_title('Falling')
        ax[2].set_title('Total')
        ax[0].set_xlabel('Period $\mu$m')
        ax[0].set_ylabel('Counts')
        ax[0].set_yscale('log')
        ax[0].set_ylim([1,10000])
        figfilename='win'+str(window_num)+'_group_periods'+str(mag)+'.png'
        fig.savefig(figfilename)
        fig.show()
        
    if stats: #print out and pickle stats
        avg=np.mean(periods)
        med=np.median(periods)
        stdv=np.std(periods)
        statfile='win'+str(window_num)+'_width_stats'+str(mag)+'.p'
        datafile='win'+str(window_num)+'_width_data'+str(mag)+'.p'
        print "-------------------STATISTICS FOR WINDOW "+str(window_num)+"---------------------"
        print '              Mean: ' + str(avg)
        print '            Median: ' + str(med)
        print 'Standard Deviation: ' + str(stdv)
        print 'Results saved in ' + statfile
        data=periods
        stdict={'mean':avg,'median':med,'stddev':stdv}
        pickle.dump(stdict,open(statfile,'wb'))
        pickle.dump(data,open(datafile,'wb'))

    return rperiods, fperiods,periods

