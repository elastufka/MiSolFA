"""
===================
correct_images.py
Erica  Lastufka 12.4.18
===================

Do the image pre-processing:

(image - dark)/ (flat - dark)

"""
import numpy as np
import os
import glob
from PIL import Image
from skimage import feature, exposure
import pickle
from matplotlib import cm
import scipy.misc

def make_flat(flatdir, nflats=100,saveim=False):
    '''Average all the flats. assume flats are after darks in sequential numbering'''
    os.chdir(flatdir)
    files=glob.glob('*.tif')
    flats=files[-nflats:]
    N=len(flats)
    h,w=Image.open(flats[0]).size
    flatarr=np.zeros((w,h),np.float)
    #now read with PIL
    for flat in flats:
        imarr=np.array(Image.open(flat),dtype=np.float)
        flatarr=flatarr+imarr/N
    if saveim:
        out=Image.fromarray(flatarr,mode='L')
        out.save('flat.tif')
        out.show()
    return flatarr

def make_dark(darkdir, ndarks=10,saveim=False):
    '''Average all the darks'''
    os.chdir(darkdir)
    files=glob.glob('*.tif')
    darks=files[:ndarks]
    N=len(darks)
    h,w=Image.open(darks[0]).size
    darkarr=np.zeros((w,h),np.float)
    #now read with PIL
    for dark in darks:
        imarr=np.array(Image.open(dark),dtype=np.float)
        darkarr=darkarr+imarr/N
    if saveim:
        out=Image.fromarray(darkarr,mode='L')
        out.save('dark.tif')
        out.show()
    return darkarr

def correct_image(im,dark,flat,saveim=False, pickleimarr=False, contrast_stretch=False,show=False):
    imc=im-dark
    flatc=flat-dark
    corrected_image=imc/flatc
    if contrast_stretch:
        #contrast stretching
        p2 = np.percentile(corrected_image, 2)
        p98 = np.percentile(corrected_image, 98)
        corrected_image = exposure.rescale_intensity(corrected_image, in_range=(p2, p98))
        savec=scipy.misc.toimage(corrected_image,high=np.max(corrected_image),low=np.min(corrected_image),mode='F')
    else:
        savec=scipy.misc.toimage(corrected_image,high=np.max(corrected_image),low=np.min(corrected_image),mode='F')
        
    if show:
        fig,ax=plt.subplots()
        ax.imshow(savec,cmap=cm.gray)
        fig.show()
    if saveim:
        savec.save(saveim)
    if pickleimarr:
        pickle.dump(corrected_image,open(pickleimarr,'wb'))
       
    #return corrected_image

def contrast_stretch_group(imlist,cropim=False,saveim=False):
    '''contrast stretch by position group - take the absolute minimum 2% of all images in group, max 98%. input is list of flatdark corrected images'''
    #open all the images
    pixlist=[]
    for im in imlist:
        imarr=np.array(Image.open(im))
        if cropim:
            imarr=imarr[cropim[0]:cropim[1],cropim[2]:cropim[3]]
        pixlist.append(imarr)
    p2 = np.percentile(pixlist, 2)
    p98 = np.percentile(pixlist, 98)
    #print p2,p98
    #now scale all the images:
    for iname,imarr in zip(imlist,pixlist):
        corrected_image = exposure.rescale_intensity(imarr, in_range=(p2, p98))
        if saveim: #should also save the percentages somewhere... write to a pickle?
            savec=scipy.misc.toimage(corrected_image,high=np.max(corrected_image),low=np.min(corrected_image),mode='F')
            imname=iname[:-12]+'groupstretch.tif' 
            savec.save(imname)

def combine_flats(prior_dir,nprior,npost,post_dir=False):
    '''combine the flats from the post and prior measurements'''
    prior_flat=make_flat(prior_dir,nflats=nprior)
    if post_dir:
        post_flat=make_flat(post_dir,nflats=npost)
        combined_flat=(prior_flat+post_flat)/2.
    else:
        combined_flat=prior_flat
        
    return combined_flat

def combine_darks(prior_dir,nprior,npost,post_dir=False):
    '''combine the dark from the post and prior measurements'''
    prior_dark=make_dark(prior_dir,ndarks=nprior)
    if post_dir:
        post_dark=make_dark(post_dir,ndarks=npost)
        combined_dark=(prior_dark+post_dark)/2.
    else:
        combined_dark=prior_dark
            
    return combined_dark

if __name__ == '__main__':
    #prior_dir='/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/BeamTest_2018-04-11/SLS_Apr2018/disk2/moireref_prior_/tif'
    #post_dir='/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/BeamTest_2018-04-11/SLS_Apr2018/disk2/moireref_post_/tif'
    #prior_dir='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018/ref_prior_Trasm/tif'
    #post_dir='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018/ref_post_Trasm/tif'
    #combined_flat=combine_flats(prior_dir,100,100)
    #combined_dark=combine_darks(prior_dir,10,10)
    ##pickle these
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018/')
    #pickle.dump(combined_flat,open('transm_combined_flat.p','wb'))
    #pickle.dump(combined_dark,open('transm_combined_dark.p','wb'))
    combined_flat=pickle.load(open('transm_combined_flat.p','rb'))
    combined_dark=pickle.load(open('transm_combined_dark.p','rb'))
    ##now get images, correct and pickle them
    os.chdir('transmQM')
    wins=[12,34]#[11,12,21,22,31,32,33,34,41,42,43,44]
    for w in wins:
        files=glob.glob('*window0'+str(w)+'*.tif')
        for f in files:
            im=np.array(Image.open(f))
            cimname=f[:-4]+'_corrected.tif'
            cim=correct_image(im,combined_dark,combined_flat,saveim=cimname,contrast_stretch=True)
            fdimname=f[:-4]+'_flatdark.tif'
            fdim=correct_image(im,combined_dark,combined_flat,saveim=fdimname,contrast_stretch=False)
    
    #combined_flat=pickle.load(open('../transm_combined_flat.p','rb'))
    #combined_dark=pickle.load(open('../transm_combined_dark.p','rb'))
    #wins=[41,42,43,44]#[11,12,21,22,31,32,33,34,41,42,43,44]
    #files=glob.glob('*flatdark.tif')
    #os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018/tra')
    #for w in wins:
    #    #wfiles=glob.glob('win'+str(w)+'*flatdark.tif') #transm
    #    wfiles=glob.glob('*window0'+str(w)+'*flatdark.tif') #moire
    #    for p in range(0,7):
    #        if len(f) < 20:
    #            #print f
    #            imlist=[f for f in wfiles if 'win'+str(w) in f]
    #            contrast_stretch_group(imlist, saveim=True)
        #contrast_stretch_group(wfiles, saveim=True)
        #im=np.array(Image.open(f))
        #imname=f[:-4]+'_flatdark.tif'
        #imname2=f[:-4]+'_corrected.tif'
        #cim=correct_image(im,combined_dark,combined_flat,saveim=imname,contrast_stretch=False)
        #cim2=correct_image(im,combined_dark,combined_flat,saveim=imname2,contrast_stretch=True)
    
    
