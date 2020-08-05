"""
===================
analyze_general.py
Erica  Lastufka 5.9.18
===================

General methods for class Analyze. To be inherited by analyze_optical and analyze_xray, etc.
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
import csv

#####################################  IMAGE PREPROCESSING ###########################################

def im2ndarray(filen):
    '''Convert .tif image to numpy.ndarray'''
    imraw = Image.open(filen) #.rotate(45)#.convert('1').convert('L')
    im=np.array(imraw)
    return im

def remove_edges(filen,im,maxx='10',maxy='13',lpix=50,rpix=50,tpix=50,bpix=50,write=True):
    '''remove edges from images with edges of the window in them'''
    if filen.endswith('.p') or filen.endswith('mosaic.tiff'):
        imshape=np.shape(im)
        im=im[50:imshape[0]-50,50:imshape[1]-50]
        return im
    else:
        cropped=False
        windex,imindex=get_index(filen)
        if imindex[0] == '01' and lpix != False: #it's got an edge on the left, trim the array (x and y are flipped for some reason)
            im=im[:,lpix:]
            cropped=True
        if imindex[0] == maxx and rpix != False: #it's got an edge on the right, trim the array
            im=im[:,:-rpix]
            cropped=True
        if imindex[1] == '01' and tpix != False: #it's got an edge on the top, so mask the edges
            im=im[tpix:,:]
            cropped=True
        if imindex[1] == maxy and bpix != False: #it's got an edge on the bottom, so mask the edges
            im=im[:bpix,:]
            cropped=True
        #print np.shape(im)
    if write and cropped:
        os.rename(filen,filen[:-4]+'_uncropped.tif') #first rename old file
        marray=Image.fromarray(im) #raw, unequalised array
        marray.save(filen)
    return im

def remove_edges_general(filen,lpix=50,rpix=50,tpix=50,bpix=50,write=True):
    '''remove edges from images with edges of the window in them'''
    im=im2ndarray(filen)
    if lpix != False: #it's got an edge on the left, trim the array (x and y are flipped for some reason)
        im=im[:,lpix:]
        cropped=True
    if rpix != False: #it's got an edge on the right, trim the array
        im=im[:,:-rpix]
        cropped=True
    if tpix != False: #it's got an edge on the top, so mask the edges
        im=im[tpix:,:]
        cropped=True
    if  bpix != False: #it's got an edge on the bottom, so mask the edges
        im=im[:bpix,:]
        cropped=True
        #print np.shape(im)
    if write and cropped:
        os.rename(filen,filen[:-4]+'_uncropped.tif') #first rename old file
        marray=Image.fromarray(im) #raw, unequalised array
        marray.save(filen)
    return im

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

def contrast_stretch(imf,saveim=True,t=98,b=2,hist=False,show=False):
    if type(imf) == str:
        im=im2ndarray(imf)
    p2 = np.percentile(im, b)
    p98 = np.percentile(im, t)
    corrected_image = exposure.rescale_intensity(im, in_range=(p2, p98))
    savec=scipy.misc.toimage(corrected_image,mode='F') #want this to be greascale!
    imname=imf[:-4]+'_corrected.tif' #why -12?
    if saveim:
        savec.save(imname)
    if hist:
        fig,ax=plt.subplots(1,2)
        h1,_,_=ax[0].hist(im.ravel(),256)
        h2,_,_=ax[1].hist(corrected_image.ravel(),256)
        ax[0].set_title('original')
        ax[1].set_title('contrast stretched')
        for a in ax:
            a.set_xlim([0,256])
            a.set_ylim([0,np.max([np.max(h1),np.max(h2)])])
        fig.show()
    if show:
        fig,ax=plt.subplots(1,2,sharex=True,sharey=True)
        ax[0].imshow(im,cmap=cm.gray)
        ax[1].imshow(corrected_image,cmap=cm.gray)
        ax[0].set_title('original')
        ax[1].set_title('contrast stretched')

        fig.show()


    return imname

def correct_image(im,dark,flat,fac=False,saveim=False, pickleimarr=False, contrast_stretch=False,show=False,ret=False, hist=False,xraydemo=False):
    if type(im) == str:
        im=im2ndarray(im)
    imc=im-dark
    flatc=flat-dark
    if not fac:
        fac=1.0
    corrected_image=imc/(flatc*fac)
    if contrast_stretch:
        #contrast stretching
        p2 = np.percentile(corrected_image, 2)
        p98 = np.percentile(corrected_image, 98)
        stretched_image = exposure.rescale_intensity(corrected_image, in_range=(p2, p98))
        savec=scipy.misc.toimage(stretched_image,mode='F')
    else:
        savec=scipy.misc.toimage(corrected_image,mode='F')

    if xraydemo: #recreat fig3 from tomcat report
        fig,ax=plt.subplots(1,4)
        ax1=ax[0]
        ax2=ax[1]
        ax3=ax[2]
        ax4=ax[3]
        ax1.imshow(im)
        ax2.imshow(corrected_image)
        ax3.imshow(flatc)
        ax4.hist(im.ravel(),256,label='original',histtype='step')
        ax4.hist(corrected_image.ravel(),256,label='corrected',histtype='step')
        ax4.hist(flatc.ravel(),256,label='flat - dark',histtype='step')
        ax4.legend(loc='upper right')
        ax1.set_title('original')
        ax2.set_title('corrected')
        ax3.set_title('flat - dark')
        plt.axis('off')
        fig.show()
        show=False
        saveim=False #"xray_demo.png"
        pickleimarr=False
        ret=False

    if show:
        fig,ax=plt.subplots()
        ax.imshow(savec,cmap=cm.gray)
        if hist:
            ax1=plt.add_subplot(1,2,1)
            ax2=plt.add_subplot(1,2,2)
            ax1.hist(im.ravel(),256)
            ax2.hist(corrected_image.ravel(),256)
            ax1.set_title('original')
            ax2.set_title('corrected')
        fig.show()
    if saveim:
        savec.save(saveim)
    if pickleimarr:
        pickle.dump(corrected_image,open(pickleimarr,'wb'))
    if ret:
        return corrected_image

def check_correction():
    ''' check that the correction works for all the datasets'''
    from skimage import exposure
    dirs=['EMmodel/SLS_Apr2018']#,'QMmodel/SLS_May2018','QMmodel/SLS_Sept2018']
    #labs=['','',']
    fig,ax=plt.subplots(1,4)#,figsize=[10,5]) #[15,10]
    for d,a in zip(dirs,ax):
        os.chdir(d)
        flat=pickle.load(open('transm'+d[:2]+'_combined_flat.p','rb'))
        dark=pickle.load(open('transm'+d[:2]+'_combined_dark.p','rb'))
        if d =='QMmodel/SLS_Sept2018':
            flat =pickle.load(open('transmQM_flat_post.p','rb'))
        dfile=glob.glob('transm'+d[:2]+'/transm*.tif')[0]
        testim=im2ndarray(dfile)
        print dfile
        ax[0].imshow(testim,origin='bottom left')
        thist,tbc=exposure.histogram(testim-dark,nbins=256)
        fhist,fbc=exposure.histogram(flat-dark,nbins=256)
        tbinmax=tbc[np.where(thist == np.max(thist[70:]))[0]]
        fbinmax=fbc[np.where(fhist == np.max(fhist[70:]))[0]]
        fac=1.0#tbinmax/fbinmax
        print fac
        #thistr=np.reshape(thist, 2048)
        #tbcr=np.reshape(tbc,2040)
        #print np.shape(thistr),np.shape(tbcr),np.shape(testim)
        #a[0].plot(tbcr,thistr,lw=2)
        ax[0].set_title('original')#(d[8:]+' original')
        cim=correct_image(testim,dark,flat,fac=fac,ret=True,saveim='testim.tif')
        ax[1].imshow(cim,origin='bottom left')#,norm=matplotlib.colors.Normalize(vmin=.9*np.min(cim),vmax=1.1*np.max(cim)))
        chist,cbc=exposure.histogram(cim*256.,nbins=256)
        #print np.shape(chist),np.shape(cbc)
        #a[1].plot(cbc,chist,lw=2)
        ax[1].set_title('corrected')#(d[8:]+' corrected')
        ax[2].imshow(flat-dark,origin='bottom left')
        print np.shape(tbc),np.shape(fbc)
        #a[2].plot(fbc,fhist,lw=2)
        ax[2].set_title('flat - dark')#(d[8:]+' flat-dark')
        #hist, bins_center = exposure.histogram(camera)
        #plt.plot(bins_center, hist, lw=2)
        ax[3].plot(range(0,256),thist,'r',label='original')
        ax[3].plot(range(0,256),chist,'b',label='corrected')
        ax[3].plot(range(0,256),fhist,'g',label='flat-dark')
        #ax[3].set_xlim([0,256])
        ax[3].set_ylim([0,np.max(cim)])
        ax[3].legend(loc='upper right')

        for aa in ax:
            aa.axis('off')
        os.chdir('../../')
    fig.show()

def show_method():
    fig,ax=plt.subplots(1,4,figsize=[12,5])
    fd=im2ndarray('win11_p3_0_flatdark.tif')
    corr=im2ndarray('win11_p3_0_corrected.tif')
    gs=im2ndarray('win11_p3_0_groupstretch.tif')
    bmean=np.mean(fd)
    imgt=np.ma.masked_greater(fd,bmean) #what happens if the value equals the mean?
    imgtfilled=np.ma.filled(imgt, 1.0)
    imlt=np.ma.masked_less(imgtfilled,bmean)
    imltfilled=np.ma.filled(imlt, 0.0)
    imarr=imltfilled

    ax[1].imshow(imarr,origin='lower left')
    ax[1].set_title('binary')
    ax[0].imshow(fd,origin='lower left')
    ax[0].set_title('corrected')
    ax[2].imshow(gs,origin='lower left')
    ax[2].set_title('group stretched')
    ax[3].imshow(corr,origin='lower left')
    ax[3].set_title('contrast stretched')
    for a in ax:
        a.axis('off')
    fig.show()

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
    newimnames=[]
    #print p2,p98
    #now scale all the images:
    for iname,imarr in zip(imlist,pixlist):
        corrected_image = exposure.rescale_intensity(imarr, in_range=(p2, p98))
        if saveim: #should also save the percentages somewhere... write to a pickle?
            savec=scipy.misc.toimage(corrected_image,high=np.max(corrected_image),low=np.min(corrected_image),mode='F')
            imname=iname[:-12]+'groupstretch.tif'
            newimnames.append(imname)
            savec.save(imname)
    return newimnames

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

def flatdark(imlist, dark, flat,gw=False):
    newimnames=[]
    for im in imlist:
        if gw:
            newimname=im[:-4]+'_corrected.tif'
        else:
            newimname=im[:-4]+'_flatdark.tif'
        newimnames.append(newimname)
        correct_image(im, dark,flat, saveim=newimname)
    return newimnames

def groupstretch(imlist):
    #first get groups of p
    imlist.sort()
    p0=[p for p in imlist if 'p0' in p]
    p1=[p for p in imlist if 'p1' in p]
    p2=[p for p in imlist if 'p2' in p]
    p3=[p for p in imlist if 'p3' in p]
    p4=[p for p in imlist if 'p4' in p]
    p5=[p for p in imlist if 'p5' in p]
    p6=[p for p in imlist if 'p6' in p]

    pgroups=[p0,p1,p2,p3,p4,p5,p6]
    newimnames=[]
    for pimlist in pgroups:
        newimnames.extend(contrast_stretch_group(pimlist,saveim=True))
    return newimnames

def EM_list(win,mag,ending='.tif'):
    import glob
    filen=glob.glob('win'+win+'*_'+mag+'X'+ending)
    newfilen=[]
    for f in filen:
        if '01' not in f and '12' not in f[6:] and '09' not in f:
            newfilen.append(f)
    return newfilen

def get_index(filen):
    '''Parse filename to get indices of image in window as well as overall window index. This helps determine whether or not it includes the edge of the window.'''
    if len(filen) > 50: #it's the whole thing
        filen=filen[filen.find('win')-1:]
    windex=filen[filen.find('win')+3:filen.find('_')]
    indices=filen[filen.find('_')+1:filen.rfind('5.0')-1]
    if indices.endswith('_'):
        indices= indices[:-1]
    index=[indices[:2],indices[3:]]
    return windex,index

def get_index_general(filen,keywin='window',keyidx='_'):
    '''Parse filename to get indices of image in window as well as overall window index. This helps determine whether or not it includes the edge of the window.'''
    if len(filen) > 50: #it's the whole thing
        filen=filen[filen.find(keywin)-1:]
    windex=filen[filen.find(keywin)+6:filen.find(keywin)+9]
    indices=filen[filen.rfind(keyidx)-7:filen.rfind(keyidx)]
    if indices.endswith('_'):
        indices= indices[:-1]
    index=[indices[:2],indices[3:]]
    return windex,index

#################################################################################################################

def mosaic(window_num, coords, all=False, plot=True,plot_coords=False, mag=5.0):
    '''Make a mosaic of the given window or all windows (individually). Input: list of window numbers, dictionary of coordinates'''
    #how to deal with overlap? make left, right, top, bottom options and compare statistically effect on edge/line fits?

    #microscope settings: mm unless pixels specified
    if mag==5.0:
        mscope={'Magnification': 5.0,
            'FOVx': 1.2512,
            'FOVy' : 0.9383,
            'ExcessX_px' : 82.0,
            'ExcessY_px' : 82.0,
            'ExcessX' : 0.1603,
            'ExcessY' : 0.1603,
            'StepX' : 1.2316,
            'StepY' : 0.9246}
    else: #15X
        mscope={'Magnification': 15.0,
            'FOVx': 0.4165,
            'FOVy' : 0.3128,
            'ExcessX_px' : 253.,
            'ExcessY_px' : 255.0,
            'ExcessX' : -0.0046,
            'ExcessY' : 0.0039,
            'StepX' : 0.4163,
            'StepY' : 0.3127}


    pixX= mscope['FOVx']#1.955 #this is too big....
    pixY=mscope['FOVy']#1.9548 #from Matej's logs since mine don't have them

    if all:
        #generate list of window numbers from dictionary
        window_num=[w['number'] for w in windows]

    for win in window_num:
        #store filenames
        fnames = [im for im in coords.keys() if int(coords[im]['number'])== win]
        fnames.sort()
        #fnames = [im for im in fnames if coords[im]['indices'][0]== '01'] #for testing on first column

        #get the coordinates of each image
        cwinx = [coords[im]['im_coords'][0] for im in fnames] #image x coordinates in mm
        cwiny = [coords[im]['im_coords'][1] for im in fnames] #image y coordinates in mm
        idxx  = [coords[im]['indices'][0] for im in fnames]
        idxy  = [coords[im]['indices'][1] for im in fnames]
        #convert mm to pixels
        cpixx=np.array(cwinx)*(640.0/pixX)
        cpixy=np.array(cwiny)*(480.0/pixY)

        #translate coords so that top left coord is at [0,0]
        cminx=np.min(cpixx)
        cmaxy=np.max(cpixy)
        print cminx,cmaxy
        cpixx=np.rint(cpixx-cminx)#cpixx-cminx
        cpixy=np.rint(cmaxy-cpixy)#cmaxy-cpixy
        #print int(np.max(cpixx))+640,int(np.max(cpixy))+480

        #make new blank image
        background=Image.new("L",[int(np.max(cpixx)),int(np.max(cpixy))], 0xff) #[0,0] is TOP left #extra 9x12 px for rounding errors
        #background=Image.new("1",[int(np.max(cpixy))+480,int(np.max(cpixx))+640], "white") #[0,0] is TOP left
        #print 'size:', int(np.max(cpixx))+649,int(np.max(cpixy))+492,
        #put things in their place
        #fnames=[fnames[2],fnames[6]]
        mx,my,fpixx,fpixy=[],[],[],[]
        #determine the actual excess pixels
        exx=640.-np.mean([c2-c1 for c1,c2 in zip(cpixx[:-13],cpixx[13:])])
        exy=480.-np.mean([c2-c1 for c1,c2 in zip(cpixy[:-1],cpixy[1:]) if c2-c1 > 0])
        print 'excess ', exx,exy
        for i,f in enumerate(fnames):#zip([2,6],fnames):#enumerate(fnames):
            im=Image.open(f)
            offsetx=0#(int(idxx[i])-1)*int(exx)#0##int(idxx[i])*6
            offsety=0#(int(idxy[i])-1)*int(exy)
            fpixx.append(int(cpixx[i])-offsetx)
            fpixy.append(int(cpixy[i])-offsety)
            # "The box argument is either a 2-tuple giving the upper left corner,
            # a 4-tuple defining the left, upper, right, and lower pixel coordinate, or None (same as (0, 0))."
            #assuming the actual coordinates are the center of the image...
            #offsetx = 640/2
            #offsety= 480/2
            box=(int(cpixx[i])-offsetx,int(cpixy[i])-offsety)#,640+int(cpixx[i]),480+int(cpixy[i]))
            mx.append(box[0])
            my.append(box[1])
            background.paste(im,box)

        mcoords=[mx,my]
        if plot: #show the mosaic
            background.show()

        if plot_coords: #just plot the coords
            fig,ax=plt.subplots()
            #ax.scatter(cwinx,cwiny)
            ax.scatter(cpixx,cpixy)
            fig.show()

        #now convert the mosaic to numpy array and save. Also save a .tiff image
        marray=np.array(background.convert('L')) #raw, unequalised array
        filename='window'+str(win)+'mosaic_'+str(mag)
        background.save(filename+'.tiff')
        pickle.dump(marray,open(filename+'.p','wb'))
        pickle.dump(mcoords, open(filename+'_coords.p','wb'))

    #return new list of filenames
    return fnames,cpixx,cpixy,idxx,idxy,fpixx,fpixy#filename+'.p'
    #return mfiles
    #return fnames,cpixx,cpixy

def numpy_mosaic(filenames, ycoords=False, pix2um=6.5,alpha=.75,xoff=-1,bpix=3,bw=False):
    '''1000 um / 6.5 um/pix =
    this still doesn't really work well :( '''
    if not ycoords:
        ycoords=range(0,len(filenames))

    ysize,xran=np.shape(im2ndarray(filenames[0]))
    #calculate y-dim of new mosaic
    ydiff=int((float(ycoords[1]-ycoords[0])*1000.)/pix2um) #mm 2 pixels
    yran=ydiff*(len(filenames)+1)
    nim=len(filenames)-1
    newarr=np.zeros([yran,xran])
    #make the np array
    arrlist=[]
    for f in filenames:
        arrlist.append(im2ndarray(f))
    for a,b,c in zip(arrlist[:-1],arrlist[1:],range(0,nim)):
        overlap_a=a[bpix:ydiff,:] #ex. a[3:154,:]
        #rest_a=a[ydiff:-bpix,:]  #ex. a[154:-3,:]
        rest_a=a[ydiff:,:]  #ex. a[154:-3,:]
        #overlap_b=b[-ydiff:-bpix,:] #ex. b[-154:-3,:]
        overlap_b=b[-ydiff:-bpix,:] #ex. b[-154:-3,:]
        rest_b=b[bpix:-ydiff,:] #ex. b[3:154,:]
        mean_overlap=(overlap_a + overlap_b)/2.
        extent_rest_a=np.shape(rest_a)[0]
        extent_rest_b=np.shape(rest_b)[0]
        extent_overlap=np.shape(mean_overlap)[0]
        c1=yran-(c+1)*extent_rest_a
        c2=yran-c*extent_rest_a
        #c3=yran-(c+3)*extent_overlap
        c4=c1-extent_overlap
        c5=c4-extent_rest_b
        #c5=yran-(c+2)*(ydiff-3)-np.shape(rest_b)[0]
        #c6=yran-(c+2)*(ydiff-3)
        newarr[c1:c2,:xran-c*xoff]=rest_a[:,c*xoff:]
        #newarr[yran-(c+1)*(ydiff-2)-1:yran-c*np.shape(rest_a)[0]-2,:]=rest_a
        #print c1,c2
        newarr[c4:c1,:xran-(c+1)*xoff]=mean_overlap[:,(c+1)*xoff:]
        #newarr[yran-(c+2)*(ydiff-2):yran-(c+1)*(ydiff-2)-1,:]=mean_overlap
        #print c4,c1
        newarr[c5:c4,:xran-(c+2)*xoff]=rest_b[:,(c+2)*xoff:]
        #print c5,c4
        #print ''
        #get mean of overlapping region
    fig,ax=plt.subplots(5,1)
    if bw:
        cmap=cm.Greys
    else:
        cmap=cm.Rainbow
    norm=matplotlib.colors.Normalize(vmin=10*np.min(newarr),vmax=.8*np.max(newarr))
    ax[0].imshow(rest_a,origin='lower left')
    ax[1].imshow(overlap_a,origin='lower left')
    ax[2].imshow(overlap_b,origin='lower left')
    ax[3].imshow(mean_overlap,origin='lower left')
    ax[4].imshow(rest_b,origin='lower left')
    print np.shape(rest_a),np.shape(mean_overlap),np.shape(rest_b)
    fig.show()

    fig,ax=plt.subplots()
    ax.set_ylim([0,yran])
    ax.set_xlim([0,xran])
    ysize=284
    for c,f in zip(ycoords,arrlist):
        #trim top and bottom 2px of each image - it's blurry
        newf=f[3:-3,:]
        ax.imshow(newf,extent=(0,xran,(nim-c)*(ydiff),(nim-c)*(ydiff)+ysize),alpha=alpha, origin='upper left',cmap=cmap,norm=norm)
        print np.shape(newf)
        print 0,xran,(nim-c)*(ydiff),(nim-c)*(ydiff)+ysize
    #ax.imshow(newarr)
    fig.show()


def angle_by_index(hough_lines, xpeak_int=False,imfilename=False,plot_idx=False, plot_angle=False,ptol=2,alpha=0.6):
    ''' Get line orientation as a function of line index. ie. layer in Gridwinder grids'''
    #determine which line goes with which xpeak to within ptol px
    theta_dict={}
    xidx,alltheta=[],[]
    if type(xpeak_int) == str: #take 3 'strips' of the images and get the xpeaks from there
        edges=pickle.load(open(xpeak_int,'rb'))
        xpeak_low = sum_peaks(clean_centers(edges[10:20,:],sigma=1.25),tol=4)
        xpeak_med= sum_peaks(clean_centers(edges[130:150,:],sigma=1.25),tol=4)
        xpeak_hi = sum_peaks(clean_centers(edges[270:280,:],sigma=1.25),tol=4)
        xl_int=[int(x) for x in xpeak_low]
        xm_int=[int(x) for x in xpeak_med]
        xh_int=[int(x) for x in xpeak_hi]
        xl_int.sort()
        xm_int.sort()
        xh_int.sort()
        #print len(xl_int),len(xm_int),len(xh_int) # if they're not the same length there's a problem....
        #xpeak_int=[xl_int,xm_int,xh_int]
    else:
        xpeak_int.sort()
    hough_lines.sort()
    for lidx,l in enumerate(hough_lines):
        ep1=l[0] #higher
        ep2=l[1] #lower
        #where are the endpoints located relative to the slices?
        #endpoint 1:
        for ep in [ep1,ep2]:
            if ep[1] >=200:
                idx,xpval=find_nearest(xh_int,ep[0])
            elif  ep[1]<=100:
                idx,xpval=find_nearest(xl_int,ep[0])
            else:
                idx,xpval=find_nearest(xm_int,ep[0])

            if np.abs(ep[0]-xpval) < ptol:
                theta=np.abs(get_angle(l))
                theta_dict[lidx]=[idx,theta,xpval,l]
                xidx.append(idx)
                alltheta.append(theta)
                break

    onex,onetheta=range(0,len(xl_int)),[]
    xiter=iter(xidx)
    thiter=iter(alltheta)
    for o in onex:
        try:
            ctheta=[thiter.next()]
            while xiter.next() == o:
                ctheta.append(thiter.next())
        except StopIteration:
            onetheta.append(np.mean(ctheta))
            break
        onetheta.append(np.mean(ctheta))
        #print o,np.mean(ctheta)

    if plot_idx and imfilename: #plot image annotated with text describing line index and label
        yidx=50
        yang=150
        fig,ax=plt.subplots()
        ax.imshow(im2ndarray(imfilename),origin='lower left',alpha=alpha)
        for k in theta_dict.keys():
            ax.text(theta_dict[k][2],yidx,k, rotation = 90.) #plot index
            ax.text(theta_dict[k][2],yang,np.round(theta_dict[k][1],2),rotation=90.) #plot index
        fig.show()

    if plot_angle: #plot angle vs line index
        fig,ax=plt.subplots()
        ax.scatter(xidx,np.abs(alltheta))
        #ax.plot(xidx,np.abs(alltheta))
        ax.scatter(onex,onetheta,c='r',marker='+')
        ax.set_xlabel('Edge index')
        ax.set_ylabel('Orientation of Edge (degrees)')
        ax.set_xlim([0,len(xl_int)])
        fig.show()
    return theta_dict,onex,onetheta

def plot_theta_by_idx(window,nmin=0,nmax=False,line=False,write=False):
    edges=glob.glob('window00'+str(window)+'*corrected_edges.p')
    hough=glob.glob('window00'+str(window)+'*hough.p')
    ad,axx,colors=[],[],[]
    c18=['k','r','g','c','m','y','b','grey','lightcoral','darkorange','yellow','lawngreen','turquoise','dodgerblue','blueviolet','orange','purple','navy']
    if not nmax:
        nmax=len(edges)
    esub=edges[nmin:nmax]
    hsub=hough[nmin:nmax]
    for e,h in zip(esub,hsub):
        ld,xidx,atheta=angle_by_index(pickle.load(open(h,'rb')),xpeak_int=e)
        ad.append(atheta)
        axx.append(xidx)

    nn=len(ad)
    allarr=np.zeros([nn,len(xidx)])
    for n in range(0,nn):
        allarr[n,:]=np.array(ad[n])
        colors.append(c18[n])#(float(n)/(384./float(nn)),(float(nn-n)/(255./float(nn))),(float(n)/(255./float(nn)))))
    #print colors
    #colors=np.random.random((len(ad),3))
    meanarr=np.mean(allarr,axis=0)
    fig,ax=plt.subplots()
    for i,x,a,c in zip(range(0,nn),axx,ad,colors):
        ax.scatter(x,a,alpha=.6,c=c,label='Image '+str(nmin+i))
        if line:
           ax.plot(x,a,c=c)
    ax.scatter(x,meanarr,marker='+',c='r',label='Mean')
    ax.plot(x,meanarr,c='r',linewidth=4)
    ax.set_xlim([0,len(x)])
    ax.set_ylim([87.5,90.1])
    ax.set_xlabel('Edge index')
    ax.set_ylabel('Angle (degrees)')
    ax.legend(loc='lower right')
    ax.set_title('Window ' + str(window))
    fig.show()

    if write: #write to a csv
        rows=[]
        row0=['Edge Index'] #header row
        for n in range(nmin,nmax):
            row0.append('Image '+str(n))
        row0.append('Mean')
        rows.append(row0)
        transarr=np.transpose(allarr)
        for i in range(0,len(meanarr)):
            foo=['Edge ' + str(i)]
            foo.extend(list(transarr[i,:]))
            foo.append(meanarr[i])
            rows.append(foo)
        with open('window00'+str(window)+'angle_idx.csv',mode='w') as outfile:
            writer=csv.writer(outfile, delimiter=',')#,quotechar='"')
            for row in rows:
                writer.writerow(row)

def anisotropic_diffusion(filepath):
    Dimension = 2
    PixelType = itk.ctype('float')
    ImageType = itk.Image[PixelType, Dimension]
    reader = itk.ImageFileReader[ImageType].New(FileName=filepath)
    reader.Update()
    image = reader.GetOutput()

    # Crop to region of interest
    roi_filter = itk.RegionOfInterestImageFilter.New(image)
    region = itk.ImageRegion[Dimension]()
    index = itk.Index[Dimension]()
    index[0] = 0
    index[1] = 0
    region.SetIndex(index)
    size = itk.Size[Dimension]()
    size[0] = 640
    size[1] = 480
    region.SetSize(size)
    roi_filter.SetRegionOfInterest(region)
    roi_filter.Update()
    roi = roi_filter.GetOutput()

    smoother = itk.CoherenceEnhancingDiffusionImageFilter.New(roi)
    # Determine good parameters here: https://insightsoftwareconsortium.github.io/ITKAnisotropicDiffusionLBR/
    smoother.SetDiffusionTime(5)
    smoother.SetLambda(0.05)
    smoother.SetEnhancement(3)
    smoother.SetNoiseScale(2)
    smoother.SetFeatureScale(10)
    smoother.SetExponent(2)
    smoother.Update()

    smootherarr=itk.GetArrayFromImage(smoother)
    return smootherarr


def Canny_edge(filen,sigma=3,mag=5.0,anisotropic=False,binary=False,gauss=False,plot=False,outfilen=False):
    #max contrast
    if filen.endswith('.p'):
        imarr=pickle.load(open(filen,'rb'))
    else:
        imarr=im2ndarray(filen)

    #im = exposure.equalize_hist(imarr)    #should probably do this AFTER removing the edges...
    #if mag==5.0:
    #    im=remove_edges(filen,imarr)
    #elif mag ==15.0:
    #    im=remove_edges(filen,imarr,maxx='12',maxy='16')
    #else:
    im=imarr

    #for the finest grids, use sigma=2
    #windex,foo=get_index(filen)
    #if windex in ['32','42','34','44']: #although I should actually use the dictionary keys to do this stuff in case we change the naming convention
    #    sigma=2
    #if windex in ['11','21']: #although I should actually use the dictionary keys to do this stuff in case we change the naming convention
    #    sigma=4
    #    gauss=3
    #if windex in ['31','41']: #although I should actually use the dictionary keys to do this stuff in case we change the naming convention
    #    sigma=3
    #    #gauss=3

    if anisotropic:
        im=anisotropic_diffusion(filen)

    #p2 = np.percentile(im, 2)
    #p98 = np.percentile(im, 98)
    #im = exposure.rescale_intensity(im, in_range=(p2, p98))

    if gauss:
        im = ndi.gaussian_filter(im, gauss)

    if binary: #classify into 1 or 0
        immin=np.min(im)
        imsig=np.std(im)
        #masked_above=np.ma.masked_greater(im,immin+3*imsig)
        oneim=np.ma.ones(np.shape(im))
        masked_below=np.ma.masked_less(im,immin+imsig) #these become zeros
        oneim.mask=masked_below.mask
        im=np.ma.filled(oneim, 0.0)
        #p2 = np.percentile(im, 2)
        #p98 = np.percentile(im, 98)
        #im = exposure.rescale_intensity(im, in_range=(p2, p98))

    edges = feature.canny(im, sigma=sigma)

    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        #ax = ax.ravel()

        ax.imshow(im, cmap=cm.gray)
        ax.imshow(np.ma.masked_where(edges == 0,edges),cmap=cm.autumn)
        ax.set_title('Input image overlaid with Canny edges')
        ax.set_axis_off()
        fig.show()

    if not outfilen:
        newfilen=filen[:-4]+'_edges.p'
    else:
        newfilen=outfilen
    pickle.dump(edges,open(newfilen,'wb'))
    return edges,newfilen

def clean_centers(edges,cfac=False,tolerance=False,plot=False,sigma=2.):
    '''clean out the Canny edge array so that pixel bins with low counts (ie probably bad edges) get deleted below a certain threshold'''
    esum=np.sum(edges,axis=0)
    cleaned_edges=np.copy(edges)
    #print sigma,np.std(esum)
    if not tolerance:
        tolerance=np.mean(esum)+sigma*np.std(esum)
    #print tolerance
    if type(cfac) !=bool:
        tolvec = sigma*np.std(esum)*cfac/np.max(cfac) #numpy 1-D array
        for i,col in enumerate(esum):
            if col < tolvec[i] and col != 0:
                #print i, col
                cleaned_edges[:,i] = False
    else:
        for i,col in enumerate(esum):
            if col < tolerance and col != 0:
                #print i, col
                cleaned_edges[:,i] = False
    if plot:
        cleaned_sum=np.sum(cleaned_edges,axis=0)
        fig,ax=plt.subplots()
        ax.plot(range(0,len(esum)),esum, label='original')
        ax.plot(range(0,len(cleaned_sum)),cleaned_sum,label='cleaned')
        if type(cfac) != bool:
            ax.plot(range(0,len(tolvec)),tolvec,label='tolerance')
        else:
            ax.axhline(tolerance,c='k',linestyle='--')
        ax.legend(loc='upper right')
        if type(cfac) != bool:
            ax.set_xlim([0,len(cfac)])
        fig.show()
    return cleaned_edges

def plot_centers_and_edges(win,p0,ang,earr=False, datafile=False,tol=2.0,cfac=False):
    fname='win'+str(win)+'_p'+str(p0)+'_'+ang+'_corrected_edges.p'
    if not datafile:
        datafile='win'+str(win)+'_width_data_p'+str(p0)+'_'+ang+'.p'
    dd=pickle.load(open(datafile,'rb'))
    if 'earr' in dd.keys():
        edges=dd['earr']
        carray=rotate(np.transpose(np.ones(dd['imshape'])),dd['rang'], reshape=True)
        cfac=np.sum(carray,axis=0)
    elif type(earr) == bool:
        edges=pickle.load(open(fname,'rb'))
    else:
        edges=earr
    cleaned_edges=clean_centers(edges,sigma=tol,cfac=cfac)
    esum=np.sum(edges,axis=0)
    tolerance=np.mean(esum)+tol*np.std(esum)
    if type(cfac) !=bool:
        tolvec = tol*np.std(esum)*cfac/np.max(cfac) #numpy 1-D array

    cleaned_sum=np.sum(cleaned_edges,axis=0)
    rising=dd['rising']
    falling=dd['falling']

    fig,ax=plt.subplots()
    ax.step(range(0,len(esum)),esum, linewidth=1,label='original',c='m')
    ax.step(range(0,len(cleaned_sum)),cleaned_sum,linewidth=2,label='cleaned',c='g')
    if type(cfac) != bool:
        ax.plot(range(0,len(tolvec)),tolvec,label='tolerance')
        ax.scatter(rising,tolerance*np.ones(len(rising)),c='r',marker='v',s=60,label='rising')
        ax.scatter(falling,tolerance*np.ones(len(falling)),c='k',marker='o',s=50,label='falling')
    else:
        ax.axhline(tolerance,c='k',linestyle='--')
        ax.scatter(rising,tolerance*np.ones(len(rising)),c='r',marker='v',s=60,label='rising')
        ax.scatter(falling,tolerance*np.ones(len(falling)),c='k',marker='o',s=50,label='falling')
    ax.legend(loc='upper right')
    ax.set_xlim([0,len(cleaned_sum)])
    ax.set_ylim([0,np.max(esum)])
    fig.show()


def sum_peaks(cleaned_edges, tol=9,irange=False):#irange=[1087,1105]):
    '''Get the locations of the peaks of a (Gaussian?) fit to the histogram of sums, return these as the width array. Tolerance is number of False allowed between Trues to still be considered a group'''
    esum=np.sum(cleaned_edges,axis=0)
    if irange:
        esum=esum[irange[0]:irange[1]]
    group_hist,group_xvec=[],[]
    #first get the groups of True values
    group=False
    for i,es in enumerate(esum[:-1]):
        #if irange:
            #print i+irange[0],es,group,group_hist,group_xvec
        if es != 0.0 and not group: #start the group for the first non-zero value
            lhist=[es]
            lxvec=[i]
            group=True
        elif es !=0.0 and group: #if aready in group, stay in group unless following tol pixels are all zeros
            lhist.append(es)
            lxvec.append(i)
            #print 'here',esum[i+1]
            #if np.mean(esum[i+1:i+tol]) == 0 :
                #print 'here'
            #    group=False
        elif es == 0.0 and group and i-lxvec[0]<tol and np.mean(esum[i:i+tol]) !=0.0: #if there is a zero but it's inside a group, keep it
            lhist.append(es)
            lxvec.append(i)
        else:
            group=False
        if not group:
            try:
                if lxvec not in group_xvec:
                    group_hist.append(lhist)
                    group_xvec.append(lxvec)
            except NameError:
                continue
    #this actually doesn't work well because there are so few points they don't make a normal distribution! So let's take the weighted average instead
    xpeak=[]
    #print group_xvec,group_hist
    for xg,yg in zip(group_xvec,group_hist):
        xpeak.append(np.average(xg, weights=yg))
        #if len(xg) > 2: #can't make a gaussian out of only 2 points...
        #    xpeak.append(gauss_peak(xg,yg))
        #else:
        #    xpeak.append(np.mean(xg))
    #if irange:
    #    print xpeak

    return xpeak#,group_xvec,group_hist

def gauss_peak(xdata,ydata,npts=50,show=False,ret_data=False):
    ''' returns peak of gaussian fit to data'''
    #x-coord does not have to be an integer!
    # Create a function which returns a Gaussian (normal) distribution.
    def gauss(x,*p0):
        if len(p0) == 1:
            a,b,c,d=p0[0]
        else:
            a,b,c,d=p0
        y = a*np.exp(-np.power((x - b), 2.)/(2. * c**2.)) + d
        return y
    #interpolate x and y so greater than 4 (or else curvefit doesn't work)
    if len(xdata) < 4:
        intfunc=interpolate.interp1d(xdata,ydata)
        xint=np.linspace(xdata[0],xdata[-1],10)
        yint=intfunc(xint)
    else:
        xint=np.array(xdata)
        yint=np.array(ydata)
    popt, pcov = curve_fit(gauss, xint, yint,p0=[0,0,0,np.mean(yint)])#, p0=p_initial, sigma=e)
    #generate high-resolution version of gaussian curve and get the maximum:
    xvec=np.linspace(xdata[0],xdata[-1],npts)
    yvec=gauss(xvec,popt) #[0]
    xpeak=xvec[np.where(yvec == np.max(yvec))[0]][0]
    if show:
        fig,ax=plt.subplots()
        ax.plot(xint,yint)
        ax.plot(xvec,yvec)
        fig.show()
    if ret_data:
        return xvec,yvec
    else:
        return xpeak #,ypeak

def get_im_rot_angle(edges,testrot=True):
    '''calculate the angle of rotation of the image by comparing the first and last row of cleaned edges'''
    row1=edges[5]
    row2=edges[-5]
    aa=np.where(row1 == True)
    bb=np.where(row2==True)
    offset=[]
    for a in aa[0]:
       idx=np.abs(bb[0]-a).argmin()
       offset.append(bb[0][idx]-a)
    print np.mean(offset)
    print np.median(offset)
    print np.std(offset)
    print offset
    #image height is 2038 (counting where the pixels are taken from)
    theta=np.arctan2(np.std(offset),2038.)
    print np.rad2deg(theta)

    if testrot: #rotate the edge array and see what happens to the offset
        rotedges=imrotate(edges,-1*np.rad2deg(theta))
        fig,ax=plt.subplots()
        ax.imshow(rotedges,origin='lower left')
        fig.show()
        #convert array to bool
        row1=rotedges[5]
        row2=rotedges[-5]
        r1mean=np.mean(row1)
        r2mean=np.mean(row2)
        row1a=[i for i,r in enumerate(row1) if r >=r1mean]
        row2a=[i for i,r in enumerate(row2) if r >=r2mean]
        aa=np.array(row1a)#np.where(row1b == True)
        bb=np.array(row2a)#np.where(row2b==True)
        offset=[]
        for a in aa:
           idx=np.abs(bb-a).argmin()
           offset.append(bb[idx]-a)
        print np.mean(offset)
        print np.median(offset)
        print np.std(offset)
        print offset
        #image height is 2038 (counting where the pixels are taken from)
        theta=np.arctan2(np.median(offset),2038.)
        print np.rad2deg(theta)
        fig,ax=plt.subplots()
        ax.step(range(0,2040),np.sum(edges,axis=0),label='original')
        ax.step(range(0,2040),np.sum(rotedges,axis=0)/256.,label='rotated')
        ax.legend(loc='upper right')
        fig.show()
        return aa,bb,rotedges

def imgrad_from_file(imfile,isum=False, imean=False):
    im=np.array(Image.open(imfile))
    imsum=np.sum(im,axis=0)
    immean=np.mean(imsum)
    imgrad=np.gradient(imsum)
    ret={"imgrad":imgrad,"imsum":False,"immean":False}
    if isum:
      ret["imsum"]=imsum
    if imean:
      ret["immean"]=immean
    return ret

# def edges_to_centers(cleaned_edges,imfile,n=9):
#     '''Additional function to be executed before sum_peaks. Finds each slit center from the edge pairs'''
#     #do I need to run rising_or_falling first then? probably....
#     xpeak_int=[int(x) for x in xpeak]
#     xpeak_int.sort()

#     imdict=imgrad_from_file(imfile,imean=True)
#     imgrad=imdict["imgrad"]

#     rising, falling,_,_,_=rising_or_falling_final(xpeak,xpeak_int,imgrad,n=n)

#     #maybe just mod rising_or_falling_final to return centers as well


def get_xpeaks(rang,rotim,rotarr,im,ef,sigma=2.,tol=9,plot=False):
    #now clean edges of the rotated arrays?
    carray=np.ones(np.shape(rotim))
    carr=np.sum(carray,axis=0)
    cleaned_edges=clean_centers(rotarr,cfac=carr, sigma=sigma,plot=plot) #make the cutoff line follow a curve determined by the local maxima? Or smooth the array ... idk
    #then group edges together using group_edges_by_idx (I think)
    xpeak=sum_peaks(cleaned_edges, tol=tol,irange=False)
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


def slit_widths_from_peaks(window_num,imfile,xpeak=False,pix2um=.65,plot=False,stats=True,gauss=False,tolerance=False,n=9,quiet=True,example=True,sigma=2.,rang=False):
    '''basically same as slit_or_slat() plus histogramming'''
    #if not ang:
    check=False
    p0ang=imfile[6:imfile.rfind('_')]
    #else:
    #    p0ang=ang
    imdict=imgrad_from_file(imfile,imean=True)
    #imsum=np.sum(im,axis=0)
    immean=imdict["immean"]
    imgrad=imdict["imgrad"]
    #imgrad2=np.gradient(imgrad)
    #gradpeaks=(np.roll(np.sign(imgrad2),1)-np.sign(imgrad2) !=0).astype(int) #locations where the gradients peak
    #gx=np.where(gradpeaks !=0)[0]
    if (example and p0ang =="p5_0") or plot:
        tplot=True
    else:
        tplot=False
    if not xpeak:
        efile=glob.glob(imfile[:-4]+'_edges.p')
        if len(efile) > 0:
            edges=pickle.load(open(efile[0],'rb'))
        else:
            edges,_=Canny_edge(imfile,sigma=3,gauss=gauss,plot=False,mag=False)
        #cleaned_edges=clean_centers(edges,sigma=tolerance)
        #xpeak=sum_peaks(cleaned_edges,tol=n,plot=plot)
    if rang:
        rotim=rotate(im2ndarray(imfile),rang,reshape=True)
        rotarr=rotate(edges,rang,reshape=True)
    else:
        rotim=im2ndarray(imfile)
        rotarr=edges
    xpeak,xpeak_int=get_xpeaks(rang,rotim,rotarr,imfile,efile,sigma=sigma,tol=n,plot=tplot)
    #xpeak_int=[int(x) for x in xpeak]
    xpeak_int.sort()

#     #insert dummy peaks based on gx. EXCLUDE these from the rising/falling analysis
#     #check for values in gx that are within 2 pixels of xpeak_int
#     xpeak_arr=np.array(xpeak_int)
#     mask,allpeaks=[],[]
#     for gxv in gx:
#         absdiff=np.abs(xpeak_arr -gxv)
#         nearest_xpeak=xpeak_int[absdiff.argmin()]
#         if absdiff.min() >n and gxv not in allpeaks:     #for those without a corresponding xpeak, insert a 'dummy' peak into the list...with a corresponding mask so we know it's a dummy
#             allpeaks.append(gxv)
#             mask.append(False)
#         elif nearest_xpeak not in allpeaks:
#             allpeaks.append(nearest_xpeak)
#             mask.append(True)

    #determine if edges are rising or falling
    #rising,falling,periodr,periodf,period,width,cIvals,cAvals

    rising, falling, periodr,periodf,period,width,cIvals,cAvals=rising_or_falling_final(xpeak,xpeak_int,imgrad,n=n,centers=True,xray=True)
    if tplot:
        from analyze_optical import plot_peaks_and_points
        plot_peaks_and_points(xpeak,rising, falling, cIvals,cAvals)
    if plot:
        #make the histogram
        fig,ax=plt.subplots()
        bins=np.arange(np.min(period),np.max(period),np.std(period)/5.)
        ax.hist(period,bins)
        #ax.set_xlim([nperiod-5,nperiod+5])
        ax.set_yscale('log')
        ax.set_ylim([1,10000])
        #figfilename='win'+str(window_num)+'_group_periods'+str(mag)+'.png'
        #fig.savefig(figfilename)
        fig.show()

    if stats: #print out and pickle stats
        avg=np.mean(period) #here period is slit width ...what if I want to plot the actual period?
        med=np.median(period)
        stdv=np.std(period)
        avgw=np.mean(width)
        medw=np.median(width)
        stdvw=np.std(width)
        if stdv >=1.:
            check=True
        statfile='win'+str(window_num)+'_width_stats_'+p0ang+'.p'
        datafile='win'+str(window_num)+'_width_data_'+p0ang+'.p'
        print "-------------------STATISTICS FOR WINDOW "+str(window_num)+"---------------------"
        print '              Mean P: ' + str(avg)
        print '            Median P: ' + str(med)
        print 'Standard Deviation P: ' + str(stdv)
        print '              Mean W: ' + str(avgw)
        print '            Median W: ' + str(medw)
        print 'Standard Deviation W: ' + str(stdvw)
        print 'Results saved in ' + statfile
        data={'period':period,'rising':rising,'falling':falling,'widths':width,'cIvals':cIvals,'cAvals':cAvals}
        stdict={'meanp':avg,'medianp':med,'stddevp':stdv,'meanw':avgw,'medianw':medw,'stddevw':stdvw}
        pickle.dump([stdict,tolerance,n],open(statfile,'wb'))
        pickle.dump(data,open(datafile,'wb'))

    return periodr,periodf,period,width,cAvals,cIvals,check

def rising_or_falling_final(xpeak,xpeak_int, imgrad,n=9, quiet=True,filter_multiples=True,filter_nominal=False, centers=False,xray=False):
    if xray: #this is reverse for x-ray and optical; for xray, light indicates slat. therefore rising => dark to light
        rising=[xpeak[i] for i in range(0,len(xpeak)-1) if np.mean(imgrad[xpeak_int[i]-1:xpeak_int[i]+2]) > 0.] #take 3 pix around xpeak
    #need xpeaks to be integers now
        falling=[xpeak[i] for i in range(0,len(xpeak)-1) if np.mean(imgrad[xpeak_int[i]-1:xpeak_int[i]+2]) < 0.]
    else:
        rising=[xpeak[i] for i in range(0,len(xpeak)-1) if np.mean(imgrad[xpeak_int[i]-1:xpeak_int[i]+2]) < 0.] #take 3 pix around xpeak
    #need xpeaks to be integers now
        falling=[xpeak[i] for i in range(0,len(xpeak)-1) if np.mean(imgrad[xpeak_int[i]-1:xpeak_int[i]+2]) > 0.]

    if not quiet:
        print "First falling: ",falling[0]
        print "First rising: ",rising[0]
        print "First 10 falling: ",falling[:10]
        print "First 10 rising: ",rising[:10]
        #print "First 10 falling masks: ",falling_mask[:10]
        #print "First 10 rising masks: ",rising_mask[:10]
        print 'f[0]-r[0]: ',falling[0]-rising[0]
        print 'f[1]-r[0]: ',falling[1]-rising[0]

    #edit to ingore falling values without a rising value between them and vice versa. This should help cut down on some errant widths at least. (should I allow these to stay for the period calculations though?)
    #if falling[0] > rising[0]:
        #not so fast....
        #width=[f-r for r,f in zip(rising,falling)]#falling-rising

    #else: #first one is falling
        #width=[f-r for r,f in zip(rising,falling[1:])]

    flist=[0 for f in falling]
    rlist=[1 for r in rising]
    rlist.extend(flist)
    masterlist=rising
    masterlist.extend(falling)
    locs, rorf = (list(t) for t in zip(*sorted(zip(masterlist, rlist)))) #list of x-coordinates and list of rising or falling codes (1 or 0)

    width,cIvals,cAvals,periodr,periodf=[],[],[],[],[]
    #now get the widths and centers
    for i,loc,code in zip(range(0,len(locs[:-2])),locs[:-2],rorf[:-2]):
        testsum=code+rorf[i+1]
        if testsum == 2 and locs[i+1]-loc >n: #two r's
            periodr.append(locs[i+1]-loc)
        elif testsum == 1: #r and f (but in what order?)
            if code ==1: #r then f
                width.append(locs[i+1]-loc)
                if centers:
                    cAvals.append((locs[i+1]-loc)/2. + loc) #slat centers, not slit centers
                if code+rorf[i+2]==2 and locs[i+2]-loc >n:
                    periodr.append(locs[i+2]-loc)
            elif code+rorf[i+2]==0 and locs[i+2]-loc >n: #f then r then f
                periodf.append(locs[i+2]-loc)
                if centers:
                   cIvals.append((locs[i+1]-loc)/2. + loc) #slit centers
        elif testsum == 0 and locs[i+1]-loc >n: #f and f
            periodf.append(locs[i+1]-loc)

    #now how to calculate the periods for the ones that have r's and f's between them anyway??

    #periodr= [rising[i+1] - rising[i] for i in range(0,len(rising)-1)]#rising-rising
    #periodf=[falling[i+1] - falling[i] for i in range(0,len(falling)-1)]#falling - falling
    period=periodr#np.mean(periodr+periodf)
    period.extend(periodf)

    if filter_multiples: #mask widths and periods that are 1.5 or more times the median value
        pmed=np.median(period)
        period=np.ma.masked_greater(period,1.5*pmed)
        period=np.ma.masked_less(period,0.5*pmed)
        wmed=np.median(width)
        width=np.ma.masked_greater(width,1.5*wmed) #should I also mask when less than?
        width=np.ma.masked_less(width,0.5*wmed) #should I also mask when less than?

    if filter_nominal != False: #mask widths and periods that are 1.5 or more times the median value
        #print filter_nominal
        period=np.ma.masked_greater(period,1.5*filter_nominal)
        period=np.ma.masked_less(period,0.5*filter_nominal)
        width=np.ma.masked_greater(width,.75*filter_nominal) #should I also mask when less than?
        width=np.ma.masked_less(width,0.25*filter_nominal) #should I also mask when less than?

    return rising,falling,periodr,periodf,period,width,cIvals,cAvals

def do_filter_nominal(period,npitch,fac=2.):
    '''periods already in um '''
    period.sort()
    try:
        lbound=np.where(np.array(period) < (1.-(fac-1.))*npitch)[0][-1]
    except IndexError:
        lbound=0
    try:
        ubound=np.where(np.array(period) > fac*npitch)[0][0]
    except IndexError:
        ubound=-1
    print(lbound, ubound, period[lbound],period[ubound])
    filtered_periods=period[lbound+1:ubound]
    multiples=period[ubound:]
    return filtered_periods, multiples

def do_filter_multiples(multiples,npitch,frac=.05):
    '''examine outliers from do_filter_nominal to see if they're multiples'''
    #take modulus?
    maxmult=int(multiples[-1]/npitch) #multiples already sorted
    mults=[m for m in multiples if m % npitch <= frac*npitch]
    #now divide them all properly?
    #print mults
    filtered_mults=[]
    for n in range(2,maxmult+1):
        for m in mults:
            #print m/n, n
            if m/n > (1.-frac)*npitch and m/n < (1.+frac)*npitch:
                filtered_mults.append(m/n)
    return filtered_mults

def slit_widths_from_centers(centers,npitch,centers2=False,pix2um=1.995,stats=False,plot=False,filter_nominal=False):
    '''basically same as slit_widths_from_peaks(). should actually streamline that one to take keyword center'''
    period=[]
    #centers is a vector of the float positions of slat centers. use these to find period now
    for n,c in enumerate(centers[:-1]):
      period.append(centers[n+1]-c)

    if type(centers2) != bool: #also have slit centers, not just slat centers
        period2=[]
        for n,c in enumerate(centers2[:-1]):
            period2.append(centers2[n+1]-c)

    if filter_nominal: #get rid of anything 2x or more of nominal -- note that this is currently only a > filter
        filtered_periods=do_filter_nominal(period,npitch,pix2um)
        if type(centers2) != bool:
            filtered_periods2=do_filter_nominal(period2,npitch,pix2um)
    else:
        filtered_periods=[]
        filtered_periods2=[]
    #flatten arrays...
    if len(filtered_periods)!=0:
        filtered_periods=np.hstack(filtered_periods)
        filtered_periods2=np.hstack(filtered_periods2)



    if plot:
        #make the histogram
        fig,ax=plt.subplots()
        bins=np.arange(np.min(period),np.max(period),np.std(period)/5.)
        ax.hist(period,bins, aplha=0.6, facecolor='g',label='Period from Slat Centers')
        if type(centers2) != bool:
            ax.hist(period2,bins, aplha=0.6, facecolor='m',label='Period from Slit Centers')
        #ax.set_xlim([nperiod-5,nperiod+5])
        ax.set_yscale('log')
        ax.set_ylim([1,len(period)])
        ax.legend()
        #figfilename='win'+str(window_num)+'_group_periods'+str(mag)+'.png'
        #fig.savefig(figfilename)
        fig.show()

    if stats: #print out and pickle stats
        avg=np.mean(period) #here period is slit width ...what if I want to plot the actual period?
        med=np.median(period)
        stdv=np.std(period)
        print "-------------------STATISTICS FOR WINDOW---------------------"
        print '              Mean: ' + str(avg)
        print '            Median: ' + str(med)
        print 'Standard Deviation: ' + str(stdv)
        if type(centers2) != bool:
            #print 'Results saved in ' + statfile
            avg=np.mean(period2) #here period is slit width ...what if I want to plot the actual period?
            med=np.median(period2)
            stdv=np.std(period2)
            print "-------------------STATISTICS FOR WINDOW---------------------"
            print '              Mean: ' + str(avg)
            print '            Median: ' + str(med)
            print 'Standard Deviation: ' + str(stdv)

    if type(centers2) != bool:
        return period, filtered_periods, period2, filtered_periods2
    else:
        return period, filtered_periods

def print_all_stats(win):
    '''what it says on the tin. need to update tho...'''
    statfiles=glob.glob('win'+str(win) + '*_stats*'+'.p')
    for st in statfiles:
        stdict=pickle.load(open(st,'rb'))
        avg=stdict['mean']
        med=stdict['median']
        stdv=stdict['stddev']
        print "-------------------STATISTICS FOR WINDOW "+str(win)+"---------------------"
        print '              Mean: ' + str(avg)
        print '            Median: ' + str(med)
        print 'Standard Deviation: ' + str(stdv)
        print 'Results saved in ' + st
        #check in tags for other stats...
        #{'mean':avg,'median':med,'stddev':stdv,'wmean':np.ma.mean(widths),'wmed':np.ma.median(widths),'wdev':np.ma.std(widths),'amean':np.mean(filtered_c),'amed':np.median(filtered_c),'adev':np.std(filtered_c),'imean':np.mean(filtered_c2),'imed':np.median(filtered_c2),'idev':np.std(filtered_c2)}
        if 'wmean' in stdict.keys():
            print '             Total Widths: ' + str(len(stdict['wmean']))
            print '              Mean Widths: ' + str(stdict['wmean'])
            print '            Median Widths: ' + str(stdict['wmed'])
            print 'Standard Deviation Widths: ' + str(stdict['wdev'])
        if 'amean' in stdict.keys():
            print '             Total Period (slat centers): ' + str(len(filtered_c))
            print '              Mean Period (slat centers): ' + str(np.mean(filtered_c))
            print '            Median Period (slat centers): ' + str(np.median(filtered_c))
            print 'Standard Deviation Period (slat centers): ' + str(np.std(filtered_c))
        if 'imean' in stdict.keys():
            print '             Total Period (slit centers): ' + str(len(filtered_c2))
            print '              Mean Period (slit centers): ' + str(np.mean(filtered_c2))
            print '            Median Period (slit centers): ' + str(np.median(filtered_c2))
            print 'Standard Deviation Period (slit centers): ' + str(np.std(filtered_c2))


def reEdge(filen,sigma=3,gauss=3,plot=True):
    '''Take Canny edge image, Gaussian blur and re-edge detect to see what happens...'''
    earr=pickle.load(open(filen,'rb')).astype(float)
    im = ndi.gaussian_filter(earr, gauss)
    edges = feature.canny(im, sigma=sigma)
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        #ax = ax.ravel()

        ax.imshow(im, cmap=cm.binary)
        ax.imshow(np.ma.masked_where(edges == 0,edges),cmap=cm.autumn)
        ax.set_title('Input image overlaid with Canny edges')
        ax.set_axis_off()
        fig.show()
    return edges

def im_peek(filen,mag=5.0,length=0.2,contrast_stretch=False):
    '''Plot Canny edges over image for a given file'''
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib_scalebar.scalebar import SI_LENGTH
    if mag==5.0:
        pix2um=(1.2512/640.)*1000.
    else:
        pix2um=0.65
    #edges=pickle.load(open(filen,'rb'))
    #imf=filen[:-8]+'.tif'
    im=im2ndarray(filen)
    if contrast_stretch:
        p2 = np.percentile(im, 2)
        p98 = np.percentile(im, 98)
        im = exposure.rescale_intensity(im, in_range=(p2, p98))

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    ax.imshow(im, cmap=cm.gray)
    #ax.imshow(np.ma.masked_where(edges == 0,edges),cmap=cm.autumn)
    scalebar = ScaleBar(pix2um,'um', SI_LENGTH,length_fraction=length) # 1 pixel = 0.2 meter
    #print scale*pix2um
    ax.add_artist(scalebar)
    #ax.set_title('Input image overlaid with Canny edges')
    ax.set_axis_off()
    fig.show()

def edge_peek(filen,mag=5.0,length=0.2):
    '''Plot Canny edges over image for a given file'''
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib_scalebar.scalebar import SI_LENGTH
    if mag==5.0:
        pix2um=(1.2512/640.)*1000.
    else:
        pix2um=0.65
    edges=pickle.load(open(filen,'rb'))
    imf=filen[:-8]+'.tif'
    im=im2ndarray(imf)
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    ax.imshow(im, cmap=cm.gray,alpha=0.4)
    ax.imshow(np.ma.masked_where(edges == 0,edges),cmap=cm.winter)
    scalebar = ScaleBar(pix2um,'um', SI_LENGTH,length_fraction=length) # 1 pixel = 0.2 meter
    #print scale*pix2um
    ax.add_artist(scalebar)
    ax.set_title('Input image overlaid with Canny edges')
    ax.set_axis_off()
    fig.show()

def hough_peek(filen,edgef=False,mag=5.0,length=0.2):
    '''Plot Hough fits over Canny edges for a given file'''
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib_scalebar.scalebar import SI_LENGTH
    lines=pickle.load(open(filen,'rb'))
    if not edgef:
        edgef=filen[:-7]+'edges.p'
    if mag==5.0:
        pix2um=(1.2512/640.)*1000.
    else:
        pix2um=0.65

    edata=pickle.load(open(edgef,'rb'))
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    ax.imshow(np.ma.masked_where(edata == 0,edata),cmap=cm.gray)
    for line in lines:
        p0, p1 = line
        ax.plot((p0[0], p1[0]), (p0[1], p1[1]), color='r')
    scalebar = ScaleBar(pix2um,'um', SI_LENGTH,length_fraction=length) # 1 pixel = 0.2 meter
    #print scale*pix2um
    ax.add_artist(scalebar)
    ax.set_title('Canny edges overlaid with probabilistic Hough')
    ax.set_axis_off()
    fig.show()


def extend_line(line,shape=[640,480],plot=False):
    """extend the line as far as it can go within the boundaries of the frame. TOP LEFT corner is origin!!"""
    start=line[0]
    end=line[1]
    dxs,dys=shape[0]-start[0],shape[1]-start[1] #offsets from origin
    deltax=np.float(end[0])-np.float(start[0])
    deltay=np.float(end[1])-np.float(start[1])
    if deltax == 0.0:
        slope = 90.
    else:
        slope = deltay/deltax #*-1 ?
    #make a line with this slope, passing through start and end, that extends over the whole frame. Get endpoints...
    #if dxs >= shape[0]/2 and dys <=shape[1]/2: #look closer to bottom right corner...assume all slopes are +-45 degrees
    xvec=np.arange(0,shape[0],1)
    #x2=np.arange(int(xvec),shape[0],1)
    y2=slope*(xvec - np.float(start[0])) +np.float(start[1])
    #else:
    #    x2=np.arange(0,int(np.float(start[0])+np.float(dxs)/np.sqrt(2.)+3),1)
    #    y2=slope*(x2 - np.float(start[0])) +np.float(start[1])

    #now get endpoints for parts of the line that are within the frame - need to re-do limit on y!
    if y2[0] < y2[-1]:
        xi=np.where(y2 >= 0.)[0][0]
        try:
            xf=np.where(y2 >=shape[1]-1)[0][0]
        except IndexError:
             xf = np.where(y2==y2[-1])[0][0]
    else:
        xf=np.where(y2 >= 0.)[0][-1]
        try:
            xi=np.where(y2 >=shape[1]-1)[0][-1]
        except IndexError:
             xi = np.where(y2==y2[0])[0][0]

    extended_line=(int(xi),int(y2[xi])),(int(xf),int(y2[xf]))
    #slopeE=float(int(y2[xf])-int(y2[xi]))/float(int(xf)-int(xi))
    #print slope,slopeE
    if plot:
        s1=extended_line[0]
        e1=extended_line[1]
        fig,ax=plt.subplots()
        ax.plot((start[0],end[0]),(start[1],end[1]))
        ax.plot((s1[0],e1[0]),(s1[1],e1[1]),'r--')
        fig.show()

    return extended_line#,xvec,y2

def same_line(lines,nang,tol=3, plot=False):
    ''' rotate all line endpoints, see which ones lie on the same vertcal axis, group, return groups'''
    theta=np.radians(nang)
    c,s=np.cos(theta),np.sin(theta)
    rmat= np.array(((c,-s),(s,c)))#rotation matrix
    lx1,lx2=[],[]
    for line in lines:
        start=line[0]
        end=line[1]
        lx1.append(rmat.dot(start)[1]) #xcoord only
        lx2.append(rmat.dot(end)[1])
    #going to assume each entry in lx1 and lx2 more or less stack up...
    #or use average
    lx1e=np.transpose(zip(lx1,range(0,len(lx1))))
    lx2e=np.transpose(zip(lx2,range(0,len(lx2))))
    #sort by lx1 or lx2
    lx1s=np.transpose(sorted(lx1e,key=lambda x:x[0])) #sorts by x coord, index carried along
    lx2s=np.transpose(sorted(lx2e,key=lambda x:x[0]))

    #now group within tol -- code copied from same_vert_line
    unique_matches=[]
    for lx1,lx2 in zip(lx1s,lx2s):
        lx1x,lidx1=lx1
        lx2x,lidx2=lx2
        if abs(lx1x-lx2x) > tol: #need a bigger tolerance if the angle is so big!
            ltol=abs(lx1x-lx2x)
        else:
            ltol=tol
        try:
            xm,midx=unique_matches[-1][0] #latest unique match
            farthest_x=np.max([np.abs(xm-lx1x),np.abs(xm-lx2x)])
            #print line, to_match
            #print farthest_x,closest_y
            if farthest_x <= tol: #must be on the same edge
                unique_matches[-1].append(lx1)
            else:
                unique_matches.append([lx1])

        except IndexError: #empty list
            unique_matches.append([lx1]) #this is just the indexes

    #now re-match line to group
    grouped_lines=[]
    for um in unique_matches:
        if um[0][0].is_integer():
            idx=0
        else:
            idx=1
        grouped_lines.append([lines[int(u[idx])] for u in um])
    return grouped_lines



# def same_line(lines,tol=3, plot=False,shape=[640,480]):
#     """Test if two or more Hough line segments lie on the same line. Returns lists of line segments that fall on the same line, within a certain tolerance at the endpoints"""
#     extended_lines,matches,singles=[],[],[]
#     duplicate=False
#     for line in lines:
#         start=line[0]
#         end=line[1]
#         extended_lines.append(extend_line(line,shape=shape)) #extend the line as far as it can go within the boundaries of the frame

#     #now test extended lines against each other
#     allxl=[(l[0][0],l[1][0]) for l in lines]
#     allyl=[(l[0][1],l[1][1]) for l in lines]
#     allx=[(exl[0][0],exl[1][0]) for exl in extended_lines]
#     #ally=[(exl[0][1],exl[1][1]) for exl in extended_lines]
#     for l,exl in zip(lines,extended_lines):
#         #pad the endpoints in x
#         ex1=exl[0]
#         ex2=exl[1]
#         if ex1[0]-tol/2 > 0:
#             padded_ex1=range(ex1[0]-tol/2,ex1[0]+tol/2,1)
#         else:
#             padded_ex1=range(ex1[0],ex1[0]+tol,1)
#         if ex2[0]+tol/2 > 0:
#             padded_ex2=range(ex2[0]-tol/2,ex2[0]+tol/2,1)
#         else:
#             padded_ex2=range(ex2[0],ex2[0]+tol,1)
#         #test if any other lines fall on the extended line, ie if any in allx are in padded_ex
#         #print padded_ex1,padded_ex2
#         for m in matches:
#             if l in m:
#                 duplicate=True
#                 #print l,m, duplicate
#             else:
#                 duplicate = False
#                 #print l,m, duplicate
#         if duplicate == False:
#             matches.append([((xl[0],yl[0]),(xl[1],yl[1])) for x,xl,yl in zip(allx,allxl,allyl) if x[0] in padded_ex1 and x[1] in padded_ex2]) #now need to append the actual line segments associated with the extended lines
#     for m in matches:
#         if len(m) == 1: #find singles
#             singles.append(m)
#             matches.remove(m)
#     #print len(singles)

#     unique_matches=[]
#     for m in matches: #now matches doesn't contain any singles
#         if m not in unique_matches:
#             unique_matches.append(m)
#             for s in singles:
#                 #print s,m
#                 if s[0] in m:
#                     singles.remove(s)
#                     #print s,m
#     #print len(singles)
#     for s in singles:
#         unique_matches.append(s)

#         #if so, ...group and return

#     if plot:
#         fig,ax=plt.subplots()
#         scale=len(unique_matches)
#         for i,mm in enumerate(unique_matches):
#             for m in mm:
#                 s1=m[0]
#                 e1=m[1]
#                 ax.plot((s1[0],e1[0]),(s1[1],e1[1]),color=cm.jet(float(i)/float(scale)))
#                 #ax.plot((s1[0],e1[0]),(s1[1],e1[1]),'r--')
#             #print i/10.,cm.jet(i/10.)
#         fig.show()

#     return unique_matches

def same_vert_line(lines,tol=2,ytol=50,plot_hist=False):
    '''specifically for vertical or near-vertical lines'''
    unique_matches=[]
    all_angles=[]
    for line in lines:
        angle=np.abs(get_angle(line))
        x0,y0=line[0]
        x1,y1=line[1]
        if x1-x0 > tol: #need a bigger tolerance if the angle is so big!
            ltol=x1-x0
        else:
            ltol=tol
        try:
            to_match=unique_matches[-1] #latest unique match
            last_ang=all_angles[-1]
            xm0,ym0=to_match[0][0]
            xm1,ym1=to_match[0][1]
            farthest_x=np.max([np.abs(xm1-x1),np.abs(xm0-x0)])
            closest_y=np.min([np.abs(ym1-y1),np.abs(ym0-y0)])
            #print line, to_match
            #print farthest_x,closest_y
            if farthest_x <= ltol and closest_y > ytol: #must be on the same edge
                unique_matches[-1].append(line)
                all_angles[-1].append(angle)
            else:
                unique_matches.append([line])
                all_angles.append([angle])

        except IndexError: #empty list
            unique_matches.append([line])
            all_angles.append([angle])

    if plot_hist:
        fig,ax=plt.subplots()
        for i,a in enumerate(all_angles):
            for eacha in a:
                ax.scatter(i,eacha)
            ax.scatter(i,np.mean(a),marker='+',c='r')
        ax.set_xlim([0,len(all_angles)])
        fig.show()
    return all_angles,unique_matches


def rising_or_falling(rotim,rotarr,immean,plot=False,imfile=False,shape=False,rerotate=False, test=False):
    """Determine if an edge is rising (dark -> light) or falling (light -> dark). Return indexed mask"""
    imsize=np.shape(rotim)
    if not shape:
        mask=np.zeros([792,792])#np.zeros([640,480])
        oshape=np.zeros([640,480])
    else:
        mask=np.zeros(shape)
        oshape=np.zeros([640,480])
    meani=np.mean(rotim, axis=0) #average along y
    meanarr=np.sum(rotarr, axis=0) # total along y (remember it's boolean)
    #for rowi,rowarr in zip(meani,meanarr):
    tv=np.where(meanarr > np.max(meanarr[0:50]))
    tvshape=np.shape(tv)[1]
    #get rising
    rising=[tv[0][i] for i in range(0,tvshape-1) if np.mean(meani[tv[0][i]:tv[0][i+1]]) > 1.25*immean and tv[0][i+1] != tv[0][i]+1 ]

    mask[:,rising]= 1
    #get falling
    falling=[tv[0][i] for i in range(0,tvshape-1) if np.mean(meani[tv[0][i]:tv[0][i+1]]) < .75*immean and tv[0][i+1] != tv[0][i]+1 ]
    mask[:,falling]= -1
    #print len(rising),len(falling)

    if test:
        print 'immean: ', immean
        print 'threshold: ',np.max(meanarr[0:50])
        print len(rising), len(falling)
        print 'first few rising: ', rising[:10]
        print 'means: ', [np.mean(meani[tv[0][i]:tv[0][i+1]]) for i in range(0,tvshape-1) if np.mean(meani[tv[0][i]:tv[0][i+1]]) > 1.25*immean and tv[0][i+1] != tv[0][i]+1 ][:10]
        print 'first few falling: ', falling[:10]
        print 'means: ',[np.mean(meani[tv[0][i]:tv[0][i+1]]) for i in range(0,tvshape-1) if np.mean(meani[tv[0][i]:tv[0][i+1]]) < .75*immean and tv[0][i+1] != tv[0][i]+1 ][:10]

    if len(rising)-len(falling) > 3.*np.max([len(rising),len(falling)]):
         #print 'check mask!'
         #plot=True
         mask=[]

    if plot:
        #import make_mask as mm
        #mm.mask_peek('win11_05_05_5.0X.tif',mask)
        if test:
            fig,ax=plt.subplots()
            ax.plot(range(0,np.shape(meanarr)[0]),meanarr)
            fig.show()

        fig,ax=plt.subplots()
        if not imfile:
            ax.imshow(rotate(np.transpose(im2ndarray('win11_05_05_5.0X.tif')),-44.79357,reshape=True),alpha=0.6,cmap=cm.gray)
        else:
            ax.imshow(rotate(np.transpose(im2ndarray(imfile)),45.03455,reshape=True),alpha=0.6,cmap=cm.gray)
        ax.imshow(mask,cmap=cm.gray,alpha=0.7)
        fig.show()

    if rerotate: #rotate it back by given angle and trim to original size
        mask=imrotate(mask,rerotate)
    #if len(rising) > 1.25* len(falling) or len(falling) > 1.25* len(rising):
        #return []
    return mask

def thin_to_binary(subarr,left):
    dim=np.shape(subarr)
    thinnedarr=np.zeros(np.shape(subarr))
    for i in range(0,dim[0]):
        row =  subarr[i,:]
        if np.sum(row) == 0:
            continue
        else:
            maxloc= np.where(row == np.max(row))
            thinnedarr[i,maxloc] =maxloc+left

            #alledges= np.where(row != 0)
            #thinnedarr[i,maxloc] =np.mean(row[alledges])

    return thinnedarr

def group_edges_by_idx(rotarr,mask,nomp,mod=3,plot=False,tolerance=.6,from_line = False):
    '''Group edges based on mask index'''
    rising=list(np.where(mask[0,:] == 1.0)[0])
    falling=list(np.where(mask[0,:] == -1.0)[0])
    rmean,fmean=[],[]
    for r in rising:
        subarr=rotarr[:,r-mod:r+mod]
        tarr = thin_to_binary(subarr,r-mod)
        #print r-mod
        if np.isfinite(np.mean(tarr)):
            total_edges=np.sum([1 for i in np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))])
            rmean.append(np.mean(tarr[np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))]))

    for f in falling:
        subarr=rotarr[:,f-mod:f+mod]
        tarr = thin_to_binary(subarr,f-mod)
        if np.isfinite(np.mean(tarr)):
            total_edges=np.sum([1 for i in np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))])
            fmean.append(np.mean(tarr[np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))]))

    rising.extend(falling)

    if plot:
        fig,ax=plt.subplots()
        ax.scatter(rising,rmean+fmean)
        fig.show()

    periodr = [rmean[j+1]-rmean[j] for j in range(0,len(rmean)-1)]
    periodr = [p for p in periodr if p > tolerance*nomp and p < (2.-tolerance)*nomp]
    #print rmean
    #print periodr
    periodf= [fmean[j+1]-fmean[j] for j in range(0,len(fmean)-1)]
    periodf = [p for p in periodf if p > tolerance*nomp and p < (2.-tolerance)*nomp]
    #print np.mean(periodr)*1.955,np.mean(periodf)*1.955#, np.mean(periodr.extend(periodf))*1.955

    if from_line: #plot and fit line for rising and falling
        fig,ax=plt.subplots(2,sharex=True)
        x=np.arange(0,.5*len(rmean),.5)
        ax[0].scatter(x,rmean)
        liner=np.polyfit(x,rmean,1)
        print(liner)
        #x=np.arange(0,2*len(rmean),2)
        yfitr = lambda x: (liner[0]*x+liner[1])
        print np.shape(x),np.shape(yfitr)
        ax[0].plot(x,yfitr(x),'k--')
        x=np.arange(0,.5*len(fmean),.5)#falling#np.arange(0,len(fmean),1)
        linef=np.polyfit(x,fmean,1)
        yfitf = lambda x: (linef[0]*x+linef[1])
        print np.shape(x),linef,np.shape(yfitf(x))
        ax[1].scatter(x,fmean)
        ax[1].plot(x,yfitf(x),'k--')
        #ax.scatter(np.arange(0,len(rising),1),rmean)
        print 'line of best fit rising: '+str(liner[0])+'*x + ' +str(liner[1])
        print 'line of best fit falling: '+str(linef[0])+'*(x) + ' +str(linef[1])
        fig.show()


    return periodr,periodf


def group_edges_by_mosaic_idx(rotarr,mask,nomp,mod=3,plot=True,tolerance=.6,from_line = True):
    '''Group edges based on mosaic mask index'''
    rising=list(np.where(mask[0,:] == 1.0)[0])
    falling=list(np.where(mask[0,:] == -1.0)[0])
    rmean,fmean=[],[]
    for r in rising:
        subarr=rotarr[:,r-mod:r+mod]
        tarr = thin_to_binary(subarr,r-mod)
        #print r-mod
        if np.isfinite(np.mean(tarr)):
            total_edges=np.sum([1 for i in np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))])
            rmean.append(np.mean(tarr[np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))]))

    for f in falling:
        subarr=rotarr[:,f-mod:f+mod]
        tarr = thin_to_binary(subarr,f-mod)
        if np.isfinite(np.mean(tarr)):
            total_edges=np.sum([1 for i in np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))])
            fmean.append(np.mean(tarr[np.where(np.logical_and(tarr != 0, np.isfinite(tarr) == True))]))

    rising.extend(falling)

    if plot:
        fig,ax=plt.subplots()
        ax.scatter(rising,rmean+fmean)
        fig.show()

    periodr = [rmean[j+1]-rmean[j] for j in range(0,len(rmean)-1)]
    periodr = [p for p in periodr if p > tolerance*nomp and p < (2.-tolerance)*nomp]
    #print rmean
    #print periodr
    periodf= [fmean[j+1]-fmean[j] for j in range(0,len(fmean)-1)]
    periodf = [p for p in periodf if p > tolerance*nomp and p < (2.-tolerance)*nomp]
    #print np.mean(periodr)*1.955,np.mean(periodf)*1.955#, np.mean(periodr.extend(periodf))*1.955

    return periodr,periodf

def test_mask(ef, im, nang):
    edges=pickle.load(open(ef,'rb'))
    ima=im2ndarray(im)
    ima=remove_edges(im,ima)
    rotim=rotate(np.transpose(ima),nang, reshape=True)
    rotarr=rotate(np.transpose(edges),nang,reshape=True)
    #print im,ef
    mask=rising_or_falling(rotim[300:950,:],rotarr[300:950,:],np.mean(rotim[300:950,100:-100]), shape=np.shape(rotarr),imfile=im,plot=True,test=True)
    print np.mean(rotim[300:950,:])

def mean_hough_angle(f_all,side=1.0,save=True):
    lines=pickle.load(open(f_all,'rb'))
    angles=[]
    for l in lines:
        angles.append(side*get_angle(l))
    meana=np.mean(angles)
    if save:
        pickle.dump(meana,open(f_all[:-2]+'_mean.p','wb'))
    return meana


def get_period_by_grouping(window_num,mag=5.0,ftags='X*',plot=False,side=1.0,pix2um=1.955,stats=True,EM=True,tolerance=0.6,mosaic=True,offset=False, fitangle=True):
    '''Put everything together for a single window'''
    #first get the files
    if EM:
        edgef=EM_list(str(window_num),str(mag),ending='_edges.p')
    else:
        edgef=glob.glob('win'+str(window_num)+'_*' + str(mag)+ftags+'_edges.p')
    if EM:
        imf=EM_list(str(window_num),str(mag),ending='.tif')
    else:
        imf=glob.glob('win'+str(window_num)+'_*' + str(mag)+'X.tif')
    if mosaic:
        edgef=['win11_mosaic_5.0_edges.p']
        imf=['window11mosaic_5.0.tiff']
    print imf[:10], edgef[:10]
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
        #print nang
        rotim=rotate(np.transpose(ima),side*nang, reshape=True)
        rotarr=rotate(np.transpose(edges),side*nang,reshape=True)
        #print im,ef
        mask=rising_or_falling(rotim[300:450,:],rotarr[300:450,:],np.mean(rotim[300:450,100:-100]), shape=np.shape(rotarr),imfile=im,plot=False)
        #group
        if np.shape(mask)[0] > 0:
            periodr,periodf=group_edges_by_idx(rotarr,mask,nperiod/pix2um,tolerance=tolerance,mod=3,plot=False)
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
        figfilename='win'+str(window_num)+'_group_periods'+str(mag)+ftags+'.png'
        fig.savefig(figfilename)
        fig.show()

    if stats: #print out and pickle stats
        avg=np.mean(periods)
        med=np.median(periods)
        stdv=np.std(periods)
        statfile='win'+str(window_num)+'_width_stats'+str(mag)+ftags+'.p'
        datafile='win'+str(window_num)+'_width_data'+str(mag)+ftags+'.p'
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

def slit_or_slat(j,row,imarr,immean,pix2um): #for j,row in enumerate(edges): #MPI this! super slow right now....
    tv=np.where(row == True)
    tvshape=np.shape(tv)[1]
    try:
        #width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1) if np.mean(imarr[tv[0][i+1]-tv[0][i],j]) < immean] #slats = space between slits
        width=[float(tv[1][i+1]-tv[1][i]) for i in range(0,tvshape-1) if np.mean(imarr[tv[1][i+1]-tv[1][i],j]) < immean] #slats = space between slits
    except ValueError:
        #width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1)]
        width=[float(tv[1][i+1]-tv[1][i]) for i in range(0,tvshape-1)]
    if width != []:
        width=filter(lambda a:a!=1.0,width) #filter out all the 1px separations
        #now convert to period. Let's start with some assumptions:
        #1) the width of the slat is about equal to the width of the slit. ideal case: width[15,13,15,13] etc
        #   Then period is simply: [width[i]+width[i+1] for i in range(0,len(width))]
        #   Of course this will probably not be the case. Worst case, there is some noise due to a bump or something:
        #   Or an edge pixel is missing so that the width is already the period
        #   width=[2,2,15,3,8,13,2,30,13]
        #   So we need some conditions. This could cost quite some computing time....we can do it iteratively in steps.
        #   we'll need to use some basic statistics to help...let's trust for now that the rows are long enough that we can do this.
        #   might not be the case with coarser grids!!!
        sigma=3*np.std(width) # 3 sigma
        mean=np.mean(width)
        period=[]
        for i in range(0,len(width)-1):
            if width[i] > mean+sigma:
                period.append(width[i])
            elif width[i] < mean-sigma or width[i+1] < mean-sigma: #discard it
                continue
            else: #add it to the previous one UNLESS the next one is a double! then just continue
                if width[i+1] < mean+sigma and  width[i+1] > mean-sigma:
                    period.append(width[i]+width[i+1])
        period=np.array(period)*pix2um # need to account for the fact that it's rotated! width is not true width!
    return period#width

def slit_or_slat_full(edges,imarr,immean,pix2um):
    '''For full image, not just row in mosaic'''
    awidths,means,sigmas=[],[],[]
    for j, row in enumerate(edges):
        tv=np.where(row == True)
        tvshape=np.shape(tv)[1]
        #try:
        #    slat_width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1) if np.mean(imarr[tv[0][i+1]-tv[0][i],j]) > immean] #slats =glold part
        #    slit_width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1) if np.mean(imarr[tv[0][i+1]-tv[0][i],j]) < immean] #slats =glold part
        #except ValueError:
            #width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1)]
        #    width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1)]
        #if slit_width != [] and slat_width !=[]:
        #    slit_width=filter(lambda a:a!=1.0,slit_width) #filter out all the 1px separations
        #    slat_width=filter(lambda a:a!=1.0,slat_width)
            #now convert to period. Let's start with some assumptions:
            #1) the width of the slat is about equal to the width of the slit. ideal case: width[15,13,15,13] etc
            #   Then period is simply: [width[i]+width[i+1] for i in range(0,len(width))]
            #   Of course this will probably not be the case. Worst case, there is some noise due to a bump or something:
            #   Or an edge pixel is missing so that the width is already the period
            #   width=[2,2,15,3,8,13,2,30,13]
            #   So we need some conditions. This could cost quite some computing time....we can do it iteratively in steps.
            #   we'll need to use some basic statistics to help...let's trust for now that the rows are long enough that we can do this.
            #   might not be the case with coarser grids!!!
            #2) Every detection of a slit/slat is followed by a detection of a slat/slit, EXCEPT when the width given is the entire period
         #   sigma=np.std(slat_width) # 1 sigma
         #   mean=np.mean(slat_width)
         #   period=[]
         #   for i in range(0,len(slat_width)-1): #let's hope these arrays are equally long....
         #       if slat_width[i] > 2*mean-4*sigma: #it's an entire period
         #           period.append(slat_width[i])
         #       elif slat_width[i] < mean-sigma or slat_width[i+1] < mean-sigma: #discard it
         #           continue
         #       else: #add it to the previous one UNLESS the next one is a double! then just continue
         #           if slat_width[i] < mean+sigma and  slat_width[i] > mean-sigma:
         #               #now we have to look at the corresponding slits:
         #               if slit_width[i] < mean+sigma and  slit_width[i] > mean-sigma:
         #                   period.append(slit_width[i]+slat_width[i+1])

        width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1) if float(tv[0][i+1]-tv[0][i]) !=1.0]
        sigma=np.std(width) # 1 sigma
        mean=np.mean(width)
        means.append(mean)
        sigmas.append(sigma)
        #print sigma,mean
        period=[]
        if sigma > mean: sigma=0.25*mean
        for i in range(0,len(width)-1): #let's hope these arrays are equally long....
            if width[i] > 2*mean-sigma and width[i] < 2*mean+sigma: #it's an entire period
                period.append(width[i])
            elif width[i+1] < mean+sigma and  width[i+1] > mean-sigma:
                    period.append(width[i]+width[i+1])

        period=np.array(period)*pix2um
        awidths.append(period)
    #flatten
    fperiod=[item for sublist in awidths for item in sublist ]
    #print np.mean(fperiod)/pix2um,np.std(fperiod)/pix2um
    return fperiod#,means,sigmas

def get_slit_width(edges,mag=5.0,im=False,window_num=False, title='',xran=[0,45],figname=False, stats=True):
    '''Currently this will only work for the finest grids, that have well-defined edges without stuff in the middle'''
    import time
    import multiprocessing as mpi
    start=time.time()
    #let's go row-by-row and get the distance between True values, and compile into a histogram:
    widths=[]
    if mag == 5.0:
        pix2um=1.954983#(1.2512/640.0)*1000.0 #from RearGrid/windowcorners.txt
        print pix2um
    else: #15X
        pix2um=.6505
        print 'pix2um',pix2um
        #if im !=False: #it's a numpy array of an image that we will use to determine if slit or slat. Let's make sure these are the same shape
    if im.any():
        imarr=im
        if np.shape(imarr) !=np.shape(edges):
            trim=raw_input("Dimensions don't match! which side or list of sides to trim? [top,bottom,left,right,all]")
            if trim == 'all':
                imarr=imarr[50:-50,50:-50]
            elif trim == 'top':
                imarr=imarr[:,:-50]
            elif trim == 'bottom':
                imarr=imarr[:,50:]
            elif trim == 'left':
                imarr=imarr[50:,:]
            elif trim == 'right':
                imarr=imarr[:-50,:]

    pool=mpi.Pool()
    print 'preparation took %.2f seconds' % (time.time()-start)
    ltime=time.time()
    #immean=np.mean(imarr) #save time and not do this every call to slit_or_slat()

    if type(edges) == list: #it's a list of edge pickles
        for earr in edges:
            edgearr=pickle.load(open(earr,'rb'))
            imf=earr[:-8]+'.tif'
            imarr=im2ndarray(imf)
            imarr=remove_edges(earr,imarr)
            p2, p98 = np.percentile(imarr, (2, 98))
            imarr = exposure.rescale_intensity(imarr, in_range=(p2, p98))

            immean=np.mean(imarr)
            #result=pool.apply_async(slit_or_slat_full,args=(edgearr,imarr,immean,pix2um))
            result=slit_or_slat_full(edgearr,imarr,immean,pix2um)
            try:
                #print np.mean(result.get())
                widths.append(result)
            except IndexError:
                continue
        m,s=divmod((time.time() - ltime),60.)
        print 'Loop took %02d:%02d ' % (m,s)
    else:
        immean=np.mean(imarr)
        for j,row in enumerate(edges): #MPI this! super slow right now....
            result=pool.apply_async(slit_or_slat,args=(j,row,imarr,immean,pix2um))
            try:
                widths.append(result.get())
            except IndexError:
                continue
            #widths.append(slit_or_slat(j,row,imarr))

        m,s=divmod((time.time() - ltime),60.)
        print 'Loop took %02d:%02d ' % (m,s)

    widths_vec=[item for sublist in widths for item in sublist ]
    m,s=divmod((time.time() - start),60.)
    print 'Processing took %02d:%02d ' % (m,s)

    if stats: #print out and pickle stats
        avg=np.mean(widths_vec)
        med=np.median(widths_vec)
        stdv=np.std(widths_vec)
        statfile='win'+str(window_num)+'_width_stats'+str(mag)+'.p'
        datafile='win'+str(window_num)+'_width_data'+str(mag)+'.p'
        print "-------------------STATISTICS FOR WINDOW "+str(window_num)+"---------------------"
        print '              Mean: ' + str(avg)
        print '            Median: ' + str(med)
        print 'Standard Deviation: ' + str(stdv)
        print 'Results saved in ' + statfile
        data=widths_vec
        stdict={'mean':avg,'median':med,'stddev':stdv}
        pickle.dump(stdict,open(statfile,'wb'))
        pickle.dump(data,open(datafile,'wb'))


    fig,ax=plt.subplots()
    bins=ax.hist(np.array(widths_vec),np.arange(np.min(widths_vec),np.max(widths_vec),pix2um))
    ax.set_xlabel('Separation of Edge Pixels ($\mu$m)')
    ax.set_ylabel('Frequency')
    ax.set_title(title)
    ax.set_yscale('log')
    #if window_num: #set the x-range by pitch
    #    nom_width=
    ax.set_xlim(xran)
    #fig.show()
    if figname:
        fig.savefig(figname)
    plt.close(fig)
    return widths,bins


def get_theta_range(nang,spread=5.,n=201):
    #define range of theta around theta_nominal
#     if type(edges) == int:
#         winnum=edges
#     elif 'win' not in edges:
#         try:
#             winnum=int(edges)
#         except ValueError:
#             edges=raw_input('What is the window number?')
#             winnum=int(edges)
#     else:
#         winnum,index=get_index(edges)
#         winnum=int(winnum)
#     if side==1.0:
#         nang=[w['nominal angle'] for w in windows if w['number'] == winnum]
#         nang=nang[0]
#     else:
#         nang=[w['nominal angle'] for w in windowsr if w['number'] == winnum]
#         nang=-1*nang[0]
    theta0= nang*(np.pi/180.)#in radians
    #tendeg2rad=np.pi/18.
    spreaddeg2rad=spread*(np.pi/180.)
    thetaran = np.linspace(theta0-spreaddeg2rad, theta0+spreaddeg2rad, num=n)#in radians
    return thetaran

def prob_hough(edges, nang,threshold=10, line_length=50, line_gap=2,retlines=False,plot=False,spread=5.,n=201, tag=False,overwrite=False):
    '''Perform probabilistic Hough fit to given set of edges'''
    import glob
    if type(edges) == str: #it's a filename
        inp=edges
        edata=pickle.load(open(edges,'rb'))
    else:
        edata=edges
        edges=raw_input('What is the window number?')
        inp=raw_input('Output file name?')
    if not overwrite:
        if tag:
            defaultn=inp[:-8]+'_hough_'+tag+'.p'
        else:
            defaultn=inp[:-8]+'_hough.p'
        names=glob.glob(defaultn)
        if names !=[]:
            return

    thetaran=get_theta_range(nang,spread=spread,n=n)
    start=time.time()
    lines = probabilistic_hough_line(edata, threshold=threshold, line_length=line_length,
                                 line_gap=line_gap, theta=thetaran)
    print 'Probabilistic Hough fitting took %.2f seconds' % (time.time() - start)
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        #ax = ax.ravel()
        ax.imshow(np.ma.masked_where(edata == 0,edata),cmap=cm.gray)
        for line in lines:
            p0, p1 = line
            ax.plot((p0[0], p1[0]), (p0[1], p1[1]), color='r')
        ax.set_title('Canny edges overlaid with probabilistic Hough')
        ax.set_axis_off()
        fig.show()

    try:
        if inp:
            newfilen=inp[:-8]+'_hough.p'
            if tag:
                newfilen=inp[:-8]+'_hough_'+tag+'.p'
            pickle.dump(lines,open(newfilen,'wb'))
            return newfilen
    except ValueError:
        pass

    if retlines:
        return lines

def straight_hough(edges,plot=False,side=1.0,spread=5.,n=201):
    '''Perform straight line Hough fit to given set of edges'''
    if type(edges) == str: #it's a filename
        inp=edges
        edata=pickle.load(open(edges,'rb'))
    else:
        edata=edges
        edges=raw_input('What is the window number?')
        inp=raw_input('Output file name?')

    thetaran=get_theta_range(edges,spread=spread,n=n,side=side)
    start=time.time()
    lines, angles,distances = hough_line(edata, theta=thetaran)
    print 'Straight line Hough fitting took %.2f seconds' % (time.time() - start)
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        #ax = ax.ravel()
        ax.imshow(np.ma.masked_where(edata == 0,edata),cmap=cm.gray)
        for line in lines:
            p0, p1 = line
            ax.plot((p0[0], p1[0]), (p0[1], p1[1]), color='r')
        ax.set_title('Canny edges overlaid with probabilistic Hough')
        ax.set_axis_off()
        fig.show()

    try:
        if inp:
            newfilen=inp[:-8]+'_hough_angles.p'
            pickle.dump(angles,open(newfilen,'wb'))
            return newfilen
    except ValueError:
        pass

def lines2data(lines):
    '''convert lines to binary numpy array'''
    data=np.zeros([640,480])
    for line in lines:
        sc,ec=line
        x1=sc[0]
        x2=ec[0]
        y1=sc[1]
        y2=ec[1]
        deltax=float(x2-x1)
        deltay=float(y2-y1)
        m=deltay/deltax #slope of line
        #print y1,y2, m, sc,ec
        for x in np.arange(min(x1,x2),max(x1,x2)):
            y=m*(x-x1)+y1
            data[x,int(y)]=1.0

    #fig,ax=plt.subplots()
    #ax.imshow(data.T,cmap=cm.binary)
    #fig.show()
    return data

def whole_window_hough_hist(window_num,path, plot=False,ymax=10000):
    '''Make histogram based off of rotated images of line plots for each separate image in complete window'''
    #this is not going to work without knowledge of the coordinates... have no way of knowing where 'origin' is supposed to be, so the pattern won't overlay correctly. I guess this means I have to mosaic first... i hate mosaicing
    #restore all the lines
    os.chdir(path)
    files=glob.glob('win'+str(window_num)+'*5.0X_bright_hough.p')
    #print files
    llist=[pickle.load(open(f,'rb')) for f in files]
    ndlist=[lines2data(lines) for lines in llist]
    #sum all arrays into one big one
    bigarray=np.sum(ndlist[0:5],axis=0)
    #rotate. This does do bilinear interpolation... not sure what kind of errors this introduces.
    for w in windows:
        if window_num == w['number']:
            angle=w['nominal angle']
    rotim=imrotate(bigarray,angle)

    #compress along the y-axis...
    bigvector=np.sum(rotim,axis=0)

    fig,ax=plt.subplots()
    #ax.plot(np.arange(640),bigvector)
    #n,bins=ax.hist(bigvector,np.linspace(0,640,640), facecolor='b')
    for nd in ndlist:
        ax.imshow(nd.T, alpha=.5,cmap=cm.gray)

    #if plot:
    #    ax.set_xlim([0,640])
    #    ax.set_xlabel('Distance across grid (px)')
    #    ax.set_ylabel('Number of edge pixels in bin')
        #ax.set_yscale('log')
        #ax.set_ylim([1,ymax])
    #fig.show()

    #return n,bins
    return ndlist,bigarray

def sg_filter(edges, order=1,window_size=641):
    '''Apply Savitzky-Golay filter to the edges to smooth them and make line detection better'''
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def dialate_and_erode(edges):
    import cv2
    blur=((3,3),1)
    erode_=(5,5)
    dilate_=(3, 3)
    earr=pickle.load(open(edges,'rb')).astype(float)*255
    im=earr*255.
    imblurred=cv2.GaussianBlur(im, blur[0], blur[1])
    imdilated=cv2.dilate(imblurred, np.ones(dilate_))*255
    imeroded=cv2.erode(imdilated, np.ones(erode_))
    ims=[imblurred,imdilated,imeroded]
    fig,axes = plt.subplots(3,sharex=True,sharey=True)
    for ax,im in zip(axes,ims):
        ax.imshow(im,cmap=cm.gray)
    fig.show()
    #imdilated=cv2.dilate(imblurred, np.ones(dilate_))*255
    return imblurred,imeroded,imdilated
    #cv2.imwrite('testim_blured.png',cv2.dilate(cv2.erode(imblurred, np.ones(erode_)), np.ones(dilate_))*255)

def cat_hough(window_num, imtags='X_', tags=['ll150','ll200','ll250','ll300'],weights=False,outtag='all',EM=True,mag=5.0):
    """Cat together the lines """
    #first get indices for the window
    if not EM:
        ofiles=glob.glob('win'+str(window_num)+'_*'+str(mag)+imtags+'edges.p')
    else:
        ofiles=EM_list(str(window_num),str(mag),ending='edges.p')
    idxs=[get_index(o)[1] for o in ofiles]
    #print ofiles[0]
    #all_lines=[]
    for idx in idxs:
        basen='win'+str(window_num)+'_'+str(idx[0])+'_'+str(idx[1])+'_'+str(mag)+imtags+'hough'
        #print basen
        all_lines=[]
        #lines=pickle.load(open(basen+'.p','rb'))
        #print len(lines)
        #all_lines=lines
        #print len(all_lines)
        for i,tag in enumerate(tags):
            try:
                lines=pickle.load(open(basen+'_'+tag+'.p','rb'))
            except IOError:
                continue
            #print len(lines)
            if not weights:
                all_lines.extend(lines)
            else:
                for j in range(0,weights[i]):
                    all_lines.extend(lines)
            #print len(all_lines)
        pickle.dump(all_lines,open(basen+'_'+outtag+'.p','wb'))

def hough_hist(lines, windownum, nang, mag=5.0,log=True,ret=False, title=False,xran=False,stats=True,figname=False,spread=5.,n=201,gaussfit=False, sameline=True,tol=3,mask45=False): #cuz of mpi
    '''Make a histogram of the line orientations returned by probabilistic Hough'''
    theta=[]
    llist=[]
    if type(lines) == list: #it's a list of filenames
        for f in lines:
            if f != None:
                llist.append(pickle.load(open(f,'rb')))
            #print f,np.shape(pickle.load(open(f,'rb')))
        #print np.shape(llist)
    for i,l in enumerate(llist):
        if i == 0:
            all_lines=l
            lines=all_lines
        if i < len(llist)-1:
            if sameline:
                grouped_lines=same_line(l,nang,tol=tol)
                for gl in grouped_lines:
                    theta.append(np.mean([get_angle(g) for g in gl]))
            else:
                all_lines=all_lines + llist[i+1]
                lines=all_lines
    if theta == []:
        for l in lines:
            #if get_length(l) >= 50.:
            #if get_angle(l) != 45.0 and get_angle(l) != -45.0: #suppress the 45.0 value as a test...
            try:
                theta.append(get_angle(l))
            except TypeError:
                continue

    print len(theta), np.min(theta),np.max(theta)
    if mask45:
        theta=np.ma.masked_equal(theta,45.0)
    thetaran=get_theta_range(nang,n=n,spread=spread)
    #make a histogram
    fig,ax=plt.subplots()

    thetaax=np.arange(np.min(theta),np.max(theta),.05)
    if theta !=[]:
        try:
            yhist,xhist=np.histogram(theta.compressed(),thetaax) #if masked array
        except AttributeError:
            yhist,xhist=np.histogram(theta,thetaax)
        foo=ax.hist(theta,thetaax)

    if gaussfit:
        #hindex=np.where(yhist > 0)[0]
        #xh = xhist[hindex] #do I really want to fit this? this is indexes not theta...
        #yh = yhist[xh]
        #print np.max(yhist), np.mean(yhist), np.std(yhist)
        if xran: #fit also within nthe given x-range
            hindexl=np.where(xhist >= xran[0])[0]
            hindexh=np.where(xhist <= xran[1])[0]
            xhist=xhist[hindexl[0]:hindexh[-1]]
            yhist=yhist[hindexl[0]:hindexh[-1]]
        if np.shape(xhist) !=np.shape(yhist):
            xhist=xhist[:-1]
        def gaussian(x, a, mean, sigma):
            return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))
        y_at_xmean=yhist[np.where(xhist>np.ma.mean(theta))[0][0]]
        popt, pcov = curve_fit(gaussian, xhist, yhist, [y_at_xmean, np.ma.mean(theta), np.ma.std(theta)])
        ax.plot(thetaax, gaussian(thetaax, *popt))
        gcurve=gaussian(thetaax, *popt)
        maxval= np.max(gcurve)#,np.min(gaussian(thetaax, *popt))
        #import matplotlib.mlab as mlab
        #gauss=mlab.normpdf(thetaax,np.mean(theta),np.var(theta))
        #ax.plot(thetaax,gauss)
        print np.where(gcurve==maxval), thetaax[np.where(gcurve==maxval)[0]]
    #if side==1.0:
    #    ang=[aa['nominal angle'] for aa in windows if aa['number']==windownum]
    #else:
    #    ang=[aa['nominal angle'] for aa in windowsr if aa['number']==windownum]
    ang=nang

    if stats: #print out and pickle stats
        avg=np.ma.mean(theta)
        med=np.ma.median(theta)
        stdv=np.ma.std(theta)
        statfile='win'+str(windownum)+'_angle_stats_'+str(mag)+'.p'
        datafile='win'+str(windownum)+'_angle_data'+str(mag)+'.p'
        print "-------------------STATISTICS FOR WINDOW "+str(windownum)+"---------------------"
        print '     Nominal Angle: ' + str(ang)
        print '              Mean: ' + str(avg)
        print '            Median: ' + str(med)
        print 'Standard Deviation: ' + str(stdv)
        print 'Results saved in ' + statfile
        stdict={'mean':avg,'median':med,'stddev':stdv}
        pickle.dump(stdict,open(statfile,'wb'))
        pickle.dump(theta,open(datafile,'wb'))


    if windownum:
        #ang=[aa['nominal angle'] for aa in windows if aa['number']==windownum]
        xran=[ang-2.,ang+2.]
        title='Slat angle distribution for window '+str(windownum) #+ ', nominal angle '+str(ang[0])[:-3] + ' degrees'
        #print title
    if not title:
        ax.set_title('Distribution of line orientations')
    else:
        ax.set_title(title)
    if xran:
        ax.set_xlim(xran)
    else:
        ax.set_xlim([thetaran[30]*(180./np.pi),thetaran[70]*(180./np.pi)])
    ax.set_xlabel('Probablistic Hough line angle (degrees)')
    if log:
        ax.set_yscale('log')
        ax.set_ylim([1,100000])
    #else:
    #    ax.set_ylim([0,1])
    fig.show()
    if figname:
        fig.savefig(figname)
    print thetaran[0]*180./np.pi,thetaran[-1]*180./np.pi
    #plt.close(fig)
    if ret:
        return theta

def straight_hough_hist(lines, windownum, mag=5.0,log=True,ret=False, title=False,xran=False,stats=True,figname=False,side=1.0,spread=5.,n=201): #cuz of mpi
    '''Make a histogram of the line orientations returned by straight line Hough'''
    theta=[]
    if type(lines) == list: #it's a list of filenames
        for f in lines:
            theta.append(pickle.load(open(f,'rb'))*180./np.pi)
    else:
        theta=pickle.load(open(f,'rb'))*180./np.pi

    #theta=np.array(theta)*180./np.pi
    thetaran=get_theta_range(windownum,side=side,n=n,spread=spread)
    #make a histogram
    fig,ax=plt.subplots()

    if theta !=[]:
        foo=ax.hist(theta,np.arange(np.min(theta),np.max(theta),.05))

    if side==1.0:
        ang=[aa['nominal angle'] for aa in windows if aa['number']==windownum]
    else:
        ang=[aa['nominal angle'] for aa in windowsr if aa['number']==windownum]

    if stats: #print out and pickle stats
        avg=np.mean(theta)
        med=np.median(theta)
        stdv=np.std(theta)
        statfile='win'+str(windownum)+'_angle_stats_'+str(mag)+'.p'
        datafile='win'+str(windownum)+'_angle_data'+str(mag)+'.p'
        print "-------------------STATISTICS FOR WINDOW "+str(windownum)+"---------------------"
        print '     Nominal Angle: ' + str(ang[0])[:-3]
        print '              Mean: ' + str(avg)
        print '            Median: ' + str(med)
        print 'Standard Deviation: ' + str(stdv)
        print 'Results saved in ' + statfile
        stdict={'mean':avg,'median':med,'stddev':stdv}
        pickle.dump(stdict,open(statfile,'wb'))
        pickle.dump(theta,open(datafile,'wb'))


    if windownum:
        #ang=[aa['nominal angle'] for aa in windows if aa['number']==windownum]
        xran=[side*ang[0]-2.,side*ang[0]+2.]
        title='Slat angle distribution for window '+str(windownum) #+ ', nominal angle '+str(ang[0])[:-3] + ' degrees'
        #print title
    if not title:
        ax.set_title('Distribution of line orientations')
    else:
        ax.set_title(title)
    if xran:
        ax.set_xlim(xran)
    else:
        ax.set_xlim([thetaran[30]*(180./np.pi),thetaran[70]*(180./np.pi)])
    ax.set_xlabel('Probablistic Hough line angle (degrees)')
    if log:
        ax.set_yscale('log')
        ax.set_ylim([1,100000])
    else:
        ax.set_ylim([0,len(lines)])
    fig.show()
    if figname:
        fig.savefig(figname)
    print thetaran[0]*180./np.pi,thetaran[-1]*180./np.pi
    #plt.close(fig)
    if ret:
        return theta

def get_length(line):
    '''Get length of line via Pythagoras'''
    deltax=line[1][0]-line[0][0]
    deltay=line[1][1]-line[0][1]
    length=(deltax**2 + deltay**2)**0.5
    return length

def get_angle(line):
    '''Get angle of line via tangent. Note that because the top left corner is 0,0 in the background image we multiply x's by -1'''
    deltax=-1.*(line[1][0]-line[0][0]) #removed -1.*
    deltay=line[1][1]-line[0][1]
    if np.abs(deltax)==0.0:
        thetadeg=90.0
    elif np.abs(deltay)==0.0:
        thetadeg=0.
    else:
        theta=np.arctan(float(deltay)/float(deltax))  #np.arctan2(float(deltay)/float(deltax))
        thetadeg=np.rad2deg(theta) #theta*180./np.pi
    return thetadeg

def get_intersect(m1,m2,b1,b2):
    '''get intersection point (x,y) between two lines'''
    x=(b2-b1)/(m1-m2)
    y=m1*x+b1
    return x,y

def find_nearest(array, value,ret_idx=True):
    ''' find closest value in array to given value'''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if ret_idx:
        return idx,array[idx]
    else:
        return array[idx]

def GW_angle_wrt_position(win,idx,ll=[50,100,150,200]):
    for l in ll:
        fn=glob.glob('window00'+str(win)+'_000'+str(idx)+'_flatdark_hough_ll'+str(l)+'.p')
        lines=pickle.load(open(fn[0],'rb'))
        theta=[get_angle(line) for line in lines]
        hough_peek(fn[0],edgef=fn[0][:fn[0].rfind('_')-5]+'edges.p')
        fig,ax=plt.subplots()
        ax.scatter(range(0,len(theta)),theta)
        ax.plot(range(0,len(theta)),theta)
        fig.show()


