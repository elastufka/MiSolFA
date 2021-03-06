"""
===================
analyze_general_tidy.py
Erica  Lastufka 5.9.18
===================

General functions for class Analyze. To be inherited by analyze_optical and analyze_xray, and called by Grating() objects. Tidied up for external users
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

#####################################  IMAGE PREPROCESSING ###########################################

def im2ndarray(filen):
    '''Convert .tif image to numpy.ndarray'''
    imraw = Image.open(filen) #.rotate(45)#.convert('1').convert('L')
    im=np.array(imraw)
    return im

def remove_edges_general(filen,lpix=50,rpix=50,tpix=50,bpix=50,write=True):
    '''crop edges from images by given number of pixels in each direction. Option to write resulting image to input filename, which adds the suffix 'uncropped' to the original filename and saves the full image there.
    Input: string name of .tif file.
    Output: numpy array of image'''
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
    '''Average all the flats. Assumes flats are after darks in sequential file numbering'''
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

def contrast_stretch(imf,tag='cstretch',saveim=True,t=98,b=2):
    '''as in http://scikit-image.org/docs/0.9.x/auto_examples/plot_equalize.html
    Input: image file name
    Output: name of new image file'''
    if type(imf) == str:
        im=im2ndarray(imf)
    p2 = np.percentile(im, b)
    p98 = np.percentile(im, t)
    corrected_image = exposure.rescale_intensity(im, in_range=(p2, p98))
    savec=scipy.misc.toimage(im,high=np.max(im),low=np.min(im),mode='F')
    imname=imf[:-12]+tag+'.tif'
    if saveim:
        savec.save(imname)
    return imname

def correct_image(im,dark,flat,fac=False,saveim=False, pickleimarr=False, contrast_stretch=False,show=False,ret=False):
    '''Dark- and flatfield correct X-ray images, option to contrast stretch the image as well. '''
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
    if ret:
        return corrected_image

def check_correction():
    ''' check that the correction works for all the datasets'''
    from skimage import exposure
    dirs=['EMmodel/SLS_Apr2018','QMmodel/SLS_May2018','QMmodel/SLS_Sept2018']
    #labs=['','',']
    fig,ax=plt.subplots(3,4,figsize=[15,10])
    for d,a in zip(dirs,ax):
        os.chdir(d)
        flat=pickle.load(open('transm'+d[:2]+'_combined_flat.p','rb'))
        dark=pickle.load(open('transm'+d[:2]+'_combined_dark.p','rb'))
        if d =='QMmodel/SLS_Sept2018':
            flat =pickle.load(open('transmQM_flat_post.p','rb'))
        dfile=glob.glob('transm'+d[:2]+'/win*_p3_0.tif')[0]
        testim=im2ndarray(dfile)
        print dfile
        a[0].imshow(testim,origin='bottom left')
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
        a[0].set_title(d[8:]+' original')
        cim=correct_image(testim,dark,flat,fac=fac,ret=True,saveim='testim.tif')
        a[1].imshow(cim,origin='bottom left')#,norm=matplotlib.colors.Normalize(vmin=.9*np.min(cim),vmax=1.1*np.max(cim)))
        chist,cbc=exposure.histogram(cim*256.,nbins=256)
        #print np.shape(chist),np.shape(cbc)
        #a[1].plot(cbc,chist,lw=2)
        a[1].set_title(d[8:]+' corrected')
        a[2].imshow(flat-dark,origin='bottom left')
        print np.shape(tbc),np.shape(fbc)
        #a[2].plot(fbc,fhist,lw=2)
        a[2].set_title(d[8:]+' flat-dark')
        #hist, bins_center = exposure.histogram(camera)
        #plt.plot(bins_center, hist, lw=2)
        a[3].plot(range(0,256),thist,'r',label='original')
        a[3].plot(range(0,256),chist,'b',label='corrected')
        a[3].plot(range(0,256),fhist,'g',label='flat-dark')
        a[3].legend(loc='upper right')

        for aa in a:
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
    ''' Simple method to dark- and flat-field correct the images. Option gw introduced for Gridwinder data.'''
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


def mosaic(window_num, coords, all=False, plot=True,plot_coords=False, mag=5.0):
    '''Make a mosaic of the given window or all windows (individually). Input: list of window numbers, dictionary of coordinates. (Probably doesn't work currently)'''
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

###################################### CREATING EDGE ARRAYS AND LINE SEGMENTS ########################################

def Canny_edge(filen,sigma=3,gauss=False,plot=False,outfilen=False):
    '''Use scikit_image's Canny edge detector algorithm (http://scikit-image.org/docs/0.9.x/auto_examples/plot_canny.html?highlight=canny%20edge) to find single-pixel edegs of a given input image. Return as ndarray'''
    #max contrast
    if filen.endswith('.p'):
        imarr=pickle.load(open(filen,'rb'))
    else:
        imarr=im2ndarray(filen)

    im=imarr

    if anisotropic:
        im=anisotropic_diffusion(filen)

    p2 = np.percentile(im, 2)
    p98 = np.percentile(im, 98)
    im = exposure.rescale_intensity(im, in_range=(p2, p98))

    if gauss:
        im = ndi.gaussian_filter(im, gauss)

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

def prob_hough(edges, nang,threshold=10, line_length=50, line_gap=2,retlines=False,plot=False,spread=5.,n=201, tag=False,overwrite=False):
    '''Perform probabilistic Hough transform (http://scikit-image.org/docs/0.9.x/auto_examples/plot_line_hough_transform.html?highlight=probabilistic%20hough) to given set of edges'''
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

def get_theta_range(nang,spread=5.,n=201):
    '''define range of theta around theta_nominal. For use with probablistic Hough fitting'''
    theta0= nang*(np.pi/180.)#in radians
    #tendeg2rad=np.pi/18.
    spreaddeg2rad=spread*(np.pi/180.)
    thetaran = np.linspace(theta0-spreaddeg2rad, theta0+spreaddeg2rad, num=n)#in radians
    return thetaran

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

def im_peek(filen,mag=5.0,length=0.2,contrast_stretch=False):
    '''Plot Canny edges over image for a given file. Pixel size defaults and filename protocols are set for MiSolFA optical/xray datasets'''
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib_scalebar.scalebar import SI_LENGTH
    if mag==5.0:
        pix2um=(1.2512/640.)*1000.
    else:
        pix2um=0.6
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
        pix2um=0.6
    edges=pickle.load(open(filen,'rb'))
    imf=filen[:-8]+'.tif'
    im=im2ndarray(imf)
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    ax.imshow(im, cmap=cm.gray)
    ax.imshow(np.ma.masked_where(edges == 0,edges),cmap=cm.autumn)
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
        pix2um=0.6

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

#################################### MANIPULATING EDGE ARRAYS ###########################################

def clean_centers(edges,cfac=False,tolerance=False,plot=False,sigma=2.):
    '''clean out the Canny edge array so that pixel bins with low counts (ie probably bad edges) get deleted below a certain threshold. Either set the threshold directly (variable tolerance) or set it as a fraction of the standard deviation (sigma, default option). Only works with nearly vertical edges.'''
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
        if cfac:
            ax.set_xlim([0,len(cfac)])
        fig.show()
    return cleaned_edges

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

def plot_centers_and_edges(win,p0,ang,earr=False, datafile=False,tol=2.0,cfac=False):
    '''Check that rising_or_falling_final has identified the rising and falling edges correctly by plotting them overlaid on sum of edge array'''
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

################################ CALCULATING PERIODS AND SLIT WIDTHS #######################################

def slit_widths_from_peaks(window_num,imfile,xpeak=False,pix2um=.65,plot=False,stats=True,gauss=False,tolerance=False,n=9,quiet=True):
    '''Primary way of getting both slit widths and periods from TOMCAT data set.'''
    #if not ang:
    check=False
    p0ang=imfile[6:imfile.rfind('_')]
    #else:
    #    p0ang=ang
    im=np.array(Image.open(imfile))
    imsum=np.sum(im,axis=0)
    immean=np.mean(imsum)
    imgrad=np.gradient(imsum)
    #imgrad2=np.gradient(imgrad)
    #gradpeaks=(np.roll(np.sign(imgrad2),1)-np.sign(imgrad2) !=0).astype(int) #locations where the gradients peak
    #gx=np.where(gradpeaks !=0)[0]
    if not xpeak:
        efile=glob.glob(imfile[:-4]+'_edges.p')
        if len(efile) > 0:
            edges=pickle.load(open(efile[0],'rb'))
        else:
            edges,_=Canny_edge(imfile,sigma=3,gauss=gauss,plot=False,mag=False)
        cleaned_edges=clean_centers(edges,sigma=tolerance)
        xpeak=sum_peaks(cleaned_edges,tol=n)

    xpeak_int=[int(x) for x in xpeak]
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

    rising, falling, periodr,periodf,periods=rising_or_falling_final(xpeak,xpeak_int,imgrad,n=n)
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
        if stdv >=1.:
            check=True
        statfile='win'+str(window_num)+'_width_stats_'+p0ang+'.p'
        datafile='win'+str(window_num)+'_width_data_'+p0ang+'.p'
        print "-------------------STATISTICS FOR WINDOW "+str(window_num)+"---------------------"
        print '              Mean: ' + str(avg)
        print '            Median: ' + str(med)
        print 'Standard Deviation: ' + str(stdv)
        print 'Results saved in ' + statfile
        data={'period':period,'rising':rising,'falling':falling,'widths':width}
        stdict={'mean':avg,'median':med,'stddev':stdv}
        pickle.dump([stdict,tolerance,n],open(statfile,'wb'))
        pickle.dump(data,open(datafile,'wb'))

    return period,width,check

def rising_or_falling_final(xpeak,xpeak_int, imgrad,n=9, quiet=True,filter_multiples=True,filter_nominal=False):
    '''For a given list of peaks, determine if the peak (in the summed edge array) is the result of a rising (dark -> light) or falling (light->dark) edge. '''
    rising=[xpeak[i] for i in range(0,len(xpeak)-1) if np.mean(imgrad[xpeak_int[i]-1:xpeak_int[i]+2]) > 0.] #take 3 pix around xpeak
    #need xpeaks to be integers now
    falling=[xpeak[i] for i in range(0,len(xpeak)-1) if np.mean(imgrad[xpeak_int[i]-1:xpeak_int[i]+2]) < 0.]

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

    width, periodr,periodf=[],[],[]
    #now get the widths
    for i,loc,code in zip(range(0,len(locs[:-2])),locs[:-2],rorf[:-2]):
        testsum=code+rorf[i+1]
        if testsum == 2 and locs[i+1]-loc >n: #two r's
            periodr.append(locs[i+1]-loc)
        elif testsum == 1: #r and f (but in what order?)
            if code ==1: #r then f
                width.append(locs[i+1]-loc)
                if code+rorf[i+2]==2 and locs[i+2]-loc >n:
                    periodr.append(locs[i+2]-loc)
            elif code+rorf[i+2]==0 and locs[i+2]-loc >n: #f then r then f
                periodf.append(locs[i+2]-loc)
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

    return rising,falling,periodr,periodf,period,width


################################### CALCULATING ANGLE DISTRIBUTION #####################################

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
                grouped_lines=same_line(l,tol=tol)
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

#####################################FUNCTIONS DEALING WITH LINES#################################################

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

def extend_line(line,shape=[640,480],plot=False):
    """extend the line as far as it can go within the boundaries of the frame. TOP LEFT corner is origin!!"""
    start=line[0]
    end=line[1]
    dxs,dys=shape[0]-start[0],shape[1]-start[1] #offsets from origin
    slope = (np.float(end[1])-np.float(start[1]))/(np.float(end[0])-np.float(start[0])) #*-1 ?
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

def same_line(lines,tol=3, plot=False):
    """Test if two or more Hough line segments lie on the same line. Returns lists of line segments that fall on the same line, within a certain tolerance at the endpoints"""
    extended_lines,matches,singles=[],[],[]
    duplicate=False
    for line in lines:
        start=line[0]
        end=line[1]
        extended_lines.append(extend_line(line)) #extend the line as far as it can go within the boundaries of the frame

    #now test extended lines against each other
    allxl=[(l[0][0],l[1][0]) for l in lines]
    allyl=[(l[0][1],l[1][1]) for l in lines]
    allx=[(exl[0][0],exl[1][0]) for exl in extended_lines]
    #ally=[(exl[0][1],exl[1][1]) for exl in extended_lines]
    for l,exl in zip(lines,extended_lines):
        #pad the endpoints in x
        ex1=exl[0]
        ex2=exl[1]
        if ex1[0]-tol/2 > 0:
            padded_ex1=range(ex1[0]-tol/2,ex1[0]+tol/2,1)
        else:
            padded_ex1=range(ex1[0],ex1[0]+tol,1)
        if ex2[0]+tol/2 > 0:
            padded_ex2=range(ex2[0]-tol/2,ex2[0]+tol/2,1)
        else:
            padded_ex2=range(ex2[0],ex2[0]+tol,1)
        #test if any other lines fall on the extended line, ie if any in allx are in padded_ex
        #print padded_ex1,padded_ex2
        for m in matches:
            if l in m:
                duplicate=True
                #print l,m, duplicate
            else:
                duplicate = False
                #print l,m, duplicate
        if duplicate == False:
            matches.append([((xl[0],yl[0]),(xl[1],yl[1])) for x,xl,yl in zip(allx,allxl,allyl) if x[0] in padded_ex1 and x[1] in padded_ex2]) #now need to append the actual line segments associated with the extended lines
    for m in matches:
        if len(m) == 1: #find singles
            singles.append(m)
            matches.remove(m)
    #print len(singles)

    unique_matches=[]
    for m in matches: #now matches doesn't contain any singles
        if m not in unique_matches:
            unique_matches.append(m)
            for s in singles:
                #print s,m
                if s[0] in m:
                    singles.remove(s)
                    #print s,m
    #print len(singles)
    for s in singles:
        unique_matches.append(s)

        #if so, ...group and return

    if plot:
        fig,ax=plt.subplots()
        scale=len(unique_matches)
        for i,mm in enumerate(unique_matches):
            for m in mm:
                s1=m[0]
                e1=m[1]
                ax.plot((s1[0],e1[0]),(s1[1],e1[1]),color=cm.jet(float(i)/float(scale)))
                #ax.plot((s1[0],e1[0]),(s1[1],e1[1]),'r--')
            #print i/10.,cm.jet(i/10.)
        fig.show()

    return unique_matches

def cat_hough(window_num, imtags='X_', tags=['ll150','ll200','ll250','ll300'],weights=False,outtag='all',EM=True,mag=5.0):
    """Concatenate all the line segments contained in given files to a single list """
    #first get indices for the window
    if not EM:
        ofiles=glob.glob('win'+str(window_num)+'_*'+str(mag)+imtags+'edges.p')
    else:
        ofiles=EM_list(str(window_num),str(mag),ending='_edges.p')
    idxs=[get_index(o)[1] for o in ofiles]
    #all_lines=[]
    for idx in idxs:
        basen='win'+str(window_num)+'_'+str(idx[0])+'_'+str(idx[1])+'_'+str(mag)+imtags+'hough'
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

