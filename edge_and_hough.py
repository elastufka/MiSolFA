"""
===================
edge_and_hough.py
Erica  Lastufka 20.9.17
===================

Module to do Canny edge detection at desired sigma level, then continue with probabilistic Hough. Methods to write/store results. Use __main__ to run the whole thing eg from a mpi interface? first maybe do some timing tests..

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

#list of dictionaries of the window numbers, pitches (mm), and nominal angles (deg)
global windows
global windowsr

windows=[{'number':11,'pitch':89.6752,  'nominal angle':-44.79357},
                {'number':21,'pitch':90.326,  'nominal angle': 45.20793},
                {'number':12,'pitch':22.4797,'nominal angle': 44.94825},
                {'number':22,'pitch':22.5203,'nominal angle':-45.05184},
                {'number':31,'pitch':45.0814, 'nominal angle':-45.10378},
                {'number':41,'pitch':44.9187, 'nominal angle': 44.8966},
                {'number':32,'pitch':18.013, 'nominal angle': 45.04146},
                {'number':42,'pitch':17.987, 'nominal angle':-44.95859},
                {'number':33,'pitch':29.9639,  'nominal angle':-44.93102},
                {'number':43,'pitch':30.0362,  'nominal angle': 45.06914},
                {'number':34,'pitch':14.991, 'nominal angle': 44.96549},
                {'number':44,'pitch':15.009, 'nominal angle':-45.03455}]

windowsr=[{'number':11,'pitch':90.326,  'nominal angle':-45.20793},
                {'number':21,'pitch':89.6752,  'nominal angle': 44.79357},
                {'number':12,'pitch':22.5203,'nominal angle': 45.05184},
                {'number':22,'pitch':22.4797,'nominal angle':-44.94825},
                {'number':31,'pitch':44.9187, 'nominal angle':-44.8966},
                {'number':41,'pitch':45.0814, 'nominal angle': 45.10378},
                {'number':32,'pitch':17.987, 'nominal angle': 44.95859},
                {'number':42,'pitch':18.013, 'nominal angle':-45.04146},
                {'number':33,'pitch':30.0362,  'nominal angle':-45.06914},
                {'number':43,'pitch':29.9639,  'nominal angle': 44.93102},
                {'number':34,'pitch':15.009, 'nominal angle': 45.03455},
                {'number':44,'pitch':14.991, 'nominal angle':-44.96549}] #dectector side angle for ease, windows swapped...


def get_coordinates(path,filen,zcoord=False):
    '''Get coordinates of all given images in list filen. If there's a random z-coordinate stuck in the middle of the filename on the log readout, deal with it'''
    #get corresponding log file
    os.chdir(path)
    logs=glob.glob('*.txt')
    print logs
    #parse log file for coordinates of image
    cdict={}
    for log in logs:
        with open(log,'rb') as f:
            lines=f.readlines()
            for i,line in enumerate(lines):
                if 'Image saved:' in line:
                    fn=line[line.find(':')+2:line.find('.tif')+4]
                    if zcoord:
                        #get rid of the z-coordinate that is in the name in the files but not in the actual files
                        fn=fn[:fn.find('Z')-1]+fn[fn.rfind('_'):]
                    if fn in filen:
                        #print fn
                        windex,index=get_index(fn)
                        info_line=lines[i-1]
                        xmm=info_line[info_line.find('X =')+4:info_line.find('mm')-1]
                        ymm=info_line[info_line.find('Y =')+4:info_line.find('Z')-4]
                        cdict[fn]={'number':windex,'indices':index,'im_coords':[float(xmm),float(ymm)]}
    #return image coordinates
    return cdict

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

def im2ndarray(filen):
    '''Convert .tif image to numpy.ndarray'''
    imraw = Image.open(filen) #.rotate(45)#.convert('1').convert('L')
    im=np.array(imraw)
    return im

def remove_edges(filen,im,maxx='10',maxy='13'):
    '''remove edges from images with edges of the window in them'''
    if filen.endswith('.p') or filen.endswith('mosaic.tiff'):
        imshape=np.shape(im)
        im=im[50:imshape[0]-50,50:imshape[1]-50]
        return im
    else:
        windex,imindex=get_index(filen)
        if imindex[0] == '01': #it's got an edge on the left, trim the array (x and y are flipped for some reason)
            im=im[:,50:]
        if imindex[0] == maxx: #it's got an edge on the right, trim the array
            im=im[:,:-50]
        if imindex[1] == '01': #it's got an edge on the top, so mask the edges
            im=im[50:,:]
        if imindex[1] == maxy: #it's got an edge on the bottom, so mask the edges
            im=im[:-50,:]
        #print np.shape(im)
    return im

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

def EM_list(win,mag,ending='.tif'):
    import glob
    filen=glob.glob('win'+win+'*_'+mag+'X'+ending)
    newfilen=[]
    for f in filen:
        if '01' not in f and '12' not in f[6:] and '09' not in f:
            newfilen.append(f)
    return newfilen

def anisotropic_diffusion_test(filepath, canny=True,sigma=3,ht=.5,lt=.2):
    Dimension = 2
    PixelType = itk.ctype('float')
    ImageType = itk.Image[PixelType, Dimension]
    reader = itk.ImageFileReader[ImageType].New(FileName=filepath)
    reader.Update()
    image = reader.GetOutput()
    #contrast stretching
    #p2 = np.percentile(image, 2)
    #p98 = np.percentile(image, 98)
    #image = exposure.rescale_intensity(image, in_range=(p2, p98))

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

    # Crop to small region of interest to experiment with smoothing parameters
    small_roi_filter = itk.RegionOfInterestImageFilter.New(roi)
    roi_region = itk.ImageRegion[Dimension]()
    index = itk.Index[Dimension]()
    index[0] = 50
    index[1] = 50
    roi_region.SetIndex(index)
    size = itk.Size[Dimension]()
    size[0] = 100
    size[1] = 100
    roi_region.SetSize(size)
    small_roi_filter.SetRegionOfInterest(roi_region)
    small_roi_filter.Update()
    small_roi = small_roi_filter.GetOutput()


    smoother = itk.CoherenceEnhancingDiffusionImageFilter.New(roi)
    # Determine good parameters here: https://insightsoftwareconsortium.github.io/ITKAnisotropicDiffusionLBR/
    smoother.SetDiffusionTime(5)
    smoother.SetLambda(0.05)
    smoother.SetEnhancement(3) #enum EnhancementType {CED, cCED, EED, cEED, Isotropic};
    smoother.SetNoiseScale(2)
    smoother.SetFeatureScale(10)
    smoother.SetExponent(2)
    smoother.Update()

    filtered_small_roi_filter = itk.RegionOfInterestImageFilter.New(smoother)
    filtered_small_roi_filter.SetRegionOfInterest(roi_region)
    filtered_small_roi_filter.Update()
    filtered_small_roi = filtered_small_roi_filter.GetOutput()

    if canny: #do the canny edge detection and overlay
        roiarr=itk.GetArrayFromImage(roi)
        smootherarr=itk.GetArrayFromImage(smoother)
        small_roiarr=itk.GetArrayFromImage(small_roi)
        filtered_small_roiarr=itk.GetArrayFromImage(filtered_small_roi)
        #contrast stretching
        def cstretch(im):
            p2 = np.percentile(im, 2)
            p98 = np.percentile(im, 98)
            im = exposure.rescale_intensity(im, in_range=(p2, p98))
            return im
        csroiarr=cstretch(roiarr)
        cssmootherarr=cstretch(smootherarr)
        cssmall_roiarr=cstretch(small_roiarr)
        csfiltered_small_roiarr=cstretch(filtered_small_roiarr)


        input_edges=feature.canny(csroiarr, sigma=sigma, high_threshold=ht, low_threshold=lt,use_quantiles=True)
        smoothed_edges=feature.canny(cssmootherarr, sigma=sigma, high_threshold=ht, low_threshold=lt,use_quantiles=True)
        roi_edges=feature.canny(cssmall_roiarr, sigma=sigma, high_threshold=ht, low_threshold=lt,use_quantiles=True)
        sroi_edges=feature.canny(csfiltered_small_roiarr, sigma=sigma, high_threshold=ht, low_threshold=lt,use_quantiles=True)

    fig, axs = plt.subplots(2, 2, figsize=(12,12))
    #print(axs)
    axs[0,0].imshow(itk.GetArrayViewFromImage(roi), cmap='gray', vmin=0)
    axs[0,0].axis('off')
    axs[0,0].set_title('Input')
    axs[0,1].imshow(itk.GetArrayViewFromImage(smoother), cmap='gray', vmin=0)
    axs[0,1].axis('off')
    axs[0,1].set_title('Smoothed')
    axs[1,0].imshow(itk.GetArrayViewFromImage(small_roi), cmap='gray', vmin=0)
    axs[1,0].axis('off')
    axs[1,0].set_title('Region of Interest')
    axs[1,1].imshow(itk.GetArrayViewFromImage(filtered_small_roi), cmap='gray', vmin=0)
    axs[1,1].axis('off')
    axs[1,1].set_title('Smoothed Region of Interest')
    if canny:
        axs[0,0].imshow(np.ma.masked_where(input_edges == 0,input_edges),cmap=cm.autumn)
        #scalebar = ScaleBar(pix2um,'um', SI_LENGTH,length_fraction=length) # 1 pixel = 0.2 meter
        #axs[0,0].add_artist(scalebar)
        axs[0,1].imshow(np.ma.masked_where(smoothed_edges == 0,smoothed_edges),cmap=cm.autumn)
        #axs[0,1].add_artist(scalebar)
        axs[1,0].imshow(np.ma.masked_where(roi_edges == 0,roi_edges),cmap=cm.autumn)
        #axs[1,0].add_artist(scalebar)
        axs[1,1].imshow(np.ma.masked_where(sroi_edges == 0,sroi_edges),cmap=cm.autumn)
        #axs[1,1].add_artist(scalebar)
    fig.show()
    fig.savefig('comparison_figure.png')

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


def Canny_edge(filen,sigma=3,mag=5.0,anisotropic=True,binary=False,gauss=False,plot=False,outfilen=False):
    #max contrast
    if filen.endswith('.p'):
        imarr=pickle.load(open(filen,'rb'))
    else:
        imarr=im2ndarray(filen)

    #im = exposure.equalize_hist(imarr)    #should probably do this AFTER removing the edges...
    if mag==5.0:
        im=remove_edges(filen,imarr)
    else:
        im=remove_edges(filen,imarr,maxx='12',maxy='16')

    #for the finest grids, use sigma=2
    windex,foo=get_index(filen)
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

    p2 = np.percentile(im, 2)
    p98 = np.percentile(im, 98)
    im = exposure.rescale_intensity(im, in_range=(p2, p98))

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
    return newfilen

def print_all_stats(win):
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


def compare_methods(win,ll=[50,100,150],tol=3,sigma=2,first=True,second=True,third=True,fourth=True,nang=-44.7935):
    '''compare 3 different pre-processing methods and their effect on the overall statistics'''
    filen=glob.glob('win'+str(win)+'*_5.0X.tif')
    filep,fileh=[],[]

    #first method - usual
    if first:
        for f in filen:
            filep.append(Canny_edge(f,mag=5.0,anisotropic=False, binary=False,outfilen=f[:-4]+"_edges.p",sigma=sigma))
        for l in ll:
            tag='ll'+str(l)
            for f in filep:
                fileh.append(prob_hough(f,line_length=l,spread=2.,tag=tag,overwrite=True)) #new list of picklefiles

        cat_hough(win,weights=[1,1,1],tags=['ll'+str(ll[0]),'ll'+str(ll[1]),'ll'+str(ll[2])],mag=5.0,EM=False)
        fileh=glob.glob('win'+str(win)+'*_5.0X_hough_all.p')
        figname='win'+str(win) + 'ang_hist.png'
        theta=hough_hist(fileh,win,figname=figname,tol=tol,side=1.,spread=2.,ret=True)
        offset=np.mean(theta)-nang
        figname='win'+str(win) + 'width_hist.png'
        if len(filen) !=0:
            w,b,f=get_period_by_grouping(win,ftags='X',mosaic=False,EM=False,tolerance=.85,offset=offset,mag=5.0,side=1,fitangle=False,plot=True)

    #second - anisotropic diffusion
    if second:
        filep,fileh=[],[]
        for f in filen:
            filep.append(Canny_edge(f,mag=5.0,anisotropic=True, binary=False,outfilen=f[:-4]+"_ani_edges.p",sigma=sigma))
        for l in ll:
            tag='ani_ll'+str(l)
            for f in filep:
                fileh.append(prob_hough(f,line_length=l,spread=2.,tag=tag,overwrite=True)) #new list of picklefiles

        cat_hough(win,imtags='X_ani_',weights=[1,1,1],tags=['ani_ll'+str(ll[0]),'ani_ll'+str(ll[1]),'ani_ll'+str(ll[2])],mag=5.0,EM=False,outtag='ani')
        fileh=glob.glob('win'+str(win)+'*_5.0X_ani_hough_ani.p')
        figname='win'+str(win) + 'ang_hist_ani.png'
        theta=hough_hist(fileh,win,figname=figname,tol=tol,side=1.,spread=2,ret=True)
        offset=np.mean(theta)-nang
        figname='win'+str(win) + 'width_hist_ani.png'
        if len(filen) !=0:
            w,b,f=get_period_by_grouping(win,ftags='X_ani*',mosaic=False,EM=False,tolerance=.85,offset=offset,mag=5.0,side=1.,fitangle=False,plot=True)

    #third - binary
    if third:
        filep,fileh=[],[]
        for f in filen:
            filep.append(Canny_edge(f,mag=5.0,anisotropic=False, binary=True,outfilen=f[:-4]+"_bin_edges.p",sigma=sigma))
        for l in ll:
            tag='bin_ll'+str(l)
            for f in filep:
                fileh.append(prob_hough(f,line_length=l,spread=2.,tag=tag,overwrite=True)) #new list of picklefiles

        cat_hough(win,imtags='X_bin_',weights=[1,1,1],tags=['bin_ll'+str(ll[0]),'bin_ll'+str(ll[1]),'bin_ll'+str(ll[2])],mag=5.0,EM=False,outtag='bin')
        fileh=glob.glob('win'+str(win)+'*_5.0X_bin_hough_bin.p')
        figname='win'+str(win) + 'ang_hist_bin.png'
        theta=hough_hist(fileh,win,figname=figname,tol=tol,side=1.,spread=2.,ret=True)
        offset=np.mean(theta)-nang
        figname='win'+str(win) + 'width_hist_bin.png'
        if len(filen) !=0:
            w,b,f=get_period_by_grouping(win,ftags='X_bin*',mosaic=False,EM=False,tolerance=.85,offset=offset,mag=5.0,side=1.,fitangle=False,plot=True)

    #fourth - both
    if fourth:
        filep,fileh=[],[]
        for f in filen:
            filep.append(Canny_edge(f,mag=5.0,anisotropic=True, binary=True,outfilen=f[:-4]+"_both_edges.p",sigma=sigma))
        for l in ll:
            tag='bin_ll'+str(l)
            for f in filep:
                fileh.append(prob_hough(f,line_length=l,spread=2.,tag=tag,overwrite=True)) #new list of picklefiles

        cat_hough(win,imtags='X_both_',weights=[1,1,1],tags=['bin_ll'+str(ll[0]),'bin_ll'+str(ll[1]),'bin_ll'+str(ll[2])],mag=5.0,EM=False,outtag='bin')
        fileh=glob.glob('win'+str(win)+'*_5.0X_both_hough_bin.p')
        figname='win'+str(win) + 'ang_hist_both.png'
        theta=hough_hist(fileh,win,figname=figname,tol=tol,side=1.,spread=2.,ret=True)
        offset=np.mean(theta)-nang
        figname='win'+str(win) + 'width_hist_both.png'
        if len(filen) !=0:
            w,b,f=get_period_by_grouping(win,ftags='X_both*',mosaic=False,EM=False,tolerance=.85,offset=offset,mag=5.0,side=1.,fitangle=False,plot=True)

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
    ax.imshow(im, cmap=cm.gray,alpha=0.6)
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

def rising_or_falling(rotim,rotarr,immean,plot=False,shape=False,rerotate=False, test=False):
    """Determine if an edge is rising (dark -> light) or falling (light -> dark). Return indexed mask"""
    imsize=np.shape(rotim)
    if not shape:
        mask=np.zeros([792,792])#np.zeros([640,480])
        mask1=np.zeros([792,792])
        oshape=np.zeros([640,480])
    else:
        mask=np.zeros(shape)
        mask1=np.zeros(shape)
        oshape=np.zeros([640,480])
    meani=np.mean(rotim, axis=0) #average along y
    imsum=np.sum(rotim,axis=0)
    imgrad=np.gradient(np.array(imsum,dtype=float))
    meanarr=np.sum(rotarr, axis=0) # total along y (remember it's boolean)
    #for rowi,rowarr in zip(meani,meanarr):
    tv=np.where(meanarr > .5*np.max(meanarr[0:50]))
    tvshape=np.shape(tv)[1]
    #get rising
    rising=[tv[0][i] for i in range(0,tvshape-1) if np.mean(meani[tv[0][i]:tv[0][i+1]]) > 1.25*immean and tv[0][i+1] != tv[0][i]+1 ]
    rising1=[tv[0][i] for i in range(0,tvshape-1) if np.mean(imgrad[tv[0][i]-1:tv[0][i]+2]) > 0.] #take 3 pix around xpeak
    mask1[:,rising1]=1
    mask[:,rising]= 1
    #get falling
    falling=[tv[0][i] for i in range(0,tvshape-1) if np.mean(meani[tv[0][i]:tv[0][i+1]]) < .75*immean and tv[0][i+1] != tv[0][i]+1 ]
    falling1=[tv[0][i] for i in range(0,tvshape-1) if np.mean(imgrad[tv[0][i]-1:tv[0][i]+2]) < 0.]
    mask[:,falling]= -1
    mask1[:,falling1]=-1



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

        fig,ax=plt.subplots(1,2)
        #if not imfile:
            #ax.imshow(rotate(np.transpose(im2ndarray('win11_05_05_5.0X.tif')),-44.79357,reshape=True),alpha=0.6,cmap=cm.gray)
        ax[0].imshow(rotim,alpha=0.6,cmap=cm.gray,origin='lower')
        ax[1].imshow(rotim,alpha=0.6,cmap=cm.gray,origin='lower')
        ax[0].plot(range(0,792),meanarr)
        ax[0].axhline(.5*np.max(meanarr[0:50]))
        ax[1].plot(range(0,792),imgrad/16.)
        #else:
            #ax.imshow(rotate(np.transpose(im2ndarray(imfile)),-45.03455,reshape=True),alpha=0.6,cmap=cm.gray)
        #    ax.imshow(imarr,alpha=0.6,cmap=cm.gray)
        ax[0].imshow(mask,cmap=cm.gray,alpha=0.7,origin='lower')
        ax[1].imshow(mask1,cmap=cm.gray,alpha=0.7,origin='lower')
        ax[0].set_title('original')
        ax[1].set_title('using gradient')
        fig.show()

    if rerotate: #rotate it back by given angle and trim to original size
        mask=imrotate(mask,rerotate)
    #if len(rising) > 1.25* len(falling) or len(falling) > 1.25* len(rising):
        #return []
    return mask1

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

def build_mosaic_mask():
    return 'foo'

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
        mask=rising_or_falling(rotim[300:450,:],rotarr[300:450,:],np.mean(rotim[300:450,100:-100]), shape=np.shape(rotarr),plot=False)
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


def get_theta_range(edges,side=1.0,spread=5.,n=201):
    #define range of theta around theta_nominal
    if type(edges) == int:
        winnum=edges
    elif 'win' not in edges:
        try:
            winnum=int(edges)
        except ValueError:
            edges=raw_input('What is the window number?')
            winnum=int(edges)
    else:
        winnum,index=get_index(edges)
        winnum=int(winnum)
    if side==1.0:
        nang=[w['nominal angle'] for w in windows if w['number'] == winnum]
        nang=nang[0]
    else:
        nang=[w['nominal angle'] for w in windowsr if w['number'] == winnum]
        nang=-1*nang[0]
    theta0= nang*(np.pi/180.)#in radians
    #tendeg2rad=np.pi/18.
    spreaddeg2rad=spread*(np.pi/180.)
    thetaran = np.linspace(theta0-spreaddeg2rad, theta0+spreaddeg2rad, num=n)#in radians
    return thetaran

def prob_hough(edges, threshold=10, line_length=50, line_gap=2,retlines=False,plot=False,side=1.0,spread=5.,n=201, tag=False,overwrite=False):
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

    thetaran=get_theta_range(edges,spread=spread,n=n,side=side)
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
            lines=pickle.load(open(basen+'_'+tag+'.p','rb'))
            #print len(lines)
            if not weights:
                all_lines.extend(lines)
            else:
                for j in range(0,weights[i]):
                    all_lines.extend(lines)
            #print len(all_lines)
        pickle.dump(all_lines,open(basen+'_'+outtag+'.p','wb'))

def hough_hist(lines, windownum, mag=5.0,log=True,ret=False, title=False,xran=False,stats=True,figname=False,side=1.0,spread=5.,n=201,gaussfit=False, sameline=True,tol=3): #cuz of mpi
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

    thetaran=get_theta_range(windownum,side=side,n=n,spread=spread)
    #make a histogram
    fig,ax=plt.subplots()

    thetaax=np.arange(np.min(theta),np.max(theta),.05)
    if theta !=[]:
        yhist,xhist=np.histogram(theta,thetaax)
        foo=ax.hist(theta,thetaax)

    if gaussfit:
        xh = np.where(yhist > 0)[0]
        yh = yhist[xh]
        #print np.max(yhist), np.mean(yhist), np.std(yhist)
        def gaussian(x, a, mean, sigma):
            return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))
        y_at_xmean=yhist[np.where(xhist>np.mean(theta))[0][0]]
        popt, pcov = curve_fit(gaussian, xh, yh, [y_at_xmean, np.mean(theta), np.std(theta)])
        ax.plot(thetaax, gaussian(thetaax, *popt))
        print np.max(gaussian(thetaax, *popt)),np.min(gaussian(thetaax, *popt))
        #import matplotlib.mlab as mlab
        #gauss=mlab.normpdf(thetaax,np.mean(theta),np.var(theta))
        #ax.plot(thetaax,gauss)
        #print np.where(gauss==np.max(gauss)), thetaax[np.where(gauss==np.max(gauss))[0]]
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
        print '     Nominal Angle: ' + str(ang[0])
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
    deltax=-1*(line[1][0]-line[0][0])
    deltay=line[1][1]-line[0][1]
    theta=np.arctan(float(deltay)/float(deltax))  #np.arctan2(float(deltay)/float(deltax))
    thetadeg=np.rad2deg(theta) #theta*180./np.pi
    return thetadeg

##plot outilers
#out1,out2=[],[]
#for l1,l2 in zip(lines1,lines2):
#    if get_angle(l1) !=-45.0:
#       #theta1.append(get_angle(l1))
#       out1.append(l1)
#    if get_angle(l2) !=-45.0:
#       #theta2.append(get_angle(l2))
#       out2.append(l2)

#fig,(ax1,ax2)=plt.subplots(1,2,figsize=(5,6), sharex=True,sharey=True)
#ax1.imshow(imraw, cmap=cm.gray)
#ax2.imshow(imraw, cmap=cm.gray)
#ax1.imshow(np.ma.masked_where(edges2 == 0,edges2),cmap=cm.cool,alpha=.7)
#3ax2.imshow(np.ma.masked_where(edges2 == 0,edges2),cmap=cm.cool,alpha=.7)

#for line in out1:
#    p0, p1 = line
#    ax1.plot((p0[0], p1[0]), (p0[1], p1[1]), color='r')
#ax1.set_xlim((0, np.shape(imraw)[0]))
#ax1.set_ylim((np.shape(imraw)[1], 0))
#ax1.set_title('Probabilistic Hough')

#ax1.set_axis_off()
#ax1.set_adjustable('box-forced')
#ax2.set_axis_off()
#ax2.set_adjustable('box-forced')

#fig.show()
