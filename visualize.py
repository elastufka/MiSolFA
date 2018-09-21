"""
===================
visualize.py
Erica  Lastufka 5.9.18
===================

Plot data passed in via Results object
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

