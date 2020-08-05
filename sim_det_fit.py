"""
===================
sim_det_fit.py
Erica  Lastufka 9.3.18
===================

Given nxm pixel detector, degrade the moire image (or 'bin') then fit sines or triangle functions to calculate the moire period and orientation

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
from scipy import optimize
from plot_divergent import plot_moire_fancy
import edge_and_hough as eh

#list of dictionaries of the window numbers, pitches (mm), and nominal angles (deg)
global moire

moire=[{'number':11,'period_px':99., 'angle': 31.,'amp':.10},
    {'number':21,'period_px':58., 'angle': -35.,'amp':.10},
    {'number':12,'period_px':43., 'angle': -49.,'amp':.10},
    {'number':22,'period_px':22., 'angle': 42,'amp':.10},
    {'number':31,'period_px':60., 'angle': 49.,'amp':.10},
    {'number':41,'period_px':40., 'angle': -50.,'amp':.10},
    {'number':32,'period_px':18., 'angle': -42.,'amp':.10},
    {'number':42,'period_px':22., 'angle': -42.,'amp':.10},
    {'number':33,'period_px':40., 'angle': 39.,'amp':.10},
    {'number':43,'period_px':35., 'angle': -39.,'amp':.10},
    {'number':34,'period_px':22., 'angle': -45.,'amp':.10},
    {'number':44,'period_px':18., 'angle': 45.,'amp':.10}]

def im2ndarray(filen):
    '''Convert .tif image to numpy.ndarray'''
    imraw = Image.open(filen) #.rotate(45)#.convert('1').convert('L')
    im=np.array(imraw)
    return im

def remove_edges(im,ncrop=50,plot=False):
    '''remove edges from images with edges of the window in them'''
    imshape=np.shape(im)
    im=im[ncrop:imshape[0]-ncrop,ncrop:imshape[1]-ncrop]
    percentage = np.float((np.shape(im)[0]*np.shape(im)[1]))/np.float((imshape[0]*imshape[1]))
    print "Percentage cropped:", (1.-percentage)*100.
    #print np.shape(im),imshape
    if plot:
        fig,ax=plt.subplots()
        ax.imshow(im,cmap=cm.gray)
        fig.show()
    return im

def contrast_stretch(im,plow=2,phigh=98):
    pl = np.percentile(im, plow)
    ph = np.percentile(im, phigh)
    im2 = exposure.rescale_intensity(im, in_range=(pl, ph))
    return im2

def det_degrade(n,m,im,pix2um,plot=True,nomsize=[110.,110.]):
    '''degrade image given a detector arrangement. nomsize is nominal size of the window and/or detector in mm'''
    imx,imy=np.size(im)
    pix2mm=pix2um/1000.
    #create bins 
    binx=[(,) for aa in range(0,)]
    binx=[(,) for aa in range(0,)]    
    #sum values of image in bins (should the whole thing be normalised first or not?)
    for nn,mm in zip(range(0,n),range(0,m)):
        degarray[nn,mm]=np.sum(im[bx[0]:bx[1],mm:])

    if plot:
        fig,ax=plt.subplots()
        ax.imshow(degarray,cmap=cm.gray)
        fig.show()
        
def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def fit_sine(im,im0,period,amp,plot=True):
    '''fit sine to gradient'''
    xax=np.arange(0,len(im))
    # Fit the first set
    fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/p[1]*(x-p[4]))+ p[2]+p[3]*x # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    ###find initial fit params####
    zl=(np.diff(np.sign(np.gradient(im))) !=0)*1
    yoff= np.where(zl != 0)[0] #location of first sign change in gradient
    #amp= np.max(im0)-np.min(im0)#(diff between min and max values of actual image?)
    slope= (np.mean(im[-10:])-np.mean(im[:10]))/float(len(xax)) #fit a line, get the slope
    xoff= yoff - (period/2.)
    
    p0 = [amp, period,yoff[0],slope,xoff[0]] # Initial guess for the parameters
    p1, success = optimize.leastsq(errfunc, p0[:], args=(xax,im))
    print p1
    fig,ax=plt.subplots()
    ax.plot(xax, im, "ro", xax, fitfunc(p1, xax), "k-") # Plot of the data and the fit
    fig.show()
    return p1

def fit_sine_interactive(win,pix2um=24.,plot=False,edges=False,hough=False,rmean=False):
    wdict = [wd for wd in moire if wd['number'] == win]
    ang= -1*wdict[0]['angle']
    period= wdict[0]['period_px']
    amp= wdict[0]['amp']
    imf='win'+str(win)+'.tif'
    im=im2ndarray(imf)
    proceed=False
    while not proceed:
        crop=raw_input('crop number of pixels: ')
        im=remove_edges(im,ncrop=int(crop),plot=True)
        pstr=raw_input('proceed? (y/n) ')
        if pstr =='y':
            proceed=True
    im=contrast_stretch(im)
    rs=rot_sum(im,ang)
    fit = True
    trim=False
    while fit:
        if hough:
            edges=False
            edata=Canny_edge(im)
            ang,im= hough2im(edata,plot=plot)
        print 'rotation angle ',ang
        rs=rot_sum(im,ang,plot=plot,edges=edges,rmean=rmean)
        if trim:
            rs=rs[trim:-trim]
        p1=fit_sine(rs,im,period,amp)
        fstr=raw_input('try again? [period, amp,angle,trim] or "n" :')
        if fstr == "n":
            fit = False
        else:
            flist=[fstr[1:-1].split(',')][0]
            period=float(flist[0])
            amp=float(flist[1])
            ang=float(flist[2])
            trim=float(flist[3])
    plotmoire=raw_input('plot moire? (y,n)')
    if plotmoire == 'y':
        periodmm=2.*p1[1]*pix2um/1000.
        print 'period (mm) ',periodmm
        plot_moire_fancy([periodmm],[90.+ang])

def calc_dist(intpoints):
    '''now calculate the distances between the intersection points'''
    dist=[]
    for i,p1 in enumerate(intpoints):
        try:
            p2=intpoints[i+1]
            dist.append(np.sqrt((p2[0]-p1[0])**2+(p2[1]-p1[1])**2))
        except IndexError:
            break
    return dist

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
