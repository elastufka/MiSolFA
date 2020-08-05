"""
===================
grating.py
Erica  Lastufka 15.5.11
===================
Given the ideal grid parameters, simulate the transmission images/profiles I should get out of them. For a SINGLE grating.

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from PIL import Image
from matplotlib import cm
import pickle
from scipy.misc import imrotate
from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit
from scipy import interpolate
import time
import glob
from fractions import Fraction

def calc_3D_grating(param_dict,source=False,size=[.05,.05]):
    '''Do the actual calculation, the other two are just wrapper methods'''
    #self.nominal={"pitch":,"orientation":,"phase":,"height":,"slit_width":]}
    #probably need the size too

    #get the pixel to um conversion
    if not source:
        pix2um=1.
    else:
        pix2um=source.pix2um

    #intialize array of zeros in the correct size
    arrx=int(size[0]*10000.*pix2um) #cm to um to pix
    arry=int(size[1]*10000.*pix2um) #cm to um to pix
    arrz=int(param_dict['height']*pix2um) #um to pix
    garr=np.zeros([arrx,arry,arrz])

    #convert quantities to pixels
    period_pix=int(pix2um*param_dict['period'])
    angle=param_dict['orientation']
    phase_pix=int(pix2um*param_dict['phase'])
    width_pix=int(pix2um*param_dict['slit_width']) #height doesn't matter...
    slat_pix=period_pix-width_pix
    #print period_pix,width_pix,slat_pix

    #calculate the number of rows you can go down before you need to shift to get the orientation
    tan=np.tan(np.deg2rad(angle))
    if tan >=1:
        angle_factor=int(np.round(tan,2))
    else:
        tan_frac=Fraction(str(np.round(tan,1))).limit_denominator(4)
        if tan_frac.numerator !=1:
            angle_factor=-1*tan_frac.numerator*tan_frac.denominator
        else:
            angle_factor=-1*tan_frac.denominator #-1 in this case indicates that it's less than 1
    #print np.round(tan,1), tan_frac, angle_factor

    #put 1's where the grid slats are. Probably easiest to do this row-by row
    #define phase to be the offset of the rising edge of the first slat w.r.t. corner (0,0,*) in the array
    for i in range(0,arry):
        nperiods=arrx/period_pix #hopefully this is integer
        #print nperiods
        if phase_pix >=0:
            garr[phase_pix:phase_pix+slat_pix,i,:]=1
        else:
            garr[0:phase_pix+slat_pix,i,:]=1
        #print phase_pix,phase_pix+slat_pix
        n=1
        while n<nperiods+2:
            start_pix=phase_pix+n*period_pix
            end_pix=start_pix+slat_pix
            #print n, start_pix,end_pix
            garr[start_pix:end_pix,i,:]=1
            #print np.sum(garr[:,i,:])
            n+=1
        #get the orientation by changing the phase pixel location
        if angle_factor >0:
            if i % angle_factor ==0:
                phase_pix+=1 #of course if this gets bigger than a slit width you have to reset to 0..
                if phase_pix >width_pix:
                    phase_pix=-width_pix
        else:
            phase_pix=phase_pix+-1*angle_factor
            #if phase_pix >width_pix:
            #    phase_pix=-width_pix

    return arrx,arry,garr

def random_gen(low, high):
    while True:
        yield random.randrange(low, high)


def test_thickness_error(arr,dev,elen,percent,slat_width=False):
    #first flatten array
    import time
    import itertools
    arr2d=np.sum(arr,axis=2)
    arrnonzero=np.where(arr2d !=0.0)
    nslatpix=len(arrnonzero[0])
    nomth=np.max(arr2d)
    tmin=nomth-dev
    tmax=nomth+dev
    npixvar=int(tmax)-int(tmin) #number of pixels in variation that is possible
    #varran=np.linspace(int(tmin),int(tmax),npixvar)
    print 'here', tmin,tmax,nslatpix
    if elen == -1: #entire array is affected
        #generate random vector of thicknesses and gaussian smooth it locally? use gaussian with same width as slat if length scale is not defined
        rvec=[]
        gen=random_gen(tmin,tmax)
        #print time.time()
        #for x in itertools.takewhile(lambda x: len(rvec) < nslatpix, gen):
        #    if len(rvec) % 1000 == 0:
        #        print len(rvec)
        #    rvec.add(x)
        while len(rvec) < nslatpix:
            rvec.append(gen.next())
    else: #use elen and percent to determine where to introduce error
        ygen=random_gen(ymin,ymax) #can i just index-ify tuples in arrnonzero?
        xgen=random_gen(xmin,xmax)
    #now assign values in rvec back to slats in arr2D
    for item,locx,locy in zip(rvec,arrnonzero[0],arrnonzero[1]):
        arr2d[locx,locy]=item
    #3D gaussian smooth to get rid of discontinuities
    #first make it a masked array so that the edges are preserved
    smoothed_arr=arr2d

    #re-inflate arr2D
    ashape=np.shape(arr)
    errarr=np.zeros([ashape[0],ashape[1],ashape[2]+dev])
    for locx,locy in zip(arrnonzero[0],arrnonzero[1]):
        lim=smoothed_arr[locx,locy]
        errarr[locx,locy,:lim]=1
    return rvec,errarr,arr2d

def zfunc(x,y,garr):
    if garr[x,y]==1:
        return 1
    else:
        return 0

def plot2D(arr):
    fig,ax=plt.subplots()
    #ax=fig.add_subplot(111,projection='3d')
    #X,Y=np.meshgrid(np.array(range(0,arrx)),np.array(range(0,arry)))
    #Z=garr[:,:,0]
    #print np.mean(Z)
    #ax.plot_trisurf(X,Y,Z)
    ax.imshow(np.transpose(arr),cmap=cm.binary,origin='bottom left')
    #ax.step(range(0,arrx),np.sum(garr,axis=1))
    #ax.set_ylim([-1,2])
    fig.show()

def plot1Dslice(arr3d,y):
    fig,ax=plt.subplots()
    #ax=fig.add_subplot(111,projection='3d')
    #X,Y=np.meshgrid(np.array(range(0,arrx)),np.array(range(0,arry)))
    #Z=garr[:,:,0]
    #print np.mean(Z)
    #ax.plot_trisurf(X,Y,Z)
    ax.plot(range(0,np.shape(arr3d)[0]),np.sum(arr3d[:,y,:],axis=1))
    #ax.step(range(0,arrx),np.sum(garr,axis=1))
    #ax.set_ylim([-1,2])
    fig.show()

def test_transm_plot():
    theta=np.array([-5,-4,-3,-2,-1,0,1,2,3,4,5])
    h=250.#250. #um
    pp=[88.,89.,90.,91.,92.] #um
    ww=[]
    for p in pp:
        ww.append((p/2.) - np.abs(h*np.tan(np.deg2rad(theta))))
    fig,ax=plt.subplots()
    for w,p in zip(ww,pp):
        ax.plot(theta,w/p,label=str(p))
    ax.legend()
    fig.show()


#if __name__ == "__main__":
    #param_dict={"period":15,"orientation":45,"phase":0,"height":250,"slit_width":7.5}
    #arrx,arry,garr=calc_3D_grating(param_dict)
    #from mpl_toolkits.mplot3d import Axes3D
    #fig=plt.figure()
    #fig,ax=plt.subplots()
    #ax=fig.add_subplot(111,projection='3d')
    #X,Y=np.meshgrid(np.array(range(0,arrx)),np.array(range(0,arry)))
    #Z=garr[:,:,0]
    #print np.mean(Z)
    #ax.plot_trisurf(X,Y,Z)
    #ax.imshow(np.transpose(garr[:,:,0]),cmap=cm.binary,origin='bottom left')
    #ax.step(range(0,arrx),np.sum(garr,axis=1))
    #ax.set_ylim([-1,2])
    #fig.show()

