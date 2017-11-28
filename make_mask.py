"""
===================
make_mask.py
Erica  Lastufka 17.11.17
===================

make mask, do convolutions

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from PIL import Image
from matplotlib import cm
import pickle
import time
import glob
from PIL import ImageChops
import multiprocessing as mpi
from scipy.misc import imrotate
    
def make_binary_mask(period,angle,pix2um=1.955,phase=0,dcycle=0.5,mosaic=True,plot=False):
    '''Given a period and orientation, make a binary mask of the ideal grating (edges?)'''
    if mosaic:
        msize=[8288,8288]
        #cropsize=[5756,5754]
        cropsize=[5656,5654]
    else:
        msize=[905,905]
        cropsize=[479,639]
    bmask=np.zeros(msize)
    #starting from pixel phase, put 1's where the slats should be...
    ppix=period/pix2um
    nedges=(2*msize[1]/ppix)/256.
    eloc=np.rint(ppix*dcycle)
    for y in range(0,msize[1]-1):
        ymod=np.mod(y-phase,np.rint(ppix))
        if y >= phase:
            if ymod == 0: #or ymod== 0:
                bmask[:,y]=1
            elif ymod== eloc: #or ymod ==0:
                bmask[:,y]=-1

    rbmask=imrotate(bmask,-1.*angle) #again because of indexing
    #crop
    #exy=(msize[0]-cropsize[0])/2 #because ndarray indexing
    #exx=(msize[1]-cropsize[1])/2
    #print exx,exy
    #rbmask=rbmask[exy:-exy,exx:-exx]
    #thin edges and make binary
    #minidx=rbmask < 100#np.max(rbmask)
    #maxidx=rbmask >= 100#np.max(rbmask)
    #rbmask[minidx] = 0
    #rbmask[maxidx] = 1
    
    if plot:
        fig,ax=plt.subplots(2)
        #ax.imshow(imrotate(bmask,angle),cmap=cm.binary)
        ax[0].imshow(bmask,cmap=cm.binary)
        ax[1].imshow(rbmask,cmap=cm.binary)
        fig.show()
    
    return rbmask

def make_index_mask(period,angle,um2pix=0.511501,phase=0,dcycle=0.5,mosaic=True,plot=False, rising=True):
    '''Given a period and orientation, make a mask of the ideal grating edges, indexing each. Rising is positive, falling negative'''
    #add width paramater to make the edges that many pixels wide... easier to index later on that way
    if mosaic:
        msize=[8288,8288]
        #cropsize=[5756,5754]
        cropsize=[5657,5654]
    else:
        msize=[905,905]
        cropsize=[479,639]
    rmask,fmask=np.zeros(msize),np.zeros(msize)
    #starting from pixel phase, put 1's where the slats should be...
    if not rising: # assume the first edge is rising, else not
        factor=-1.
    else:
        factor = 1.
    idx=0
    ppix=period*um2pix
    nedges=(2*msize[1]/ppix)/256.
    eloc=np.rint(ppix*dcycle)
    for y in range(0,msize[1]-1):
        ymod=np.mod(y-phase,np.rint(ppix))
        if y >= phase:
            if ymod == 0: #or ymod== 0:
                rmask[:,y]=factor*float(idx)/float(nedges)
                #print factor*float(idx)/float(nedges)
                idx+=1
            elif ymod== eloc: #or ymod ==0:
                fmask[:,y]=-1*factor*float(idx)/float(nedges)
                idx+=1
                #print -1*factor*float(idx)/float(nedges)
    #print np.max(rmask),np.min(fmask),idx
    rrmask=imrotate(rmask,-1*angle,interp='nearest') 
    rfmask=-1*imrotate(fmask,-1*angle,interp='nearest') #again because of indexing
    #crop
    exy=(msize[0]-cropsize[0])/2 #because ndarray indexing
    exx=(msize[1]-cropsize[1])/2
    #print exx,exy
    rrmask=rrmask[exy:-exy,exx:-exx]
    rfmask=rfmask[exy:-exy,exx:-exx]
    #get rid of interpolation of zeros
    rmedidx=np.where(rrmask ==np.median(rrmask))
    fmedidx=np.where(rfmask ==np.median(rfmask))
    rrmask[rmedidx]=0
    rfmask[fmedidx]=0
    #now re-adjust so that indexing of edges proceeds linearly....
    rrmax=np.max(rrmask)
    rrmin=np.min(rrmask[np.where(rrmask != 0)])
    rfmax=np.min(rfmask)
    rfmin=np.max(rfmask[np.where(rfmask != 0)])
    #print rfmax,rfmin,np.mean(rfmask)
    rrstep=(rrmax-rrmin)/(nedges*128.)
    rfstep=(rfmax-rfmin)/(nedges*128.)
    if rrstep <=2:
        rrmask=rrmask*(2./rrstep)
        rrmax=np.max(rrmask)
        rfmin=np.max(rfmask[np.where(rfmask != 0)])
        rrstep=(rrmax-rrmin)/(nedges*128.)
    if rfstep <=2:
        rfmask=rfmask*(2./np.abs(rfstep))
        #print np.max(rfmask),np.min(rfmask),np.mean(rfmask)
        rfmax=np.min(rfmask)
        rfmin=np.max(rfmask[np.where(rfmask != 0)])
        rfstep=(rfmax-rfmin)/(nedges*128.)

    print rrmin,rrmax,rrstep
    print rfmin,rfmax,rfstep
    
    for i,j in zip(np.arange(rrmin,rrmax+rrstep,rrstep),np.arange(0,int(nedges*128.)+2,1)): 
        rrmask[np.where(np.logical_and(rrmask <= int(i)+int(rrstep/2.),rrmask > int(i)-int(rrstep/2.)))] = j
        #print int(i)+int(rrstep/2.),int(i)-int(rrstep/2.),j, np.where(np.logical_and(rrmask < int(i)+int(rrstep/2.),rrmask > int(i)-int(rrstep/2.)))
    for i,j in zip(np.arange(-1*rfmax,-1*(rfmin+rfstep),rfstep),-1*np.arange(0,int(nedges*128.)+1,1)): 
        rfmask[np.where(np.logical_and(rfmask >= -1*(int(i)-int(rfstep/2.)),rfmask <= -1*(int(i)+int(rfstep/2.))))] = j
        #print -1*int(i)+int(rfstep/2.),-1*int(i)-int(rfstep/2.),j#,np.where(np.logical_and(rfmask > -1*(int(i)-int(rfstep/2.)),rfmask <= -1*(int(i)+int(rfstep/2.))))
    #print rrmin,rrmax,rrstep
    #print rfmin,rfmax,rfstep
    
    
    if plot:
        fig,ax=plt.subplots(2)
        #ax.imshow(imrotate(bmask,angle),cmap=cm.binary)
        f1=ax[0].imshow(rmask+fmask)
        cbar=fig.colorbar(f1)
        f2=ax[1].imshow(np.transpose(rrmask+rfmask))
        cbar=fig.colorbar(f2)
        fig.show()

    print np.mean(rrmask),np.mean(rfmask)
    return rmask,fmask,rrmask,rfmask,nedges

def chop_index_mask(mask,coords,filen,peek=False):
    '''Slice up full mosaic index mask into pieces corresponding to individual images. Coordinates are top left corner of image. Save index masks.'''
    #get mosaic coordinates
    if type(coords) == str:
        coords=pickle.load(open(coords,'rb'))
    imsize=[640,480]

    for f,x,y in zip(filen,coords[0],coords[1]):
        m=mask[x:x+imsize[0],y-imsize[1]:y]
        filename=f[:-4]+'_imask.p'
        pickle.dump(m,open(filename,'wb'))

    if peek:
        mask_peek(f,m)

    return m

def mask_peek(f,m,pix2um=1.955,length=.2):
    import edge_and_hough as eh
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib_scalebar.scalebar import SI_LENGTH
   
    if type(m) == str:
        m=pickle.load(open(m,'rb'))
    fig,ax=plt.subplots()
    ax.imshow(eh.im2ndarray(f),alpha=.6)
    ax.imshow(np.transpose(m),alpha=.6)
    scalebar = ScaleBar(pix2um,'um', SI_LENGTH,length_fraction=length)
    ax.add_artist(scalebar)
    #ax[1].add_artist(scalebar)
    fig.show()
    print np.shape(m)
    
def compare_with_edges(edges,mask,plot=False):
    '''Compare array of Canny edges with the rotated mask'''
    #first check that the arrays are the same size
    if np.shape(edges) != np.shape(mask):
        print np.shape(edges)
        print np.shape(mask)
        print "Array sizes don't match! Exiting... "
        return 'foo'
    result=edges*mask
    if plot:
        fig,ax=plt.subplots(2,sharex=True,sharey=True)
        ax[0].imshow(edges,cmap=cm.binary)
        ax[0].imshow(mask,cmap=cm.Reds,alpha=0.6)
        #ax.imshow(imrotate(bmask,angle),cmap=cm.binary)
        ax[1].imshow(result,cmap=cm.binary)
        fig.show()
        
    return np.sum(result)

def fit_window(window_num,nperiod=False,nangle=False,prange=[-2.,2.],arange=[-.1,.1],pstep=.1,astep=0.05,plot_hist=True):
    ''' Run fit_mask for each image in window. make histogram of periods, angles and duty cycles'''
    import edge_and_hough as eh #should bring in the global windows
    fnames = glob.glob('win'+str(window_num)+'_*5.0X_edges.p')
    if not nperiod:
        nperiod=nang=[w['pitch'] for w in eh.windows if w['number'] == window_num][0]
    if not nangle:
        nangle=nang=[w['nominal angle'] for w in eh.windows if w['number'] == window_num][0]
    stats=[]
    for f in fnames:
        edges=pickle.load(open(f,'rb'))[1:,1:]
        stats.append([f,fit_mask(edges, nperiod,nangle,prange=prange,arange=arange,pstep=pstep,astep=astep,plot=False)])

    period=[w[0] for s in stats for w in s[1][1]]  #[w[0] for s in ss for w in s[2]]
    angle=[w[1] for s in stats for w in s[1][1]]
    means=[s[1][0] for s in stats]
    
    if plot_hist:
        fig,ax=plt.subplots(2)
        ax[0].hist(period,np.arange(nperiod+prange[0],nperiod+prange[1],pstep)) #don't log scale there's only 108 points.... 
        ax[0].set_xlim([nperiod+prange[0],nperiod+prange[1]])
        ax[1].hist(angle,np.arange(nangle+arange[0],nangle+arange[1],astep))
        ax[1].set_xlim([nangle+arange[0],nangle+arange[1]])
        fig.show()

    print np.mean(means)
    return period,angle

def fit_mask(edges,nperiod,nangle,prange=[-2.,.5],arange=[-.1,.1],pstep=.1,astep=.05,phaserange=[0,10],plot=True,mosaic=False,dcycle=False):
    '''Find the best mask to fit the data. Do this on a small image first and then fit the phase on the mosaic!'''
    means,params=[],[]
    for p in np.arange(prange[0],prange[1],pstep):
        for a in np.arange(arange[0],arange[1],astep):
            if dcycle:
                for f in np.arange(phaserange[0],phaserange[1],1):
                    for d in np.arange(0.3,0.7,.1):
                        cmask=make_binary_mask(nperiod+p,nangle+a,phase=f,dcycle=d,mosaic=mosaic)
                        means.append(compare_with_edges(edges,cmask))
                        params.append([nperiod+p,nangle+a,f,d])                        
            else:
                for f in np.arange(phaserange[0],phaserange[1],1):
                    cmask=make_binary_mask(nperiod+p,nangle+a,phase=f,mosaic=mosaic)
                    means.append(compare_with_edges(edges,cmask))
                    params.append([nperiod+p,nangle+a,f])
                
    idx=np.where(means == np.max(means))[0]#[0] #should probabl return all of these...
    print np.max(means),params[idx[0]]
    if plot:
        dmask=make_binary_mask(params[idx[0]][0],params[idx[0]][1],phase=params[idx[0]][2],mosaic=mosaic)
        compare_with_edges(edges,dmask,plot=True)
        print np.max(means)/np.sum(dmask)
    return np.max(means),[params[i] for i in idx]
        
