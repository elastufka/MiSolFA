"""
===================
make_mosaic.py
Erica  Lastufka 16.11.17
===================

make mosaic

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
    index=[indices[:2],indices[3:]]
    return windex,index

def mosaic(window_num, coords=False, all=False, plot=True,plot_coords=False, mag=5.0,zcoord=False,offsetx=0,offsety=0):
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
            'FOVx': 0.4163,
            'FOVy' : 0.3127,
            'ExcessX_px' : -7.0,
            'ExcessY_px' : 6.0,
            'ExcessX' : -0.0046,
            'ExcessY' : 0.0039,
            'StepX' : 0.4163,
            'StepY' : 0.3127}
        

    pixX= mscope['FOVx']#1.955 #this is too big....
    pixY=mscope['FOVy']#1.9548 #from Matej's logs since mine don't have them
    
    if all:
        #generate list of window numbers from dictionary 
        window_num=[11,21,12,22,31,41,32,42,33,43,34,44]#[w['number'] for w in windows]

    for win in window_num:
        if not coords:
            filen=glob.glob('win'+str(win)+'*'+str(mag)+'X.tif')
            coords=get_coordinates('.',filen,zcoord=zcoord)
        #store filenames
        fnames = [im for im in coords.keys() if int(coords[im]['number'])== win]
        #fnames = [im for im in fnames if coords[im]['indices'][0]== '01'] #for testing on first column
        
        #get the coordinates of each image
        cwinx = [coords[im]['im_coords'][0] for im in fnames] #image x coordinates in mm
        cwiny = [coords[im]['im_coords'][1] for im in fnames] #image y coordinates in mm
        idxx  = [coords[im]['indices'][1] for im in fnames]
        idxy  = [coords[im]['indices'][0] for im in fnames]
        #convert mm to pixels
        cpixx=np.array(cwinx)*(640.0/pixX) 
        cpixy=np.array(cwiny)*(480.0/pixY)
        coords=False
        #translate coords so that top left coord is at [0,0]
        cminx=np.min(cpixx)
        cmaxy=np.max(cpixy)
        #print cminx,cmaxy
        cpixx=cpixx-cminx#np.rint(cpixx-cminx)
        cpixy=cmaxy-cpixy#np.rint(cmaxy-cpixy)
        #print int(np.max(cpixx))+640,int(np.max(cpixy))+480

        #make new blank image
        background=Image.new("L",[int(np.max(cpixx))+714,int(np.max(cpixy))+554], 0xff) #[0,0] is TOP left #extra 9x12 px for rounding errors
        #background=Image.new("1",[int(np.max(cpixy))+480,int(np.max(cpixx))+640], "white") #[0,0] is TOP left
        print 'size:', int(np.max(cpixx))+649,int(np.max(cpixy))+492, 
        #put things in their place
        #fnames=[fnames[2],fnames[6]]
        for i,f in enumerate(fnames):#zip([2,6],fnames):#enumerate(fnames):
            im=Image.open(f)
            #offsetx=int(idxx[i])*6
            #offsety=int(idxy[i])*9
            box=(int(cpixx[i])+offsetx,int(cpixy[i])-offsety)#,640+int(cpixx[i]),480+int(cpixy[i]))
            background.paste(im,box)

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
        
    #return new list of filenames
    #return fnames,cpixx,cpixy#filename+'.p'
    #return mfiles
    #return fnames,cpixx,cpixy

def align(mosaic1,mosaic2,xran=[0,50],yran=[0,50]):
    if type(mosaic1) == str:
        #    mosaic1=pickle.load(open(mosaic1,'rb'))
        mosaic1=Image.open(mosaic1)
    if type(mosaic2) == str: #assume it's a pickle
        #    mosaic2=pickle.load(open(mosaic2,'rb'))
        mosaic2=Image.open(mosaic2)

    means,results=[],[]#np.empty([3])
    start=time.time()
    #pool=mpi.Pool()
    print range(xran[0],xran[1]),range(yran[0],yran[1])
    for x in range(xran[0],xran[1]): #should mpi
        for y in range(yran[0],yran[1]):
            #dargs=[mosaic1,mosaic2,x,y]
            #results.append(pool.apply_async(difference_image,args=(dargs,)))
            #means=np.append(means,result.get(),axis=0) #does this work?
            di=difference_image(mosaic1,mosaic2,x,y)
            means.append(di) #does this work?
            #print di
            #print result.get()

    #for r in results:
    #    print r.get(timeout=30)
    #    means.append(r.get())
    #find minimum
    xx=[means[i][0] for i in range(0,len(means)-1)]
    yy=[means[i][1] for i in range(0,len(means)-1)]
    mm=[means[i][2] for i in range(0,len(means)-1)]
    idx=np.where(mm == np.min(mm))[0][0]
    print np.min(mm),xx[idx],yy[idx]
    print 'Processing took %.2f seconds' % (time.time() - start)
    difference_image(mosaic1,mosaic2,xx[idx],yy[idx],plot=True)    

def difference_image(mosaic1,mosaic2,xalign,yalign,rotate=0.,plot=False):
    '''Make a difference image of 2 mosaics'''
    #mosaic1=dargs[0]
    #mosaic2=dargs[1]
    #xalign=dargs[2]
    #yalign=dargs[3]
    if type(mosaic1) == str:
        #    mosaic1=pickle.load(open(mosaic1,'rb'))
            mosaic1=Image.open(mosaic1)
    if type(mosaic2) == str: #assume it's a pickle
        #    mosaic2=pickle.load(open(mosaic2,'rb'))
        mosaic2=Image.open(mosaic2)

    mosaic2=ImageChops.offset(mosaic2,xalign,yalign)
    difference=ImageChops.difference(mosaic1,mosaic2)
    if plot:
        fig,ax=plt.subplots()
        cax=ax.imshow(difference)
        cbar=fig.colorbar(cax)
        fig.show()

    mdiff=np.mean(difference.crop([500,500,5195,5195]))
    return [xalign,yalign,mdiff]
