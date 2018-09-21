"""
===================
moire_mosaic.py
Erica  Lastufka 12.4.18
===================

Read the logfiles, make the mosaic

"""
import numpy as np
import os
import glob
from PIL import Image
from PIL import ImageDraw
from skimage import feature, exposure
import pickle
from matplotlib import cm

def get_coordinates(logfile,tag='groupstretch'):
    '''Get coordinates of all given images in list filen. If there's a random z-coordinate stuck in the middle of the filename on the log readout, deal with it'''
    #get corresponding log file
    cdict={}
    with open(logfile,'rb') as f:
        lines=f.readlines()
        for i,line in enumerate(lines):
            if 'moireQM_step' in line:
                stepn=int(line[12:15])#int(line[10:13])
                win=int(line[22:25])#int(line[20:23])
            else:
                llist=line.split('\t')
                idx=int(llist[0])+1
                xmm=float(llist[1])
                ymm=float(llist[2])
                if tag:
                    fn='moireQM_step'+'{0:03d}'.format(stepn)+'_window0'+str(win)+'_'+'{0:04d}'.format(idx)+ '_'+tag +'.tif'
                else:
                    fn='moireQM_step'+'{0:03d}'.format(stepn)+'_window0'+str(win)+'_'+'{0:04d}'.format(idx) +'.tif'
                    
                cdict[fn]={'number':win,'step':stepn,'im_coords':[float(xmm),float(ymm)]}
    #return image coordinates
    return cdict

def make_mosaic(win, coords,pix2um=6.5,plot=True,plot_coords=True):
    '''Make a mosaic of the given window or all windows (individually). Input: list of window numbers, dictionary of coordinates'''            
    #store filenames
    fnames = [im for im in coords.keys() if int(coords[im]['number'])== win]
    fnames.sort(reverse=True)
        
    #get the coordinates of each image - x's don't matter here since they are all the same
    #cwinx = [coords[im]['im_coords'][0] for im in fnames] #image x coordinates in mm
    cwiny = [coords[im]['im_coords'][1] for im in fnames] #image y coordinates in mm
    #idxx  = [int(im[-6:-4]) for im in fnames]
    #idxx  = [int(im[-6:-4]) for im in fnames]

    #convert mm to pixels -here, image size is 330x1620
    #cpixx=np.array(cwinx)*(1620./(pix2um*1000.)) 
    cpixy=(np.array(cwiny)*1000.)/(pix2um)

    #translate coords so that top left coord is at [0,0]
    #cminx=np.min(cpixx)
    cmin=np.min(cpixy)
    #print cminx,cmaxy
    #cpixx=np.rint(cpixx-cminx)#cpixx-cminx
    cpixy=np.rint(cpixy-cmin)#cmaxy-cpixy
    #print int(np.max(cpixx))+640,int(np.max(cpixy))+480

    #make new blank image
    background=Image.new("F",[2040,1700], 0xff) #[0,0] is TOP left #extra 9x12 px for rounding errors
    print np.shape(background)
    mx,my,fpixx,fpixy=[],[],[],[]
    #print 'excess ', exx,exy
    for i,f in enumerate(fnames):#zip([2,6],fnames):#enumerate(fnames):
        im1=Image.open(f)#.convert("RGB")
        print np.min(im1),np.max(im1),np.mean(im1)
        
        #trim top and bottom 10 px and correct cpixy
        #im2=im1.crop([0,5,1620,250])
        #print np.min(im2),np.max(im2),np.mean(im2)
        #print cpixy[i]
        #cpixy[i]=cpixy[i]-10
        #print cpixy[i]
        
        #im2=Image.open(fnames[i+1]).convert("RGBA")
        #yoffset=int(cpixy[i]-cpixy[i+1]) #20 #trim off the top and bottom 20 pix because they are dark
        #cropped_im=im.crop((0,yoffset,1620,330-yoffset))
        #cropped_im.show()

        #add alpha channel for the top/bottom 20 pix
        #transparent_bottom=(0,0,1620,100)
        #transparent_top=(0,300,1620,300)
        #mask=Image.new('L',im.size,color=255)
        #draw=ImageDraw.Draw(mask)
        #draw.rectangle(transparent_bottom,fill=0)
        #draw.rectangle(transparent_top,fill=0)
        #im.putalpha(mask)
        #im.convert("RGBA")
        #print f,cpixy[i],im.mode
        # "The box argument is either a 2-tuple giving the upper left corner,
        # a 4-tuple defining the left, upper, right, and lower pixel coordinate, or None (same as (0, 0))."
        #assuming the actual coordinates are the center of the image...
        #offsetx = 640/2
        #offsety= 480/2
        #yoffset=20 #trim off the top and bottom 20 pix because they are dark
        box1=(0,int(cpixy[i]))#,640+int(cpixx[i]),480+int(cpixy[i]))
        #box2=(0,int(cpixy[i+1]))#,640+int(cpixx[i]),480+int(cpixy[i]))
        #newimg1=Image.new('RGBA',size=(1620,330+yoffset),color=(0,0,0,0))
        #newimg1.paste(im2,(0,yoffset))
        #newimg1.paste(im1,(0,0))
        #newimg2=Image.new('RGBA',size=(1620,330+yoffset),color=(0,0,0,0))
        #newimg2.paste(im1,(0,0))
        #newimg2.paste(im2,(0,yoffset))
        #blended_im=Image.blend(newimg1,newimg2,alpha=0.5)
        
        background.paste(im1,box1)

    if plot: #show the mosaic
        #newimg1.show()
        #newimg2.show()
        #blended_im.show()
        print np.mean(background),np.min(background),np.max(background)
        background.show()
        #im.putalpha(mask)
        #im.show()
    if plot_coords: #just plot the coords
        fig,ax=plt.subplots()
        #ax.scatter(cwinx,cwiny)
        ax.scatter(np.zeros(len(cpixy)),cpixy)
        fig.show()
            
    #now convert the mosaic to numpy array and save. Also save a .tiff image
    #marray=np.array(background.convert('L')) #raw, unequalised array
    filename='window'+str(win)+'mosaic'
    savec=scipy.misc.toimage(background,high=np.max(background),low=np.min(background),mode='F')
    savec.save(filename+'.tiff')
    #pickle.dump(marray,open(filename+'.p','wb'))
    #pickle.dump(mcoords, open(filename+'_coords.p','wb'))
        
    #return new list of filenames
    #return fnames,cpixx,cpixy,idxx,idxy,fpixx,fpixy#filename+'.p'
    return fnames,cpixy,background
