"""
===================
smart_mosaic.py
Erica  Lastufka 19.12.17
===================

Make a mosaic where the edges line up to within 1 pixel, starting with the original coordinates. Based off mosaic() in edge_and_hough.py

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
    
def get_coordinates(path,filen,zcoord=False):
    '''Get coordinates of all given images in list filen. If there's a random z-coordinate stuck in the middle of the filename on the log readout, deal with it'''
    #get corresponding log file
    os.chdir(path)
    logs=glob.glob('*.txt')
    #print logs
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

def pad(earray, indices,ofname, plot=False):
    '''Pad the array with 0's if it had its edges trimmed'''
    idxx,idxy=indices[0],indices[1]
    eshape=np.shape(earray)
    ishape=np.zeros([480,640])
    #print np.shape(ishape)
    if np.shape(ishape)[1] != eshape[1] and np.shape(ishape)[0] != eshape[0]: #trimmed from both
        if idxx < 5 and idxy <5: #it was trimmed from the top left
            ishape[480-eshape[0]:,640-eshape[1]:] = earray
            #print np.shape(ishape[:480-eshape[0],640-eshape[1]:])
        else:
            ishape[:-(480-eshape[0]),:-(640-eshape[1])] = earray
    elif np.shape(ishape)[1] != eshape[1]:
        if idxx < 5: #it was trimmed from the left
            ishape[:,640-eshape[1]:] = earray
            print np.shape(ishape[:,640-eshape[1]:])
        else:
            ishape[:,:-(640-eshape[1])] = earray
    elif np.shape(ishape)[0] != eshape[0]:
        #if idxy < 5: # it was trimmed from the top
            #ishape[480-eshape[0]:,:] = earray
        ishape[:-(480-eshape[0]),:] = earray
        #else: # it was trimmed from the bottom
            #ishape[:-(480-eshape[0]),:] = earray
        #    ishape[480-eshape[0]:,:] = earray

    if plot: #check
        fig,ax=plt.subplots()
        ax.imshow(ishape,cmap=cm.binary)
        fig.show()
    #re-pickle padded array
    pickle.dump(ishape, open(ofname,'wb'))
    
def smart_align(coords1,coords2,p1,p2,order='t',tol=10, fixed_c=False, plot=False):
    '''align images to each other within a certain tolerance. Return the relative coordinates in pixels of the aligned images'''
    #print coords2
    oldscore1,oldscore2=0.,0.
    bestcoords1=([0,0],[0,0])
    if not fixed_c: #form: [id, [cx,cy]]
        permutations = [[coords1[0]+a,coords1[1]+b,coords2[0]+c,coords2[1]+d] for a in range(0,tol) for b in range(0,tol) for c in range(0,tol) for d in range(0,tol)] #make the list of permutations
        bestcoords1=(coords1,coords2)
    elif fixed_c[0] == 0: #p1 is fixed
        dx=fixed_c[1][0]-coords2[0]
        permutations = [[fixed_c[1][0],fixed_c[1][1],coords2[0]+c+dx,coords2[1]+d] for c in range(0,tol) for d in range(0,tol)]
        bestcoords1=(fixed_c[1],coords2)
        coords1=fixed_c[1]
    else:
        dx=fixed_c[1][0]-coords2[0]
        permutations = [[coords1[0]+a+dx,coords1[1]+b,fixed_c[1][0],fixed_c[1][1]] for a in range(0,tol) for b in range(0,tol)]
        bestcoords1=(coords1,fixed_coords[1])
        coords2=fixed_coords[1]
        
    #print 'dy' ,(480-abs(coords1[1]-coords2[1]))
    #print 'permutations ' + str(len(permutations))
    for x1,y1,x2,y2 in permutations:
    #get adjacent rows to compare
        dy = (480-abs(y1-y2))
        if order == 't': #p1 is on top of p2
            arow1=p1[-dy:,:]
            #print np.sum(arow1),np.shape(arow1)
            arow2=p2[:dy,:]
        elif order == 'l': 
            arow1=p1[:,590:]
            arow2=p2[:,:50]
        score1=align([x1,y1], [x2,y2],arow1,arow2,order=order) #different if fixed_coords....
        #print score1, [[coords1[0]+x1,coords1[1]+y1],[coords2[0]+x2,coords2[1]+y2]]
        if score1 > oldscore1:
            oldscore1=score1
            bestcoords1=([x1,y1],[x2,y2])
            #print dy,score1,bestcoords1#,([coords1[0]+x1,coords1[1]+y1],[coords2[0]+x2,coords2[1]+y2])

    #account for any initial conditions 
    if plot:
        im1=Image.open(im2n)
        im2=Image.open(im1n)
        #fig,ax=plt.subplots()
        background=Image.new("L",[680,5737], 0xff) #[0,0] is TOP left #extra 9x12 px for rounding errors
        mx,my=[],[]
        box1=(int(bestcoords1[0][0]),int(bestcoords1[0][1]))
        box2=(int(bestcoords1[1][0]),int(bestcoords1[1][1]))
        #mx.append(box1[0])
        #my.append(box1[1])
        #mx.append(box2[0])
        #my.append(box2[1])
        background.paste(im2,box2) 
        background.paste(im1,box1)
        #background.paste(im2,box2) 
        background.show()
    #print 'score: ',oldscore1
    #print 'coords: ',bestcoords1
    offset1=bestcoords1[0]
    offset2=bestcoords1[1]
    #aligned_coords1=[coords1[0]+bestcoords1[0][0],coords1[1]+bestcoords1[0][1]]
    #aligned_coords2=[coords2[0]+bestcoords1[1][0],coords2[1]+bestcoords1[1][1]]
    return offset1,offset2#oldscore1, #aligned_coords1,aligned_coords2,box1,box2#,cpixx,cpixy#aligned_pixx,aligned_pixy

def align(c1,c2,a1,a2, order = 't'):
    '''return the quality of the alignment'''
    #if order == 't' or 'b': #x coordinate doesn't need an offset
    #keep first imaged 'fixed' and move 2nd image relative to its corodinates
    dx=c1[0]-c2[0]
    dy=((480-abs(c1[1]-c2[1]))) #difference between y-coordinates
    #print dy
    #rmatrix,smatrix=np.transpose(np.zeros([abs(dx)+640,abs(dy)+50])),np.transpose(np.zeros([abs(dx)+640,abs(dy)+50]))
    #mshape=np.shape(rmatrix)
    if order == 't' or order == 'b':
        xlim,ylim=640,dy
        rmatrix,smatrix=np.transpose(np.zeros([abs(dx)+640,abs(dy)])),np.transpose(np.zeros([abs(dx)+640,abs(dy)]))
    else:
        xlim,ylim=480,82
        rmatrix,smatrix=np.transpose(np.zeros([abs(dx)+50,abs(dy)+480])),np.transpose(np.zeros([abs(dx)+50,abs(dy)+480]))
        
    if dx == 0:
        x1i,x1j,x2i,x2j=0,xlim,0,xlim
    elif dx >0:
        x1i,x1j,x2i,x2j=dx,xlim+dx,0,xlim
    else:
        x1i,x1j,x2i,x2j=0,xlim,abs(dx),xlim+abs(dx) #basically treat 2nd image as reference in this case

    y1i,y1j,y2i,y2j=0,ylim,0,ylim

    #ydiff = 50-dy
    #if ydiff < 0:
    #    ydiff=dy-50
    #else:
    #    y1i,y1j,y2i,y2j=0,ylim,abs(dy),ylim+abs(dy)

    rmatrix[:,x1i:x1j] = a1
    #if ydiff != 0:
    #    smatrix[y2i:y2j,x2i:x2j] = a2[:-ydiff,:]
    #else:
    smatrix[:,x2i:x2j] = a2
        
    #print np.shape(rmatrix),np.shape(smatrix)
    #print rmatrix[-20:-10:,-18:],np.sum(rmatrix[50:110,:])
    #print smatrix[-20:-10:,-18:], np.sum(smatrix[50:110,:])
    #print np.sum(rmatrix),np.sum(smatrix)
    qmatrix=rmatrix*smatrix #now multiply by the second matrix in the correct position and return the sum
    q = np.sum(qmatrix)
    return q

#def align_shifts(ap_evens,ap_odds):
#    '''go from relative shifts to absolute'''
#    #if there's one with 0 shift, use that as the reference

#    #otherwise use the first one

    

def smart_column(win, colnum, filen = False, plot=True, mag=5.0,tol=10,path='.',zcoord=False):
    '''Make a mosaic of the given window or all windows (individually). Input: list of window numbers, dictionary of coordinates'''
    #how to deal with overlap? make left, right, top, bottom options and compare statistically effect on edge/line fits?
    if not filen:
        filen=glob.glob('win'+str(win)+'_0'+str(colnum)+'*'+str(mag)+'X.tif')
    #filen.sort()
    #print len(filen),filen
    coords=get_coordinates(path,filen,zcoord=zcoord)
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
        window_num=[w['number'] for w in windows]

    #for win in window_num:
        #store filenames
    fnames = [im for im in coords.keys() if int(coords[im]['number'])== win]
    pnames = [im[:-4]+'_edges.p' for im in coords.keys() if int(coords[im]['number'])== win]
    fnames.sort()
    pnames.sort()
    #pnames = [im[:-4]+'_edges.p' for im in fnames if coords[im]['indices'][0]== '01'] #for testing on first column
        
    #get the coordinates of each image
    cwinx = [coords[im]['im_coords'][0] for im in fnames] #image x coordinates in mm
    cwiny = [coords[im]['im_coords'][1] for im in fnames] #image y coordinates in mm
    idxx  = [coords[im]['indices'][1] for im in fnames]
    idxy  = [coords[im]['indices'][0] for im in fnames]
        #convert mm to pixels
    cpixx=np.array(cwinx)*(640.0/pixX) 
    cpixy=np.array(cwiny)*(480.0/pixY)

    #translate coords so that top left coord is at [0,0]
    cminx=np.min(cpixx)
    cmaxy=np.max(cpixy)
    #print cminx,cmaxy
    cpixx=cpixx-cminx#np.rint(cpixx-cminx)
    cpixy=cmaxy-cpixy#np.rint(cmaxy-cpixy)
    #print int(np.max(cpixx))+640,int(np.max(cpixy))+480

    #pad the arrays if necessary
    for i,p in enumerate(pnames):
        edgearr=pickle.load(open(p,'rb'))
        if np.shape(edgearr) != (480,640):
            pad(edgearr,[idxx,idxy],p)

    pair=[[p,pnames[i+1]] for i,p in enumerate(pnames[:-1])]
    cpairs=[[[int(cpx),int(cpy)],[int(cpixx[i+1]),int(cpixy[i+1])]] for i,cpx,cpy in zip(range(0,len(cpixx)-1),cpixx,cpixy)]
    newcoords={}
    newcoords[pair[0][0]]=[0,0]
    
    for pp,cp in zip(pair,cpairs):
        p1=pickle.load(open(pp[0],'rb'))
        p2=pickle.load(open(pp[1],'rb'))
        order='t'
        if pp[0] in newcoords.keys():
            fixed_c = [0,[newcoords[pp[0]][0],newcoords[pp[0]][1]]]
            #print 'fixed_c for ' + pp[0], fixed_c
        else:
            fixed_c = False
        #print 'aligning '+pp[0]+' and '+pp[1]
        aligned_pix1,aligned_pix2=smart_align(cp[0],cp[1],p1,p2,tol=tol,order=order,fixed_c = fixed_c) #do the alignment
        if pp[0] not in newcoords.keys():
            newcoords[pp[0]]=aligned_pix1#+cp[0]
        if pp[1] not in newcoords.keys():
            newcoords[pp[1]]=aligned_pix2#+cp[1]

    aligned_pixx=[newcoords[k][0] for k in sorted(newcoords.keys())]
    aligned_pixy=[newcoords[k][1] for k in sorted(newcoords.keys())]
        
    if plot:
        #make new blank image
        background=Image.new("L",[640+tol*8,5737], 0xff) #[0,0] is TOP left #extra 9x12 px for rounding errors
        for i,f in enumerate(fnames):#zip([2,6],fnames):#enumerate(fnames):
            im=Image.open(f) #does this work for pickles?
            box=(int(aligned_pixx[i]),int(aligned_pixy[i]))
            background.paste(im,box)
        background.show()
                    
    #now convert the mosaic to numpy array and save. Also save a .tiff image
    #marray=np.array(background.convert('L')) #raw, unequalised array
    #filename='window'+str(win)+'mosaic_'+str(mag)
    #background.save(filename+'.tiff')
    #pickle.dump(marray,open(filename+'.p','wb'))
    #pickle.dump(mcoords, open(filename+'_coords.p','wb'))
    
    return newcoords#,aligned_pixx,aligned_pixy#filename+'.p'

def smart_mosaic(win, tol=10,mag=5.0):
    '''Sew the aligned columns together, using FIXED dx'''
    fnames=glob.glob('win'+str(win)+'*'+str(mag)+'X.tif')
    fnames.sort()
    cols = range(1,10)
    ncs=[]

    xcdict,ycdict={},{}
    fac=630
    dy=10
    #adjust the dictionary accordingly

    for c in cols:
        colcoords=smart_column(win, c,plot=False)
        for k in colcoords.keys():
            newkey=k[:-8]+'.tif'
            ycdict[newkey]=colcoords[k][1]+dy*(9-c)
            xcdict[newkey]=colcoords[k][0]+c*fac
        
    background=Image.new("L",[8000,8000], 0xff) #[0,0] is TOP left #extra 9x12 px for rounding errors
    for f in fnames:#zip([2,6],fnames):#enumerate(fnames):
        im=Image.open(f) #does this work for pickles?
        box=((xcdict[f],ycdict[f]))
        background.paste(im,box)
    background.show()

    return fnames,xcdict,ycdict
 
