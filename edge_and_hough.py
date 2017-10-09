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
import pickle
from scipy.misc import imrotate
import time

#list of dictionaries of the window numbers, pitches (mm), and nominal angles (deg)
global windows

windows=[{'number':11,'pitch':0.09,  'nominal angle':-44.79357},
                {'number':21,'pitch':0.09,  'nominal angle': 45.20793},
                {'number':12,'pitch':0.0225,'nominal angle': 44.94825},
                {'number':22,'pitch':0.0225,'nominal angle':-45.05184},
                {'number':31,'pitch':0.045, 'nominal angle':-45.10378},
                {'number':41,'pitch':0.045, 'nominal angle': 44.89660},
                {'number':32,'pitch':0.018, 'nominal angle': 45.04146},
                {'number':42,'pitch':0.018, 'nominal angle':-44.95859},
                {'number':33,'pitch':0.03,  'nominal angle':-44.93102},
                {'number':43,'pitch':0.03,  'nominal angle': 45.06914},
                {'number':34,'pitch':0.015, 'nominal angle': 44.96549},
                {'number':44,'pitch':0.015, 'nominal angle':-45.03455}]

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

def im2ndarray(filen):
    '''Convert .tif image to numpy.ndarray'''
    imraw = Image.open(filen) #.rotate(45)#.convert('1').convert('L')
    im=np.array(imraw)
    return im

def remove_edges(filen,im):
    '''remove edges from images with edges of the window in them'''
    if filen.endswith('.p') or filen.endswith('mosaic.tiff'):
        imshape=np.shape(im)
        im=im[50:imshape[0]-50,50:imshape[1]-50]
        return im
    else:
        windex,imindex=get_index(filen)
        #print imindex[0],imindex[1]
        if imindex[0] == '01': #it's got an edge on the left, trim the array (x and y are flipped for some reason)
            im[:,:50]=0#im[:,50:]
        if imindex[0] == '09': #it's got an edge on the right, trim the array
            im[:,590:]=0#im[:,:590]
        if imindex[1] == '01': #it's got an edge on the top, so mask the edges
            im[:50,:]=0#im[50:,:]
        if imindex[1] == '12': #it's got an edge on the bottom, so mask the edges
            im[430:,:]=0#im[:430,:]
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
            'FOVx': 0.4163,
            'FOVy' : 0.3127,
            'ExcessX_px' : -7.0,
            'ExcessY_px' : 6.0,
            'ExcessX' : -0.0046,
            'ExcessY' : 0.0039,
            'StepX' : 0.4163,
            'StepY' : 0.3127}
        

    pixX=mscope['FOVx']#1.955 #this is too big....
    pixY=mscope['FOVy']#1.9548 #from Matej's logs since mine don't have them
    
    if all:
        #generate list of window numbers from dictionary 
        window_num=[w['number'] for w in windows]

    for win in window_num:
        #store filenames
        fnames = [im for im in coords.keys() if int(coords[im]['number'])== win]
        #fnames = [im for im in fnames if coords[im]['indices'][0]== '01'] #for testing on first column
        
        #get the coordinates of each image
        cwinx = [coords[im]['im_coords'][0] for im in fnames] #image x coordinates in mm
        cwiny = [coords[im]['im_coords'][1] for im in fnames] #image y coordinates in mm
        
        #convert to pixels
        cpixx=np.array(cwinx)*(640.0/pixX) 
        cpixy=np.array(cwiny)*(480.0/pixY)

        #translate coords so that top left coord is at [0,0]
        cminx=np.min(cpixx)
        cmaxy=np.max(cpixy)
        #print cminx,cmaxy
        cpixx=np.rint(cpixx-cminx)
        cpixy=np.rint(cmaxy-cpixy)
        #print int(np.max(cpixx))+640,int(np.max(cpixy))+480

        #make new blank image
        background=Image.new("L",[int(np.max(cpixx))+649,int(np.max(cpixy))+492], 0xff) #[0,0] is TOP left #extra 9x12 px for rounding errors
        #background=Image.new("1",[int(np.max(cpixy))+480,int(np.max(cpixx))+640], "white") #[0,0] is TOP left

        #put things in their place
        #fnames=[fnames[2],fnames[6]]
        for i,f in enumerate(fnames):#zip([2,6],fnames):#enumerate(fnames):
            im=Image.open(f)
            box=(int(cpixx[i]),int(cpixy[i]))#,640+int(cpixx[i]),480+int(cpixy[i]))
            background.paste(im,box)

        if plot: #show the mosaic
            #fig,ax=plt.subplots()
            #ax.scatter(cwinx,cwiny)
            #ax.imshow(background, cmap=cm.gray)
            #fig.show()
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
    return filename+'.p'
    #return mfiles
    #return fnames,cpixx,cpixy


def Canny_edge(filen,sigma=3,gauss=False,plot=False):
    #max contrast
    if filen.endswith('.p'):
        imarr=pickle.load(open(filen,'rb'))
    else:
        imarr=im2ndarray(filen)

    im = exposure.equalize_hist(imarr)    
    im=remove_edges(filen,im)

    # Equalization - another method?
    #selem = disk(30)
    #img_eq = rank.equalize(img, selem=selem)

    if type(gauss) == int:
        im = ndi.gaussian_filter(im, gauss)

    #for the finest grids, use sigma=2
    windex,foo=get_index(filen)
    if windex in ['32','42','34','44']: #although I should actually use the dictionary keys to do this stuff in case we change the naming convention
        sigma=2
    if windex in ['11','21']: #although I should actually use the dictionary keys to do this stuff in case we change the naming convention
        sigma=6
        
    edges = feature.canny(im, sigma=sigma)

    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        #ax = ax.ravel()

        ax.imshow(im, cmap=cm.gray)
        ax.imshow(np.ma.masked_where(edges == 0,edges),cmap=cm.autumn)
        ax.set_title('Input image overlaid with Canny edges')
        ax.set_axis_off()
        fig.show()

    newfilen=filen[:-4]+'_edges.p'
    pickle.dump(edges,open(newfilen,'wb'))
    return newfilen

def is_neighbor(pix1,pix2):
    '''Determine if 2 pixels are touching, given their coordinates'''
    if pix1[0]==pix2[0] and (pix1[1]+1 == pix2[1] or pix1[1]-1 == pix2[1]):
        return pix2
    elif pix1[1]==pix2[1] and (pix1[0]+1 == pix2[0] or pix1[0]-1 == pix2[0]):
        return pix2
    elif (pix1[1]+1 == pix2[1] or pix1[1]-1 == pix2[1]) and (pix1[0]+1 == pix2[0] or pix1[0]-1 == pix2[0]): #diagonal
        return pix2

def group_adjacent(coords):
    '''Return list(s) of adjacent coordinates'''
    import copy
    all_neighbors=[]
    for i,c in enumerate(coords):
        #print 'coord ',c
        c2=copy.deepcopy(coords)
        c2.remove(c) #now this is all the coordinates except the one we're looking at
        #print coords,c2
        neighbors=[is_neighbor(c,c1) for c1 in c2 if is_neighbor(c,c1) is not None] #this is the list of all neighbors of c
        if neighbors != []:
            all_neighbors.append([[c]+neighbors]) #now we have a list of all coords c and their neighbors [c1,c2,...cn]
    #which neighbors have neighbors? return groups of neighboring pixels
    groups=[]

    #first let's sort the neighbors list to make it easier to search and group:
    all_neighbors_sorted=[sorted(all_neighbors[j][0]) for j in range(0,len(all_neighbors)-1)]
    backup=copy.deepcopyall_neighbors_sorted()
    #now we can start element-by-element grouping:
    for row in all_neighbors_sorted:
        all_neighbors_sorted.remove(row) #all other rows
        for coord in row:
            for orow in all_neighbors_sorted:
                if coord in orow:
                    #append new ones
                    for orow_element in orow:
                        if orow_element not in row: #my God this is ugly
                            row.append(orow_element) 
    
    return all_neighbors_sorted,'ha'

def dir_adjacent(coords):
    '''Find the direction that adjacent pixels are oriented in'''
    #first sort by x,y

def cleanup_slits(edges,n=10):
    '''Remove false edges from inside the slits'''
    #if not divisible, set back to default:
    if np.mod(480,n) !=0 or np.mod(640,n) !=0:
        print 'Not divisible by ' + str(n) + '! resetting to default n=10'
        n=10
    #Let's treat NxN chunks of pixels at a time
    xran=range(0,480,n)
    yran=range(0,640,n)

    #first get chunks
    chunks=[]
    for i,xin in enumerate(xran):
        for j,yin in enumerate(yran):
            if i == len(xran)-1: xmax=479
            else: xmax=xran[i+1]
            if j == len(yran)-1: ymax=639
            else: ymax=yran[j+1]
            chunks.append(edges[xin:xmax,yin:ymax])

    zipcoords10=[[i,j] for i in range(0,n) for j in range(0,n)]
    zipcoordsy9=[[i,j] for i in range(0,n) for j in range(0,n-1)]
    zipcoordsx9=[[i,j] for i in range(0,n-1) for j in range(0,n)]
    zipcoords9=[[i,j] for i in range(0,n-1) for j in range(0,n-1)]
    #now see which values in each chunk are True and get their indexes in the chunk    
    for chunk in chunks[29:30]:
        try:
            tc=[zc for zc in zipcoords10 if chunk[zc[0],zc[1]]==True]
        except IndexError:
            try:
                tc=[zc for zc in zipcoordsy9 if chunk[zc[0],zc[1]]==True]
            except IndexError:
                try:
                    tc=[zc for zc in zipcoordsx9 if chunk[zc[0],zc[1]]==True]
                except IndexError:
                    tc=[zc for zc in zipcoords9 if chunk[zc[0],zc[1]]==True]

        #print tc            
        fig,ax=plt.subplots()
        ax.imshow(chunk,cmap=cm.binary)
        fig.show()
        #now see which values in the list are adjacent
        alln,group=group_adjacent(tc)

        #now get the directions and lengths of the shapes
        #for group in groups:
        #    grouplen=len(group[1])
        #    print grouplen
        #    groupdir=dir_adjacent(group)
        #    #if they are not in the 'right direction', or too short, mask/delete them
        #    if grouplen < 5:
        #        chunk[group[0]]=False
        #        for coord in group[1]:
        #            chunk[coord]=False
        #fig,ax=plt.subplots()
        #ax.imshow(chunk,cmap=cm.binary)
        #fig.show()
    
    return group,chunk,tc,alln


def slit_or_slat(j,row,imarr): #for j,row in enumerate(edges): #MPI this! super slow right now....
    tv=np.where(row == True)
    tvshape=np.shape(tv)[1]
    immean=np.mean(imarr)
    try:
        width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1) if np.mean(imarr[tv[0][i+1]-tv[0][i],j]) < immean] #slats = space between slits
    except ValueError:
        width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1)]
        if width != []:
            width=filter(lambda a:a!=1.0,width) #filter out all the 1px separations
            width=np.array(width)*pix2num
    return width

def get_slit_width(edges,im=False,mag=5.0,window_num=False, title=''):
    '''Currently this will only work for the finest grids, that have well-defined edges without stuff in the middle'''
    import time
    import multiprocessing as mpi
    start=time.time()
    #let's go row-by-row and get the distance between True values, and compile into a histogram:
    widths=[]
    if mag == 5.0:
        pix2um=(1.954983/640.0)*1000.0 #from RearGrid/windowcorners.txt
    else: #15X
        pix2um= (1.2512/640.0)*1000.0
    #pixX=mscope['FOVx']#1.955 #this is too big....
    #pixY=mscope['FOVy']#1.9548 #from Matej's logs since mine don't have them
    #if im !=False: #it's a numpy array of an image that we will use to determine if slit or slat. Let's make sure these are the same shape
    try:
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
    except ValueError:
        pass

    pool=mpi.Pool()
    
    #def slit_or_slat(j,row,imarr): #for j,row in enumerate(edges): #MPI this! super slow right now....
    #    tv=np.where(row == True)
    #    tvshape=np.shape(tv)[1]
    #    immean=np.mean(imarr)
    #    try:
    #        width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1) if np.mean(imarr[tv[0][i+1]-tv[0][i],j]) < immean] #slats = space between slits
    #    except ValueError:
    #        width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1)]
    #        if width != []:
    #            width=filter(lambda a:a!=1.0,width) #filter out all the 1px separations
    #            width=np.array(width)*pix2num
    #    return width
    
    if type(edges) == list: #it's a list of edge arrays
        for earr in edges:
            for j,row in enumerate(edges): #MPI this! super slow right now....
                result=pool.apply_async(slit_or_slat,args=(j,row,imarr))
                widths.append(result.get())

    else:
        for j,row in enumerate(edges): #MPI this! super slow right now....
            result=pool.apply_async(slit_or_slat,args=(j,row,imarr))
            widths.append(result.get())
            #widths.append(slit_or_slat(j,row,imarr))

    widths_vec=[item for sublist in widths for item in sublist ]
    print 'Processing took %.2f seconds' % (time.time() - start)

    fig,ax=plt.subplots()
    bins=ax.hist(np.array(widths_vec),np.linspace(0,40,40))
    ax.set_xlabel('Separation of Edge Pixels ($\mu$m)')
    ax.set_ylabel('Frequency')
    ax.set_title(title)
    ax.set_yscale('log')
    #if window_num: #set the x-range by pitch
    #    nom_width=
    #    ax.set_xlim([0,])
    fig.show()
    return widths,bins
        
    
def get_theta_range(edges):
    #define range of theta around theta_nominal
    if 'win' not in edges:
        try:
            winnum=int(edges)
        except ValueError:
            edges=raw_input('What is the window number?')
            winnum=int(edges)
    else:
        if edges.startswith('window'):
            winnum=int(edges[6:8])
        else:
            winnum=int(edges[3:5])

    nang=[w['nominal angle'] for w in windows if w['number'] == winnum]
    theta0= nang[0]*(np.pi/180.)#in radians
    tendeg2rad=np.pi/18.
    thetaran = np.linspace(theta0-tendeg2rad, theta0+tendeg2rad, num=101)#in radians
    return thetaran

def prob_hough(edges, threshold=10, line_length=50, line_gap=2,retlines=False,plot=False):
    '''Perform probabilistic Hough fit to given set of edges'''
    if type(edges) ==str: #it's a filename
        inp=edges
        edata=pickle.load(open(edges,'rb'))
    else:
        edata=edges
        edges=raw_input('What is the window number?')

    thetaran=get_theta_range(edges)
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
            pickle.dump(lines,open(newfilen,'wb'))
            return newfilen
    except ValueError:
        pass

    if retlines:
        return lines
        
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
    fig.show()
        
    #return n,bins
    return ndlist,bigarray
    
def hough_hist(lines, windownum, log=False,ret=False, title=False):
    '''Make a histogram of the line orientations returned by Hough'''
    theta=[]
    llist=[]
    if type(lines) == list: #it's a list of filenames
        for f in lines:
            llist.append(pickle.load(open(f,'rb')))
        print np.shape(llist)
        for i,l in enumerate(llist):
            if i == 0:
                all_lines=l
            if i < len(llist)-1:
                all_lines=all_lines + llist[i+1]
        lines=all_lines
    for l in lines:
        if get_length(l) >= 50. :
            theta.append(-1.*get_angle(l))

    thetaran=get_theta_range(windownum)
    #make a histogram
    fig,ax=plt.subplots()

    foo=ax.hist(theta,np.linspace(thetaran[0]*(180./np.pi),thetaran[50]*(180./np.pi),101), facecolor='b')
    if not title:
        ax.set_title('Distribution of line orientations')
    else:
        ax.set_title(title)        
    ax.set_xlim([thetaran[0]*(180./np.pi),thetaran[100]*(180./np.pi)])
    ax.set_xlabel('Probablistic Hough line angle (degrees)')
    if log:
        ax.set_yscale('log')
        ax.set_ylim([1,10000])
    else:
        ax.set_ylim([0,len(lines)])   
    fig.show()
    if ret:
        return theta

def show_Hough_outliers():
    print 'f'

def get_length(line):
    '''Get length of line via Pythagoras'''
    deltax=line[1][0]-line[0][0]
    deltay=line[1][1]-line[0][1]
    length=(deltax**2 + deltay**2)**0.5
    return length

def get_angle(line):
    '''Get angle of line via tangent'''
    deltax=line[1][0]-line[0][0]
    deltay=line[1][1]-line[0][1]
    theta=np.arctan(float(deltay)/float(deltax))
    thetadeg=theta*180./np.pi
    return thetadeg

def plot_all():
    print 'oo'

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
