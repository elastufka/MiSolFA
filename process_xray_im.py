"""
===================
process_xray_im.py
Erica  Lastufka 1.2.18
===================

Find the period and orientation of the moire patterns from the X-ray lamp test images

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


def im2ndarray(filen):
    '''Convert .tif image to numpy.ndarray'''
    imraw = Image.open(filen) #.rotate(45)#.convert('1').convert('L')
    im=np.array(imraw)
    return im

def remove_edges(im,ncrop=50):
    '''remove edges from images with edges of the window in them'''
    imshape=np.shape(im)
    im=im[ncrop:imshape[0]-ncrop,ncrop:imshape[1]-ncrop]
    percentage = np.float((np.shape(im)[0]*np.shape(im)[1]))/np.float((imshape[0]*imshape[1]))
    print "Percentage cropped:", (1.-percentage)*100.
    #print np.shape(im),imshape
    return im

def contrast_stretch(im,plow=2,phigh=98):
    pl = np.percentile(im, plow)
    ph = np.percentile(im, phigh)
    im2 = exposure.rescale_intensity(im, in_range=(pl, ph))
    return im2

def Canny_edge(im,sigma=3,mag=5.0,gauss=False,plot=False,outfilen=False):
    #max contrast
    if gauss:
        im = ndi.gaussian_filter(im, gauss)

    #contrast stretching
    #im=contrast_stretch(im)

    edges = feature.canny(im, sigma=sigma)

    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        #ax = ax.ravel()

        ax.imshow(im, cmap=cm.gray)
        ax.imshow(np.ma.masked_where(edges == 0,edges),cmap=cm.autumn)
        ax.set_title('Input image overlaid with Canny edges')
        ax.set_axis_off()
        fig.show()

    #if not outfilen:
    #    newfilen=filen[:-4]+'_edges.p'
    #else:
    #    newfilen=outfilen
    #pickle.dump(edges,open(newfilen,'wb'))
    return edges#newfilen

def interactive_lines(im): #im is contrast-streched
    global pcoords
    pcoords=[]
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    ax.imshow(im, cmap=cm.gray)
    connection_id = fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect('pick_event', onpick)
    cid = fig.canvas.mpl_connect('key_press_event', on_key)
    fig.show()
    while pcoords==[]:
        plt.pause(10)
    return pcoords

# function to draw lines - from matplotlib examples.  Note you don't need
# to keep a reference to the lines drawn, so I've removed the class as it
# is overkill for your purposes
def draw_line(startx,starty):
        ax = plt.gca()
        xy = plt.ginput(1)
        x = [startx,xy[0][0]]
        y = [starty,xy[0][1]]
        line = ax.plot(x,y, picker=5) # note that picker=5 means a click within 5 pixels will "pick" the Line2D object
        ax.figure.canvas.draw()        

def onclick(event):
    """
    This implements click functionality.  If it's a double click do something,
    else ignore.
    Once in the double click block, if its a left click, wait for a further 
    click and draw a line between the double click co-ordinates and that click
    (using ginput(1) - the 1 means wait for one mouse input - a higher number
    is used to get multiple clicks to define a polyline)
    If the double click was a right click, draw the fixed radius circle

    """
    if event.dblclick:
        if event.button == 1:
            # Draw line
            draw_line(event.xdata,event.ydata) # here you click on the plot
        elif event.button == 3:
            # Write to figure
            plt.figtext(3, 8, 'boxed italics text in data coords', style='italic', bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
            circ = plt.Circle((event.xdata, event.ydata), radius=0.07, color='g', picker = True)
            ax.add_patch(circ)
            ax.figure.canvas.draw()
        else:
            pass # Do nothing

def onpick(event):    
    """
    Handles the pick event - if an object has been picked, store a
    reference to it.  We do this by simply adding a reference to it
    named 'stored_pick' to the axes object.  Note that in python we
    can dynamically add an attribute variable (stored_pick) to an 
    existing object - even one that is produced by a library as in this
    case
    """
    thisline = event.artist #the picked object is available as event.artist
    # print(this_artist) #For debug just to show you which object is picked
    xdata = thisline.get_xdata()
    ydata = thisline.get_ydata()
    ind = [0]#event.ind
    print ind, len(xdata),len(ydata)
    points = [[xdata[ind][0], ydata[ind][0]],[xdata[[-1][0]], ydata[[-1][0]]]]
    pcoords.append(points)
    print('onpick points:', points)

    #plt.gca().picked_object = this_artist

def on_key(event):
    """
    Function to be bound to the key press event
    If the key pressed is delete and there is a picked object,
    remove that object from the canvas
    """
    if event.key == u'delete':
        ax = plt.gca()
        if ax.picked_object:
            ax.picked_object.remove()
            ax.picked_object = None
            ax.figure.canvas.draw()

def get_nline(line,image,flip=True,ngen=3): #line is [x0,y0],[x1,y1]
    meanth=np.arctan((line[1][1]-line[0][1])/(line[0][1]-line[0][0])) #tan theta= dy/dx
    #x=line[1][0]-line[0][0]
    thetan=-1./np.tan(meanth) #slope of perpendicular line is -1/slope = -1/tan(theta)
    print "angle of original line: ", np.rad2deg(meanth)
    print "angle of perpendicular line: ", np.rad2deg(thetan)
    y0n = (np.cos(thetan)) / np.sin(thetan) #line endpoints are (0,0), (x_max,y1)
    y1n = (image.shape[1] * np.cos(thetan)) / np.sin(thetan) #need to get this into a useful format...
    #print y0n,y1n
    if y0n <0 and y1n <0: #shift into the frame
        flip=-1
    else:
        flip=1
    nline=[[(0,y1n+np.max([flip*y0n,flip*y1n])), (image.shape[1],y0n+np.max([flip*y0n,flip*y1n]))]]
    ng=0
    while ng < ngen:
        fac=((ng+.5)*flip*+np.max([flip*y0n,flip*y1n]))/ngen
        nline.append([(0,y1n+fac),(image.shape[1],y0n+fac)])
        ng+=1

    print 'nline',nline[0]
    return nline,np.rad2deg(meanth)

def straight_hough(image,edges,plot=False,side=1.0,spread=5.,n=201,nlines=2):
    '''Perform straight line Hough fit to given set of edges'''
    lines=[]
    h,theta,d=hough_line(edges)
    hlp= hough_line_peaks(h, theta, d)
    meanth = np.mean(hlp[1]) #use this to re-compute the hough to finer resolution in angle
    spreaddeg2rad=spread*(np.pi/180.)  
    thetaran = np.linspace(meanth-spreaddeg2rad, meanth+spreaddeg2rad, num=n)#in radians
    h,theta,d=hough_line(edges, theta=thetaran)
    hlp= hough_line_peaks(h, theta, d)
    morientation=np.rad2deg(np.mean(hlp[1]))
    print 'mean orientation (degrees): ', morientation
    thetan=-1./np.tan(meanth) #slope of perpendicular line is -1/slope = -1/tan(theta)
    print "angle of perpendicular line: ", np.rad2deg(thetan)
    y0n = (np.cos(thetan)) / np.sin(thetan) #line endpoints are (0,0), (x_max,y1)
    y1n = (image.shape[1] * np.cos(thetan)) / np.sin(thetan) #need to get this into a useful format...
    #print y0n,y1n
    if y0n <0 and y1n <0: #shift into the frame
        flip=-1
    else:
        flip=1
    nline=[(0,y1n+np.max([flip*y0n,flip*y1n])), (image.shape[1],y0n+np.max([flip*y0n,flip*y1n]))]

    if nlines > 1: #generate other normal lines
        nline=[nline]
        n=0
        while n < nlines:
            fac=((n+.5)*flip*+np.max([flip*y0n,flip*y1n]))/nlines
            nline.append([(0,y1n+fac),(image.shape[1],y0n+fac)])
            n+=1
        
    for _, angle, dist in zip(*hlp):
        y0 = (dist - 0 * np.cos(angle)) / np.sin(angle) #line endpoints are (0,y0), (x_max,y1)
        y1 = (dist - image.shape[1] * np.cos(angle)) / np.sin(angle) #need to get this into a useful format... 
        #print 0, image.shape[1],y0,y1
        lines.append(((0, y0),(image.shape[1],y1)))
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        ax.imshow(image, cmap=cm.gray)
        for l in lines:
            ax.plot((l[0][0], l[1][0]), (l[0][1], l[1][1]), '-r')
            #print l
        #ax.plot((0,image.shape[1]), (y0n, y1n), '-g') #intersection line
        ax.plot((nline[0][0],nline[1][0]), (nline[0][1],nline[1][1]), '-g')
        #ax.set_xlim((0, image.shape[1]))
        #ax.set_ylim((image.shape[0], 0))
        ax.set_axis_off()
        ax.set_title('Detected lines')
        fig.show()
    return hlp,lines,nline,morientation

def get_intersection_point(line, nline):
    '''Get the intersection points between the Hough lines and the line normal to them. Use this to calculate the period '''
    #lines in format: (x0,y0), (x1,y1)
    x0,y0,x1,y1=float(line[0][0]),float(line[0][1]),float(line[1][0]),float(line[1][1]) 
    x0n,y0n,x1n,y1n=float(nline[0][0]),float(nline[0][1]),float(nline[1][0]),float(nline[1][1])

    xdiff = (x0 - x1, x0n -x1n)
    ydiff = (y0-y1,y0n-y1n) 

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')

    d = (det(*line), det(*nline))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return (x, y)

    #m1=(y1-y0)/(x1-x0)
    #m2=(y1n-y0n)/(x1n-x0n)
    #print m1,m2
    #xint = (m1*x1-m2*x1n+y1n-y1)/(m1-m2)
    #yint = m1*(xint-x1)+y1
    #return (xint,yint)

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

def group_int_points(intpoints,threshold):
    '''Group intersection points with a certain line, return periods calculated'''
    mgroups,sgroups,ggroups=[],[],[] #list of coordinate pair INDICES that are closer to each other than a certain threshold
    for j,i in enumerate(intpoints[:-1]):
        cdist=calc_dist([i,intpoints[j+1]])
        if cdist[0] <= threshold:
            #print cdist[0]
            mgroups.append(j)#((j,j+1))
        else:
            sgroups.append(j)
    
    #group again
    import more_itertools as mit
    ggroups=[list(group) for group in mit.consecutive_groups(mgroups)]

    #print 'ggroups',ggroups
    #average over group:
    newintpoints=[]
    for s in sgroups:
        newintpoints.append((intpoints[s][0],intpoints[s][1]))
    for g in ggroups:
        #this is a list of indices that can be used in intpoints
        xav=np.mean([intpoints[x][0] for x in g])
        yav=np.mean([intpoints[y][1] for y in g])
        #print xav,yav
        newintpoints.append((xav,yav))
    if newintpoints !=[]:
        intpoints=newintpoints
    #group into rising/falling
    odd = [d for i,d in enumerate(intpoints) if i % 2 == 1]
    even = [d for i,d in enumerate(intpoints) if i % 2 == 0]
    podd=calc_dist(odd)
    peven=calc_dist(even)
    print "mean period: ",np.mean([np.mean(podd),np.mean(peven)])
    return odd,even,podd,peven

def calc_period_from_lines(lines,nline,plot=False,image=False,threshold=2,figname='figure'):
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    lines.sort()
    intpoints,allpodd,allpeven=[],[],[]
    for k,nl in enumerate(nline):
        for l in lines:
            #print get_intersection_point(l,nline)
            intpoints.append(get_intersection_point(l,nl))

        #group intersection points if necessary:
        odd,even,podd,peven=group_int_points(intpoints,threshold)
        allpodd.append(podd)
        allpeven.append(peven)
        mperiod=np.mean([np.median(podd),np.median(peven)])
        bins=np.linspace(mperiod-(mperiod/5.),mperiod+(mperiod/5.),20)
        ax[1].hist(podd,bins,facecolor=(.1,float(k)/float(len(nline)),.2),alpha=0.6,label='line '+str(k)+ ' odd')
        ax[1].hist(peven,bins,facecolor=(.2,float(k)/float(len(nline)),.1),alpha=0.6,label='line '+str(k)+' even')

    #calculate mean period
    modd=[a for o in allpodd for a in o]
    meven=[a for o in allpeven for a in o]
    mperiod=np.mean([np.median(modd),np.median(meven)])
    #mdev=np.min([np.std(modd),np.std(meven)])
    #print mperiod, mdev
    
    ax[0].imshow(image, cmap=cm.gray,alpha=0.5)
    for l in lines:
        ax[0].plot((l[0][0], l[1][0]), (l[0][1], l[1][1]), '-r')
    for nl in nline:
        ax[0].plot((nl[0][0],nl[1][0]), (nl[0][1],nl[1][1]), '-g')
    ipx=[x[0] for x in intpoints]
    ipy=[x[1] for x in intpoints]
    ax[0].scatter(ipx,ipy)
    ax[0].set_xlim((0, image.shape[1]+50))
    ax[0].set_ylim((image.shape[0]+50, 0))
    ax[0].set_axis_off()
    ax[0].set_title('Detected lines and intersection points')
    #period histogram
    #n,bins=ax[1].hist(p1odd,,label='')
    ax[1].set_xlim((mperiod-(mperiod/5.),mperiod+(mperiod/5.)))
    #ax[1].set_ylim((image.shape[0]+50, 0))
    ax[1].set_title('Calculated Period (pixels)')
    ax[1].legend()
    if plot:
        fig.show()
    else:
        plt.savefig(figname+'.png')
   
    return intpoints, allpodd,allpeven
        
def edge_and_plot(filen,ncrop=10,sigma=3,show=True,gauss=False,fname=False,nlines=2,threshold=2,interactive=False):
    print filen
    im=im2ndarray(filen)
    im=remove_edges(im,ncrop=ncrop)
    #if show:
    #    im_peek(im)
    im=contrast_stretch(im)
    if show:
        im_peek(im)
    if interactive:
        lines=interactive_lines(im)
        #lines=[[t[0],tpls[i+1][0]] for i,t in enumerate(tpls[:-1])]
        nline,morientation=get_nline(lines[0],im)
        #nline=[nline
    else:
        edges=Canny_edge(im,sigma=sigma,gauss=gauss, plot=show)
        hlp,lines,nline,morientation=straight_hough(im,edges,nlines=nlines,plot=False)
    ips, odd, even=calc_period_from_lines(lines, nline,plot=show, image=im,figname=fname,threshold=threshold)
    return ips, odd, even,nline,morientation

def calc_all_windows(data_dir, pix2um,infile=False,sigma=3,gauss=False,nlines=2,json=False,plot=True,filename=False,tag='plot',tol=2,interactive=False):
    '''Segment/images need to first be cropped into individual window images'''
    os.chdir(data_dir)
    #get files
    if not infile:
        imfiles=glob.glob('win*.tif') #should probably throw an error if they're not all there
    else:
        imfiles=[infile]
    data_dict={"windows":[],"data_dir":data_dir}
    for fn in imfiles: #w in [11,21,12,22,31,41,32,42,33,43,34,44]: #get period for each window
        ips,odd,even,nline,morientation=edge_and_plot(fn,show=plot,sigma=sigma,gauss=gauss,fname=fn[:-3]+tag,threshold=tol,interactive=interactive)
        w=fn[3:fn.find('.')]
        #calculate mean period
        modd=[a for o in odd for a in o]
        meven=[a for o in even for a in o]
        d={"number":int(w),"mperiod":np.mean([np.mean(modd),np.mean(meven)]),"morientation":morientation,"nlines":nline,"odd":odd,"even":even,"pix2um":pix2um,"image":fn}
        data_dict["windows"].append(d) #might need to modify this
        
    #save the dictionary
    if not filename:
        filename=data_dir[data_dir.rfind('/')+1:]
    pickle.dump(data_dict,open(filename+'.p','wb'))
    if json:
        import json
        with open(filename+'.json','w') as outfile:
            json_string=json.dump(data_dict,outfile)
    return data_dict

def rcalc(Mp, Mo,d=7):
    '''reverse calcualte periods and orientations using formulas in Diego's notes'''
    #does condition 3 mean that we can't use that here? he is presuming a certain orientation of the pattern
    #or I can 're-orient' the pattern... 
    fac=d/(2.*Mp) #d is in mm, Mp should also be in mm!
    alpha=np.rad2deg(np.arctan2(1./(1-fac)))
    beta=np.rad2deg(np.arctan2(1./(1+fac)))
    p=q
    q=p
    return p,alpha,q,beta

def im_peek(im,mag=5.0,length=0.2):
    '''Plot Canny edges over image for a given file'''
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib_scalebar.scalebar import SI_LENGTH
    if mag==5.0:
        pix2um=(1.2512/640.)*1000.
    else:
        pix2um=0.6
    #edges=pickle.load(open(filen,'rb'))
    #imf=filen[:-8]+'.tif'
    #im=im2ndarray(filen)
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    ax.imshow(im, cmap=cm.gray)
    #ax.imshow(np.ma.masked_where(edges == 0,edges),cmap=cm.autumn)
    scalebar = ScaleBar(pix2um,'um', SI_LENGTH,length_fraction=length) # 1 pixel = 0.2 meter
    #print scale*pix2um
    ax.add_artist(scalebar)
    #ax.set_title('Input image overlaid with Canny edges')
    ax.set_axis_off()
    fig.show()
    
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


