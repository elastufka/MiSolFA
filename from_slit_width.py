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
from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit
from scipy import interpolate
import time
import glob

#list of dictionaries of the window numbers, pitches (mm), and nominal angles (deg)
global windows
global windowsr



def im2ndarray(filen):
    '''Convert .tif image to numpy.ndarray'''
    imraw = Image.open(filen)#.convert('B') #.rotate(45)#.convert('1').convert('L')
    im=np.array(imraw)
    return im

def get_intersect(m1,m2,b1,b2):
    '''get intersection point (x,y) between two lines'''
    x=(b2-b1)/(m1-m2)
    y=m1*x+b1
    return x,y

def scat_allP(win,widths=True,periods=False, bsum=False,meanfit=True,figname=False,theta_ran=False,print_height=True,legend_on=True):
    sumvecs=[]
    errvecs=[]
    fitsplus=[]
    fitsminus=[]
    allfitsplus,allfitsminus=[],[]
    for p in range(0,7):
        if bsum:
            files=glob.glob('win'+win+'_p'+str(p)+'*'+'_groupstretch.tif')
            xv,sv=basic_sum(files,groupstretch=True,show=False)
            #fix order for plotting connecting line
            #xvmid=xv.index(0.0)
            #xv.sort()
            xvmid=xv.index(0.0)
            #svsort=(sv[:xvmid][::-1]+sv[xvmid:])/np.max(sv)
            aa=zip(xv,sv/np.max(sv))
            aa.sort()
            xv.sort()
            svsort=[b[1] for b in aa]
            sumvecs.append(svsort)
        if widths or periods:
            if widths:
                files=glob.glob('win'+win+'_width_stats_p'+str(p)+'*.p')
            else:
                files=glob.glob('win'+win+'_period_stats_p'+str(p)+'*.p')
                
            files.sort()
            #print files
            if theta_ran:
                badf=[]
                for f in files:
                    theta=float(f[f.rfind('_')+1:-2])
                    if theta < theta_ran[0] or theta > theta_ran[1]:
                        badf.append(f)
                for b in badf:
                    files.remove(b)
            svsort,std,xv,xaxlabels=[],[],[],[]
            for f in files:
                sdict=pickle.load(open(f,'rb'))
                svsort.append(sdict['mean'])
                std.append(sdict['stddev'])
                xv.append(np.tan(np.deg2rad(float(f[f.rfind('_')+1:-2])))) #these should already be sorted correctly
                xaxlabels.append(float(f[f.rfind('_')+1:-2]))
            aa=zip(xv,svsort,std)
            aa.sort()
            xv.sort()
            svsort2=[b[1] for b in aa]
            sumvecs.append(svsort2)
            stdsort=[b[2] for b in aa]
            errvecs.append(stdsort)
            xvmid=xv.index(0.0)
            allfitsplus.append(np.polyfit(xv[:xvmid-1],svsort2[:xvmid-1],1))
            allfitsminus.append(np.polyfit(xv[xvmid+1:],svsort2[xvmid+1:],1))
            #print np.shape(allfitsplus),allfitsplus
    meanvec=np.mean(sumvecs,axis=0)
    fplus=np.polyfit(xv[:xvmid-1],meanvec[:xvmid-1],1)
    fminus=np.polyfit(xv[xvmid+1:],meanvec[xvmid+1:],1)
    #print np.min(svsort), np.max(svsort),fplus,fminus
    fitsplus.append(fplus)
    fitsminus.append(fminus)

    labels=['p0','p1','p2','p3','p4','p5','p6']
    cols=['b','c','m','k','r','g','y']
    fig,ax=plt.subplots(2,sharex=True)
    #ax=fig.add_subplot(112)
    for sv,c,lab in zip(sumvecs,cols,labels):
        ax[0].plot(xv,sv,'-o',color=c,label=lab)
    ax[0].axvline(color='k',linestyle='dashed')
    if meanfit:
        meanplusslope=np.mean([fp[0] for fp in fitsplus])
        meanplusintercept=np.mean([fp[1] for fp in fitsplus])
        meanplusline=meanplusslope*np.array(xv[:xvmid+1]) + meanplusintercept
        #print meanplusline
        meanminusslope=np.mean([fm[0] for fm in fitsminus])
        meanminusintercept=np.mean([fm[1] for fm in fitsminus])
        meanminusline=meanminusslope*np.array(xv[xvmid:]) + meanminusintercept
        lineintx,lineinty=get_intersect(meanplusslope,meanminusslope,meanplusintercept,meanminusintercept)
        #print meanminusline
        ax[0].plot(xv[:xvmid+1],meanplusline,'k--')
        ax[0].plot(xv[xvmid:],meanminusline,'k--')
        yran=ax[0].get_ylim()
        ymax=float(yran[1])
        xran=ax[0].get_xlim()
        xmin=float(xran[0])
        xmax=float(xran[1])
        ax[0].text(.75*xmin,.2*ymax,'avg. slope: ' + str(np.round(meanplusslope,3)))
        ax[0].text(.75*xmin,.1*ymax,'avg. intercept: ' + str(np.round(meanplusintercept,3)))
        ax[0].text(.5*xmax,.2*ymax,'avg. slope: ' + str(np.round(meanminusslope,3)))
        ax[0].text(.5*xmax,.1*ymax,'avg. intercept: ' + str(np.round(meanminusintercept,3)))
        ax[0].text(.15*xmin,.3*ymax,'adj. peak: (' + str(np.round(lineintx,3))+ ',' + str(np.round(lineinty,3))+')')
        
    #ax.set_ylim([0,1.1])
    #ax.set_xlim([xv[0],xv[-1]])
    ax[1].set_xlabel('tan(Angle)') #should I set the lables to the actual angle values? xaxlabels... but they would be in the wrong place
    if widths:
       ax[0].set_title('Window '+ win + ' slit width as a function of angle') 
       ax[0].set_ylabel('Slit width ($\mu$m)')
    elif periods:
       ax[0].set_title('Window '+ win + ' period as a function of angle') 
       ax[0].set_ylabel('Period ($\mu$m)')
        
    else:
        ax[0].set_title('Window '+ win + ' pixel sum as a function of angle')
        ax[0].set_ylim([0,1.1])
        ax[0].set_ylabel('Percent of maximum transmission')
    if legend_on:
        ax[0].legend(loc='upper right')

    #make the lower subplot showing the error bars
    #ax2=fig.add_subplot(121)
    #generate the points for the fit lines
    meanpluspoints=(meanplusslope*np.array(xv[:xvmid+1]) + meanplusintercept).tolist()
    meanminuspoints=(meanminusslope*np.array(xv[xvmid+1:]) + meanminusintercept).tolist()
    meanpluspoints.extend(meanminuspoints)
    #should plot the error bars here
    for sv,c,lab,err in zip(sumvecs,cols,labels,errvecs):
        yv=np.array(sv)-np.array(meanpluspoints)
        ax[1].plot(xv,yv,'o',color=c,label=lab)
        ax[1].errorbar(xv,yv,yerr=err)
    ax[1].axhline(color='k',linestyle='dashed')

    #re-adjust plot sizes
    box1=ax[0].get_position()
    box2=ax[1].get_position()
    ax[0].set_position([box1.x0,box1.y0*.7,box1.width,box1.height*1.5])
    ax[1].set_position([box2.x0,box2.y0,box2.width,box2.height*.7])
    if figname:
        plt.savefig(figname)
        data=[xv,sumvecs,fitsplus,fitsminus,meanplusslope,meanplusintercept,meanminusslope,meanminusintercept,lineintx,lineinty]
        pickle.dump(data,open(figname[:-4]+'.p','wb'))
    else:
        fig.show()
    if print_height:
        print allfitsplus,allfitsminus
        for fp,fm,i in zip(allfitsplus,allfitsminus,range(0,len(allfitsplus))):
            heights,meanheight=calc_height_from_fit(xv,fp[0],fp[1],fm[0],fm[1])
            print 'p',i,' mean height: ', meanheight
    return xv,sumvecs,errvecs

def calc_height_from_fit(xv,l1m,l1b,l2m,l2b):
    ''' Input is line parameters: slopes, intercepts, xvec'''
    #first find w_a, or width at theta=0 using the intercept point of the two lines
    theta_off,w_a=get_intersect(l1m,l2m,l1b,l2b) #theta is offset by the coordinate x of the intercept, and w_a is the y-value
    xvmid=xv.index(0.0) #assume xv is tan(theta) and the line params have been calculated with tan(theta)
    w_mvec=(l1m*np.array(xv[:xvmid+1]) + l1b).tolist()
    meanminuspoints=(l2m*np.array(xv[xvmid+1:]) + l2b).tolist()
    w_mvec.extend(meanminuspoints)
    heights=[np.abs((w_a - w_m)/tantheta) for tantheta,w_m in zip(xv,w_mvec) if tantheta !=0.0] #should do this per p
    #if theta is changing uniformly and so is w_m, shouldn't all heights be the same?
    meanheight=np.mean(heights)
    return heights,meanheight

def calc_height_from_data(xv,sumvecs):
    '''Calculate height from actual data points'''
    #first find w_a, or width at theta=0 using the intercept point of the two lines
    theta_off,w_a=get_intersect(l1m,l2m,l1b,l2b) #theta is offset by the coordinate x of the intercept, and w_a is the y-value
    xvmid=xv.index(0.0) #assume xv is tan(theta) and the line params have been calculated with tan(theta)
    w_mvec=(l1m*np.array(xv[:xvmid+1]) + l1b).tolist()
    meanminuspoints=(l2m*np.array(xv[xvmid+1:]) + l2b).tolist()
    w_mvec.extend(meanminuspoints)
    heights=[np.abs((w_a - w_m)/tantheta) for tantheta,w_m in zip(xv,w_mvec) if tantheta !=0.0] #should do this per p
    #if theta is changing uniformly and so is w_m, shouldn't all heights be the same?
    meanheight=np.mean(heights)
    return heights,meanheight

def Canny_edge(filen,sigma=3,gauss=False,plot=False,outfilen=False):
    #max contrast
    from skimage import io
    if filen.endswith('.p'):
        imarr=pickle.load(open(filen,'rb'))
    else:
        imarr=io.imread(filen)#im2ndarray(filen)
    im=imarr#.astype(int)
    print np.min(im),np.max(im)
    #convert to binary?

    if gauss:
        im = ndi.gaussian_filter(im, gauss)

    edges = feature.canny(im, sigma=sigma)

    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        #ax = ax.ravel()

        ax.imshow(im, cmap=cm.gray,alpha=.4)
        ax.imshow(np.ma.masked_where(edges == 0,edges),cmap=cm.autumn)
        ax.set_title('Input image overlaid with Canny edges')
        ax.set_axis_off()
        fig.show()

    if not outfilen:
        newfilen=filen[:-4]+'_edges.p'
    else:
        newfilen=outfilen
    pickle.dump(edges,open(newfilen,'wb'))
    return edges

def clean_centers(edges,tolerance=False,plot=False,sigma=2.):
    '''clean out the Canny edge array so that pixel bins with low counts (ie probably bad edges) get deleted below a certain threshold'''
    esum=np.sum(edges,axis=0)
    cleaned_edges=np.copy(edges)
    if not tolerance:
        tolerance=np.mean(esum)+sigma*np.std(esum)
    #print tolerance
    for i,col in enumerate(esum):
        if col < tolerance and col != 0:
            #print i, col
            cleaned_edges[:,i] = False
    if plot:
        cleaned_sum=np.sum(cleaned_edges,axis=0)
        fig,ax=plt.subplots()
        ax.plot(range(0,len(cleaned_sum)),cleaned_sum)
        fig.show()
    return cleaned_edges

def sum_peaks(cleaned_edges, tol=2,irange=False):#irange=[1087,1105]):
    '''Get the locations of the peaks of a (Gaussian?) fit to the histogram of sums, return these as the width array. Tolerance is number of False allowed between Trues to still be considered a group'''
    esum=np.sum(cleaned_edges,axis=0)
    if irange:
        esum=esum[irange[0]:irange[1]]
    group_hist,group_xvec=[],[]
    #first get the groups of True values
    group=False
    for i,es in enumerate(esum[:-1]):
        if irange:
            print i+irange[0],es,group
        if es != 0.0 and not group:
            lhist=[es]
            lxvec=[i]
            group=True
        elif es !=0.0 and group:
            lhist.append(es)
            lxvec.append(i)
            #print 'here',esum[i+1]
            if np.mean(esum[i+1:i+tol]) == 0 :
                #print 'here'
                group=False
        if not group:
            try:
                if lhist not in group_hist:
                    group_hist.append(lhist)
                    group_xvec.append(lxvec)
            except NameError:
                continue
    #now build little histograms out of the groups and fit gaussians to them. return the new x-vector (float!) of peak locations
    #this actually doesn't work well because there are so few points they don't make a normal distribution! So let's take the weighted average instead
    xpeak=[]
    #print group_xvec,group_hist
    for xg,yg in zip(group_xvec,group_hist):
        xpeak.append(np.average(xg, weights=yg))
        #if len(xg) > 2: #can't make a gaussian out of only 2 points...
        #    xpeak.append(gauss_peak(xg,yg))
        #else:
        #    xpeak.append(np.mean(xg))
        
    return xpeak

def gauss_peak(xdata,ydata,npts=50,show=False,ret_data=False):
    ''' returns peak of gaussian fit to data'''
    #x-coord does not have to be an integer!
    # Create a function which returns a Gaussian (normal) distribution.
    def gauss(x,*p0):
        if len(p0) == 1:
            a,b,c,d=p0[0]            
        else:
            a,b,c,d=p0
        y = a*np.exp(-np.power((x - b), 2.)/(2. * c**2.)) + d
        return y
    #interpolate x and y so greater than 4 (or else curvefit doesn't work)
    if len(xdata) < 4:
        intfunc=interpolate.interp1d(xdata,ydata)
        xint=np.linspace(xdata[0],xdata[-1],10)
        yint=intfunc(xint)
    else:
        xint=np.array(xdata)
        yint=np.array(ydata)
    popt, pcov = curve_fit(gauss, xint, yint,p0=[0,0,0,np.mean(yint)])#, p0=p_initial, sigma=e)
    #generate high-resolution version of gaussian curve and get the maximum:
    xvec=np.linspace(xdata[0],xdata[-1],npts)
    yvec=gauss(xvec,popt) #[0]
    xpeak=xvec[np.where(yvec == np.max(yvec))[0]][0]
    if show:
        fig,ax=plt.subplots()
        ax.plot(xint,yint)
        ax.plot(xvec,yvec)
        fig.show()
    if ret_data:
        return xvec,yvec
    else:
        return xpeak #,ypeak

def find_normal(edges, thetaran=[-.1,.1],crotate=False,ctol=.05):
    '''do a staight line Hough transform to find the mean angle of rotation. Rotate the image by the correctionif requested.'''
    h,theta,d=hough_line(edges,np.linspace(thetaran[0],thetaran[1],100))
    lps=hough_line_peaks(h, theta, d)
    meanang=np.mean(lps[1])
    print meanang
    #fig,ax=plt.subplots()
    #ax.imshow(edges,cmap=cm.gray)
    #for _, angle, dist in zip(*hough_line_peaks(h, theta, d)):
    #    y0 = (dist - 0 * np.cos(angle)) / np.sin(angle)
    #    y1 = (dist - np.shape(edges)[1] * np.cos(angle)) / np.sin(angle)
    #    ax.plot((0, np.shape(edges)[1]), (y0, y1), '-r')
    #ax.set_xlim([0,2048])
    #ax.set_ylim([0,2040])
    #fig.show()
    if crotate:
        if np.abs(meanang) > ctol: #rotate the image and return. else return False
            rotim=np.imrotate(e,meanang)
            return rotim
        else:
            return False
    else:
        return False

def slit_widths_from_peaks(window_num,imfile,xpeak=False,pix2um=.65,plotwidth=False,plotperiod=False,stats=True,gauss=False,tolerance=False,n=5):
    '''basically same as slit_or_slat() plus histogramming'''
    p0ang=imfile[6:imfile.rfind('_')]
    im=np.array(Image.open(imfile))
    imsum=np.sum(im,axis=0)
    immean=np.mean(imsum)
    
    if not xpeak:
        efile=glob.glob(imfile[:-4]+'_edges.p')
        if len(efile) > 0:
            edges=pickle.load(open(efile[0],'rb'))
        else:
            edges=Canny_edge(imfile,sigma=3,gauss=gauss,plot=False)
        cleaned_edges=clean_centers(edges,sigma=tolerance)
        xpeak=sum_peaks(cleaned_edges,tol=n)
    
    xpeak_int=[int(x) for x in xpeak]
    
    width=[float(xpeak[i+1]-xpeak[i]) for i in range(0,len(xpeak)-1) if np.mean(imsum[xpeak_int[i]:xpeak_int[i+1]]) > immean] #need xpeaks to be integers now
    peven=[float(xpeak[i+2]-xpeak[i]) for i in range(0,len(xpeak)-2) if i % 2 == 0]
    podd=[float(xpeak[i+2]-xpeak[i]) for i in range(0,len(xpeak)-2) if i % 2 != 0]
    if width != []:
        width=filter(lambda a:a > n,width) #filter out all the 0- n px separations
        median=np.median(width) #assume it gets it right more often than not. might not be the case for the finer grids...
        #print median
        widths=[]
        for i in range(0,len(width)-1):
            #which is smaller, width % median, width % n*median?
            factor=float(int(width[i])/int(median))
            if factor == 0:
                factor = 1.
            #print width[i],factor
            widths.append(width[i]/factor)
        widths=np.array(widths)*pix2um # need to account for the fact that it's rotated! width is not true width!
    else:
        widths =[]

    #do the same for podd and peven
    if podd != []:
        podd=filter(lambda a:a > n,podd) #filter out all the 0- n px separations
        peven=filter(lambda a:a > n,peven) #filter out all the 0- n px separations
        odd_median=np.median(podd)
        even_median=np.median(peven)
        #print podd,peven
        podds,pevens=[],[]
        for i in range(0,len(peven)-1):
            #which is smaller, width % median, width % n*median?
            factor=float(int(peven[i])/int(even_median))
            if factor == 0:
                factor = 1.
            #print width[i],factor
            pevens.append(pix2um*(peven[i]/factor))
        for i in range(0,len(podd)-1):
            #which is smaller, width % median, width % n*median?
            factor=float(int(podd[i])/int(odd_median))
            if factor == 0:
                factor = 1.
            #print width[i],factor
            podds.append(pix2um*(podd[i]/factor))
        #pevens=np.array(pevens)*pix2um # need to account for the fact that it's rotated! width is not true width!
        #podds=np.array(podds)*pix2um # need to account for the fact that it's rotated! width is not true width!
    else:
        widths =[]
    if plotwidth:
        #make the histogram 
        fig,ax=plt.subplots()
        bins=np.arange(np.min(widths),np.max(widths),np.std(widths)/5.)
        ax.hist(widths,bins)
        #ax.set_xlim([nperiod-5,nperiod+5])
        ax.set_yscale('log')
        ax.set_ylim([1,10000])
        #figfilename='win'+str(window_num)+'_group_periods'+str(mag)+'.png'
        #fig.savefig(figfilename)
        fig.show()

    if plotperiod:
        #make the histogram 
        fig,ax=plt.subplots()
        bins=np.arange(np.min(pevens),np.max(pevens),np.std(pevens)/5.)
        ax.hist(pevens,bins,facecolor='g',alpha=0.6)
        ax.hist(podds,bins,facecolor='m',alpha=0.6)
        ax.set_yscale('log')
        ax.set_ylim([1,10000])
        #figfilename='win'+str(window_num)+'_group_periods'+str(mag)+'.png'
        #fig.savefig(figfilename)
        fig.show()
        
    if stats: #print out and pickle stats
        avgw=np.mean(widths)
        medw=np.median(widths)
        stdvw=np.std(widths)
        pevens.extend(podds)#=pevens+podds
        avgp=np.mean(pevens)
        medp=np.median(pevens)
        stdvp=np.std(pevens)
        wstatfile='win'+str(window_num)+'_width_stats_'+p0ang+'.p'
        wdatafile='win'+str(window_num)+'_width_data_'+p0ang+'.p'
        pstatfile='win'+str(window_num)+'_period_stats_'+p0ang+'.p'
        pdatafile='win'+str(window_num)+'_period_data_'+p0ang+'.p'
        print "-------------------STATISTICS FOR WINDOW "+str(window_num)+"---------------------"
        print '              Mean width: ' + str(avgw)
        print '            Median width: ' + str(medw)
        print 'Standard Deviation width: ' + str(stdvw)
        print '              Mean period: ' + str(avgp)
        print '            Median period: ' + str(medp)
        print 'Standard Deviation period: ' + str(stdvp)
        print 'Results saved in ' + wstatfile + ', '+ pstatfile
        wstdict={'mean':avgw,'median':medw,'stddev':stdvw}
        pstdict={'mean':avgp,'median':medp,'stddev':stdvp}
        pickle.dump(wstdict,open(wstatfile,'wb'))
        pickle.dump(widths,open(wdatafile,'wb'))
        pickle.dump(pstdict,open(pstatfile,'wb'))
        pickle.dump(np.array(pevens),open(pdatafile,'wb'))

    return pevens,widths


def combine_width_data(window_num, ang, period=False,stats=True,save=False,plot=False,n=20):
    if not period:
        dfiles=glob.glob('win'+str(window_num)+'_width_data*_'+str(ang)+'.p')
    else:
        dfiles=glob.glob('win'+str(window_num)+'_period_data*_'+str(ang)+'.p')
        
    allwidths=[]
    for d in dfiles:
        allwidths.extend(pickle.load(open(d,'rb')).tolist())
    if stats:
        avg=np.mean(allwidths)
        med=np.median(allwidths)
        stdv=np.std(allwidths)
        print "-------------------STATISTICS FOR WINDOW "+str(window_num)+", ANGLE " + str(ang)+ "---------------------"
        print '              Mean: ' + str(avg)
        print '            Median: ' + str(med)
        print 'Standard Deviation: ' + str(stdv)
    if save:
        #print 'Results saved in ' + statfile
        data=allwidths
        stdict={'mean':avg,'median':med,'stddev':stdv}
        if not period:
            statfile='win'+str(window_num)+'_allwidth_stats_'+str(ang)+'.p'
            datafile='win'+str(window_num)+'_allwidth_data_'+str(ang)+'.p'
        else:
            statfile='win'+str(window_num)+'_allperiod_stats_'+str(ang)+'.p'
            datafile='win'+str(window_num)+'_allperiod_data_'+str(ang)+'.p'
            
        pickle.dump(stdict,open(statfile,'wb'))
        pickle.dump(data,open(datafile,'wb'))
        print 'Results saved in ' + statfile
    if plot:
        fit_gaussian_to_hist(allwidths,win=window_num,ang=ang,n=n)

    return allwidths

def batch_window(win,angvec,tolvec=False,nvec=False):
    for j,ang in enumerate(angvec):
        ifiles=glob.glob('win'+str(win)+'_*'+str(ang)+'_corrected.tif')
        if tolvec:
            tol=tol[j]
        else:
            tol=False
        if nvec:
            n=nvec[j]
        else:
            n=5
        for i in ifiles:
            pp,bb=slit_widths_from_peaks(win,i,tolerance=tol,n=n)
        combine_width_data(win,ang,period=True,save=True)
    scat_allP(str(win),widths=False,periods=True, bsum=False,meanfit=True,figname=False,theta_ran=False)
        
def fit_gaussian_to_hist(widths,win=False, ang=False,bran=False,n=20):
    '''overplot gaussian on histogram limited by xran. Make a big histogram first of all the positions....'''
    import matplotlib.mlab as mlab
    import math
    #first filter periods. don't need to alter the bins
    if bran:
        bins=np.linspace(bran[0],bran[1],n)
    else:
        bins=n#np.linspace(np.min(widths),np.max(widths),n)
    
    mu=np.mean(widths)
    variance=np.std(widths)
    gauss=mlab.normpdf(bins,mu,variance)
    #plot
    fig,ax=plt.subplots()
    nn,bn,patches=ax.hist(widths,bins,normed=1)
    y=mlab.normpdf(bn,mu,variance)
    ax.plot(bn,y,'r--')
    if win:
        if ang:
            ax.set_title('Window '+str(win) + '  width histogram at angle '+str(ang))
        else:
            ax.set_title('Window '+str(win) + ' width histogram')
    #ax.plot(xdata,gauss)
    #if xran:
    ax.set_xlim([np.min(widths),np.max(widths)])
    #ax.set_yscale('log')
    #ax.set_ylim([1,10000])
    fig.show()
    #return xdata,ydata,gauss

def im_peek(filen,pix2um=0.65,length=0.2):
    '''Plot Canny edges over image for a given file'''
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib_scalebar.scalebar import SI_LENGTH
    #edges=pickle.load(open(filen,'rb'))
    #imf=filen[:-8]+'.tif'
    im=im2ndarray(filen)
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    ax.imshow(im, cmap=cm.gray)
    #ax.imshow(np.ma.masked_where(edges == 0,edges),cmap=cm.autumn)
    scalebar = ScaleBar(pix2um,'um', SI_LENGTH,length_fraction=length) # 1 pixel = 0.2 meter
    #print scale*pix2um
    ax.add_artist(scalebar)
    #ax.set_title('Input image overlaid with Canny edges')
    ax.set_axis_off()
    fig.show()

def edge_peek(filen,pix2um=.65,length=0.2):
    '''Plot Canny edges over image for a given file'''
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib_scalebar.scalebar import SI_LENGTH
    edges=pickle.load(open(filen,'rb'))
    imf=filen[:-8]+'.tif'
    im=im2ndarray(imf)
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    ax.imshow(im, cmap=cm.gray,alpha=.4)
    ax.imshow(np.ma.masked_where(edges == 0,edges),cmap=cm.autumn)
    scalebar = ScaleBar(pix2um,'um', SI_LENGTH,length_fraction=length) # 1 pixel = 0.2 meter
    #print scale*pix2um
    ax.add_artist(scalebar)
    ax.set_title('Input image overlaid with Canny edges')
    ax.set_axis_off()
    fig.show()
    
def hough_peek(filen,edgef=False,pix2um=0.65,length=0.2):
    '''Plot Hough fits over Canny edges for a given file'''
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib_scalebar.scalebar import SI_LENGTH
    lines=pickle.load(open(filen,'rb'))
    if not edgef:
        edgef=filen[:-7]+'edges.p'
    pix2um=0.65
        
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
    
def rising_or_falling(rotim,rotarr,immean,plot=False,imfile=False,shape=False,rerotate=False, test=False):
    """Determine if an edge is rising (dark -> light) or falling (light -> dark). Return indexed mask"""
    imsize=np.shape(rotim)
    if not shape:
        mask=np.zeros([792,792])#np.zeros([640,480])
        oshape=np.zeros([640,480])
    else:
        mask=np.zeros(shape)
        oshape=np.zeros([640,480])
    meani=np.mean(rotim, axis=0) #average along y
    meanarr=np.sum(rotarr, axis=0) # total along y (remember it's boolean)
    #for rowi,rowarr in zip(meani,meanarr):
    tv=np.where(meanarr > 100.)#np.max(meanarr[0:50]))
    tvshape=np.shape(tv)[1]
    #get rising
    rising=[tv[0][i] for i in range(0,tvshape-1) if np.mean(meani[tv[0][i]:tv[0][i+1]]) > 1.25*immean and tv[0][i+1] != tv[0][i]+1 ] 
    
    mask[:,rising]= 1
    #get falling
    falling=[tv[0][i] for i in range(0,tvshape-1) if np.mean(meani[tv[0][i]:tv[0][i+1]]) < .75*immean and tv[0][i+1] != tv[0][i]+1 ]
    mask[:,falling]= -1
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

        fig,ax=plt.subplots()
        #if not imfile:
        #    ax.imshow(im2ndarray('win11_05_05_5.0X.tif'),alpha=0.6,cmap=cm.gray)
        #else:
        #    ax.imshow(im2ndarray(imfile),alpha=0.4,cmap=cm.gray)
        ax.imshow(mask,cmap=cm.binary)
        fig.show()

    if rerotate: #rotate it back by given angle and trim to original size
        mask=imrotate(mask,rerotate)
    #if len(rising) > 1.25* len(falling) or len(falling) > 1.25* len(rising):
        #return []
    return mask

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
        ax.scatter(rising[:-len(falling)],rmean,c='b')
        ax.scatter(falling, fmean,c='r')
        fig.show()

    periodr = [rmean[j+1]-rmean[j] for j in range(0,len(rmean)-1)]
    #periodr = [p for p in periodr if p > tolerance*nomp and p < (2.-tolerance)*nomp]
    #print rmean
    #print periodr
    periodf= [fmean[j+1]-fmean[j] for j in range(0,len(fmean)-1)]
    #periodf = [p for p in periodf if p > tolerance*nomp and p < (2.-tolerance)*nomp]
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

'''
def test_mask(ef, im, nang):
    edges=pickle.load(open(ef,'rb'))
    ima=im2ndarray(im)
    ima=remove_edges(im,ima)
    rotim=rotate(np.transpose(ima),nang, reshape=True)
    rotarr=rotate(np.transpose(edges),nang,reshape=True)
    #print im,ef
    mask=rising_or_falling(rotim[300:950,:],rotarr[300:950,:],np.mean(rotim[300:950,100:-100]), shape=np.shape(rotarr),imfile=im,plot=True,test=True)
    print np.mean(rotim[300:950,:])
'''

def mean_hough_angle(f_all,side=1.0,save=True):
    lines=pickle.load(open(f_all,'rb'))
    angles=[]
    for l in lines:
        angles.append(side*get_angle(l))
    meana=np.mean(angles)
    if save:
        pickle.dump(meana,open(f_all[:-2]+'_mean.p','wb'))
    return meana


def get_period_by_grouping(window_num,ang,plot=True,nperiod=90.,pix2um=0.65,stats=True,EM=True,tolerance=0.6):
    '''Put everything together for a single window'''
    #first get the files
    imf=glob.glob('win'+str(window_num)+'*_'+ang+'*_corrected.tif')[0]
    edgef=glob.glob('win'+str(window_num)+'*_'+ang+'*_corrected_edges.p') #presume these are the 'cleaned' edges
    print imf,edgef
    rperiods,fperiods=[],[]
    badmasks=0
    for im,ef in zip([imf],edgef):
        #print im,ef
        edges=pickle.load(open(ef,'rb'))
        ima=im2ndarray(im)
        #now make mask
        mask=rising_or_falling(ima,edges,np.mean(ima), shape=np.shape(edges),imfile=im,plot=plot)
        #group
        if np.shape(mask)[0] > 0:
            periodr,periodf=group_edges_by_idx(edges,mask,nperiod/pix2um,tolerance=tolerance,mod=3,plot=plot)
            rperiods.extend(periodr)
            fperiods.extend(periodf)
        else:
            badmasks+=1
            print im

    print 'bad masks ', badmasks
            
    #periods=rperiods #need to copy this ... right now it extends rperiods too
    periods=rperiods+fperiods

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
        #figfilename='win'+str(window_num)+'_group_periods'+str(mag)+'.png'
        #fig.savefig(figfilename)
        fig.show()
        
    if stats: #print out and pickle stats
        avg=np.mean(periods)
        med=np.median(periods)
        stdv=np.std(periods)
        #statfile='win'+str(window_num)+'_width_stats'+str(mag)+'.p'
        #datafile='win'+str(window_num)+'_width_data'+str(mag)+'.p'
        print "-------------------STATISTICS FOR WINDOW "+str(window_num)+"---------------------"
        print '              Mean: ' + str(avg)
        print '            Median: ' + str(med)
        print 'Standard Deviation: ' + str(stdv)
        #print 'Results saved in ' + statfile
        data=periods
        stdict={'mean':avg,'median':med,'stddev':stdv}
        #pickle.dump(stdict,open(statfile,'wb'))
        #pickle.dump(data,open(datafile,'wb'))

    return rperiods, fperiods,periods

def find_gaps(tv):
    '''find the gaps in the slats, ie where there are long groups of True values in the true vector. Keep only the first and last value
    Deprecated by use of clean_centers'''
    tv_fixed=[]
    for i,t in enumerate(tv[:-1]):
        if i !=0:
            if tv[i-1] == False and tv[i] == True:
                tv_fixed.append(i)
            elif tv[i] == True and tv[i+1] == False:
                tv_fixed.append(i)                
            #else:
                #tv_fixed.append(False)
        #else:
        #    tv_fixed.append(i)
    return tv_fixed

def slit_or_slat(j,row,imrow,immean,pix2um,n): #for j,row in enumerate(edges): #MPI this! super slow right now....
    tv=np.where(row == True)[0]
    #tv=find_gaps(row)
    tvlen=len(tv)
    #print tvlen
    #try:
    #width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1) if np.mean(imarr[tv[0][i+1]-tv[0][i],j]) < immean] #slats = space between slits
    width=[float(tv[i+1]-tv[i]) for i in range(0,tvlen-1) if np.mean(imrow[tv[i]:tv[i+1]]) > immean] #slats = space between slits . white = 255, black=0
    #width=[i for i in range(0,tvlen-1) ] #slats = space between slits
    #except ValueError:
    #width=[float(tv[0][i+1]-tv[0][i]) for i in range(0,tvshape-1)]
    #    width=[float(tv[1][i+1]-tv[1][i]) for i in range(0,tvshape-1)]
    #print len(width)
    if width != []:
        width=filter(lambda a:a > n,width) #filter out all the 0- n px separations
        #now convert to period. Let's start with some assumptions:
        #1) the width of the slat is about equal to the width of the slit. ideal case: width[15,13,15,13] etc
        #   Then period is simply: [width[i]+width[i+1] for i in range(0,len(width))]
        #   Of course this will probably not be the case. Worst case, there is some noise due to a bump or something:
        #   Or an edge pixel is missing so that the width is already the period
        #   width=[2,2,15,3,8,13,2,30,13]
        #   So we need some conditions. This could cost quite some computing time....we can do it iteratively in steps.
        #   we'll need to use some basic statistics to help...let's trust for now that the rows are long enough that we can do this.
        #   might not be the case with coarser grids!!!
        #sigma=3*np.std(width) # 3 sigma
        median=np.median(width) #assume it gets it right more often than not. might not be the case for the finer grids...
        print median
        period=[]
        for i in range(0,len(width)-1):
            #which is smaller, width % median, width % n*median?
            factor=float(int(width[i])/int(median))
            if factor == 0:
                factor = 1.
            print width[i],factor
            period.append(width[i]/factor)
        period=np.array(period)*pix2um # need to account for the fact that it's rotated! width is not true width!
    else:
        period =[]
    return period,width,tv

def get_slit_width(edges,p0ang=5.0,pix2um=.65,im=False,window_num=False, title='',xran=[80,120],figname=False, stats=True,n=5,show=False):
    '''Currently this will only work for the finest grids, that have well-defined edges without stuff in the middle'''
    import time
    import multiprocessing as mpi
    start=time.time()
    #let's go row-by-row and get the distance between True values, and compile into a histogram:
    widths=[]
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
        for j,(row,imrow) in enumerate(zip(edges,imarr)): #MPI this! super slow right now....
            #imrow=imarr[:,j]
            result=pool.apply_async(slit_or_slat,args=(j,row,imrow,immean,pix2um,n))
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
        statfile='win'+str(window_num)+'_'+p0ang +'_width_stats.p'
        datafile='win'+str(window_num)+'_'+p0ang+'_width_data.p'
        print "-------------------STATISTICS FOR WINDOW "+str(window_num)+"---------------------"
        print '              Mean: ' + str(avg)
        print '            Median: ' + str(med)
        print 'Standard Deviation: ' + str(stdv)
        print 'Results saved in ' + statfile
        data=widths_vec
        stdict={'mean':avg,'median':med,'stddev':stdv}
        pickle.dump(stdict,open(statfile,'wb'))
        pickle.dump(data,open(datafile,'wb'))

    if show:    
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
        #plt.close(fig)
        fig.show()
    return data, stdict #avg

def test_tilt():
    testf=glob.glob('win44_p0*corrected.tif')
    #testf=glob.glob('win11_p0*width_stats.p')
    avgvec=[]
    xvec=[]
    for t in testf:#[:2]+testf[5:8]:
        xvec.append(float(t[9:-14]))
        edges=Canny_edge(t)
        #edges=t[:-4]+'_edges.p'#Canny_edge(t,gauss=3)
        ee=pickle.load(open(edges,'rb'))
        p0ang=t[t.find('_')+1:t.rfind('_')]
        d,st=get_slit_width(ee,window_num=44,p0ang=p0ang,im=np.array(Image.open(t)),xran=[70,110],n=2,title=edges,show=False)
        avgvec.append(st['mean'])
        #avgvec.append(avg)
    #data=zip(xvec,yvec)
    #data.sort(0)
    #for t in testf:
    #    xvec.append(float(t[9:-14]))    
    #    stats=pickle.load(open(t,'rb'))
    #    avgvec.append(stats['mean'])
    data=zip(xvec,avgvec)
    data.sort()
    xdat=[tp[0] for tp in data]
    ydat=[tp[1] for tp in data]
    fig,ax=plt.subplots()
    ax.plot(xdat,ydat,'o-')
    fig.show()
    return xdat,ydat
        
def get_theta_range(nang,side=1.0,spread=5.,n=201):
    #define range of theta around theta_nominal
    theta0= nang*(np.pi/180.)#in radians
    #tendeg2rad=np.pi/18.
    spreaddeg2rad=spread*(np.pi/180.)  
    thetaran = np.linspace(theta0-spreaddeg2rad, theta0+spreaddeg2rad, num=n)#in radians
    return thetaran

def prob_hough(edges, threshold=10, line_length=50, line_gap=2,retlines=False,plot=False,side=1.0,spread=5.,n=201, tag=False,overwrite=False):
    '''Perform probabilistic Hough fit to given set of edges'''
    nang=0
    import glob
    if type(edges) == str: #it's a filename
        inp=edges
        edata=pickle.load(open(edges,'rb'))
    else:
        edata=edges
        #edges=raw_input('What is the window number?')
        inp=raw_input('Output file name?')
    if not overwrite:
        if tag:
            defaultn=inp[:-8]+'_hough_'+tag+'.p'
        else:
            defaultn=inp[:-8]+'_hough.p'
        names=glob.glob(defaultn)
        if names !=[]:
            return

    thetaran=get_theta_range(nang,spread=spread,n=n,side=side)
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
    nang=0.
    if type(edges) == str: #it's a filename
        inp=edges
        edata=pickle.load(open(edges,'rb'))
    else:
        edata=edges
        edges=raw_input('What is the window number?')
        inp=raw_input('Output file name?')

    thetaran=get_theta_range(nang,spread=spread,n=n,side=side)
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

def cat_hough(window_num, tags=['ll150','ll200','ll250','ll300'],weights=False,outtag='all',EM=True,mag=5.0):
    """Cat together the lines """
    #first get indices for the window
    if not EM:
        ofiles=glob.glob('win'+str(window_num)+'_*'+str(mag)+'X*_edges.p')
    else:
        ofiles=EM_list(str(window_num),str(mag),ending='_edges.p')
    idxs=[get_index(o)[1] for o in ofiles]
    #all_lines=[]
    for idx in idxs:
        basen='win'+str(window_num)+'_'+str(idx[0])+'_'+str(idx[1])+'_'+str(mag)+'X_hough'
        lines=pickle.load(open(basen+'.p','rb'))
        #print len(lines)
        all_lines=lines
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
    ang=[0.]
    if type(lines) == list: #it's a list of filenames
        for f in lines:
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
        foo=ax.hist(theta,thetaax)#,normed=True)

    if gaussfit:
        import matplotlib.mlab as mlab
        gauss=mlab.normpdf(thetaax,np.mean(theta),np.sqrt(np.var(theta)))
        ax.plot(thetaax,gauss)
        print np.where(gauss==np.max(gauss)), thetaax[np.where(gauss==np.max(gauss))[0]]
    #if side==1.0:
    #    ang=[aa['nominal angle'] for aa in windows if aa['number']==windownum]
    #else:
    #    ang=[aa['nominal angle'] for aa in windowsr if aa['number']==windownum]            

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
    else:
        ax.set_ylim([0,10])   
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
    if deltax == 0 or deltay == 0:
        #theta=np.arctan(float(deltay)/float(deltax))  #np.arctan2(float(deltay)/float(deltax))
        thetadeg=0.0#np.rad2deg(theta) #theta*180./np.pi
    else:   
        theta=np.arctan(float(deltay)/float(deltax))  #np.arctan2(float(deltay)/float(deltax))
        thetadeg=np.rad2deg(theta) #theta*180./np.pi
    if thetadeg > 89.:
        thetadeg = 90.-thetadeg
    return thetadeg

