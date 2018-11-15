"""
===================
analyze_xray.py
Erica  Lastufka 5.9.18
===================

Implementations of general Analyze methods specific to X-ray lamp or beam data
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
import analyze_general as ag
import results

def analyze_window(win,model,data_dict, method,tolvec=False,sigma=False):
    '''Run the whole analysis for a window. Return a Results object'''
    #prep: check for processing level of the files. should be at least level 1 (renamed)
    updates=[]
    pkeys=data_dict.keys() #i is the number of folders - for ex. there might have been a correction directory
    if 'level02' not in pkeys: #need to flatdark correct
        flat=pickle.load(open('../transm'+model+'_combined_flat.p','rb'))
        dark=pickle.load(open('../transm'+model+'_combined_dark.p','rb'))
        newdd={}
        for k in data_dict['level01'].keys():
            newdd[k]=ag.flatdark(data_dict['level01'][k],dark,flat)
        updates.append({'key':'level02', 'data':newdd})
        data_dict['level02']=newdd

    if 'level03' not in pkeys: #need to use flatdark corrected
        cstretch={}
        for k in data_dict['level02'].keys():
            cl=[]
            for im in data_dict['level02'][k]:
                cl.append(ag.contrast_stretch(im))
            cstretch[k]=cl
        updates.append({'key':'level03', 'data':cstretch})
        data_dict['level03']=cstretch

    if 'level04' not in pkeys: #need to groupstretch
        gs={}
        for k in data_dict['level02'].keys():
            gs[k]=ag.contrast_stretch_group(data_dict['level02'][k],saveim=True)
        updates.append({'key':'level04','data':gs}) #should use the flatdark correct ones for this
        data_dict['level04']=gs
#         except KeyError: #have to use the 'updates' dict
#             gs = ag.groupstretch(updates[0]['data'])
#             updates.append({'key':'level04','data':ag.groupstretch(updates[0]['data'])}) #should use the flatdark correct ones for this

    if method =='widths': #need edges and such
        if 'level05' not in pkeys: #need to edge
            nvec={}
            nlist=[]
            for k in data_dict['level03'].keys():
                for im in data_dict['level03'][k]:
                    _,nfn=ag.Canny_edge(im,mag=False)
                    nlist.append(nfn)
                nvec[k]=nlist
            updates.append({'key':'level05', 'data':nvec}) #shouldn't this be sorted by p though?
            data_dict['level05']=nvec
        #elif 'level06' not in pkeys: #need to do line fits
        #    try:
        #        updates.append({'key':'level06', 'data':ag.prob_hough(data_dict['level05'])}) #needs more args
        #    except KeyError:
        #        updates.append({'key':'level06', 'data':ag.prob_hough(updates[-1]['data'])})


    #now execute the chosen method
    if method == 'sum':
        xv,sv=basic_sum_allP(data_dict['level04'])
        #xvec,sumvec,errvec=scat_allP(win,xv,sv)
        #    def __init__(self,kind,method, errors, data,filenames=False,stats=False):

        res=results.Results(win,'transm','sum',False,{'xdata':xv,'ydata':sv})
    elif method == 'binary':
        xv,sv=basic_sum_allP(data_dict['level02'],binary=True)
        #xvec,sumvec,errvec=scat_allP(win,xv,sv)
        res=results.Results(win,'transm','binary',False,{'xdata':xv,'ydata':sv})
    elif method == 'widths':
        #run batch_window
        if tolvec and sigma:
            xvec,sumvecp,sumveca,sumpp,sumpa=calc_widths_allP(win,data_dict['level03'],tolvec=tolvec,sigma=sigma)
        else:
            xvec,sumvecp,sumveca,sumpp,sumpa=calc_widths_allP(win,data_dict['level03'])

        meanvecp,meanpp=[],[]
        maskvecp,maskpp=False,False
        try:
            flags=data_dict['level03']['flagged']
            maskvecp,maskpp=[],[]

        except KeyError:
            flags=[]
        for p in range(0,7):
            meanvecp.append([np.mean(sumvecp[i,p]) for i in range(0,len(xvec))])
            meanpp.append([np.mean(sumpp[i,p]) for i in range(0,len(xvec))])
            #build masks too
            if flags !=[]:
                lvecp,lpp=[],[]
                fullimnames=['win'+str(win)+'_p'+str(p)+'_'+str('{:g}'.format(x))+'_corrected.tif' for x in xvec]
                for im in fullimnames:
                    if im in flags:
                        lvecp.append(True)
                        lpp.append(True)
                    else:
                        lvecp.append(False)
                        lpp.append(False)
                maskvecp.append(lvecp)
                maskpp.append(lpp)


        meanvecp=np.ma.array(meanvecp,mask=maskvecp)
        meanveca=np.transpose(meanvecp)#np.mean(sumveca,axis=1)
        #meanpp=np.array([[np.mean(svp) for svp in sumvp] for sumvp in sumpp])
        #for a in range(0,len(xvec)):
        #    meanpp.append([np.mean(sumpp[i,a]) for i in range(0,7)])

        meanpp=np.ma.array(meanpp,mask=maskpp)
        meanpa=np.transpose(meanpp)#np.mean(sumveca,axis=1)
        ydata=meanvecp/meanpp #duty cycle, grouped by p
        ydataa=meanveca/meanpa #duty cycle, grouped by angle
        res=results.Results(win,'transm','widths',False,{'xdata':xvec,'ydata':ydata,'mean_ang':ydataa, \
                                                     'raw_widths_P':sumvecp,'raw_widths_A':sumveca, \
                                                    'raw_periods_P':sumpp,'raw_periods_A':sumpa, \
                                                    'mean_widths_P':meanvecp,'mean_widths_A':meanveca, \
                                                    'mean_periods_P':meanpp,'mean_periods_A':meanpa},params={'tolvec':tolvec, 'sigma':sigma})

    return updates,res

def calc_widths_allP(win,data_dict,tolvec=False,sigma=False):
    '''should I deal with flagged data here? send the arrays back as masked arrays? '''
    widths,allwidths,allperiods,check=[],[],[],[]
    files=data_dict['p0']
    #flags=data_dict['flagged']
    xvec=[float(im[im.find("_")+4:-14]) for im in files]
    xvec.sort()
    xvecstr=[str(xv) for xv in xvec]
    for k,xv in enumerate(xvecstr):
        if xv.endswith('.0'):
            xvecstr[k]=xv[:-2]
    #allwidthsp,allwidthsa=[],[]
    #print xvecstr
    for i,p in enumerate(range(0,7)):
        for j,ang in enumerate(xvecstr):
            fnames=data_dict['p'+str(p)]
            im=[fn for fn in fnames if fn=='win'+str(win)+'_p'+str(p)+'_'+ang+'_corrected.tif'][0]
            if tolvec:
                tol=tolvec[j]
            else:
                tol=6
            #print im
            pp,ww,ck=ag.slit_widths_from_peaks(win,im,plot=False,tolerance=sigma,n=tol) #do even if flagged
            if ck: #and im not in flags:
                check.append(im)
            #if im in flags: #mask?
            #    mask.append() # empty array of same shape?
            allwidths.append(np.array(ww))
            allperiods.append(np.array(pp))
    #print np.shape(allwidths),np.shape(allperiods)
    #sort into grouped by P and grouped by angle
    #allwidthsp=allwidths
    allwidthsp=[allwidths[i::len(xvec)] for i in range (0,len(xvec))]#[allwidths[i::7] for i in range(0,7)]
    allwidthsa=[allwidths[i::7] for i in range (0,7)]
    allperiodsp=[allperiods[i::len(xvec)] for i in range (0,len(xvec))]
    allperiodsa=[allperiods[i::7] for i in range (0,7)]
    print "Check these images: "
    print check
    return xvec,np.array(allwidthsp),np.array(allwidthsa),np.array(allperiodsp),np.array(allperiodsa)#,errvec

def check_single_image(imname,tolerance=False,sigma=2.):
    win=int(imname[3:5])
    p0=int(imname[7])
    ang=imname[9:-14]
    print win,p0,ang
    edges=pickle.load(open(imname[:-4]+'_edges.p','rb'))
    ce=ag.clean_centers(edges,plot=False,sigma=sigma)
    ag.plot_centers_and_edges(win,p0,ang,tol=sigma)
    recalc =True
    while recalc: #do center cleaning and slit or slat again
        tol=int(raw_input('Tolarance? (int) '))
        sig=float(raw_input('Sigma? (float) '))
        ce=ag.clean_centers(edges,plot=False,sigma=sig)
        xp=ag.sum_peaks(ce,tol=tol)
        pp,ww,ck=ag.slit_widths_from_peaks(win, imname, xpeak=xp,plot=True, tolerance=sig,n=tol,quiet=False)
        ag.plot_centers_and_edges(win,p0,ang,tol=sig)
        pp.sort()
        ww.sort()
        print "10 smallest periods: ", pp[:10]
        print "10 largest periods: ", pp[-10:]
        print "10 smallest widths: ", ww[:10]
        print "10 largest widths: ", ww[-10:]
        rc =(raw_input('Recalculate?'))
        if rc == 1 or rc == 'True' or rc == 'true' or rc == 'Y':
            recalc=True
        else:
            recalc=False

def rename_files(win,flist,logfile=False):
    '''rename from default TOMCAT names to something more precise'''
    #read the log file until the correct window is found
    import os
    if not logfile:
        logfile=glob.glob('*.log')[0]
    with open(logfile) as f:
        lines=f.readlines()
    for l in lines:
        if 'window0'+str(win) in l: #start paying attention after this
            pre=l[:l.rfind('_')] #ex: transmQM_step000_win001
            idx=lines.index(l)

    newnames=[]
    nim=len(flist)/7
    flist.sort()
    for i,f in enumerate(flist):
        step,x,_,_,ang=lines[idx+i+1].split('\t')
        newname='win'+str(win)+'_p'+str(int(step)/nim)+'_'+ang.strip()[:-1]+'.tif' #need to round up! or do this post facto
        #print f,newname
        newnames.append(newname)
        os.rename(f,newname)
    return flist,newnames

def basic_sum(imlist,flaglist=False,binary=False,groupstretch=True,show=True):
    '''make a plot of the sum of the image as a function of angle'''
    sumvec=[]
    xvec=[]
    mask=[]
    for im in imlist:
        imarr=ag.im2ndarray(im)
        if binary: #end tag is flatdark
            xvec.append(float(im[im.find("_")+4:-13]))
            #convert each image to binary representation before summing. Use the mean as the cutoff
            bmean=np.mean(imarr)
            imgt=np.ma.masked_greater(imarr,bmean) #what happens if the value equals the mean?
            imgtfilled=np.ma.filled(imgt, 1.0)
            imlt=np.ma.masked_less(imgtfilled,bmean)
            imltfilled=np.ma.filled(imlt, 0.0)
            imarr=imltfilled
        elif groupstretch:
            xvec.append(float(im[im.find("_")+4:-17]))
        sumvec.append(float(np.sum(imarr)))
        if flaglist:
            if im in flaglist:
                mask.append(True)
            else:
                mask.append(False)
        else:
            mask.append(False)

    #if flaglist:
    #    sumvec=np.ma.array(sumvec,mask=mask)
    #else:
    #    sumvec=np.ma.array(sumvec) #nothing is masked

    if show:
        fig,ax=plt.subplots()
        ax.scatter(xvec,sumvec/np.max(sumvec),s=25)
        #ax.set_ylim([0,1])
        ax.set_xlabel('Angle')
        fig.show()
    return xvec,sumvec,mask

def basic_sum_allP(data_dict,binary=False):
    sumvecs,masks=[],[]
    for p in range(0,7):
        files=data_dict['p'+str(p)]#glob.glob('win'+win+'_p'+str(p)+'*'+'_groupstretch.tif')
        try:
            flagged_files=data_dict['flagged']
        except KeyError:
            flagged_files=False
        #first fix order so can avoid that later...
        if binary:
            xv=[float(im[im.find("_")+4:-13]) for im in files]
        else:
            xv=[float(im[im.find("_")+4:-17]) for im in files]
        xvmid=xv.index(0.0)
        aa=zip(xv,files)
        aa.sort()
        xv.sort()
        fsort=[b[1] for b in aa]

        xvec,sv,mask=basic_sum(fsort,flagged_files,binary=binary,show=False) #they should be in the right order now
        #fix order for plotting connecting line
        sumvecs.append(sv/np.max(sv))
        masks.append(mask)
        #sumvecs.append(svsort) #want this to be a masked array...

    svarray=np.ma.array(np.array(sumvecs),mask=np.array(masks))
    return xvec,svarray

def widths_allP(win,theta_ran=False):
    '''operates on stat files ... '''
    sumvecs=[]
    errvecs=[]
    allfitsplus,allfitsminus=[],[]
    for p in range(0,7):
        files=glob.glob('win'+win+'_width_stats_p'+str(p)+'*.p') # don't be dumb this is the sum of all p's ...
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
    return xv,sumvecs,errvecs,allfitsplus,allfitsminus

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

def combine_width_data(window_num, ang, stats=True,save=False,plot=True,n=20):
    dfiles=glob.glob('win'+str(window_num)+'_width_data*_'+str(ang)+'.p')
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
        statfile='win'+str(window_num)+'_allwidth_stats_'+str(ang)+'.p'
        datafile='win'+str(window_num)+'_allwidth_data_'+str(ang)+'.p'
        pickle.dump(stdict,open(statfile,'wb'))
        pickle.dump(data,open(datafile,'wb'))
        print 'Results saved in ' + statfile
    if plot:
        fit_gaussian_to_hist(allwidths,win=window_num,ang=ang,n=n)

    return allwidths

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


def compare_methods(win,binary,sums,widths,yran=False,figname=False):
    '''compare 3 different analysis methods'''
    #binary,sums and widths are Results objects
    xvec=binary.data['xdata']
    xvmid=xvec.index(0.0)
    binary_ydata=np.array(binary.data['ydata'])
    sum_ydata=np.array(sums.data['ydata'])
    width_ydata=np.array(widths.data['ydata']) #have to multiply by 2

    #get averages
    binary_av=np.mean(binary_ydata,axis=0)
    sum_av=np.mean(sum_ydata,axis=0)
    wav=np.mean(width_ydata,axis=0)
    fac=1./np.max(wav)
    width_av=fac*wav#np.mean(width_ydata,axis=0)
    width_ydata=fac*width_ydata

    #get the total average
    tot_av=(binary_av+sum_av+width_av)/3.
    fitsplus=tot_av[:len(tot_av)/2 +1]
    fitsminus=tot_av[len(tot_av)/2:]

    #get upper bounds at each p
    binary_max=np.ndarray.max(binary_ydata,axis=0)
    sum_max=np.ndarray.max(sum_ydata,axis=0)
    width_max=np.ndarray.max(width_ydata,axis=0)

    #get lower bounds at each p
    binary_min=np.ndarray.min(binary_ydata,axis=0)
    sum_min=np.ndarray.min(sum_ydata,axis=0)
    width_min=np.ndarray.min(width_ydata,axis=0)

    fig,ax=plt.subplots(2,sharex=True,figsize=[10,8])
    #ax=fig.add_subplot(112)
    ax[0].fill_between(xvec,binary_min,binary_max,alpha=.6,color='c',)
    ax[0].fill_between(xvec,sum_min,sum_max,alpha=.6,color='m',)
    ax[0].fill_between(xvec,width_min,width_max,alpha=.6,color='y',)

    ax[0].plot(xvec,binary_av,'-o',color='c',label='binary')
    ax[0].plot(xvec,sum_av,'-o',color='m',label='sum')
    ax[0].plot(xvec,width_av,'-o',color='y',label='widths')
    #ax[0].plot(xvec,tot_av,'--',color='g',label='widths')

    ax[0].axvline(color='k',linestyle='dashed')
    if not yran:
        ax[0].set_ylim([.1,1.1])
    else:
        ax[0].set_ylim(yran)

    #meanplusslope=np.mean([fp[0] for fp in fitsplus])
    #meanplusintercept=np.mean([fp[1] for fp in fitsplus])
    meanplusslope,meanplusintercept=np.polyfit(np.array(xvec[:xvmid+1]),fitsplus,1)#meanplusslope*np.array(xv[:xvmid+1]) + meanplusintercept
    meanplusline=meanplusslope*np.array(xvec[:xvmid+1])+meanplusintercept
    #print meanplusline
    #meanminusslope=np.mean([fm[0] for fm in fitsminus])
    #meanminusintercept=np.mean([fm[1] for fm in fitsminus])
    meanminusslope,meanminusintercept=np.polyfit(np.array(xvec[xvmid:]),fitsminus,1)
    meanminusline=meanminusslope*np.array(xvec[xvmid:]) + meanminusintercept
    lineintx,lineinty=ag.get_intersect(meanplusslope,meanminusslope,meanplusintercept,meanminusintercept)
    #print meanminusline
    ax[0].plot(xvec[:xvmid+1],meanplusline,'k--')
    ax[0].plot(xvec[xvmid:],meanminusline,'k--')
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

    #now plot the differences from the average
    mp=list(meanplusline)[:-1]
    mp.extend(list(meanminusline))
    mpa=np.array(mp)
    ax[1].fill_between(xvec,binary_min-mpa,binary_max-mpa,alpha=.6,color='c',)
    ax[1].fill_between(xvec,sum_min-mpa,sum_max-mpa,alpha=.6,color='m',)
    ax[1].fill_between(xvec,width_min-mpa,width_max-mpa,alpha=.6,color='y',)

    ax[1].plot(xvec,binary_av-mpa,'-o',color='c')
    ax[1].plot(xvec,sum_av-mpa,'-o',color='m')
    ax[1].plot(xvec,width_av-mpa,'-o',color='y')
    ax[1].axhline(color='k',linestyle='dashed')

    #ax[1].set_ylim([np.min(errvecs),np.max(errvecs)])

    #re-adjust plot sizes
    box1=ax[0].get_position()
    box2=ax[1].get_position()
    ax[0].set_position([box1.x0,box1.y0*.7,box1.width,box1.height*1.5])
    ax[1].set_position([box2.x0,box2.y0,box2.width,box2.height*.7])

    ax[1].set_xlabel('tan(Angle)') #should I set the lables to the actual angle values? xaxlabels... but they would be in the wrong place
    #     if widths:
    #        ax[0].set_title('Window '+ win + ' slit width as a function of angle')
    #        ax[0].set_ylabel('Slit width ($\mu$m)')
    #    else:
    ax[0].set_title('Window '+ str(win))
    ax[0].set_ylabel('Percent of maximum transmission')
    ax[0].legend(loc='upper right')
    if figname:
        plt.savefig(figname)
    else:
        fig.show()


