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
import os
from grating import cbp

def analyze_window(grating,j, method,tolvec=False,sigma=False,angle_only=False,psigma=2.):
    '''Run the whole analysis for a window. Return a Results object'''
    #prep: check for processing level of the files. should be at least level 1 (renamed)
    win=grating.win
    model=grating.model
    folder=grating.data.Xdata.path
    side=grating.side
    nominal=grating.nominal
    data_dict=grating.data.Xdata.transm[j]
    flags=data_dict['flags']
    #os.chdir(data_dict['folder'])

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

#     if 'level04' not in pkeys: #need to groupstretch IF p0-7 exists
#         gs={}
#         for k in data_dict['level02'].keys():
#             gs[k]=ag.contrast_stretch_group(data_dict['level02'][k],saveim=True)
#         updates.append({'key':'level04','data':gs}) #should use the flatdark correct ones for this
#         data_dict['level04']=gs
#         except KeyError: #have to use the 'updates' dict
#             gs = ag.groupstretch(updates[0]['data'])
#             updates.append({'key':'level04','data':ag.groupstretch(updates[0]['data'])}) #should use the flatdark correct ones for this

    if angle_only: #just calculate the relative orientations of line fits to images
        dkeys=[k for k in data_dict['level03'].keys() if k != 'flagged']
        if 'level05' not in pkeys: #need to get edges
            nvec={}
            nlist=[]
            for k in dkeys:
                for im in data_dict['level03'][k]:
                    _,nfn=ag.Canny_edge(im,mag=False)
                    nlist.append(nfn)
                nvec[k]=nlist
            updates.append({'key':'level05', 'data':nvec}) #shouldn't this be sorted by p though?
            data_dict['level05']=nvec
        if 'level06' not in pkeys: #need to do Hough line fits
            nvec={}
            nlist=[]
            for k in dkeys:
                for im in data_dict['level05'][k]:
                    ag.prob_hough(im,0.0,line_length=200,spread=2.,tag='ll200',overwrite=True)
                    nfn=im[:-7]+'hough_ll200.p'#,'rb'
                    nlist.append(nfn)
                nvec[k]=nlist
            updates.append({'key':'level06', 'data':nvec}) #shouldn't this be sorted by p though?
            data_dict['level06']=nvec
        #need to do this one file at a time
        #hdata=np.mean(data_dict['level06']['data'])
        figname='win'+str(win)+'_hough_hist.png'
        hfiles=[]
        for k in dkeys:
            hfiles.extend(data_dict['level06'][k])
        theta=ag.hough_hist(hfiles,win,0.0,figname=figname,spread=2.,mask45=False,ret=True) #make histogram

        #add theta to what's in the logfile...
        theta0=read_logfile(data_dict['folder']+'.log', win)

        try:
            grating.results['widths'].angle_stats={'theta0':theta0,'theta':theta,'theta_std':np.nanstd(theta)}
            return [], grating.results['widths']

        except AttributeError or KeyError:
            pass
        #xpeaks=[]#clean centers, sum peaks
        #theta_dict=ag.angle_by_index(hdata,xpeaks)#check angle vs index
        #res=results.Results(win,'transm','angles',False,{'theta':theta,'theta_dict':theta_dict},params={'tolvec':tolvec, 'sigma':sigma})
        #return updates, res

    kind ={'type':'transm','model':model,'folder':folder,'side':side}

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
        xv,sv=basic_sum_allP(data_dict['level04'],flags=flags) #IF NOT FLAGGED
        #xvec,sumvec,errvec=scat_allP(win,xv,sv)
        tdata={'xdata':xv,'ydata':sv}
        res=results.Results(win,kind,flags,method='sum',data=tdata,nominal=nominal)

    elif method == 'binary':
        xv,sv=basic_sum_allP(data_dict['level02'],binary=True)
        #xvec,sumvec,errvec=scat_allP(win,xv,sv)
        res=results.Results(win,'transm','binary',False,{'xdata':xv,'ydata':sv})
    elif method == 'widths':
        #run batch_window
        xvec,sumvecp,sumveca,sumpp,sumpa=calc_widths_allP(win,data_dict['level03'],flags,tolvec=tolvec,sigma=sigma,rotang=False,psigma=psigma)

        meanvecp,meanpp=[],[]
        maskvecp,maskpp=False,False
        #try:
        #    flags=data_dict['level03']['flagged']
        #    maskvecp,maskpp=[],[]

        #except KeyError:
        #    flags=[]
        for p in range(0,7):
            meanvecp.append([np.mean(sumvecp[i,p]) for i in range(0,len(xvec))])
            meanpp.append([np.mean(sumpp[i,p]) for i in range(0,len(xvec))])
            #build masks too
            #if flags !=[]:
            #    lvecp,lpp=[],[]
            #    fullimnames=['win'+str(win)+'_p'+str(p)+'_'+str('{:g}'.format(x))+'_corrected.tif' for x in xvec]
            #    for im in fullimnames:
            #        if im in flags:
            #            lvecp.append(True)
            #            lpp.append(True)
            #        else:
            #            lvecp.append(False)
            #            lpp.append(False)
            #    maskvecp.append(lvecp)
            #    maskpp.append(lpp)


        meanvecp=np.ma.array(meanvecp,mask=maskvecp)
        meanveca=np.transpose(meanvecp)#np.mean(sumveca,axis=1)
        #meanpp=np.array([[np.mean(svp) for svp in sumvp] for sumvp in sumpp])
        #for a in range(0,len(xvec)):
        #    meanpp.append([np.mean(sumpp[i,a]) for i in range(0,7)])

        meanpp=np.ma.array(meanpp,mask=maskpp)
        meanpa=np.transpose(meanpp)#np.mean(sumveca,axis=1)
        ydata=meanvecp/meanpp #duty cycle, grouped by p
        ydataa=meanveca/meanpa #duty cycle, grouped by angle
        #win,kind,method=False,errors=False,data=False,nominal=False,filenames=False,stats=False,params=False):

        if angle_only:
            tdata={'xdata':xvec,'ydata':ydata,'mean_ang':ydataa, 'raw_widths_P':sumvecp,'raw_widths_A':sumveca,'raw_periods_P':sumpp,'raw_periods_A':sumpa,'mean_widths_P':meanvecp,'mean_widths_A':meanveca,'mean_periods_P':meanpp,'mean_periods_A':meanpa,'theta0':theta0,'theta':theta, 'theta_std':np.nanstd(theta)}
        else:
            tdata={'xdata':xvec,'ydata':ydata,'mean_ang':ydataa, 'raw_widths_P':sumvecp,'raw_widths_A':sumveca,'raw_periods_P':sumpp,'raw_periods_A':sumpa,'mean_widths_P':meanvecp,'mean_widths_A':meanveca,'mean_periods_P':meanpp,'mean_periods_A':meanpa}
        params={'tolvec':tolvec, 'sigma':sigma}
        res=results.Results(win,kind,flags,method='widths',data=tdata,nominal=nominal,params=params)

    return updates,res

def calc_widths_allP(win,data_dict,flags,tolvec=False,sigma=False,rotang=False,psigma=2.0):
    '''should I deal with flagged data here? send the arrays back as masked arrays? '''
    widths,allwidths,allperiods,check=[],[],[],[]
    files=data_dict['p0']
    if "/" in files[0]:
        filest=[f[f.rfind("/"):] for f in files]
    else:
        filest=files
    #flags=data_dict['flagged']
    xvec=[float(im[im.find("_")+4:-14]) for im in filest] #18 if filest is corrected_edges
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
            #im=[fn for fn in fnames if
            im='win'+str(win)+'_p'+str(p)+'_'+ang+'_corrected.tif'#][0]
            if tolvec:
                tol=tolvec[j]
            else:
                tol=6
            #print im
            if rotang: #get mean angle from hough lines
                lines=pickle.load(open('win'+str(win)+'_p'+str(p)+'_'+ang+'_corrected_hough_ll200.p','rb'))
                grouped_lines=ag.same_line(lines,0.0,tol=2)
                for gl in grouped_lines:
                    theta=np.mean([ag.get_angle(g) for g in gl])
                rang=90.-np.mean(theta)
                if rang < 0.25:
                    rang = False
            else:
                rang = False

            pf,pr,pp,ww,ca,ci,ck=ag.slit_widths_from_peaks(win,im,plot=False,tolerance=sigma,n=tol,rang=rang,sigma=psigma) #do even if flagged
            #make sure from_centers is true
            if ck and im[6:-14] not in flags: #and im not in flags:
                check.append(im)
            #if im in flags: #mask?
            #    mask.append() # empty array of same shape?
            allwidths.append(np.array(ww))
            periods_from_centers=[]
            for n,c in enumerate(ca[:-1]):
                periods_from_centers.append(ca[n+1]-c)
            allperiods.append(np.array(periods_from_centers)) #use period calculated from slat centers (this is just the positions of slat centers tho...
            if p ==0:
                print ang,np.mean(periods_from_centers),np.mean(ww)
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

def read_logfile(filename,win):
    ''' transm_step011_window012_     -28.178	 9.157	314.500	 0.000 '''
    with open(filename) as f:
        lines=f.readlines()
    wstr='window0'+str(win)
    for l in lines:
        if l.startswith('transm') and wstr in l:
            ls=l.split('\t')
            theta0=ls[2]
    return float(theta0)

def GWcalc_widths_allP(files,tolvec=False,sigma=False):
    '''should I deal with flagged data here? send the arrays back as masked arrays? '''
    widths,allwidths,allperiods,check=[],[],[],[]

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

    if plot:
        #make the histogram
        if type(periods) == np.ma.core.MaskedArray:
            print 'masked arrays'
            periods=periods.compressed()
            periodr=periodr.compressed()
            periodf=periodf.compressed()
        fig,ax=plt.subplots(1,3,sharex=True,sharey=True)
        bins=np.arange(np.min(periods),np.max(periods),np.std(periods)/2.)
        ax[0].hist(periodr,bins,facecolor='g')
        ax[1].hist(periodf,bins,facecolor='r')
        ax[2].hist(periods,bins)
        ax[0].set_xlim([grating.nominal['pitch']-5,grating.nominal['pitch']+5])
        ax[0].set_title('Rising')
        ax[1].set_title('Falling')
        ax[2].set_title('Total')
        ax[0].set_xlabel('Period $\mu$m')
        ax[0].set_ylabel('Counts')
        ax[0].set_yscale('log')
        ax[0].set_ylim([1,10000])
        figfilename='win'+str(win)+'_group_periods5.0X'+ftags+'.png'
        fig.savefig(figfilename)
        fig.show()

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
        sumvec.append(float(np.nanmean(imarr))) #use nanmean just in case the image has been cropped...none have as far as I know
        if flaglist:
            if im[6:im.rfind('_')] in flaglist:
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
    #print mask
    return xvec,sumvec,mask

def basic_sum_allP(data_dict,binary=False,flags=False):
    sumvecs,masks=[],[]
    for p in range(0,7):
        files=data_dict['p'+str(p)]#glob.glob('win'+win+'_p'+str(p)+'*'+'_groupstretch.tif')
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

        xvec,sv,mask=basic_sum(fsort,flaglist=flags,binary=binary,show=False) #they should be in the right order now
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


def compare_methods(win,sums,widths,wnom,yran=False,figname=False,title=False,pix2um=.65,xvmid=False,excl=False,ex_mid='L'):
    '''compare 3 different analysis methods'''
    #binary,sums and widths are Results objects
    xvec=sums.data['xdata']
    if not xvmid:
        xvmid=xvec.index(0.0)
    xvec=np.array(xvec)
    #binary_ydata=np.array(binary.data['ydata'])
    sum_ydata=sums.data['ydata']
    width_ydata=np.array(widths.data['mean_widths_P'])*pix2um #have to multiply by 2
    #adjust ydata... take averages along ax=0?

    #get averages
    #binary_av=np.mean(binary_ydata,axis=0)
    sum_av=np.nanmean(sum_ydata,axis=0)
    wav=np.nanmean(width_ydata,axis=0)
    #fac=1./np.max(wav)
    width_av=wav#fac*wav#np.mean(width_ydata,axis=0)
    width_ydata=width_ydata

    #get upper bounds at each p
    #binary_max=np.ndarray.max(binary_ydata,axis=0)
    sydt=np.transpose(sum_ydata)
    wydt=np.transpose(width_ydata)
    sum_max=np.array([np.nanmax(y) for y in sydt])
    width_max=[np.nanmax(y) for y in wydt]

    #get lower bounds at each p
    #binary_min=np.ndarray.min(binary_ydata,axis=0)
    sum_min=np.array([np.nanmin(y) for y in sydt])
    width_min=[np.nanmin(y) for y in wydt]

    #exclude points if needed. excl is a mask!
    if type(excl)==list:
        sum_av=np.ma.masked_array(sum_av,mask=excl)
        width_av=np.ma.masked_array(width_av,mask=excl)
        sum_max=np.ma.masked_array(sum_max,mask=excl)
        sum_min=np.ma.masked_array(sum_min,mask=excl)
        width_max=np.ma.masked_array(width_max,mask=excl)
        width_min=np.ma.masked_array(width_min,mask=excl)
        xvec=np.ma.masked_array(xvec,mask=excl)
        #apply mask to xvec too

    #get the total average
    #tot_av=(sum_av+width_av)/2.
    if ex_mid=='L':
        xfplus=xvec[:xvmid]
        xfminus=xvec[xvmid:]
        fitsplus=width_av[:xvmid]
        fitsminus=width_av[xvmid:]

    elif ex_mid=='R': #exclude midpoint from line fit on left
        xfplus=xvec[:xvmid+1]
        xfminus=xvec[xvmid+1:]
        fitsplus=width_av[:xvmid+1]
        fitsminus=width_av[xvmid+1:]

    else: #fit both
        xfplus=xvec[:xvmid]
        xfminus=xvec[xvmid-1:]
        fitsplus=width_av[:xvmid]
        fitsminus=width_av[xvmid-1:]

#     elif ex_mid=='R':
#         pend=xvmid
#         mstart=xvmid+1
#     else: #default exclude from both
#         pend=xvmid-1
#         mstart=xvmid+1



    (meanplusslope,meanplusintercept),perr,_,_,_=np.ma.polyfit(xfplus,fitsplus,1,full=True)
    (meanminusslope,meanminusintercept),merr,_,_,_=np.ma.polyfit(xfminus,fitsminus,1,full=True)
    lineintx,lineinty=ag.get_intersect(meanplusslope,meanminusslope,meanplusintercept,meanminusintercept)

    #now plot the differences from the average
    xfminus=np.insert(xfminus,0,xvec[xvmid-1],axis=0)
    if ex_mid =='R': #ad xvmid to this one
        xfminus=np.insert(xfminus,1,xvec[xvmid],axis=0)
    if ex_mid == 'L':
        xfplus=np.insert(xfplus,len(xfplus),xvec[xvmid],axis=0)
    if ex_mid == 'B':
        xfminus=np.insert(xfminus,1,xvec[xvmid],axis=0)
        xfplus=np.insert(xfplus,len(xfplus),xvec[xvmid],axis=0)
    xfplus=np.insert(xfplus,len(xfplus),xvec[xvmid+1],axis=0)
    meanminusline=meanminusslope*xfminus + meanminusintercept
    meanplusline=meanplusslope*xfplus+meanplusintercept
    #print xvec
    #print xfplus,xfminus

    if ex_mid =="R":
        mp=list(meanplusline[:-1])
        mp.extend(list(meanminusline[2:]))
    elif ex_mid =="B": #include midpoint in both
        mp=list(meanplusline[:-3])
        mp.extend(list(meanminusline[2:]))
    else:
        mp=list(meanplusline[:-2])
        mp.extend(list(meanminusline[1:]))
    mpa=np.array(mp)

    ddict={'xvec':xvec,'xfplus':xfplus,'xfminus':xfminus,'width_av':width_av,'width_min':width_min,'width_max':width_max,'fitsplus':fitsplus,'fitsminus':fitsminus,'sum_av':sum_av,'sum_min':sum_min,'sum_max':sum_max,'mpa':mpa}
    fdict={'plus_slope':meanplusslope,'minus_slope':meanminusslope,'plus_line':meanplusline,'minus_line':meanminusline,"int_x":lineintx,'int_y':lineinty,'xvmid':xvmid,'excl':excl,'xvmid':xvmid,'plus_err':perr,'minus_err':merr}

    compare_methods_plot(win,ddict,fdict,yran=yran,title=title,figname=figname)

    return ddict,fdict

def compare_methods_plot(win,ddict,fdict,yran=False,title=False,figname=False):

    fig,ax=plt.subplots(2,sharex=True,figsize=[10,8])
    ax_sum = ax[0].twinx()
    ax_sum.set_zorder(0)
    ax[0].set_zorder(1)
    ax_diff = ax[1].twinx()

    #ax=fig.add_subplot(112)
    #ax[0].fill_between(xvec,binary_min,binary_max,alpha=.6,color='c',)
    ax_sum.fill_between(ddict['xvec'],ddict['sum_min'],ddict['sum_max'],alpha=.6,color=cbp[1],label='sum')
    ax[0].fill_between(ddict['xvec'],ddict['sum_min']*fdict['int_y'],ddict['sum_max']*fdict['int_y'],alpha=.6,color=cbp[1],label='sum')
 #ax[0].fill_between(xvec[:1],width_av[:1],width_av[:1],alpha=.6,color=cbp[1],label='sum')
    #ax[0].fill_between(xvec,width_min,width_max,alpha=.6,color='y',)

    #ax[0].plot(xvec,binary_av,'-o',color='c',label='binary')
    #ax[0].plot(xvec,sum_av,'-o',color='m',label='sum')
    ax[0].plot(ddict['xvec'],ddict['width_av'],'-o',color=cbp[0],linewidth=2,label='widths')
    #ax[0].plot(xvec,tot_av,'--',color='g',label='widths')

    ax[0].axvline(color='k',linestyle='dashed')
    if not yran:
        ax[0].set_ylim([np.min(ddict['width_min']),np.max(ddict['width_max'])])
        ax_sum.set_ylim([np.min(ddict['width_min'])/fdict['int_y'],np.max(ddict['width_max'])/fdict['int_y']])
        #ax_sum.set_ylim([.1,1.1])
    else:
        ax[0].set_ylim(yran)
        ax_sum.set_ylim([yran[0]*100./fdict['int_y'],yran[1]*100./fdict['int_y']])
        #what does this mean for the sum?

    #meanplusslope=np.mean([fp[0] for fp in fitsplus])
    #meanplusintercept=np.mean([fp[1] for fp in fitsplus])
        #print meanminusline
    ax[0].plot(ddict['xfplus'],fdict['plus_line'],'k--')
    ax[0].plot(ddict['xfminus'],fdict['minus_line'],'k--')
    ymin=np.ma.min(ddict['width_av'])
    ymax=np.ma.max(ddict['width_av'])
    xran=ax[0].get_xlim()
    xmin=float(xran[0])
    xmax=float(xran[1])
    #ax[0].text(.75*xmin,.9*ymax,'avg. slope: ' + str(np.round(fdict['plus_slope'],3)))
    ##ax[0].text(.75*xmin,.1*ymax,'avg. intercept: ' + str(np.round(meanplusintercept,3)))
    #ax[0].text(.5*xmax,.9*ymax,'avg. slope: ' + str(np.round(fdict['minus_slope'],3)))
    ##ax[0].text(.5*xmax,.1*ymax,'avg. intercept: ' + str(np.round(meanminusintercept,3)))
    #ax[0].text(.15*xmin,ymin,'adj. peak: (' + str(np.round(fdict['int_x'],3))+ ',' + str(np.round(fdict['int_y'],3))+')')


    #ax[1].fill_between(xvec,binary_min-mpa,binary_max-mpa,alpha=.6,color='c',)
    ax[1].fill_between(ddict['xvec'],(ddict['sum_min']*fdict['int_y'])-ddict['mpa'],(ddict['sum_max']*fdict['int_y'])-ddict['mpa'],alpha=.6,color=cbp[1],)
    ax[1].fill_between(ddict['xvec'],ddict['width_min']-ddict['mpa'],ddict['width_max']-ddict['mpa'],alpha=.6,color=cbp[0],)

    #ax[1].plot(xvec,binary_av-mpa,'-o',color='c')
    ax[1].plot(ddict['xvec'],(ddict['sum_av']*fdict['int_y'])-ddict['mpa'],'-o',color=cbp[1],linewidth=2)
    ax[1].plot(ddict['xvec'],ddict['width_av']-ddict['mpa'],'-o',color=cbp[0],linewidth=2)
    ax[1].axhline(color='k',linestyle='dashed')

    #ax[1].set_ylim([np.min(errvecs),np.max(errvecs)])

    #re-adjust plot sizes
    box1=ax[0].get_position()
    box2=ax[1].get_position()
    ax[0].set_position([box1.x0,box1.y0*.7,box1.width,box1.height*1.5])
    ax_sum.set_position([box1.x0,box1.y0*.7,box1.width,box1.height*1.5])
    ax[1].set_position([box2.x0,box2.y0,box2.width,box2.height*.7])
    ax_diff.set_position([box2.x0,box2.y0,box2.width,box2.height*.7])

    ax[1].set_xlabel('Angle (degrees)') #should I set the lables to the actual angle values? xaxlabels... but they would be in the wrong place
    #     if widths:
    #        ax[0].set_title('Window '+ win + ' slit width as a function of angle')
    #        ax[0].set_ylabel('Slit width ($\mu$m)')
    #    else:
    ax[1].set_ylabel('Difference from fit ($\mu$m)')
    ax_diff.set_ylabel('Difference from fit (%)')
    yrand=ax[1].get_ylim()
    ax_diff.set_ylim([yrand[0]*100./fdict['int_y'],yrand[1]*100./fdict['int_y']])
    if title:
        ax[0].set_title(title+' Window '+ str(win))
    else:
        ax[0].set_title('Window '+ str(win))
    ax_sum.set_ylabel('Percent of maximum transmission')
    ax[0].set_ylabel('Slit Width ($\mu$m)')
    ax[0].legend(loc='upper right')

#     fig.canvas.draw()
#     sticks=[float(t.get_text) for t in ax_sum.get_yticklabels()]
#     new_sticks=sticks/lineinty
#     ax_sum.set_yticklabels(new_sticks)

    if figname:
        plt.savefig(figname)
    else:
        fig.show()
