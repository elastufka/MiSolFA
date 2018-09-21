"""
===================
analyze_optical.py
Erica  Lastufka 20.9.17
===================

Implementations of general Analyze methods specific to optical data. Based on get_all_angles_and_widths.py

Input: Grating object
Output: Result object

"""
import analyze_general as ag
import glob
import time
import multiprocessing as mpi
import os

def prep_all_images(grating,folder=False):
    '''Edge trimming, contrast optimizing, etc. of images '''
    if folder:
        os.chdir(folder)
    else:
        os.chdir(grating.Folder.optical)
    win=grating.win
    filen=glob.glob('win'+str(win)+'*_5.0X.tif')
    return filen,filep,fileh

def get_Canny_edges(fn,mag=5.0,sigma=3,tag='_edges'):
    '''Wrapper for running Canny_edge for the optical case'''
    edges=ag.Canny_edge(fn, mag=mag, outfilen=f[:-4]+tag+".p",sigma=sigma)
    return edges

def run_prob_hough(fp,ll=[50,100,150,200,250]):
    '''Wrapper for running prob_hough for the optical case'''


def cat_hough():
    return hough_all

def write_stats()

def write_results():

def analyze_optical(grating, folder=False, ll=[50,100,150,200,250]):
    '''main method - run everything'''
    filen,filep,fileh=prep_all_images(grating,folder=folder)

    if filep==[] or len(filep) != len(filen): #need to generate the pickles by doing the edges
        filep=[get_Canny_edges(fn) for fn in filen]#.append(eh.Canny_edge(f,mag=5.0,outfilen=f[:-4]+"_ani_edges.p")) #new list of picklefiles
        #if fileh==[]or len(fileh) != len(filen): #need to generate the pickles by doing the edges
        #fileh=[run_prob_hough(fp) for fp in filep]
        #fileh.append(eh.prob_hough(f,line_length=ll[0],spread=2.,overwrite=False)) #new list of picklefiles
        for l in ll:
            tag='ll'+str(l)
            for f in filep:
                fileh.append(eh.prob_hough(f,line_length=l,spread=2.,tag=tag,overwrite=False)) #new list of picklefiles

    hough_all=cat_hough(grating.win,weights=[1,1,1,1,1],tags=['ll50','ll100','ll150','ll200','ll250'],mag=5.0,EM=False)
    #fileh=glob.glob('win'+str(win)+'*_5.0X_hough_all.p')
        #if win in [11,21,32,42] : #use Hough fits to get the widths for coarse windows since Canny edges have false detections
        #    continue
        #    filep=fileh #are these really 1 px thick though? need to convert them to boolean arrays

    #    titlec='Slit angle distribution for window '+str(win) #+ ', nominal pitch '+pitch[k] + ' $\mu$m at 5X magnification'
        #figname='win'+str(win) + 'width_hist.png'
        #if len(filen) !=0:
        #    w,b,f=eh.get_period_by_grouping(win,mosaic=False,EM=False,tolerance=.85,offset=0.,pix2um=.65047,mag=15.0,side=1.)

    figname='win'+str(win) + 'ang_hist.png'
     #   print figname
        #results.append(pool.apply_async(eh.hough_hist, args=(fileh,win,figname)))
    eh.hough_hist(hough_all,win,figname=figname,tol=tol[k],side=1.,spread=2.)

    #moire patterns
    #import process_xray_im as px
    #os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/2018-02-02-Xraytest')
    #for k,win in enumerate(wins):
    #    filen=glob.glob('win'+str(win)+'.tif')
    #    filep,fileh=[],[]
    #    for f in filen:
    #        im=px.im2ndarray(f)
    #        im=px.remove_edges(im,ncrop=20)
    #        im=px.contrast_stretch(im)
    #        edges=px.Canny_edge(im,sigma=4,gauss=4, plot=True)
    #        pname='win'+str(win)+'_edges.p'
    #        pickle.dump(edges,open(pname,'wb'))
    #        filep.append(pname)
    #    for l in ll[1:]:
    #        tag='ll'+str(l)
    #        for f in filep:
    #            fileh.append(eh.prob_hough(f,line_length=l,spread=2.,tag=tag,overwrite=False)) #new list of picklefiles

        #eh.cat_hough(win,weights=[1,2,3,4,5],mag=15.0,EM=False)
    #    fileh=glob.glob('win'+str(win)+'*_hough_all.p')

    #    titlec='Slit angle distribution for window '+str(win) #+ ', nominal pitch '+pitch[k] + ' $\mu$m at 5X magnification'
        #figname='win'+str(win) + 'width_hist.png'
        #if len(filen) !=0:
        #    w,b,f=eh.get_period_by_grouping(win,mosaic=False,EM=False,tolerance=.85,offset=0.,pix2um=.65047,mag=15.0,side=1.)

     #   figname='win'+str(win) + 'ang_hist.png'
     #   print figname
        #results.append(pool.apply_async(eh.hough_hist, args=(fileh,win,figname)))
     #   eh.hough_hist(fileh,win,figname=figname,tol=tol[k],side=-1.,spread=2.)


    #rear window, EM 5.0
    #os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/EMassembly_2017_10_03_RearGrid')
    #for k,win in enumerate(windows):
        #filen=eh.EM_list(str(win),'5.0',ending='.tif')
    #    filep=glob.glob('win'+str(win)+'*_5.0X_bright_edges.p')#[]#eh.EM_list(str(win),'5.0',ending='_edges.p')
    #    fileh=[]#eh.EM_list(str(win),'5.0',ending='_hough.p')

        #if filep==[] or len(filep) != len(filen): #need to generate the pickles by doing the edges
        #for f in filen:
        #    filep.append(eh.Canny_edge(f)) #new list of picklefiles
        #if fileh==[]or len(fileh) != len(filen): #need to generate the pickles by doing the edges
    #    for f in filep:
    #        fileh.append(eh.prob_hough(f,line_length=ll[0],spread=2.)) #new list of picklefiles
    #    for l in ll[1:]:
    #        tag='ll'+str(l)
    #        for f in filep:
    #            fileh.append(eh.prob_hough(f,line_length=l,spread=2.,tag=tag)) #new list of picklefiles

    #    cat_hough(win,weights=[1,2,3,4,5],EM=False)
    #    fileh=glob.glob('win'+str(win)+'*_5.0X_*hough_all.p')
        #if win in [11,21,32,42] : #use Hough fits to get the widths for coarse windows since Canny edges have false detections
        #    continue
        #    filep=fileh #are these really 1 px thick though? need to convert them to boolean arrays

   #     titlec='Slit angle distribution for window '+str(win) #+ ', nominal pitch '+pitch[k] + ' $\mu$m at 5X magnification'
        #figname='win'+str(win) + 'width_hist.png'
        #if len(filen) !=0:
        #    w,b=eh.get_slit_width(filep,mag=5.0,window_num=win,title=titlec,xran=[float(pitch[k])-float(pitch[k]),float(pitch[k])+float(pitch[k])],figname=figname)

    #    figname='win'+str(win) + 'ang_hist.png'
    #    print figname
        #results.append(pool.apply_async(eh.hough_hist, args=(fileh,win,figname)))
    #    eh.hough_hist(fileh,win,figname=figname,tol=tol[k])


    # #front window, EM 15.0
    # os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/EMassembly_2017_10_03_FrontGrid_highmag')
    # for k,win in enumerate(windows):
    #     filen=glob.glob('win'+str(win)+'*_15.0X.tif')
    #     filep=glob.glob('win'+str(win)+'*_15.0X_edges.p')
    #     fileh=glob.glob('win'+str(win)+'*_15.0X_hough.p')

    #     if filep==[] or len(filep) != len(filen): #need to generate the pickles by doing the edges
    #         for f in filen:
    #             filep.append(eh.Canny_edge(f)) #new list of picklefiles
    #     if fileh==[]or len(fileh) != len(filen): #need to generate the pickles by doing the edges
    #         for f in filep:
    #             fileh.append(eh.prob_hough(f)) #new list of picklefiles
    #     #if win in [11,21,32,42] : #use Hough fits to get the widths for coarse windows since Canny edges have false detections
    #     #    continue
    #     #    filep=fileh #are these really 1 px thick though? need to convert them to boolean arrays

    #     titlec='Slit width distribution for window '+str(win) + ', nominal pitch '+pitch[k] + ' $\mu$m at 15X magnification'
    #     figname='win'+str(win) + 'width_hist_15.png'
    #     if len(filen) !=0:
    #         w,b=eh.get_slit_width(filep,mag=15.0,window_num=win,title=titlec,xran=[float(pitch[k])-float(pitch[k]),float(pitch[k])+float(pitch[k])],figname=figname)

    #         figname='win'+str(win) + 'ang_hist_15.png'
    #         print figname
    #     #results.append(pool.apply_async(eh.hough_hist, args=(fileh,win,figname)))
    #         eh.hough_hist(fileh,win,figname=figname,mag=15.0)

    #rear window, EM 15.0
    # os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/EMassembly_2017_10_03_RearGrid')
    # for k,win in enumerate(windows):
    #     filen=glob.glob('win'+str(win)+'*_15.0X.tif')
    #     filep=glob.glob('win'+str(win)+'*_15.0X_edges.p')
    #     fileh=glob.glob('win'+str(win)+'*_15.0X_hough.p')

    #     if filep==[] or len(filep) != len(filen): #need to generate the pickles by doing the edges
    #         for f in filen:
    #             filep.append(eh.Canny_edge(f)) #new list of picklefiles
    #     if fileh==[]or len(fileh) != len(filen): #need to generate the pickles by doing the edges
    #         for f in filep:
    #             fileh.append(eh.prob_hough(f)) #new list of picklefiles
    #     #if win in [11,21,32,42] : #use Hough fits to get the widths for coarse windows since Canny edges have false detections
    #     #    continue
    #     #    filep=fileh #are these really 1 px thick though? need to convert them to boolean arrays

    #     #titlec='Slit width distribution for window '+str(win) + ', nominal pitch '+pitch[k] + ' $\mu$m at 15X magnification'
    #     #figname='win'+str(win) + 'width_hist_15.png'
    #     if len(filen) !=0:
    #         #w,b=eh.get_slit_width(filep,mag=15.0,window_num=win,title=titlec,xran=[float(pitch[k])-float(pitch[k]),float(pitch[k])+float(pitch[k])],figname=figname)

    #         figname='win'+str(win) + 'ang_hist_15.png'
    #         print figname
    #     #results.append(pool.apply_async(eh.hough_hist, args=(fileh,win,figname)))
    #         eh.hough_hist(fileh,win,mag=15.0,figname=figname)

