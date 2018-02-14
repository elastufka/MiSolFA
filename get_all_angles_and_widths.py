"""
===================
get_all_angles_and_widths.py
Erica  Lastufka 13.10.17
===================

What the title says

"""
import edge_and_hough as eh
import glob
import time
import multiprocessing as mpi
import os

if __name__ == '__main__':
    wins=[11,22,22,31,41,32,42,33,43,34,44]
    ll=[100,150,200,250,300]
    tol=[16,16,4,4,8,8,3,3,5,5,3,3]
    #front window, EM 5.0:
    #print 'here'
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/EMassembly_2017_10_03_RearGrid')
    #os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/mw469sub2634_2017_12_11')
    for k,win in enumerate(wins):
        filen=glob.glob('win'+str(win)+'*_5.0X.tif')#eh.EM_list(str(win),'5.0',ending='.tif')
        filep=eh.EM_list(str(win),'5.0',ending='edges.p')#glob.glob('win'+str(win)+'*_5.0X_edges.p')
        fileh=[]#eh.EM_list(str(win),'5.0',ending='_hough.p')

        #if filep==[] or len(filep) != len(filen): #need to generate the pickles by doing the edges
        #    for f in filen:
        #        filep.append(eh.Canny_edge(f,mag=15.0)) #new list of picklefiles
        #if fileh==[]or len(fileh) != len(filen): #need to generate the pickles by doing the edges
        #for f in filep:
        #    fileh.append(eh.prob_hough(f,line_length=ll[0],spread=2.,overwrite=False)) #new list of picklefiles
        for l in ll[1:]:
            tag='ll'+str(l)
            for f in filep:
                fileh.append(eh.prob_hough(f,line_length=l,spread=2.,tag=tag,overwrite=False)) #new list of picklefiles

        #eh.cat_hough(win,weights=[1,2,3,4,5],mag=15.0,EM=False)
        fileh=glob.glob('win'+str(win)+'*_5.0X_hough_all.p')
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
        eh.hough_hist(fileh,win,figname=figname,tol=tol[k],side=-1.,spread=2.)
            
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
