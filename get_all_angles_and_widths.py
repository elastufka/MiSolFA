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
    windows=[34,44]#,31,41]#[11,12,21,22,31,32,33,34,41,42,43,44]
    pitch=['20','20']#,'60','60']#['90','22.5','90','22.5','45','18','30','15','45','18','30','15']
    ang=['45.10','-45.03']#,'-45.138','44.862']#['-44.79','44.95','45.21','-45.10','-45.10','45.10','-44.93','44.97','44.90','-44.96','45.10','-45.03']
    #pool=mpi.Pool()
    #results=[]
    #front window, prototypes::
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2737_2017_06_01')
    for k,win in enumerate(windows):
        filen=glob.glob('win'+str(win)+'*_bright.tif')
        filep=[]#glob.glob('win'+str(win)+'*_bright_edges.p')
        fileh=[]#glob.glob('win'+str(win)+'*_bright_hough.p')
        #filen=eh.highres_list(str(win),ending='.tif')
        #filep=eh.highres_list(str(win),ending='_edges.p')
        #fileh=eh.highres_list(str(win),ending='_hough.p')

        if filep==[] or len(filep) != len(filen): #need to generate the pickles by doing the edges
            for f in filen:
                filep.append(eh.Canny_edge(f)) #new list of picklefiles
        if fileh==[]or len(fileh) != len(filen): #need to generate the pickles by doing the edges
            for f in filep:
                fileh.append(eh.prob_hough(f)) #new list of picklefiles
        #if win in [11,21,32,42] : #use Hough fits to get the widths for coarse windows since Canny edges have false detections
        #    continue
        #    filep=fileh #are these really 1 px thick though? need to convert them to boolean arrays

        titlec='Slit width distribution for window '+str(win) + ', nominal pitch '+pitch[k] + ' $\mu$m at 5X magnification'
        figname='win'+str(win) + 'width_hist.png'
        if len(filen) !=0:
            w,b=eh.get_slit_width(filep,mag=5.0,window_num=win,title=titlec,xran=[float(pitch[k])-float(pitch[k]),float(pitch[k])+float(pitch[k])],figname=figname)
        
            figname='win'+str(win) + 'ang_hist.png'
            print figname
        #results.append(pool.apply_async(eh.hough_hist, args=(fileh,win,figname)))
            eh.hough_hist(fileh,win,figname=figname)
        #results[0].get()

    #rear window, prototypes:
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02')
    for k,win in enumerate(windows):
        filen=glob.glob('win'+str(win)+'*_bright.tif')
        filep=[]#glob.glob('win'+str(win)+'*_bright_edges.p')
        fileh=[]#glob.glob('win'+str(win)+'*_bright_hough.p')
        #filen=eh.highres_list(str(win),ending='.tif')
        #filep=eh.highres_list(str(win),ending='_edges.p')
        #fileh=eh.highres_list(str(win),ending='_hough.p')

        #if filep==[] or len(filep) != len(filen): #need to generate the pickles by doing the edges
        for f in filen:
            filep.append(eh.Canny_edge(f)) #new list of picklefiles
        #if fileh==[]or len(fileh) != len(filen): #need to generate the pickles by doing the edges
        for f in filep:
            fileh.append(eh.prob_hough(f)) #new list of picklefiles
        #if win in [11,21,32,42] : #use Hough fits to get the widths for coarse windows since Canny edges have false detections
        #    continue
        #    filep=fileh #are these really 1 px thick though? need to convert them to boolean arrays

        titlec='Slit width distribution for window '+str(win) + ', nominal pitch '+pitch[k] + ' $\mu$m at 5X magnification'
        figname='win'+str(win) + 'width_hist.png'
        if len(filen) !=0:
            w,b=eh.get_slit_width(filep,mag=5.0,window_num=win,title=titlec,xran=[float(pitch[k])-float(pitch[k]),float(pitch[k])+float(pitch[k])],figname=figname)
        
            figname='win'+str(win) + 'ang_hist.png'
            print figname
        #results.append(pool.apply_async(eh.hough_hist, args=(fileh,win,figname)))
            eh.hough_hist(fileh,win,figname=figname)

    #front window, EM 5.0:
    #print 'here'
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/EMassembly_2017_10_02_FrontGrid')
    for k,win in enumerate(windows):
        filen=eh.EM_list(str(win),'5.0',ending='.tif')
        filep=[]#eh.EM_list(str(win),'5.0',ending='_edges.p')
        fileh=[]#eh.EM_list(str(win),'5.0',ending='_hough.p')

        #if filep==[] or len(filep) != len(filen): #need to generate the pickles by doing the edges
        for f in filen:
            filep.append(eh.Canny_edge(f)) #new list of picklefiles
        #if fileh==[]or len(fileh) != len(filen): #need to generate the pickles by doing the edges
        for f in filep:
            fileh.append(eh.prob_hough(f)) #new list of picklefiles
        #if win in [11,21,32,42] : #use Hough fits to get the widths for coarse windows since Canny edges have false detections
        #    continue
        #    filep=fileh #are these really 1 px thick though? need to convert them to boolean arrays

        titlec='Slit width distribution for window '+str(win) + ', nominal pitch '+pitch[k] + ' $\mu$m at 5X magnification'
        figname='win'+str(win) + 'width_hist.png'
        if len(filen) !=0:
            w,b=eh.get_slit_width(filep,mag=5.0,window_num=win,title=titlec,xran=[float(pitch[k])-float(pitch[k]),float(pitch[k])+float(pitch[k])],figname=figname)
        
            figname='win'+str(win) + 'ang_hist.png'
            print figname
        #results.append(pool.apply_async(eh.hough_hist, args=(fileh,win,figname)))
            eh.hough_hist(fileh,win,figname=figname)
            
    #rear window, EM 5.0
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/EMassembly_2017_10_03_RearGrid')
    for k,win in enumerate(windows):
        filen=eh.EM_list(str(win),'5.0',ending='.tif')
        filep=[]#eh.EM_list(str(win),'5.0',ending='_edges.p')
        fileh=[]#eh.EM_list(str(win),'5.0',ending='_hough.p')
        #fileh=[]
        #if filep==[] or len(filep) != len(filen): #need to generate the pickles by doing the edges
        for f in filen:
            filep.append(eh.Canny_edge(f)) #new list of picklefiles
        #if fileh==[]or len(fileh) != len(filen): #need to generate the pickles by doing the edges
        for f in filep:
            fileh.append(eh.prob_hough(f,side=-1.0)) #new list of picklefiles
        
        #if win in [11,21,32,42] : #use Hough fits to get the widths for coarse windows since Canny edges have false detections
        #    continue
        #    filep=fileh #are these really 1 px thick though? need to convert them to boolean arrays

        titlec='Slit width distribution for window '+str(win) + ', nominal pitch '+pitch[k] + ' $\mu$m at 5X magnification'
        figname='win'+str(win) + 'width_hist.png'
        if len(filen) !=0:
            w,b=eh.get_slit_width(filep,mag=5.0,window_num=win,title=titlec,xran=[float(pitch[k])-float(pitch[k]),float(pitch[k])+float(pitch[k])],figname=figname)
        
            figname='win'+str(win) + 'ang_hist.png'
            print figname
        #results.append(pool.apply_async(eh.hough_hist, args=(fileh,win,figname)))
            eh.hough_hist(fileh,win,figname=figname,side=-1.0)
            
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
