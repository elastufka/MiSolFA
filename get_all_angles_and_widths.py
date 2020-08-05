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
    wins=[11,21,12,22,31,41,32,42,33,43,34,44]
    ll=[20,30,50,100]#[50,100,150,200,250]

    #front window EM - need to get the additional images for win 21 from a previous run...
#     tols=[10,12,4,4,7,7,3,3,2,2,1,1]
#     sigma=[2.,2.,1.5,1.5,2.,2.,1.5,1.5,.75,.75,1.0,1.0]
#     sides=[1,1,1,1,1,1,1,1,1,1,1,1]
#     cparams=[{'tpix':350,'lpix':400,'rpix':False,'bpix':False},{'tpix':False,'lpix': 300,'rpix':False,'bpix':400},{'tpix':400,'rpix':400,'lpix':False,'bpix':False},{'lpix':False,'tpix':False,'rpix':400,'bpix':400},{'tpix':350,'lpix':400,'rpix':False,'bpix':False},{'tpix':False,'lpix': 300,'rpix':False,'bpix':400},{'tpix':400,'rpix':400,'lpix':False,'bpix':False},{'lpix':False,'tpix':False,'rpix':400,'bpix':400},{'tpix':350,'lpix':400,'rpix':False,'bpix':False},{'tpix':False,'lpix': 300,'rpix':False,'bpix':400},{'tpix':400,'rpix':400,'lpix':False,'bpix':False},{'lpix':False,'tpix':False,'rpix':400,'bpix':400}]

    #rear window EM
#     tols=[12,12,4,4,7,7,3,3,5,5,1,1]
#     sigma=[2.,2.,1.5,1.5,2.,2.,1.5,1.5,1.75,1.75,1.0,1.0]
#     sides=[0,0,0,0,0,0,0,0,0,0,0,0]
#     cparams=[{'lpix':250,'tpix':100,'rpix':100,'bpix':240},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':50,'rpix':False,'bpix':300},{'lpix':250,'tpix':50,'rpix':False,'bpix':300},{'lpix':250,'tpix':50,'rpix':False,'bpix':300},{'lpix':250,'tpix':50,'rpix':False,'bpix':300}]

#     #rear window QM
    tols=[7,7,7,7,5,5,4,4,2,2,1,1]
    sigma=[2.,2.,2.,2.,2.,2.,1.75,1.75,1.5,1.5,1.25,1.25]
    sides=[0,0,0,0,0,0,0,0,0,0,0,0]
    cparams=[{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':30,'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'tpix':60,'rpix':False,'bpix':False}]

    #front window QM. Sept=True,May=False
#     tols=[7,7,7,7,5,5,4,4,2,2,1,1]
#     sigma=[2.,2.,2.,2.,2.,2.,1.75,1.75,1.5,1.5,1.25,1.25]
#     sides=[1,1,1,1,1,1,1,1,1,1,1,1]
#     cparams=[{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':30,'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'tpix':60,'rpix':False,'bpix':False},{'tpix':60,'rpix':False,'bpix':False}]

    n=0
    for w,t,sig,s,cp in zip(wins[n:],tols[n:],sigma[n:],sides[n:],cparams[n:]):
        #tw=Grating(w,May=True,EM=True,side=s)
        #tw=Grating(w,May=True,EM=False,side=s) #QM front grid
        #tw=Grating(w,Apr=True,May=False) #EM rear window - the real one...did I fix that?
        tw=Grating(w, May=True, EM=False,side=s)
        #cropping
        #tw.trim_optical_edges(**cp)
        tw.__get_and_store_data__()
        tw.parameterize_optical(coverwrite=True,hoverwrite=True, tol=t,sigma=sig)

    #plot derived vs. nominal parameters, with error bars
    means,noms,yerr=[],[],[]
    for w in win:
        dummy=Grating(w,May=False,Sept=True,EM=False, side=0)
        os.chdir(dummy.data.Odata.path)
        stats=pickle.load(open('win'+str(w)+'_width_stats5.0Xa.p','rb'))
        means.append(stats['mean'])
        noms.append(dummy.nominal['pitch'])
        yerr.append([-1*stats['stddev'],stats['stddev']])
    fig,ax=plt.subplots()
    ax.errorbar(wins,np.array(means)-np.array(noms),yerr=np.transpose(np.array(yerr)))
    fig.suptitle('Sub 2765 - EM Rear Grid')
    ax.set_xlabel('Window Number')
    ax.set_ylabel('Measured Period Difference ($\mu$m)')
    fig.show()


    #front window, EM 5.0:
    #print 'here'
#     os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/EMassembly_2017_11_15_FrontGrid')
#     #os.chdir('/Users/wheatley/Documents/Solar/MiSolFA//mw469sub2634_2017_12_11')
#     for k,win in enumerate(wins):
#         filen=glob.glob('win'+str(win)+'*_5.0X.tif')#eh.EM_list(str(win),'5.0',ending='.tif')
#         filep=glob.glob('win'+str(win)+'*_5.0X_edges.p')#[]#eh.EM_list(str(win),'5.0',ending='edges.p')#glob.glob('win'+str(win)+'*_5.0X_edges.p')
#         fileh=glob.glob('win'+str(win)+'*_5.0X_hough_all.p')#[]#eh.EM_list(str(win),'5.0',ending='_hough.p')

#         if filep==[] or len(filep) != len(filen): #need to generate the pickles by doing the edges
#             for f in filen:
#                 filep.append(eh.Canny_edge(f,mag=5.0,outfilen=f[:-4]+"_edges.p")) #new list of picklefiles
#         #if fileh==[]or len(fileh) != len(filen): #need to generate the pickles by doing the edges
#         #for f in filep:
#         #    fileh.append(eh.prob_hough(f,line_length=ll[0],spread=2.,overwrite=False)) #new list of picklefiles
#         #for l in ll:
#         #    tag='ll'+str(l)
#         #    for f in filep:
#         #        fileh.append(eh.prob_hough(f,line_length=l,spread=2.,tag=tag,overwrite=False)) #new list of picklefiles

#         #eh.cat_hough(win,weights=[1,1,1,1,1],tags=['ll50','ll100','ll150','ll200','ll250'],mag=5.0,EM=False)
#         #fileh=glob.glob('win'+str(win)+'*_5.0X_hough_all.p')
#         #if win in [11,21,32,42] : #use Hough fits to get the widths for coarse windows since Canny edges have false detections
#         #    continue
#         #    filep=fileh #are these really 1 px thick though? need to convert them to boolean arrays

#         titlec='Slit angle distribution for window '+str(win) #+ ', nominal pitch '+pitch[k] + ' $\mu$m at 5X magnification'
#         figname='win'+str(win) + 'width_hist.png'
#         #if len(filen) !=0:
#         w,b,f=eh.get_period_by_grouping(win,mosaic=False,EM=True,plot=True,tolerance=.85,offset=0.,mag=5.0,side=1.)

#         figname='win'+str(win) + 'ang_hist.png'
#      #   print figname
#         #results.append(pool.apply_async(eh.hough_hist, args=(fileh,win,figname)))
#         eh.hough_hist(fileh,win,figname=figname,tol=tol[k],side=1.,spread=2.)

    #rear window, EM 5.0
#     os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw469sub2765_2018_01_31')
#     for k,win in enumerate(wins):
#         filen=eh.EM_list(str(win),'5.0',ending='.tif')
#         filep=glob.glob('win'+str(win)+'*_5.0X_edges.p')#[]#eh.EM_list(str(win),'5.0',ending='_edges.p')
#         fileh=[]#eh.EM_list(str(win),'5.0',ending='_hough.p')

#         if filep==[] or len(filep) != len(filen): #need to generate the pickles by doing the edges
#             for f in filen:
#                 filep.append(eh.Canny_edge(f,mag=5.0,outfilen=f[:-4]+"_edges.p")) #new list of picklefiles
#         if fileh==[]or len(fileh) != len(filen): #need to generate the pickles by doing the edges
#             for f in filep:
#                 fileh.append(eh.prob_hough(f,line_length=ll[0],spread=2.,overwrite=True,side=-1)) #new list of picklefiles
#         for l in ll[1:]:
#             tag='ll'+str(l)
#             for f in filep:
#                fileh.append(eh.prob_hough(f,line_length=l,spread=2.,tag=tag,overwrite=True,side=-1)) #new list of picklefiles

#         eh.cat_hough(win,weights=[1,1,1,1,1],tags=['ll100','ll150','ll200','ll250'],EM=False)
#         fileh=glob.glob('win'+str(win)+'*_5.0X_*hough_all.p')
        #if win in [11,21,32,42] : #use Hough fits to get the widths for coarse windows since Canny edges have false detections
        #    continue
        #    filep=fileh #are these really 1 px thick though? need to convert them to boolean arrays

#         titlec='Slit angle distribution for window '+str(win) #+ ', nominal pitch '+pitch[k] + ' $\mu$m at 5X magnification'
#         figname='win'+str(win) + 'width_hist.png'
#         if len(filen) !=0:
#            #w,b=eh.get_slit_width(filep,mag=5.0,window_num=win,title=titlec,xran=[float(pitch[k])-float(pitch[k]),float(pitch[k])+float(pitch[k])],figname=figname)
#             w,b,f=eh.get_period_by_grouping(win,mosaic=False,EM=True,plot=True,tolerance=.85,offset=0.,mag=5.0,side=-1.,fitangle=False)
#         figname='win'+str(win) + 'ang_hist.png'
        #print figname
        #results.append(pool.apply_async(eh.hough_hist, args=(fileh,win,figname)))
        #if win != 11:
        #    eh.hough_hist(fileh,win,figname=figname,tol=tol[k],spread=2.)


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
