"""
===================
fit_all_transm.py
Erica  Lastufka 11.11.18
===================

What the title says

"""
import edge_and_hough as eh
import glob
import time
import multiprocessing as mpi
import os

if __name__ == '__main__':
    wins=[11,21,12,22,31,41,32,42,33,43,34,44] #no win12 in this data set! it's fake...

    #front window EM - need to get the additional images for win 21 from a previous run...
    sub='sub2737'
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmEM')
    hbounds=[254.9,255.1]#[245,265] #from manufacturer
    wbounds=[[41.3456,48.0712],[42.7218,47.5184],[10.3721,12.1817],[10.39175,12.18535],[21.4977,23.5777],[21.36925,23.57585],[8.265,9.757],[8.2697,9.7079],[13.87375,16.11395],[13.9257,16.1361],[6.80025,8.29385],[6.7704,8.2992]]
    #wbounds=[[89.4168/2.,89.6752/2.],[90.2402/2.,90.326/2.],[22.4797/2.,22.5538/2.],[22.5203/2.,22.5771/2.],[45.0754/2.,45.0814/2.],[44.9187/2.,44.9451/2.],[18.013/2.,18.022/2.],[17.9776/2.,17.987/2.],[29.9639/2.,29.9877/2.],[30.0362/2.,30.0618/2.],[14.991/2.,15.0941/2.],[15.009/2.,15.07/2.]]
    fbounds=[[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[2,-2],[1,-1],[1,-1]]

    #rear window EM
    #sub='sub2765'
    #hbounds=[279.9,280.1]#[255,305] #from manufacturer
    #os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/SLS_Apr2018/transmEM')
    #wbounds=[[42.34495,47.58595],[43.4142,46.4464],[10.40045,12.21945],[10.31165,12.24685],[21.55615,23.39435],[21.5916,23.526],[8.2406,9.8016],[8.20185,9.85905],[14.06385,16.01765],[13.98065,16.02005],[6.759,8.3018],[6.7719,8.2731]]
    #wbounds=[[89.9309/2.,90.326/2.],[89.6752/2.,89.8606/2.],[22.5203/2.,22.6199/2.],[22.4797/2.,22.5585/2.],[44.9187/2.,44.9505/2.],[45.0814/2.,45.1176/2.],[17.987/2.,18.0422/2.],[18.013/2.,18.0609/2.],[30.0362/2.,30.0815/2.],[29.9639/2.,30.0007/2.],[15.009/2.,15.0608/2.],[14.991/2.,15.05/2.]]
    ybounds=False#[-.0001,.0001]
    dc_bounds=False#[.4999,.5001]
    #fbounds=[[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[2,-2],[1,-1],[1,-1],[1,-2],[1,-1],[1,-1]]

#     tols=[12,12,4,4,7,7,3,3,5,5,1,1]
#     sigma=[2.,2.,1.5,1.5,2.,2.,1.5,1.5,1.75,1.75,1.0,1.0]
#     sides=[0,0,0,0,0,0,0,0,0,0,0,0]
#     cparams=[{'lpix':250,'tpix':100,'rpix':100,'bpix':240},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':100,'rpix':100,'bpix':300},{'lpix':250,'tpix':50,'rpix':False,'bpix':300},{'lpix':250,'tpix':50,'rpix':False,'bpix':300},{'lpix':250,'tpix':50,'rpix':False,'bpix':300},{'lpix':250,'tpix':50,'rpix':False,'bpix':300}]

#     #front window QM - sub 3501
    #sub='sub3501'
    #wbounds=[[87.2244/2.,91.9756/2.],[87.4751/2.,93.1249/2.],[43.3418/2.,46.4582/2.],[43.638/2.,46.562/2.],[29.145/2.,30.855/2.],[29.1235/2.,31.0765/2.],[21.6905/2.,23.3095/2.],[21.6713/2.,23.3287/2.],[17.3539/2.,18.6461/2.],[17.3559/2.,18.6441/2.],[14.2387/2.,15.7613/2.],[14.0687/2.,15.9313/2.]]
    #wbounds=[[89.5116,89.6],[90.2881,90.3],[44.9,44.9411],[45.0883,45.1],[30,30.0389],[30.1,30.1934],[22.5,22.5928],[22.5,22.5505],[18,18.0438],[18,18.0804],[15,15.0655],[15,15.0394]] #use the same from 3501
    #fbounds=[[1,-2],[1,-1],[0,-1],[0,-1],[1,-4],[3,-2],[1,-1],[1,-1],[3,-1],[0,-4],[1,-1],[2,-1]]
    #xvmids=[False, 5, 2, 2,2,3,3,False,2,2,False,1]
    #os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018/transmQM')

#     tols=[7,7,7,7,5,5,4,4,2,2,1,1]
#     sigma=[2.,2.,2.,2.,2.,2.,1.75,1.75,1.5,1.5,1.25,1.25]
#     sides=[0,0,0,0,0,0,0,0,0,0,0,0]
#     cparams=[{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':30,'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'tpix':60,'rpix':False,'bpix':False}]
    #rear window QM. Sept=True,May=False
    sub='sub3437'
    hbounds=[225,325]#[284.9,285.1]#[225,325] #nominally 285, but they say possible less than 250
    wbounds=[[87.9244/2.,92.6756/2.],[86.7751/2.,92.4249/2.],[43.5418/2.,46.6582/2.],[43.438/2.,46.362/2.],[29.245/2.,30.955/2.],[29.0235/2.,30.9765/2.],[21.6905/2.,23.3095/2.],[21.6713/2.,23.3287/2.],[17.3539/2.,18.6461/2.],[17.3559/2.,18.6441/2.],[14.2387/2.,15.7613/2.],[14.0687/2.,15.9313/2.]]
    fbounds=[[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[1,-1],[2,-2],[1,-1],[1,-1]]
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmQM')
    dc_bounds=False#[.4999,.5001]
#     tols=[7,7,7,7,5,5,4,4,2,2,1,1]
#     sigma=[2.,2.,2.,2.,2.,2.,1.75,1.75,1.5,1.5,1.25,1.25]
#     sides=[1,1,1,1,1,1,1,1,1,1,1,1]
#     cparams=[{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':30,'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'rpix':False,'bpix':False},{'lpix':60,'tpix':60,'rpix':False,'bpix':False},{'tpix':60,'rpix':False,'bpix':False},{'tpix':60,'rpix':False,'bpix':False}]
    tag='_all'
    n=0
    vdicts=[]
    for w,wb,fb,xvm in zip(wins[n:],wbounds[n:],fbounds[n:],xvmids[n:]):
        dcb=[.48,.52]
        #print w,np.array(wb)*(1./np.array(dcb))
        bounds=[wb,np.array(wb)*(1./np.array(dcb)),hbounds,[-1,1],dcb,[.05,.2]]
        #tw=Grating(w,May=True,EM=False,side=s) #QM front grid
        try:
            tw=pickle.load(open(sub+'_win'+str(w)+ '.p','rb'))#restore pickles with the transm['widths'] results
        except IOError:
            os.chdir('../transmQM')
            tw=pickle.load(open(sub+'_win'+str(w)+ '.p','rb'))
        #Grating(w,Apr=True,May=False) #EM rear window - the real one...did I fix that?
        #cropping
        res=tw.results['widths']
        res.nominal=tw.nominal
        for p in range(0,7): #pindex
            figname=sub+'_win'+str(w)+'_p'+str(p)+'_fit'+tag
            res.fit_profile_to_data_optimize(p,tw.nominal['pitch'],tw.nominal['height'],plot=True ,ret=False,height_bounds=hbounds, width_bounds=wb,dc_bounds=dc_bounds,yoff_bounds=ybounds,fit_bounds=fb,fix_dc=False,show=figname,xvmid=False)

        res.calc_nominal_transm()
        res.scat_allP(meanfit=False,figname=sub+'_win'+str(w)+'_fit_profiles'+tag+'.png',theta_ran=False,print_height=False,fit='both',yran=[.2,.7])
        #store
        #print wb
        vd=res.fit_plots(bounds=bounds,ret=True)
        vdicts.append(vd)
        pickle.dump(vd,open('../'+sub+'_win'+str(w)+'fit_vdict.p','wb'))
        print "==============STATISTICS FOR FIT TO WIN "+str(w)+"=================="
        print "          Average Width: " ,np.mean(vd['width_yvals'])
        print "         Nominal Period: ",tw.nominal['pitch']
        print "         Average Period: " , np.mean(vd['period_yvals'])
        print "         Average Height: ",np.mean(vd['height_yvals'])
        print "     Average Duty Cycle: ",np.mean(vd['dc'])
        print "       Average Y-Offset: ",np.mean(vd['yoff'])

        pickle.dump(tw,open(sub+'_win'+str(w)+ '_fit.p','wb'))
        os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmQM')


    #plot for the whole grating
    xvals = range(0,7)
    colors=['k','r','m','c','g','b','y','pink','brown','gold','orange','darkred'] #12 colors
    labels=['11','21','12','22','31','41','32','42','33','43','34','44']
    pitch=[v['pitch_nom'] for v in vdicts]
    print pitch
    fig,ax=plt.subplots(2,3, figsize=[12,12],sharex=True)
    width_yvals=[v['width_yvals'] for v in vdicts]
    for wy,p,c in zip(width_yvals,pitch,colors):
        ax[0][0].plot(xvals,np.array(wy)-p/2.,'o-', color=c)
    ax[0][0].set_title('fit_width -nominal (um)')
    yran=ax[0][0].get_ylim()
    ymax=float(yran[-1])
    #ax[0][0].text(1,.8*ymax,'mean width: ' + str(np.round(np.mean(width_yvals),3))+ ' $\mu$m')

    period_yvals=[v['period_yvals'] for v in vdicts]
    for py,p,c in zip(period_yvals,pitch,colors):
        ax[0][1].plot(xvals,np.array(py)-p,'o-',color=c)
    ax[0][1].set_title('fit period -nominal (um)')

    height_yvals=[v['height_yvals'] for v in vdicts]
    #hcov=[self.fit_both['p'+str(k)]['error'][1] for k in xvals]
    #print hcov
    #ax[0][2].errorbar(xvals,height_yvals,yerr=np.array(hcov))
    for hy,c in zip(height_yvals,colors):
        ax[0][2].plot(xvals,hy,'o-',color=c)
    ax[0][2].set_title('fit_height (um)')

    th_yvals=[v['peak_theta'] for v in vdicts]
    #th_yvals2=[self.fit_height['p'+str(k)]['peak_theta'] for k in xvals]
    for th,c,lab in zip(th_yvals,colors,labels):
        ax[1][0].plot(xvals,th,'o-',color=c,label=lab)
    #ax[2].plot(xvals,th_yvals2,'o-',label='fixed width')
    ax[1][0].legend(loc='upper left')
    ax[1][0].set_title('dtheta (degrees)')

    dc=[v['dc'] for v in vdicts]
    for d,c in zip(dc,colors):
        ax[1][1].plot(xvals,d,'o-',color=c)
    ax[1][1].set_title('Duty Cycle')

    #ycov=[self.fit_both['p'+str(k)]['error'][3] for k in xvals]
    yoff=[v['yoff'] for v in vdicts]
    for yo,c in zip(yoff,colors):
        ax[1][2].plot(xvals,yo,'o-',color=c)
    #ax[2].plot(xvals,th_yvals2,'o-',label='fixed width')
    ax[1][2].set_title('Y-offset')

#     if bounds: #set all the boundaries
#         ax[0][0].set_ylim(np.array(bounds[0])-pitch/2.)
#         #ax[0][1].set_ylim(np.array(bounds[1])-pitch)
#         ax[0][2].set_ylim(bounds[2])
#         ax[1][0].set_ylim(bounds[3])
#         ax[1][1].set_ylim(bounds[4])
#         ax[1][2].set_ylim(bounds[5])

    #ax[2].legend()
    for a in ax[1]:
        a.set_xlabel('p')
    for a in ax[0]:
        a.set_xlabel('p')
    fig.suptitle('sub3437')
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
