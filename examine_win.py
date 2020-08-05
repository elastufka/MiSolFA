#get analyze optical results and generate tables for entire window
import os

def examine_single_grating(win,fac=1.5,percent=2,plots=False):
    g=pickle.load(open(win,'rb'))
    rr=g.results['optical']
    rr._get_stats()
    print rr.period_stats
    rr.recalc_sc()
    print rr.period_stats
    rr.filter_optical_widths('periods_from_slat_centers',fac=fac,percent=percent)
    rr.filter_optical_widths('periods_from_slit_centers',fac=fac,percent=percent)
    print rr.period_stats
    #print tex line
    #rr.opt_stats_tex()
    if plots:
        rr.plot_angle_hist()
        rr.plot_period_hist(nbins=50)
    #pickle object
    os.chdir(g.data.Odata.path)
    pickle.dump(g,open(win,'wb'))

def examine_grating(wins, model='EM',side=1.0,Apr=True,May=False,Sept=False,rerun=False,plots=False,fac=1.5,percent=2.,coverwrite=False, hoverwrite=False,example=False):
    #wins=[11,12,21,22,31,41,32,42,33,43,34,44]
    #wins=[11,21,22,31,41,32,42,33,43,34,44]
    if model == 'EM':
        EM=True
    else:
        EM=False
    for w in wins:
        g=Grating(w,EM=EM,side=side,Apr=Apr,May=May,Sept=Sept)
        if rerun: #rerun analysis from Hough
            g.parameterize_optical(coverwrite=coverwrite,hoverwrite=hoverwrite,example=example,filter_nominal=False) #results goes into self.results['optical']
        if g.results == {} or 'period_stats' not in g.results['optical'].__dict__.keys(): #load the results object the other way
            import results
            kind={'type':'optical','model':g.model,'folder':g.data.Odata.path,'side':g.side}
            g.results['optical']=results.Results(w,kind,nominal=g.nominal)
        #do filter_nominal
        rr=g.results['optical']
        rr._get_stats()
        print rr.period_stats
        rr.recalc_sc()
        print rr.period_stats
        rr.filter_optical_widths('periods_from_slat_centers',fac=fac,percent=percent)
        rr.filter_optical_widths('periods_from_slit_centers',fac=fac,percent=percent)
        #print tex line
        rr.opt_stats_tex()
        rr.opt_params_tex()
        if plots:
            rr.plot_angle_hist()
            rr.plot_period_hist(nbins=50)
        #pickle object
        oname='win'+str(w)+'_grating.p'
        os.chdir(g.data.Odata.path)
        pickle.dump(g,open(oname,'wb'))

#method for just getting tables out of grating objects...

def all_opt_stats_tex(folder=False):
    wins=[11,21,12,22,31,41,32,42,33,43,34,44]
    if folder:
        os.chdir(folder)
    for w in wins:
        oname='win'+str(w)+'_grating.p'
        g=pickle.load(open(oname,'rb'))
        g.results['optical'].opt_stats_tex()

def all_xray_stats_tex(folder=False, wins=False,offset=1):
    if not wins:
        wins=[11,21,12,22,31,41,32,42,33,43,34,44]
    if folder:
        os.chdir(folder)
    def rv(value):
        return str(np.round([value],decimals=4)[0])

    for w in wins:
        oname='win'+str(w)+'_grating.p'
        g=pickle.load(open(oname,'rb'))
            #do the same for X-ray results
        win=str(g.win)
        side=g.side
        npitch=g.nominal['pitch']
        try:
            nang=g.nominal['nominal angle']
        except KeyError:
            nang=g.nominal['orientation']
        #req = self.requirements
        pxstats=g.results['widths'].period_stats
        axstats=g.results['widths'].angle_stats
        postats=g.results['optical'].period_stats
        aostats=g.results['optical'].angle_stats
        th0=axstats['theta0']
        axdev=axstats['theta_std']
        axdiff=axstats['delta_nom']
        axmean=axdiff+nang
        aomean=aostats['mean']
        aodev=aostats['stddev']
        aodiff=aomean-nang
        pxmean=pxstats['mean']
        pxdev=pxstats['stddev']
        pxdiff=pxmean-npitch
        pomean=postats['mean']
        podev=postats['stddev']
        podiff=pomean-npitch

        try:
            aooff=aostats['th_off']
            axoff=axstats['th_off']
        except KeyError:
            aooff=0.
            axoff=0.
        #for the results table. win | npitch | xpmeas | pxdiff | podiff | pxdev | podev | nang | axmean | axdiff | aodiff | axdev | aodev
        line=win + ' & ' + str(npitch) +' & '+rv(pxmean)+' & '+rv(pxdiff)+' & '+rv(podiff)+ ' & '+ rv(pxdev) +' & '+rv(podev) +' & '+str(nang)+' & '+rv(axmean)+' & '+rv(axdiff-axoff)+' & '+rv(aodiff-aooff)+ ' & '+ rv(axdev) +' & '+rv(aodev)+' \\\\'
        #self.tex['stats']=line
        print line


def all_opt_params_tex(folder=False):
    wins=[11,21,12,22,31,41,32,42,33,43,34,44]
    if folder:
        os.chdir(folder)
    for w in wins:
        oname='win'+str(w)+'_grating.p'
        g=pickle.load(open(oname,'rb'))
        g.results['optical'].opt_params_tex()

def xray_fits_tex(EM=False,wins=False):
    ''' for the results table. win | front width | front chisq | front dc | front height | front height dev | rear width | rear chisq | rear dc | rear height | rear height dev '''
    def rv(value):
        return str(np.round([value],decimals=4)[0])
    if not wins:
        wins=[11,21,12,22,31,41,32,42,33,43,34,44]
    if EM and not wins:
        wins=[11,21,31,41,33,43,12,22,32,42,34,44]

    if EM:
        front_folder='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/EMassembly_2017_11_15_FrontGrid'
        rear_folder='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw469sub2765_2018_01_31'
    else:
        front_folder='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw561sub3501_2018_06_26'
        rear_folder='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw501sub3437_2018_05_04'


    for w in wins:
        oname='win'+str(w)+'_grating.p'
        os.chdir(front_folder)
        fg=pickle.load(open(oname,'rb'))
        os.chdir(rear_folder)
        rg=pickle.load(open(oname,'rb'))
        win=str(fg.win)

        if 'fdict' not in fg.results.keys():
            #it's win 42 or 43, no results for fits
            fwidth='-'
            fdc='-'
            fchisq='-'
        else:
            #front grid stuff
            fwidth=fg.results['fdict']['int_y']
            fpitch=fg.results['widths'].period_stats['mean']
            fdc=fwidth/fpitch
            try:
                if len(fg.results['fdict']['plus_err'])==0:
                    fchisq=rv(fg.results['fdict']['minus_err'][0])
                elif len(fg.results['fdict']['minus_err'])==0:
                    fchisq=rv(fg.results['fdict']['plus_err'][0])
                else:
                    fchisq=rv(np.mean([fg.results['fdict']['plus_err'],fg.results['fdict']['minus_err']]))
            except KeyError:
                fchisq='-'
            fheights=[fg.results['sum'].fit_both[k]['fit_height'] for k in fg.results['sum'].fit_both if k !='bounds']
            fheight=np.mean(fheights)
            fheight_dev=np.std(fheights)
            ferror=np.mean([fg.results['sum'].fit_both[k]['error'] for k in fg.results['sum'].fit_both if k !='bounds'])

        #front grid stuff
        rwidth=rg.results['fdict']['int_y']
        rpitch=rg.results['widths'].period_stats['mean']
        rdc=rwidth/rpitch
        try:
            if len(rg.results['fdict']['plus_err'])==0 and len(rg.results['fdict']['minus_err'])!=0:
                rchisq=rv(rg.results['fdict']['minus_err'][0])
            elif len(rg.results['fdict']['minus_err'])==0 and len(rg.results['fdict']['plus_err'])!=0:
                rchisq=rv(rg.results['fdict']['plus_err'][0])
            else:
                rchisq=rv(np.mean([rg.results['fdict']['plus_err'],rg.results['fdict']['minus_err']]))
        except KeyError:
            rchisq='-'
        rheights=[rg.results['sum'].fit_both[k]['fit_height'] for k in rg.results['sum'].fit_both if k !='bounds']
        rheight=np.mean(rheights)
        rheight_dev=np.std(rheights)
        rerror=np.mean([rg.results['sum'].fit_both[k]['error'] for k in rg.results['sum'].fit_both if k !='bounds'])

        if 'fdict' in fg.results.keys():
            line=win +' & '+rv(fwidth)+' & '+fchisq+' & '+rv(fdc)+' & '+ rv(fheight)+ ' & '+ rv(ferror)+' & ' +rv(rwidth)+' & '+rchisq+' & '+rv(rdc)+' & '+ rv(rheight)+ ' & '+ rv(rerror)+' \\\\'
        else:
            line=win +' & '+fwidth+' & '+fchisq+' & '+fdc+' & '+ rv(fheight)+ ' & '+ rv(ferror)+' & ' +rv(rwidth)+' & '+rchisq+' & '+rv(rdc)+' & '+ rv(rheight)+ ' & '+ rv(rerror)+' \\\\'

        print line

def plot_diff(folder=False,title='title',offset=1,ayran=[-.5,.5],xray=False,legend=True,lloc='upper right',EM=False):
    wins=[11,21,12,22,31,41,32,42,33,43,34,44]
    if EM:
        wins=[11,21,31,41,33,43,12,22,32,42,34,44]
    if folder:
        os.chdir(folder)
    pm,am,pdev,adev,pnom,anom,preq,areq=[],[],[],[],[],[],[],[]
    if xray:
        xpm,xam,xpdev,xadev=[],[],[],[]
    for w in wins:
        oname='win'+str(w)+'_grating.p'
        g=pickle.load(open(oname,'rb'))
        pm.append(g.results['optical'].period_stats['amean'])
        am.append(g.results['optical'].angle_stats['mean'])
        pdev.append(g.results['optical'].period_stats['adev'])
        adev.append(g.results['optical'].angle_stats['stddev'])
        pnom.append(g.nominal['pitch'])
        try:
            anom.append(g.nominal['nominal angle'])
        except KeyError:
            anom.append(g.nominal['orientation'])
        preq.append(g.results['optical'].requirements['pitch'])
        areq.append(g.results['optical'].requirements['angle'])
        if xray:
            xpm.append(g.results['widths'].period_stats['mean'])
            xam.append(np.nanmean(g.results['widths'].angle_stats['delta_nom']))
            xpdev.append(g.results['widths'].period_stats['stddev'])
            xadev.append(g.results['widths'].angle_stats['theta_std'])

    adiff=np.array(am)-np.array(anom)
    preq.append(preq[-1])

    if offset == 1 and EM: #by segment
        segBadiff=[adiff[0],adiff[1],adiff[6],adiff[7]]
        segCDadiff=[adiff[2],adiff[3],adiff[4],adiff[5],adiff[8],adiff[9],adiff[10],adiff[11]]
        amB=np.mean(segBadiff)
        amCD=np.mean(segCDadiff)
        print 'seg B offset: ',amB
        print 'seg CD offset: ',amCD
        for i in [0,1,6,7]:
            adiff[i]=adiff[i]-amB
        for j in [2,3,4,5,8,9,10,11]:
            adiff[j]=adiff[j]-amCD
        if xray: #calc offsets here too
            xsegBadiff=[xam[0],xam[1],xam[6],xam[7]]
            xsegCDadiff=[xam[2],xam[3],xam[4],xam[5],xam[8],xam[9],xam[10],xam[11]]
            XamB=np.mean(xsegBadiff)
            XamCD=np.mean(xsegCDadiff)
            print 'seg B offset (Xray): ',XamB
            print 'seg CD offset (Xray): ',XamCD
            for i in [0,1,6,7]:
                xam[i]=xam[i]-XamB
            for j in [2,3,4,5,8,9,10,11]:
                xam[j]=xam[j]-XamCD


    elif offset == 2: #whole thing.
        aoff=np.mean(adiff)
        adiff=adiff-aoff
        print 'offset: ',aoff
        if xray: #calc offsets here too
            Xaoff=np.mean(xam)
            print 'offset (Xray): ',Xaoff
            xam=xam-Xaoff

    fig,ax=plt.subplots()
    xvals=range(0,len(wins))
    xrvals=range(0,13)
    if xray:
        xrayvals=np.array(xvals)+.5
        ax.scatter(xrayvals,np.array(xpm)-np.array(pnom),c='r',label='Measured-Nominal, Xray')
        #ax.errorbar(xrayvals,np.zeros(len(wins)),yerr=np.array(xpdev),elinewidth=2,fmt=None,ecolor='b')
        ax.errorbar(xrayvals,np.array(xpm)-np.array(pnom),yerr=np.array(xpdev),ecolor='r',fmt=None,label="Std. Dev (Xray)")

    ax.scatter(xvals,np.array(pm)-np.array(pnom),c='k',label='Measured-Nominal, opt.')
    ax.fill_between(xrvals,-1*np.array(preq),np.array(preq),alpha=.6,color='b',label="Requirements")
    ax.errorbar(xvals,np.array(pm)-np.array(pnom),yerr=np.array(pdev),fmt=None,ecolor='k',label="Std. Dev, opt.")

    ax.set_xticks(xvals)
    ax.set_xticklabels(wins) #not if xray...
    ax.set_xlabel('Window')
    ax.set_ylabel('Measured Pitch - Nominal Pitch ($\mu$m)')
    ax.set_xlim(-1,len(wins))
    ax.set_ylim(-1,1)
    ax.set_title(title + ' Results - Pitch')
    if legend:
        ax.legend(loc=lloc)
    fig.show()

    fig,ax=plt.subplots()
    if xray:
        ax.scatter(xrayvals,np.array(xam),c='r',label='Measured-Nominal, Xray')
        #ax.errorbar(xrayvals,np.zeros(len(wins)),yerr=np.array(xadev),elinewidth=2,fmt=None,ecolor='b')
        ax.errorbar(xrayvals,np.array(xam),yerr=np.array(xadev),ecolor='r',fmt=None,label="Std. Dev (Xray)")
    ax.scatter(xvals,adiff,c='k',label='Measured-Nominal, opt.')
    ax.fill_between(xrvals,-1*np.array(preq),np.array(preq),alpha=.6,color='b',label="Requirements")
    ax.errorbar(xvals,adiff,yerr=np.array(adev),ecolor='k',fmt=None,label="Std. Dev, opt.")
    #ax.set_xticks(ticks=range(0,len(wins)),labels=wins)
    ax.set_xticks(range(0,len(wins)))
    ax.set_xticklabels(wins)
    ax.set_xlim(-1,len(wins))
    ax.set_ylim([ayran[0],ayran[1]])
    ax.set_xlabel('Window')
    ax.set_ylabel('Measured Angle - Nominal Angle (degrees)')
    ax.set_title(title + ' Results - Angle')
    if legend:
        ax.legend(loc=lloc)
    fig.show()


def get_diff(folder, wins, offset=2,xray=False):
    os.chdir(folder)
    pm,am,pdev,adev,pnom,anom,preq,areq=[],[],[],[],[],[],[],[]
    if xray:
        xpm,xam,xpdev,xadev=[],[],[],[]
    for w in wins:
        oname='win'+str(w)+'_grating.p'
        g=pickle.load(open(oname,'rb'))
        pm.append(g.results['optical'].period_stats['amean'])
        am.append(g.results['optical'].angle_stats['mean'])
        pdev.append(g.results['optical'].period_stats['adev'])
        adev.append(g.results['optical'].angle_stats['stddev'])
        pnom.append(g.nominal['pitch'])
        try:
            anom.append(g.nominal['nominal angle'])
        except KeyError:
            anom.append(g.nominal['orientation'])
        preq.append(g.results['optical'].requirements['pitch'])
        areq.append(g.results['optical'].requirements['angle'])
        if xray:
            xpm.append(g.results['widths'].period_stats['mean'])
            xam.append(np.nanmean(g.results['widths'].angle_stats['delta_nom']))
            xpdev.append(g.results['widths'].period_stats['stddev'])
            xadev.append(g.results['widths'].angle_stats['theta_std'])

    adiff=np.array(am)-np.array(anom)
    #preq.append(preq[-1])

    if offset == 1: #by segment
        segBadiff=[adiff[0],adiff[1],adiff[6],adiff[7]]
        segCDadiff=[adiff[2],adiff[3],adiff[4],adiff[5],adiff[8],adiff[9],adiff[10],adiff[11]]
        amB=np.mean(segBadiff)
        amCD=np.mean(segCDadiff)
        print 'seg B offset: ',amB
        print 'seg CD offset: ',amCD
        for i in [0,1,6,7]:
            adiff[i]=adiff[i]-amB
        for j in [2,3,4,5,8,9,10,11]:
            adiff[j]=adiff[j]-amCD
        if xray: #calc offsets here too
            xsegBadiff=[xam[0],xam[1],xam[6],xam[7]]
            xsegCDadiff=[xam[2],xam[3],xam[4],xam[5],xam[8],xam[9],xam[10],xam[11]]
            XamB=np.mean(xsegBadiff)
            XamCD=np.mean(xsegCDadiff)
            print 'seg B offset (Xray): ',XamB
            print 'seg CD offset (Xray): ',XamCD
            for i in [0,1,6,7]:
                xam[i]=xam[i]-XamB
            for j in [2,3,4,5,8,9,10,11]:
                xam[j]=xam[j]-XamCD


    elif offset == 2: #whole thing.
        aoff=np.mean(adiff)
        adiff=adiff-aoff
        print 'offset: ',aoff
        if xray: #calc offsets here too
            Xaoff=np.mean(xam)
            print 'offset (Xray): ',Xaoff
            xam=xam-Xaoff
    if xray:
        return pm,am,pdev,adev,pnom,anom,preq,areq,adiff,xpm,xam,xpdev,xadev
    else:
        return pm,am,pdev,adev,pnom,anom,preq,areq,adiff

def plot_diff_together(title='title',ayran=[-.5,.5],xray=False,legend=True,lloc='upper right',EM=False,offset=2):
    wins=[11,21,12,22,31,41,32,42,33,43,34,44]
    ff='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw561sub3501_2018_06_26'
    rf='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw501sub3437_2018_05_04'
    if EM:
        wins=[11,21,31,41,33,43,12,22,32,42,34,44]
        ff='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/EMassembly_2017_11_15_FrontGrid'
        rf='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw469sub2765_2018_01_31'

    print ff
    print rf
    if xray:
        fpm,fam,fpdev,fadev,fpnom,fanom,fpreq,fareq,fadiff,fxpm,fxam,fxpdev,fxadev=get_diff(ff,wins,xray=True,offset=offset)
        rpm,ram,rpdev,radev,rpnom,ranom,rpreq,rareq,radiff,rxpm,rxam,rxpdev,rxadev=get_diff(rf,wins,xray=True,offset=offset)
    else:
        fpm,fam,fpdev,fadev,fpnom,fanom,fpreq,fareq,fadiff=get_diff(ff,wins,offset=offset)
        rpm,ram,rpdev,radev,rpnom,ranom,rpreq,rareq,radiff=get_diff(rf,wins,offset=offset)

    fig,ax=plt.subplots()
    xfvals=range(0,len(wins))
    xrvals=range(12,2*len(wins))

    bgv=np.array(range(-2,12))+.5
    gb=3.+bgv

    ax.fill_between(bgv,-1*gb,gb,alpha=.2,color='gray')

    if xray:
        xrayfvals=np.array(xfvals)+.5
        xrayrvals=np.array(xrvals)+.5
        ax.scatter(xrayfvals,np.array(fxpm)-np.array(fpnom),c='r',label='Measured-Nominal, Xray')
        ax.scatter(xrayrvals,np.array(rxpm)-np.array(rpnom),c='r')
        #ax.errorbar(xrayvals,np.zeros(len(wins)),yerr=np.array(xpdev),elinewidth=2,fmt=None,ecolor='b')
        ax.errorbar(xrayfvals,np.array(fxpm)-np.array(fpnom),yerr=np.array(fxpdev),ecolor='r',fmt=None,label="Std. Dev (Xray)")
        ax.errorbar(xrayrvals,np.array(rxpm)-np.array(rpnom),yerr=np.array(rxpdev),ecolor='r',fmt=None)
        ax.fill_between(xrayfvals[-2:],-1*np.array(fpreq[-2:]),np.array(fpreq[-2:]),alpha=.6,color='b')
        ax.fill_between(xrayrvals[-2:],-1*np.array(rpreq[-2:]),np.array(rpreq[-2:]),alpha=.6,color='b')

    ax.scatter(xfvals,np.array(fpm)-np.array(fpnom),c='k',label='Measured-Nominal, opt.')
    ax.scatter(xrvals,np.array(rpm)-np.array(rpnom),c='k')
    ax.fill_between(xfvals,-1*np.array(fpreq),np.array(fpreq),alpha=.6,color='b')
    ax.fill_between(xrvals,-1*np.array(rpreq),np.array(rpreq),alpha=.6,color='b',label="Requirements")
    ax.errorbar(xfvals,np.array(fpm)-np.array(fpnom),yerr=np.array(fpdev),fmt=None,ecolor='k',label="Std. Dev, opt.")

    ax.errorbar(xrvals,np.array(rpm)-np.array(rpnom),yerr=np.array(rpdev),fmt=None,ecolor='k')

    ax.set_xticks(xfvals+xrvals)
    ax.set_xticklabels(wins+wins) #not if xray...
    ax.set_xlabel('Window')
    ax.set_ylabel('Measured Pitch - Nominal Pitch ($\mu$m)')
    ax.set_xlim(-1,25)
    ax.set_ylim(-1.5,1.5)
    ax.set_title(title + ' Results - Pitch')
    if legend:
        ax.legend(loc=lloc)
    fig.show()

    fig,ax=plt.subplots()
    ax.fill_between(bgv,-1*gb,gb,alpha=.2,color='gray')
    if xray:
        if not EM:
            mask=[0,0,0,0,0,0,0,0,0,1,1,0]
            fxam=np.ma.masked_array(fxam,mask=mask)
            fxadev=np.ma.masked_array(fxadev,mask=mask)
            print fxam
        else:
            fxam=np.array(fxam)
            fxadev=np.array(fxadev)

        ax.scatter(xrayfvals,fxam,c='r',label='Measured-Nominal, Xray')
        ax.scatter(xrayrvals,np.array(rxam),c='r')
        #ax.errorbar(xrayvals,np.zeros(len(wins)),yerr=np.array(xadev),elinewidth=2,fmt=None,ecolor='b')
        ax.errorbar(xrayfvals,fxam,yerr=fxadev,ecolor='r',fmt=None,label="Std. Dev (Xray)")
        ax.errorbar(xrayrvals,np.array(rxam),yerr=np.array(rxadev),ecolor='r',fmt=None)
        ax.fill_between(xrayfvals[-2:],-1*np.array(fareq[-2:]),np.array(fareq[-2:]),alpha=.6,color='b')
        ax.fill_between(xrayrvals[-2:],-1*np.array(rareq[-2:]),np.array(rareq[-2:]),alpha=.6,color='b')

    ax.scatter(xfvals,fadiff,c='k',label='Measured-Nominal, opt.')
    ax.scatter(xrvals,radiff,c='k')
    ax.fill_between(xfvals,-1*np.array(rareq),np.array(rareq),alpha=.6,color='b')
    ax.fill_between(xrvals,-1*np.array(rareq),np.array(rareq),alpha=.6,color='b',label="Requirements")
    ax.errorbar(xfvals,fadiff,yerr=np.array(fadev),ecolor='k',fmt=None,label="Std. Dev, opt.")
    ax.errorbar(xrvals,radiff,yerr=np.array(radev),ecolor='k',fmt=None)
    #ax.set_xticks(ticks=range(0,len(wins)),labels=wins)
    ax.set_xticks(range(0,24))
    ax.set_xticklabels(wins+wins)
    ax.set_xlim(-1,25)
    ax.set_ylim([ayran[0],ayran[1]])
    ax.set_xlabel('Window')
    ax.set_ylabel('Measured Angle - Nominal Angle (degrees)')
    ax.set_title(title + ' Results - Angle')
    if legend:
        ax.legend(loc=lloc)
    fig.show()

def plot_diff_widths(idict=False,legend=True,lloc='upper right',EM=False,dump_dict=False):
    wins=[11,21,12,22,31,41,32,42,33,43,34,44]
    if EM:
        wins=[11,21,31,41,33,43,12,22,32,42,34,44]
        g='EM1'
    if EM:
        front_folder='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/EMassembly_2017_11_15_FrontGrid'
        rear_folder='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw469sub2765_2018_01_31'
    else:
        front_folder='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw561sub3501_2018_06_26'
        rear_folder='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw501sub3437_2018_05_04'
        g='EM2'

    fwdiff,fdcdiff,fwerr,fdcerr=[],[],[],[]
    rwdiff,rdcdiff,rwerr,rdcerr=[],[],[],[]

    if idict:
        fwdiff=idict['fwdiff']
        fdcdiff=idict['fdcdiff']
        fwerr=idict['fwerr']
        fdcerr=idict['fdcerr']

        rwdiff=idict['rwdiff']
        rdcdiff=idict['rdcdiff']
        rwerr=idict['rwerr']
        rdcerr=idict['rdcerr']

    else:
        for w in wins:
            oname='win'+str(w)+'_grating.p'
            os.chdir(front_folder)
            fg=pickle.load(open(oname,'rb'))
            os.chdir(rear_folder)
            rg=pickle.load(open(oname,'rb'))
            win=str(fg.win)

            if 'fdict' not in fg.results.keys():
                #it's win 42 or 43, no results for fits
                fwidth=np.nan
                fdc=np.nan
                fchisq=np.nan
            else:
                #front grid stuff
                fwidth=fg.results['fdict']['int_y']
                fpitch=fg.results['widths'].period_stats['mean']
                fdc=fwidth/fpitch
                try:
                    if len(fg.results['fdict']['plus_err'])==0:
                        fchisq=fg.results['fdict']['minus_err'][0]
                    elif len(fg.results['fdict']['minus_err'])==0:
                        fchisq=fg.results['fdict']['plus_err'][0]
                    else:
                        fchisq=np.mean([fg.results['fdict']['plus_err'],fg.results['fdict']['minus_err']])

                except KeyError:
                    fchisq=np.nan
            #front grid stuff
            fwdiff.append(fwidth-fg.nominal['pitch']/2.)
            fdcdiff.append(fdc-0.5)
            fwerr.append(fchisq)

            #front grid stuff
            rwidth=rg.results['fdict']['int_y']
            rpitch=rg.results['widths'].period_stats['mean']
            rdc=rwidth/rpitch
            try:
                if len(rg.results['fdict']['plus_err'])==0 and len(rg.results['fdict']['minus_err'])!=0:
                    rchisq=rg.results['fdict']['minus_err'][0]
                elif len(rg.results['fdict']['minus_err'])==0 and len(rg.results['fdict']['plus_err'])!=0:
                    rchisq=rg.results['fdict']['plus_err'][0]
                else:
                    rchisq=np.mean([rg.results['fdict']['plus_err'],rg.results['fdict']['minus_err']])
            except KeyError:
                rchisq=np.nan
            rwdiff.append(rwidth-rg.nominal['pitch']/2.)
            rdcdiff.append(rdc-0.5)
            rwerr.append(rchisq)

    fig,ax=plt.subplots(2,sharex=True,sharey=True)
    #plt.ylabel('Measured Width - Nominal Width ($\mu$m)')
    xvals=range(0,12)
    ax[0].scatter(xvals,fwdiff,c='c',marker='v',label='Duty Cycle')
    ax[0].scatter(xvals,fwdiff,c='k',label='Slit Width')
    #ax.fill_between(xrvals,-1*np.array(preq),np.array(preq),alpha=.6,color='b',label="Requirements")
    ax[0].errorbar(xvals,fwdiff,yerr=np.array(fwerr),fmt=None,ecolor='k',label="Fit Error")
    dcax0=ax[0].twinx()
    dcax0.scatter(xvals,fdcdiff,c='c',marker='v',s=40,label='Measured-Nominal Duty Cycle')
    dcax0.set_title(g+' Front Grid')

    ax[1].set_xticks(xvals)
    ax[1].set_xticklabels(wins) #not if xray...
    ax[1].set_xlabel('Window')
    #plt.ylabel('Measured Width - Nominal Width ($\mu$m)')
    ax[1].set_xlim(-1,len(wins))
    #ax[0].set_ylim(-1,1)
    #ax.set_title(title + ' Results - Pitch')
    if legend:
        ax[0].legend(loc=lloc)

    ax[1].scatter(xvals,rwdiff,c='k',label='Measured-Nominal Slit Width')
    #ax.fill_between(xrvals,-1*np.array(preq),np.array(preq),alpha=.6,color='b',label="Requirements")
    ax[1].errorbar(xvals,rwdiff,yerr=np.array(rwerr),fmt=None,ecolor='k',label="Fit Error")
    dcax1=ax[1].twinx()
    dcax1.scatter(xvals,rdcdiff,c='c',marker='v',s=40,label='Measured-Nominal Duty Cycle')
    fig.text(0.06, 0.5, 'Measured Width - Nominal Width ($\mu$m)', ha='center', va='center', rotation='vertical')
    fig.text(0.97, 0.5, 'Measured Duty Cycle - Nominal Duty Cycle', ha='center', va='center', rotation='vertical')
    #plt.ylabel('Measured Duty Cycle - Nominal Duty Cycle')
    dcax1.set_title(g+' Rear Grid')

    #ax.set_xticks(ticks=range(0,len(wins)),labels=wins)
    #ax.set_xticks(range(0,len(wins)))
    #ax.set_xticklabels(wins)
    #ax.set_xlim(-1,len(wins))
    #ax.set_ylim([ayran[0],ayran[1]])
    #ax.set_xlabel('Window')
    #ax.set_ylabel('Measured Angle - Nominal Angle (degrees)')
    #ax.set_title(title + ' Results - Angle')

    fig.show()

    if dump_dict:
        dd={'fwdiff':fwdiff,'fdcdiff':fdcdiff,'fwerr':fwerr,'fdcerr':fdcerr,'rwdiff':rwdiff,'rdcdiff':rdcdiff,'rwerr':rwerr,'rdcerr':rdcerr}
        return dd


