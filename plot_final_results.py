def plot_diffs():
    wins=[11,21,12,22,31,41,32,42,33,43,34,44]
    means,noms,yerr=[],[],[]
    x=np.linspace(1,12,12)
    for w in wins:
        dummy=Grating(w,Sept=True,EM=False, side=1)
        os.chdir(dummy.data.Odata.path)
        stats=pickle.load(open('win'+str(w)+'_width_stats5.0Xa.p','rb'))
        #stats=pickle.load(open('win'+str(w)+'_angle_stats_5.0.p','rb'))
        means.append(stats['mean'])
        noms.append(dummy.nominal['pitch'])
        #noms.append(dummy.nominal['orientation'])
        yerr.append(stats['stddev'])
    fig,ax=plt.subplots()
    ax.errorbar(x, np.array(means)-np.array(noms), yerr=np.transpose(np.array(yerr)),
            fmt='o', ecolor='g', capthick=2)
    fig.suptitle('Sub 3501 - QM Front Grid')
    ax.set_xticks(x)
    ax.set_xticklabels(wins)
    ax.set_xlabel('Window Number')
    ax.set_ylabel('Measured Period Difference ($\mu$m)')
    ax.set_xlim([0,13])
    #ax.set_ylim([-.5,.5])
    fig.show()

def print_transm_period_stats(pix2um=0.65):
    wins=[11,21,12,22,31,41,32,42,33,43,34,44]
    #EM front
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/SLS_Apr2018/transmEM')
    print '--------EM Front Window----------'
    sub='sub2765'
    for w in wins:
        if w !=12:
            try:
                tw=pickle.load(open(sub+'_win'+str(w)+ '.p','rb'))#restore pickles with the transm['widths'] results
            except IOError:
                os.chdir('../transmEM1')
                tw=pickle.load(open(sub+'_win'+str(w)+ '.p','rb'))

            mpp=tw.results['widths'].data['mean_periods_P']
            print w, np.mean(mpp)*pix2um,np.std(mpp)*pix2um
            os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/SLS_Apr2018/transmEM')
    #EM rear
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmEM')
    print '--------EM Rear Window----------'
    sub='sub2737'
    for w in wins:
        try:
            tw=pickle.load(open(sub+'_win'+str(w)+ '.p','rb'))#restore pickles with the transm['widths'] results
        except IOError:
            os.chdir('../transmEM2')
            tw=pickle.load(open(sub+'_win'+str(w)+ '.p','rb'))

        mpp=tw.results['widths'].data['mean_periods_P']
        print w, np.mean(mpp)*pix2um,np.std(mpp)*pix2um
        os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmEM')
    #QM Front
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018/transmQM')
    sub='sub3501'
    print '--------QM Front Window----------'
    for w in wins:
        try:
            tw=pickle.load(open(sub+'_win'+str(w)+ '.p','rb'))#restore pickles with the transm['widths'] results
        except IOError:
            os.chdir('../transmQMa')
            tw=pickle.load(open(sub+'_win'+str(w)+ '.p','rb'))

        mpp=tw.results['widths'].data['mean_periods_P']
        print w, np.mean(mpp)*pix2um,np.std(mpp)*pix2um
        os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018/transmQM')

    #QM Rear
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmQM')
    sub='sub3437'
    print '--------QM Rear Window----------'
    for w in wins:
        try:
            tw=pickle.load(open(sub+'_win'+str(w)+ '.p','rb'))#restore pickles with the transm['widths'] results
        except IOError:
            os.chdir('../transmQMa')
            tw=pickle.load(open(sub+'_win'+str(w)+ '.p','rb'))

        mpp=tw.results['widths'].data['mean_periods_P']
        print w, np.mean(mpp)*pix2um,np.std(mpp)*pix2um
        os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmQM')

