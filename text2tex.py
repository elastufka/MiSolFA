"""
===================
text2tex.py
Erica  Lastufka 16.10.17
===================

make it a latex table

"""
import glob
import pickle
import os

def rv(value):
    return str(np.round([value],decimals=4)[0])

def make_full_table():
    frontdir='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw561sub3501_2018_06_26'
    reardir='/Users/wheatley/Documents/Solar/MiSolFA/OpticalAnalysis/mw501sub3437_2018_05_04'
    fwlines,falines,rwlines,ralines=[],[],[],[]
    os.chdir(frontdir)
    statfiles=glob.glob('*stats5.0Xa.p')
    statfiles.extend(glob.glob('*_angle_stats_5.0.p'))
    statfiles.sort()
    print statfiles
    for sf in statfiles:
        win=sf[3:5]
        sdict=pickle.load(open(sf,'rb'))
        if 'width' in sf:
            wline=win+' & & '+rv(sdict['mean'])+' & '+rv(sdict['median'])+' & '+rv(sdict['stddev'])+' & '
            #print wline
            fwlines.append(wline)
        else:
            aline=rv(sdict['mean'])+' & '+rv(sdict['median'])+' & '+rv(sdict['stddev'])+' & '
            #print aline
            falines.append(aline)
    os.chdir(reardir)
    statfiles=glob.glob('*stats5.0Xa.p')
    statfiles.extend(glob.glob('*_angle_stats_5.0.p'))
    statfiles.sort()
    print statfiles
    for sf in statfiles:
        win=sf[3:5]
        sdict=pickle.load(open(sf,'rb'))
        if 'width' in sf:
            wline=rv(sdict['mean'])+' & '+rv(sdict['median'])+' & '+rv(sdict['stddev'])+' & '
            #print wline
            rwlines.append(wline)
        else:
            aline=rv(sdict['mean'])+' & '+rv(sdict['median'])+' & '+rv(sdict['stddev'])+' \\\n'
            ralines.append(aline)
            #print aline

    #let's hope things are ordered:

    #means
    #meanw=[np.mean(float(fwlines))]
    #meana=
    lines=[fwlines[i]+falines[i]+rwlines[i]+ralines[i] for i in range(0,len(fwlines))]
    ll=' '.join(lines)
    os.chdir('../')
    with open('QMresults_table.txt','wb') as f:
        f.write(ll)

def fit_table():
    frontdir='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018'
    reardir='/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018'
    fwlines,falines,rwlines,ralines=[],[],[],[]
    os.chdir(frontdir)
    statfiles=glob.glob('sub3501_win*fit_vdict.p')
    #statfiles.extend(glob.glob('*_angle_stats_5.0.p'))
    #statfiles.sort()
    print statfiles
    for sf in statfiles:
        win=sf[11:13]
        vd=pickle.load(open(sf,'rb'))
        wline=win+' & '+rv(np.mean(vd['width_yvals']))+' & '+rv(np.mean(vd['period_yvals']))+' & '+rv(np.mean(vd['period_yvals'])-np.mean(vd['pitch_nom']))+' & '+rv(np.mean(vd['height_yvals']))+' & '+rv(np.mean(vd['dc']))+' & '+rv(np.mean(vd['yoff']))
        fwlines.append(wline)
    os.chdir(reardir)
    statfiles=glob.glob('sub3437_win*fit_vdict.p')
    #statfiles.extend(glob.glob('*_angle_stats_5.0.p'))
    #statfiles.sort()
    print statfiles
    for sf in statfiles:
        win=sf[11:13]
        vd=pickle.load(open(sf,'rb'))
        wline=' & '+rv(np.mean(vd['width_yvals']))+' & '+rv(np.mean(vd['period_yvals']))+' & '+rv(np.mean(vd['period_yvals'])-np.mean(vd['pitch_nom']))+' & '+rv(np.mean(vd['height_yvals']))+' & '+rv(np.mean(vd['dc']))+' & '+rv(np.mean(vd['yoff']))+'\\\n'
            #print wline
        rwlines.append(wline)

    #let's hope things are ordered:
    #means
    #meanw=[np.mean(float(fwlines))]
    #meana=
    lines=[fwlines[i]+rwlines[i] for i in range(0,len(fwlines))]
    ll=' '.join(lines)
    os.chdir('../')
    with open('fit_results_table_QMall.txt','wb') as f:
        f.write(ll)

