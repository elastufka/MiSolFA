"""
===================
plot_params.py
Erica  Lastufka 16.10.17
===================

Get data to dictionaries. plot.

"""
import glob
import pickle
import os
#from edge_and_hough import windows as windowsf
#from edge_and_hough import windowsr

#print windowsf[0]['pitch'],windowsf[0]['nominal angle']
def load_dict(front_dir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2737_2017_06_01',rear_dir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02',mag=5.0):
    windowsf=[{'number':11,'pitch':89.6752,  'nominal angle':-44.79357},
                {'number':21,'pitch':90.326,  'nominal angle': 45.20793},
                {'number':12,'pitch':22.4797,'nominal angle': 44.94825},
                {'number':22,'pitch':22.5203,'nominal angle':-45.05184},
                {'number':31,'pitch':45.0814, 'nominal angle':-45.10378},
                {'number':41,'pitch':44.9187, 'nominal angle': 44.8966},
                {'number':32,'pitch':18.013, 'nominal angle': 45.04146},
                {'number':42,'pitch':17.987, 'nominal angle':-44.95859},
                {'number':33,'pitch':29.9639,  'nominal angle':-44.93102},
                {'number':43,'pitch':30.0362,  'nominal angle': 45.06914},
                {'number':34,'pitch':14.991, 'nominal angle': 44.96549},
                {'number':44,'pitch':15.009, 'nominal angle':-45.03455}]

    windowsr=[{'number':11,'pitch':90.326,  'nominal angle':-45.20793},
                {'number':21,'pitch':89.6752,  'nominal angle': 44.79357},
                {'number':12,'pitch':22.5203,'nominal angle': 45.05184},
                {'number':22,'pitch':22.4797,'nominal angle':-44.94825},
                {'number':31,'pitch':44.9187, 'nominal angle':-44.8966},
                {'number':41,'pitch':45.0814, 'nominal angle': 45.10378},
                {'number':32,'pitch':17.987, 'nominal angle': 44.95859},
                {'number':42,'pitch':18.013, 'nominal angle':-45.04146},
                {'number':33,'pitch':30.0362,  'nominal angle':-45.06914},
                {'number':43,'pitch':29.9639,  'nominal angle': 44.93102},
                {'number':34,'pitch':15.009, 'nominal angle': 45.03455},
                {'number':44,'pitch':14.991, 'nominal angle':-44.96549}] #dectector side angle for ease, windows swapped...

    
    windows=[{'side':'f','number':11,'pitch':windowsf[0]['pitch'],  'nangle':windowsf[0]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':21,'pitch':windowsf[1]['pitch'],  'nangle': windowsf[1]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':12,'pitch':windowsf[2]['pitch'],'nangle': windowsf[2]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':22,'pitch':windowsf[3]['pitch'],'nangle':windowsf[3]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':31,'pitch':windowsf[4]['pitch'], 'nangle':windowsf[4]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':41,'pitch':windowsf[5]['pitch'], 'nangle': windowsf[5]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':32,'pitch':windowsf[6]['pitch'], 'nangle':windowsf[6]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':42,'pitch':windowsf[7]['pitch'], 'nangle':windowsf[7]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':33,'pitch':windowsf[8]['pitch'],  'nangle':windowsf[8]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':43,'pitch':windowsf[9]['pitch'],  'nangle': windowsf[9]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':34,'pitch':windowsf[10]['pitch'], 'nangle': windowsf[10]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':44,'pitch':windowsf[11]['pitch'], 'nangle':windowsf[11]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':11,'pitch':windowsr[0]['pitch'],  'nangle':windowsr[0]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':21,'pitch':windowsr[1]['pitch'],  'nangle': windowsr[1]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':12,'pitch':windowsr[2]['pitch'],'nangle':windowsr[2]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':22,'pitch':windowsr[3]['pitch'],'nangle':windowsr[3]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':31,'pitch':windowsr[4]['pitch'], 'nangle':windowsr[4]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':41,'pitch':windowsr[5]['pitch'], 'nangle':windowsr[5]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':32,'pitch':windowsr[6]['pitch'], 'nangle':windowsr[6]['nominal angle'] ,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':42,'pitch':windowsr[7]['pitch'], 'nangle':windowsr[7]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':33,'pitch':windowsr[8]['pitch'],  'nangle':windowsr[8]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':43,'pitch':windowsr[9]['pitch'],  'nangle':windowsr[9]['nominal angle'] ,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':34,'pitch':windowsr[10]['pitch'], 'nangle': windowsr[10]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':44,'pitch':windowsr[11]['pitch'], 'nangle':windowsr[1]['nominal angle'],'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0}]

    
    #frontdir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2737_2017_06_01'
    #reardir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02'

    os.chdir(front_dir)
    statfiles=glob.glob('*stats*'+str(mag)+'.p')
    #statfiles=glob.glob('*stats.p')
    for sf in statfiles:
        win=sf[3:5]
        idx=[i for i in range(0,12) if windows[i]['number']==int(win)][0]
        sdict=pickle.load(open(sf,'rb'))
        if 'width' in sf:
            windows[idx]['pmean']=sdict['mean']
            windows[idx]['pmed']=sdict['median']
            windows[idx]['pvar']=sdict['stddev']
        else:
            windows[idx]['amean']=sdict['mean']
            windows[idx]['amed']=sdict['median']
            windows[idx]['avar']=sdict['stddev']
            
    os.chdir(rear_dir)
    statfiles=glob.glob('*stats*'+str(mag)+'.p')
    #statfiles=glob.glob('*stats.p')
    for sf in statfiles:
        win=sf[3:5]
        idx=[i for i in range(12,24) if windows[i]['number']==int(win)][0]
        sdict=pickle.load(open(sf,'rb'))
        if 'width' in sf:
            windows[idx]['pmean']=sdict['mean']
            windows[idx]['pmed']=sdict['median']
            windows[idx]['pvar']=sdict['stddev']
        else:
            windows[idx]['amean']=sdict['mean']
            windows[idx]['amed']=sdict['median']
            windows[idx]['avar']=sdict['stddev']

    os.chdir('../')
    pickle.dump(windows,open('windows'+str(mag)+'.p','wb'))
    return windows

def repickle():
    frontdir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2737_2017_06_01'
    reardir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02'

    os.chdir(frontdir)
    statfiles=glob.glob('*stats.p')
    for sf in statfiles:
        sdict=pickle.load(open(sf,'rb'))
        data=sdict['data']
        stats={'mean':sdict['mean'],'median':sdict['median'],'stddev':sdict['stddev']}
        dname=sf[:-7]+'data.p'
        pickle.dump(data,open(dname,'wb'))
        pickle.dump(stats,open(sf,'wb'))
        
    os.chdir(reardir)
    statfiles=glob.glob('*stats.p')
    for sf in statfiles:
        sdict=pickle.load(open(sf,'rb'))
        data=sdict['data']
        stats={'mean':sdict['mean'],'median':sdict['median'],'stddev':sdict['stddev']}
        dname=sf[:-7]+'data.p'
        pickle.dump(data,open(dname,'wb'))
        pickle.dump(stats,open(sf,'wb'))

    os.chdir('../')

def plot_params(window,title1='Front Assembly',title2='Rear Assembly',ptype=False):
    error=np.zeros(12)+1.995

    #front window - pitch vs. expected values
    wnums=[window[i]['number'] for i in range(0,12)]
    nomsa=[window[i]['nangle'] for i in range(0,12)]
    pmeans=[window[i]['pmean'] for i in range(0,12)]
    print pmeans
    perrs=[window[i]['pvar'] for i in range(0,12)]
    noms=[]
    for i in range(0,12):
        noms.append(window[i]['pitch'])
    x=np.linspace(1,12,12)
    fig,ax=plt.subplots()
    ax.errorbar(x, np.array(pmeans)-np.array(noms), yerr=np.array(perrs),
            fmt='o', ecolor='g', capthick=2)
    #pl.plot(x, y, 'k', color='#CC4F1B')
    #ax.fill_between(x,-error, error,alpha=0.5, facecolor='#FF9848')
    ax.set_xticks(x)
    ax.set_xticklabels(wnums)
    ax.set_title(title1)
    ax.set_xlabel('Window Number')
    ax.set_xlim([0,13])
    ax.set_ylim([-2,2])
    ax.set_ylabel('Mean Pitch - Nominal Pitch, $\mu$m')
    fig.show()
    #return noms
    
    #front window - angle vs. expected values
    #wnums=[window[i]['number'] for i in range(0,12)]
    ameds=[np.abs(window[i]['amean']) for i in range(0,12)]
    aerrs=[window[i]['avar'] for i in range(0,12)]
    noms=[window[i]['nangle'] for i in range(0,12)]
    x=np.linspace(1,12,12)
    fig,ax=plt.subplots()
    ax.errorbar(x, ameds, yerr=aerrs,
            fmt='o', ecolor='g', capthick=2,label='Mean Value')
    ax.scatter(x,np.abs(np.array(noms)),color='r',marker='v',s=40, label='Nominal Value')
    #pl.plot(x, y, 'k', color='#CC4F1B')
    #ax.fill_between(x,noms-error, noms+error,alpha=0.5, facecolor='#FF9848')
    #ax.set_xtickinterval(x[:-1])
    ax.set_xticks(x)
    ax.set_xticklabels(wnums)
    ax.set_title(title1)
    ax.set_xlabel('Window Number')
    ax.set_xlim([0,13])
    ax.set_ylim([44.25,45.75])
    ax.set_ylabel('Absolute Value of Slat Orientation in degrees')
    ax.legend()
    fig.show()
    
    #rear window - pitch vs. expected values
    wnums=[window[i]['number'] for i in range(12,24)]
    if ptype: #shuffle the rear windows to correspond to the correct numbers
        real21=window[13]
        print real21['amean']
        real11=window[12]
        print real11['amean']       
        real22=window[14]
        real12=window[15]
        real41=window[16]
        real31=window[17]
        real42=window[18]
        real32=window[19]
        real43=window[20]
        real33=window[21]
        real44=window[22]
        real34=window[23]
        window[12]['number']=21
        window[12]['nangle']=-44.72509
        window[13]['number']=11
        window[13]['nangle']=45.27757
        window[14]['number']=22
        window[14]['nangle']=44.93102
        window[15]['number']=12
        window[15]['nangle']=-45.06914
        window[16]['number']=41
        window[16]['nangle']=-45.13845    
        window[17]['number']=31
        window[17]['nangle']=44.86222
        window[18]['number']=42
        window[18]['nangle']=45.0553
        window[19]['number']=32
        window[19]['nangle']=-44.94481
        window[20]['number']=43
        window[20]['nangle']=-44.90807
        window[21]['number']=33
        window[21]['nangle']=45.09223
        window[22]['number']=44
        window[22]['nangle']= 44.954 
        window[23]['number']=34
        window[23]['nangle']=-45.04608
        nwins=[window[0],window[1],window[2],window[3],window[4],window[5],window[6],window[7],window[8],window[9],window[10],window[11],window[13],window[12],window[15],window[14],window[17],window[16],window[19],window[18],window[21],window[20],window[23],window[22]]
        window=nwins

    nomsa=[window[i]['nangle'] for i in range(12,24)]        
    pmeans=[window[i]['pmean'] for i in range(12,24)]
    perrs=[window[i]['pvar'] for i in range(12,24)]
    noms=[]
    for i in range(0,12):
        noms.append(window[i]['pitch'])
    x=np.linspace(1,12,12)
    fig,ax=plt.subplots()
    ax.errorbar(x, np.array(pmeans)-np.array(noms), yerr=np.array(perrs),
            fmt='o', ecolor='g', capthick=2)
    #pl.plot(x, y, 'k', color='#CC4F1B')
    ax.fill_between(x,-error, error,alpha=0.5, facecolor='#FF9848')
    #ax.set_xtickinterval(x[:-1])
    ax.set_xticks(x)
    ax.set_xticklabels(wnums)
    ax.set_title(title2)
    ax.set_xlabel('Window Number')
    ax.set_xlim([0,13])
    ax.set_ylim([-15,15])
    ax.set_ylabel('Mean Width - Nominal Width, $\mu$m')
    fig.show()
    #return noms
    
    #rear window - angle vs. expected values
    #wnums=[window[i]['number'] for i in range(0,12)]
    ameds=[np.abs(window[i]['amean']) for i in range(12,24)]
    aerrs=[window[i]['avar'] for i in range(12,24)]
    noms=[window[i]['nangle'] for i in range(12,24)]
    x=np.linspace(1,12,12)
    fig,ax=plt.subplots()
    ax.errorbar(x, ameds, yerr=aerrs,
            fmt='o', ecolor='g', capthick=2,label='Mean Value')
    ax.scatter(x,np.abs(np.array(noms)),color='r',marker='v',s=40,label='Nominal Value')
    #pl.plot(x, y, 'k', color='#CC4F1B')
    #ax.fill_between(x,noms-error, noms+error,alpha=0.5, facecolor='#FF9848')
    #ax.set_xtickinterval(x[:-1])
    ax.set_xticks(x)
    ax.set_xticklabels(wnums)
    ax.set_title(title2)
    ax.set_xlabel('(Rear) Window Number')
    ax.set_xlim([0,13])
    ax.set_ylim([44.25,45.75])#    ax.set_ylim([44.25,46.75])
    ax.set_ylabel('Absolute Value of Slat Orientation in degrees')
    #ax.legend()
    fig.show()
    return window

def rv(value):
    return str(np.round([value],decimals=4)[0])

def dict2tex(window,filename=False):
    front=window[:12]
    rear=window[12:]
    fwlines,falines,rwlines,ralines=[],[],[],[]
    window_numbers=['11','21','12','22','31','41','32','42','33','43','34','44']
    for i,j,win in zip(front,rear,window_numbers):
        fwline=win+' & & '+rv(i['pmean'])+' & '+rv(i['pmed'])+' & '+rv(i['pvar'])+' & '
        faline=rv(i['amean'])+' & '+rv(i['amed'])+' & '+rv(i['avar'])+' & '
        rwline=rv(j['pmean'])+' & '+rv(j['pmed'])+' & '+rv(j['pvar'])+' & '
        raline=rv(j['amean'])+' & '+rv(j['amed'])+' & '+rv(j['avar'])+' \\'
        fwlines.append(fwline)
        falines.append(faline)
        rwlines.append(rwline)
        ralines.append(raline)
           
    lines=[fwlines[i]+falines[i]+rwlines[i]+ralines[i] for i in range(0,len(fwlines))]
    os.chdir('../')
    if not filename:
        filename='results_table.txt'
    with open(filename,'wb') as f:
        f.write('\n'.join(lines))
    #print os.get_cwd()

