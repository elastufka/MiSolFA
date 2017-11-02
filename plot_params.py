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

def load_dict(front_dir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2737_2017_06_01',rear_dir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02',mag=5.0):
    windows=[{'side':'f','number':11,'npitch':0.12,  'nangle':-44.72509,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':21,'pitch':0.12,  'nangle': 45.27757,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':12,'pitch':0.03,'nangle': 44.93102,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':22,'pitch':0.03,'nangle':-45.06914,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':31,'pitch':0.06, 'nangle':-45.13845,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':41,'pitch':0.06, 'nangle': 44.86222,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':32,'pitch':0.025, 'nangle': 45.05530,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':42,'pitch':0.025, 'nangle':-44.94481,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':33,'pitch':0.04,  'nangle':-44.90807,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':43,'pitch':0.04,  'nangle': 45.09223,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':34,'pitch':0.020, 'nangle': 44.954,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'f','number':44,'pitch':0.02, 'nangle':-45.04608,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':11,'npitch':0.12,  'nangle':45.27757,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':21,'pitch':0.12,  'nangle': -44.72509,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':12,'pitch':0.03,'nangle': -45.06914,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':22,'pitch':0.03,'nangle':44.93102,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':31,'pitch':0.06, 'nangle':44.86222,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':41,'pitch':0.06, 'nangle': -45.13845,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':32,'pitch':0.025, 'nangle': -44.94481,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':42,'pitch':0.025, 'nangle':45.0553,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':33,'pitch':0.04,  'nangle':45.09223,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':43,'pitch':0.04,  'nangle': -44.90807,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':34,'pitch':0.02, 'nangle': -45.04608,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0},
                {'side':'r','number':44,'pitch':0.02, 'nangle':44.954,'pmean':0.0,'pmed':0.0,'pvar':0.0,'amean':0.0,'amed':0.0,'avar':0.0}]

    
    #frontdir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2737_2017_06_01'
    #reardir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02'

    os.chdir(front_dir)
    #statfiles=glob.glob('*stats*'+str(mag)+'.p')
    statfiles=glob.glob('*stats.p')
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
    #statfiles=glob.glob('*stats*'+str(mag)+'.p')
    statfiles=glob.glob('*stats.p')
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
    pickle.dump(windows,open('windows.p','wb'))
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
    pmeans=[window[i]['pmean'] for i in range(0,12)]
    perrs=[window[i]['pvar'] for i in range(0,12)]
    noms=[]
    for i in range(0,12):
        try:
            noms.append(window[i]['npitch']*1000)
        except KeyError:
            noms.append(window[i]['pitch']*1000)
    x=np.linspace(1,12,12)
    fig,ax=plt.subplots()
    ax.errorbar(x, np.array(pmeans)-np.array(noms), yerr=np.array(perrs),
            fmt='o', ecolor='g', capthick=2)
    #pl.plot(x, y, 'k', color='#CC4F1B')
    ax.fill_between(x,-error, error,alpha=0.5, facecolor='#FF9848')
    ax.set_xticks(x)
    ax.set_xticklabels(wnums)
    ax.set_title(title1)
    ax.set_xlabel('Window Number')
    ax.set_xlim([0,13])
    ax.set_ylim([-15,15])
    ax.set_ylabel('Mean Width - Nominal Width, $\mu$m')
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

        
    pmeans=[window[i]['pmean'] for i in range(12,24)]
    perrs=[window[i]['pvar'] for i in range(12,24)]
    noms=[]
    for i in range(0,12):
        try:
            noms.append(window[i]['npitch']*1000)
        except KeyError:
            noms.append(window[i]['pitch']*1000)
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
    return str(np.round([value],decimals=2)[0])

def dict2tex(window):
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
    with open('results_table.txt','wb') as f:
        f.write('\n'.join(lines))
    #print os.get_cwd()

