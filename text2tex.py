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
    return str(np.round([value],decimals=2)[0])

if __name__ == '__main__':
    frontdir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2737_2017_06_01'
    reardir='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02'
    fwlines,falines,rwlines,ralines=[],[],[],[]
    os.chdir(frontdir)
    statfiles=glob.glob('*stats.p')
    for sf in statfiles:
        win=sf[3:5]
        sdict=pickle.load(open(sf,'rb'))
        if 'width' in sf:
            wline=win+' & & '+rv(sdict['mean'])+' & '+rv(sdict['median'])+' & '+rv(sdict['stddev'])+' & '
            print wline
            fwlines.append(wline)
        else:
            aline=rv(sdict['mean'])+' & '+rv(sdict['median'])+' & '+rv(sdict['stddev'])+' & '
            print aline
            falines.append(aline)
    os.chdir(reardir)
    statfiles=glob.glob('*stats.p')
    for sf in statfiles:
        win=sf[3:5]
        sdict=pickle.load(open(sf,'rb'))
        if 'width' in sf:
            wline=rv(sdict['mean'])+' & '+rv(sdict['median'])+' & '+rv(sdict['stddev'])+' & '
            print wline
            rwlines.append(wline)
        else:
            aline=rv(sdict['mean'])+' & '+rv(sdict['median'])+' & '+rv(sdict['stddev'])+' \\\n'
            ralines.append(aline)
            print aline
            
    #let's hope things are ordered:

    #means
    #meanw=[np.mean(float(fwlines))]
    #meana=
    lines=[fwlines[i]+falines[i]+rwlines[i]+ralines[i] for i in range(0,len(fwlines))]
    os.chdir('../')
    with open('results_table.txt','wb') as f:
        f.write(lines)
