"""
===================
rename_beamfiles.py
Erica  Lastufka 12.4.18
===================

Rename the files output from the beamline program to something more useable for me

original name: transm_step000_window011_0001.tif

step is pretty much meaningless
final number=position + angle tilt
eg. 11 angular steps * 7 positions = 77 total steps

"""
import numpy as np
import os
import glob


seq1=[-8.75,-7,-5.25,-3.5,-1.75,0,1.75,3.5,5.25,7,8.75]
seq2=[-5,-4,-3,-2,-1,0,1,2,3,4,5]
seq3=[-3,-2.25,-1.5,-0.75,0,0.75,1.5,2.25,3]
seq4=[-2.25,-1.5,-0.75,0,0.75,1.5,2.25]
seq5=[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]
seq6=[-1.5,-1,-0.5,0,0.5,1,1.5]

#pdict0={'11':seq1,'22':seq1,'31':seq2,'42':seq2,'34':seq3,'43':seq3,'21':seq4,'12':seq4,'32':seq5,'41':seq5,'44':seq6,'33':seq6}
pdictEM={'11':seq1,'21':seq1,'31':seq2,'41':seq2,'33':seq3,'43':seq3,'22':seq4,'12':seq4,'32':seq5,'42':seq5,'44':seq6,'34':seq6}
pdict1={'22':seq4,'41':seq2,'21':seq1,'33':seq3,'42':seq6,'34':seq6} #10/5 42 should have seq5 .. check...apparently not...

def parse_tail(tailval,seq):
    '''figure out what the tail number actually means'''
    p=tailval/(len(seq))
    ang=seq[tailval % len(seq)]
    return p,ang    

if __name__ == '__main__':
    #os.chdir('SLS_Apr2018/disk2/transm')
    os.chdir('/Volumes/LaCie/MiSolFA/transm1')
    #os.chdir('SLS_May2018/transmEM2')
    #os.chdir('transm')
    files=glob.glob('*.tif')
    names=''
    for f in files:#[100:120]:
        win=f[-11:-9]
        tailval=int(f[-8:-4]) #DO I NEED TO SUBTRACT 1 FROM THIS??
        #p,ang=parse_tail(tailval-1,pdictEM[win])
        p,ang=parse_tail(tailval-1,pdict1[win])
        #newname='../SLS_Apr2018/disk2/transm/win'+win+'_p'+str(p)+'_'+str(ang)+'.tif'
        newname=f+':'+'win'+win+'_p'+str(p)+'_'+str(ang)+'.tif,'
        #print f,p,ang,newname
        #os.rename(f,newname)
        names=names+newname
    os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/BeamTest_2018-04-11/SLS_Apr2018/disk2/transm')
    with open('filenames.csv','w') as f:
        f.write(names)
        
