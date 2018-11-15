"""
===================
transm_all.py
Erica  Lastufka 18.9.18
===================
Run everything
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from PIL import Image
from matplotlib import cm
import pickle
from scipy.misc import imrotate
from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit
from scipy import interpolate
import time
import glob
import random
import itertools
import analyze_general as ag

wins=[11,12,21,22,31,41,32,42,33,43,34,44]#[11,12,21,22,31,#
# Check these images:
# ['win22_p1_-1.5_corrected.tif', 'win22_p3_-2.25_corrected.tif', 'win22_p3_-1.5_corrected.tif', 'win22_p3_-0.75_corrected.tif', 'win22_p3_0.75_corrected.tif', 'win22_p3_1.5_corrected.tif', 'win22_p3_2.25_corrected.tif', 'win22_p6_-2.25_corrected.tif', 'win22_p6_-1.5_corrected.tif', 'win22_p6_-0.75_corrected.tif', 'win22_p6_0_corrected.tif', 'win22_p6_0.75_corrected.tif', 'win22_p6_1.5_corrected.tif', 'win22_p6_2.25_corrected.tif']

# Check these images:
# ['win33_p1_-0.75_corrected.tif', 'win33_p3_-1.5_corrected.tif', 'win33_p4_-3_corrected.tif', 'win33_p4_-1.5_corrected.tif', 'win33_p4_-0.75_corrected.tif', 'win33_p5_1.5_corrected.tif', 'win33_p6_-0.75_corrected.tif', 'win33_p6_0_corrected.tif', 'win33_p6_2.25_corrected.tif']

def readflags(flagfile):
    flaglist=[[],[],[],[],[],[],[],[],[],[],[],[]]
    with open(flagfile) as f:
       lines=f.readlines()
    for l in lines:
       if l != '\n':
           ll=l[:-1].split(',')
           #ll[2]=ll[2][:-2]
           windex=wins.index(int(ll[0]))
           flaglist[windex].append(ll[1:]) #need int if just p
    return flaglist


wins=[11,12,21,22,31,41,32,42,33,34,43,44] #41 wasn't run?
#tol=[5,5,2,2,5,5,] #vector that changes or is constant based on the tilt angle
sigmas=[2.,1.25,2.,.85,1.5,1.5,.75,0.8,1.5,1.2,.75,.75]#,1.0,1.2,1.,1.]#[2.,2.,1.5,1.5,2.,2.,1.5,1.5,2.,2.,1.,1.] #one per each window
tolvecs=[[9,9,7,7,7,5,7,7,7,9,9],[5,5,5,6,6,6,5,5,5],[9,9,7,7,7,5,7,7,7,9,9],[5,5,5,6,6,6,5,5,5],[5,5,5,7,7,9,7,7,7,5,5,5],[5,5,5,5,5,5,5,5,5,5,5,5],[5,5,5,6,6,6,5,5,5],[5,5,5,6,6,6,5,5,5],[5,5,5,6,6,6,5,5,5],[5,5,5,6,6,6,5,5,5],[4,4,4,5,5,5,4,4,4],[4,4,4,5,5,5,4,4,4]]#,[5,5,5,6,6,6,5,5,5],[5,5,5,6,6,6,5,5,5],[],[]]
Aflags=readflags('/Users/wheatley/Documents/Solar/MiSolFA/code/MiSolFA/Aprilflags.txt')
#April
# os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/SLS_Apr2018/transmEM')

# #for w,f,s in zip(wins,Aflags,sigmas):
# for w,s,t,f in zip(wins,sigmas,tolvecs,Aflags):
# #w=44
# #i=11
# #s=sigmas[i]
# #f=Aflags[i]
#   #t=tolvecs[i]
#   #   o=Grating(w,Apr=True,May=False)#pickle.load(open('win'+str(w)+'.p','rb'))#Grating(w,May=True)
#     base='sub2765_win'+str(w)
#   #   bfigname=base+'_binary.png'
# #   sfigname=base+'_sum.png'
#     wfigname=base+'_widths.png'
#     ffigname=base+'_widths_filtered.png'
#     if f != []:
#         for fl in f:
#             if fl[0] =='p':
#                o.flag_data(p=int(fl[1]))
#             else:
#                o.flag_data(p0ang=fl[1])
# #   o.parameterize_transmission('binary',plot=False)
# #   o.results['binary'].scat_allP(figname=bfigname)
# #   o.parameterize_transmission('sum',plot=False)
# #   o.results['sum'].scat_allP(figname=sfigname)
# #   o.results['widths'].scat_allP(figname=wfigname,yran=[.1,.6])
#     #o.results['widths'].filter_outliers_widths()
#     #o.results['widths'].scat_allP(figname=ffigname,yran=[.1,.6])
# #   o.compare_transm_methods(figname=base+'_compare.png')
# #now fit each individual transmission profile, compile and save the results
# #def fit_profile_to_data(self,pindex,period,height,plot=True,fixed='height'):
#     try:
#        o=pickle.load(open(base+'.p','rb'))
#     except IOError:
#        os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/SLS_Apr2018/transmEM1')
#        o=pickle.load(open(base+'.p','rb'))

#     o.parameterize_transmission('widths',plot=False,sigma=s,tolvec=t) #add sigmas here
#     o.results['widths'].scat_allP(figname=wfigname,yran=[.1,.6])
#     o.results['widths'].filter_outliers_widths()
#     o.results['widths'].scat_allP(figname=ffigname,yran=[.1,.6])
#        #for pindex in range(0,7):
#        #    o.results['widths'].fit_profile_to_data_optimize(pindex,o.nominal['pitch'],o.nominal['height'],plot=False, fix_dc=False)
#     pickle.dump(o,open(base+'.p','wb'))
#        #o.results['widths'].fit_plots(figname=base+'_fits.png')
#     os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/EMmodel/SLS_Apr2018/transmEM')

#May QM and Sept - different window order
sigmas=[2.,1.5,2.,1.5,1.25,1.25,1.25,1.25,.75,.75,.75,.75]#,1.0,1.2,1.,1.]#[2.,2.,1.5,1.5,2.,2.,1.5,1.5,2.,2.,1.,1.] #one per each window
tolvecs=[[9,9,7,7,7,5,7,7,7,9,9],[6,5,5,5,6,6,6,5,5,5,6],[9,9,7,7,7,5,7,7,7,9,9],[6,5,5,5,6,6,6,5,5,5,6],[5,5,5,7,7,9,7,7,7,5,5,5],[5,5,5,5,5,5,5,5,5,5,5,5],[5,5,5,6,6,6,5,5,5],[5,5,5,6,6,6,5,5,5],[5,5,5,6,6,6,5,5,5],[5,5,5,6,6,6,5,5,5],[4,4,4,5,5,5,4,4,4],[4,4,4,5,5,5,4,4,4]]#,[5,5,5,6,6,6,5,5,5],[5,5,5,6,6,6,5,5,5],[],[]]


#wins=[42,33,43,34,44]
# 21: P0: 0 (-1.75)? -> p0 is out of sequence!
# P1: -1.75 ghost. Out of sequence! 1.75,3.5,5.25 ,7 bad. Flag p1!
# P2: -8.75 bad. Maybe the whole window is out of sequence?
# p5 :8.75 (7 is probably ok)
# 11:P3,p4 : -8.5 partially blocked
# P5: 5.25,7,8.5 partially blocked (last 2 the worst)
# Mflags=[[['p',0],['p',1],['p0_ang','p2_-8.75']],[['p0ang','p3_-8.5'],['p0ang','p4_-8.5'],['p0ang','p5_-7'],['p0ang','p5_-8.5']],False,False,False,False,False,False,False,False,False]
# #May
Mflags=readflags('/Users/wheatley/Documents/Solar/MiSolFA/code/MiSolFA/MayQMflags.txt')
# os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmEM')
# for w,s,t,f in zip(wins,sigmas,tolvecs,Mflags):
#   #o=Grating(w,May=True, EM=True)#pickle.load(open('win'+str(w)+'.p','rb'))#Grating(w,May=True)
#   if w == 12:
#       continue
#   base='sub2737_win'+str(w)
#   try:
#     o=pickle.load(open(base+'.p','rb'))
#   except IOError:
#     os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmEM2')
#     o=pickle.load(open(base+'.p','rb'))

#   base='sub2737_win'+str(w)
#   bfigname=base+'_binary.png'
#   sfigname=base+'_sum.png'
#   wfigname=base+'_widths.png'
#   ffigname=base+'_widths_filtered.png'
#   if f:
#      for fl in f:
#         if fl[0] =='p':
#             o.flag_data(p=fl[1])
#         else:
#             o.flag_data(p0ang=fl[1])
#   #o.parameterize_transmission('binary',plot=False)
#   #o.results['binary'].scat_allP(figname=bfigname)
#   #o.parameterize_transmission('sum',plot=False)
#   #o.results['sum'].scat_allP(figname=sfigname)
#   o.parameterize_transmission('widths',plot=False,sigma=s,tolvec=t)
#   o.results['widths'].scat_allP(figname=wfigname,yran=[.1,.6])
#   o.results['widths'].filter_outliers_widths()
#   o.results['widths'].scat_allP(figname=ffigname,yran=[.1,.6])

# #   o.results['widths'].filter_outliers_widths()
# #   o.results['widths'].scat_allP(figname=ffigname,yran=[.1,.6])
# #   o.compare_transm_methods(figname=base+'_compare.png')
#   pickle.dump(o,open(base+'.p','wb'))
#   os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmEM')

os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmQM')
for w,s,t,f in zip(wins,sigmas,tolvecs,Mflags):
  #o=Grating(w,May=True, EM=True)#pickle.load(open('win'+str(w)+'.p','rb'))#Grating(w,May=True)
  base='sub3437_win'+str(w)
  o=pickle.load(open(base+'.p','rb'))
  o.__get_and_store_data__()
  bfigname=base+'_binary.png'
  sfigname=base+'_sum.png'
  wfigname=base+'_widths.png'
  ffigname=base+'_widths_filtered.png'
  if f !=[]:
     for fl in f:
        if fl[0] =='p':
            o.flag_data(p=fl[1])
        else:
            o.flag_data(p0ang=fl[1])
  #o.parameterize_transmission('binary',plot=False)
  #o.results['binary'].scat_allP(figname=bfigname)
  #o.parameterize_transmission('sum',plot=False)
  #o.results['sum'].scat_allP(figname=sfigname)
  o.parameterize_transmission('widths',plot=False,sigma=s,tolvec=t)
  o.results['widths'].scat_allP(figname=wfigname,yran=[.1,.6])
  o.results['widths'].filter_outliers_widths()
  o.results['widths'].scat_allP(figname=ffigname,yran=[.1,.6])

# #   o.results['widths'].filter_outliers_widths()
# #   o.results['widths'].scat_allP(figname=ffigname,yran=[.1,.6])
# #   o.compare_transm_methods(figname=base+'_compare.png')
  pickle.dump(o,open(base+'.p','wb'))
  os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_May2018/transmQM')


# for w in wins:
#   o=Grating(w,May=True, EM=False)#pickle.load(open('win'+str(w)+'.p','rb'))#Grating(w,May=True)
#   base='sub3437_win'+str(w)
#   bfigname=base+'_binary.png'
#   sfigname=base+'_sum.png'
#   wfigname=base+'_widths.png'
#   ffigname=base+'_widths_filtered.png'
#   o.parameterize_transmission('binary',plot=False)
#   o.results['binary'].scat_allP(figname=bfigname)
#   o.parameterize_transmission('sum',plot=False)
#   o.results['sum'].scat_allP(figname=sfigname)
#   o.parameterize_transmission('widths',plot=False)
#   o.results['widths'].scat_allP(figname=wfigname,yran=[.1,.6])
#   o.results['widths'].filter_outliers_widths()
#   o.results['widths'].scat_allP(figname=ffigname,yran=[.1,.6])
#   o.compare_transm_methods(figname=base+'_compare.png')
#   pickle.dump(o,open(base+'.p','wb'))

# wins=[33,43,34,44]

# Sflags=readflags('/Users/wheatley/Documents/Solar/MiSolFA/code/MiSolFA/Septflags.txt')
# os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018/transmQM')
# for w,s,t,f in zip(wins,sigmas,tolvecs,Sflags):
#   #o=Grating(w,May=True, EM=True)#pickle.load(open('win'+str(w)+'.p','rb'))#Grating(w,May=True)
#   base='sub3501_win'+str(w)
#   try:
#     o=pickle.load(open(base+'.p','rb'))
#   except IOError:
#     os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018/transmQMa')
#     o=pickle.load(open(base+'.p','rb'))

#   bfigname=base+'_binary.png'
#   sfigname=base+'_sum.png'
#   wfigname=base+'_widths.png'
#   ffigname=base+'_widths_filtered.png'
#   if f !=[]:
#      for fl in f:
#         if fl[0] =='p':
#             o.flag_data(p=fl[1])
#         else:
#             o.flag_data(p0ang=fl[1])
#   #o.parameterize_transmission('binary',plot=False)
#   #o.results['binary'].scat_allP(figname=bfigname)
#   #o.parameterize_transmission('sum',plot=False)
#   #o.results['sum'].scat_allP(figname=sfigname)
#   o.parameterize_transmission('widths',plot=False,sigma=s,tolvec=t)
#   o.results['widths'].scat_allP(figname=wfigname,yran=[.1,.6])
#   o.results['widths'].filter_outliers_widths()
#   o.results['widths'].scat_allP(figname=ffigname,yran=[.1,.6])

# #   o.results['widths'].filter_outliers_widths()
# #   o.results['widths'].scat_allP(figname=ffigname,yran=[.1,.6])
# #   o.compare_transm_methods(figname=base+'_compare.png')
#   pickle.dump(o,open(base+'.p','wb'))
#   os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/QMmodel/SLS_Sept2018/transmQM')
