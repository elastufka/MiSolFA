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

# wins=[11,21,12,22,31,41,32,42,33,43,34,44]
# Aflags=[]
# #April
# for w in wins:
#   o=Grating(w,Apr=True,May=False)#pickle.load(open('win'+str(w)+'.p','rb'))#Grating(w,May=True)
#   base='win'+str(w)
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


wins=[42,33,43,34,44]
# 21: P0: 0 (-1.75)? -> p0 is out of sequence!
# P1: -1.75 ghost. Out of sequence! 1.75,3.5,5.25 ,7 bad. Flag p1!
# P2: -8.75 bad. Maybe the whole window is out of sequence?
# p5 :8.75 (7 is probably ok)
# 11:P3,p4 : -8.5 partially blocked
# P5: 5.25,7,8.5 partially blocked (last 2 the worst)
# Mflags=[[['p',0],['p',1],['p0_ang','p2_-8.75']],[['p0ang','p3_-8.5'],['p0ang','p4_-8.5'],['p0ang','p5_-7'],['p0ang','p5_-8.5']],False,False,False,False,False,False,False,False,False]
# #May
# for w,flags in zip(wins,Mflags):
#   o=Grating(w,May=True, EM=True)#pickle.load(open('win'+str(w)+'.p','rb'))#Grating(w,May=True)
#   base='sub2737_win'+str(w)
#   bfigname=base+'_binary.png'
#   sfigname=base+'_sum.png'
#   wfigname=base+'_widths.png'
#   ffigname=base+'_widths_filtered.png'
#   if flags:
#      for fl in flags:
#         if fl[0] =='p':
#             o.flag_data(p=fl[1])
#         else:
#             o.flag_data(p0ang=fl[1])
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


for w in wins:
  o=Grating(w,May=False,Sept=True, EM=False)#pickle.load(open('win'+str(w)+'.p','rb'))#Grating(w,May=True)
  base='sub3501_win'+str(w)
  bfigname=base+'_binary.png'
  sfigname=base+'_sum.png'
  wfigname=base+'_widths.png'
  ffigname=base+'_widths_filtered.png'
  o.parameterize_transmission('binary',plot=False)
  o.results['binary'].scat_allP(figname=bfigname)
  o.parameterize_transmission('sum',plot=False)
  o.results['sum'].scat_allP(figname=sfigname)
  o.parameterize_transmission('widths',plot=False)
  o.results['widths'].scat_allP(figname=wfigname,yran=[.1,.6])
  o.results['widths'].filter_outliers_widths()
  o.results['widths'].scat_allP(figname=ffigname,yran=[.1,.6])
  o.compare_transm_methods(figname=base+'_compare.png')
  pickle.dump(o,open(base+'.p','wb'))
