"""
===================
mpi_im2bin.py
Erica  Lastufka 20.9.17
===================

Use multiprocessing to turn complete set of MiSolFA window images into binary representations ie Canny edges. Option to continue to 2nd level, performing probabilistic Hough transform on Canny edges.

"""
import multiprocessing as mpi
import glob
import os
import time
import numpy as np
import edge_and_hough as eh

#start the clock
start=time.time()

#where's the data at
#fdir = '/Users/wheatley/Documents/Solar/MiSolFA/prototypes/'
#os.chdir(fdir)
#folders=[fdir+x for x in glob.glob('*/') if x.startswith('mw')]

#finest windows - let's use 2sigma for those! build this into eh.Canny_edge(), easier to control there

#create the task list
tasks=[]
#for f in folders:
    #files=glob.glob(f+'/win*bright.tif')
files=eh.EM_list('44','5.0',ending='.tif')#glob.glob('win44*5.0X.tif')
print files[0]
#files=glob.glob('win*5.0X.tif')
print 'files: ',len(files)
    #for fi in files:
        #if 'win11' in fi:
#tasks.append(files)
    #print len(tasks), len(files)
flat_list=files
#flat_list = [item for sublist in tasks for item in sublist] #to start with let's only test on 20
#print flat_list

#create the pool, distribute & process
#with pool as mpi.Pool: #this way I don't have to worry about cleanup
pool=mpi.Pool()
#distribute the work
results,results_hough=[],[]
print 'flat_list: ', len(flat_list)
print 'results: ', len(results)
for t in flat_list:
    results.append( pool.apply_async(eh.Canny_edge, args=(t,)) )

print 'Canny edge detection took %.2f seconds' % (time.time() - start)
print 'results: ', len(results)
for r1 in results:
    tt=r1.get()
    results_hough.append( pool.apply_async(eh.prob_hough, args=(tt,)) )

    #for r1,r2 in zip(results,results_hough):
    #    rr=r1.get()
    #    tt=r2.get()
    #    print 'new files ' + rr[rr.rfind('/'):] +' and '+ tt[tt.rfind('/'):]+' created'
    
filep=eh.EM_list('44','5.0',ending='_edges.p')#glob.glob('win44*_edges.p')
fileh=eh.EM_list('44','5.0',ending='_hough.p')#glob.glob('win44*_hough.p')

win=44
pitch='20'
titlec='Slit width distribution for window '+str(win) + ', nominal pitch '+pitch + ' $\mu$m at 5X magnification'
figname='win'+str(win) + 'width_hist.png'

w,b=eh.get_slit_width(filep,mag=5.0,window_num=str(win),title=titlec,xran=[float(pitch)-float(pitch),float(pitch)+float(pitch)],figname=figname)
        
figname='win'+str(win) + 'ang_hist.png'
print figname
eh.hough_hist(fileh,win,figname=figname,side=1.0)

#stop the clock, print time:
print 'Processing took %.2f seconds' % (time.time() - start)

##just to compare...

#start=time.time()
#for f in flat_list:
#    slow=eh.Canny_edge(f)
#    print 'new file ' + slow[slow.rfind('/'):] +' created'
#print 'Processing took %.2f seconds' % (time.time() - start)

#2X slower for 10 images
