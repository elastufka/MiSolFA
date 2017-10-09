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
fdir = '/Users/wheatley/Documents/Solar/MiSolFA/prototypes/'
os.chdir(fdir)
folders=[fdir+x for x in glob.glob('*/') if x.startswith('mw')]

#finest windows - let's use 2sigma for those! build this into eh.Canny_edge(), easier to control there

#create the task list
tasks=[]
for f in folders:
    files=glob.glob(f+'/win*bright.tif')
    for fi in files:
        if 'win11' in fi:
            tasks.append(files)
    #print len(tasks), len(files)

flat_list = [item for sublist in tasks for item in sublist][0:120] #to start with let's only test on 20
#print flat_list

#create the pool, distribute & process
#with pool as mpi.Pool: #this way I don't have to worry about cleanup
pool=mpi.Pool()
#distribute the work
results,results_hough=[],[]
for t in flat_list:
    results.append( pool.apply_async(eh.Canny_edge, args=(t,)) )

for r1 in results:
    tt=r1.get()
    results_hough.append( pool.apply_async(eh.prob_hough, args=(tt,)) )

    #for r1,r2 in zip(results,results_hough):
    #    rr=r1.get()
    #    tt=r2.get()
    #    print 'new files ' + rr[rr.rfind('/'):] +' and '+ tt[tt.rfind('/'):]+' created'
    
#stop the clock, print time:
print 'Processing took %.2f seconds' % (time.time() - start)

##just to compare...

#start=time.time()
#for f in flat_list:
#    slow=eh.Canny_edge(f)
#    print 'new file ' + slow[slow.rfind('/'):] +' created'
#print 'Processing took %.2f seconds' % (time.time() - start)

#2X slower for 10 images
