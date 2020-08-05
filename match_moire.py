
 #######################################
#match_moire.py
# Erica Lastufka 08/2/18
# Find the best match of the observed and predicted moire patterns within a certain given parameter range
#######################################

#######################################
# Usage:

#
######################################

from numpy import arange, sin, pi
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pickle
import os
import numpy as np
from scipy.ndimage import rotate
import itertools
global windows
global windowsr
import plot_divergent as pd
import multiprocessing as mpi
import time
from scipy import optimize

windows=[{'number':11,'pitch':89.6752,  'nominal angle':-44.79357,'offset':[-27,27]},
                {'number':21,'pitch':90.326,  'nominal angle': 45.20793,'offset':[-27,13]},
                {'number':12,'pitch':22.4797,'nominal angle': 44.94825, 'offset':[-13,27]},
                {'number':22,'pitch':22.5203,'nominal angle':-45.05184 ,'offset':[-13,13]},
                {'number':31,'pitch':45.0814, 'nominal angle':-45.10378,'offset':[-27,-13]},
                {'number':41,'pitch':44.9187, 'nominal angle': 44.8966,'offset':[-13,-13]},
                {'number':32,'pitch':18.013, 'nominal angle': 45.04146,'offset':[-27,-27]},
                {'number':42,'pitch':17.987, 'nominal angle':-44.95859,'offset':[-13,-27]},
                {'number':33,'pitch':29.9639,  'nominal angle':-44.93102,'offset':[13,-13]},
                {'number':43,'pitch':30.0362,  'nominal angle': 45.06914,'offset':[27,-13]},
                {'number':34,'pitch':14.991, 'nominal angle': 44.96549,'offset':[13,-27]},
                {'number':44,'pitch':15.009, 'nominal angle':-45.03455,'offset':[27,-27]}]
#include offsets
windowsr=[{'number':11,'pitch':90.326,  'nominal angle':-45.20793,'offset':[27,27]},
                {'number':21,'pitch':89.6752,  'nominal angle': 44.79357,'offset':[27,13]},
                {'number':12,'pitch':22.5203,'nominal angle': 45.05184,'offset':[13,27]},
                {'number':22,'pitch':22.4797,'nominal angle':-44.94825,'offset':[13,13]},
                {'number':31,'pitch':44.9187, 'nominal angle':-44.8966,'offset':[27,-13]},
                {'number':41,'pitch':45.0814, 'nominal angle': 45.10378,'offset':[13,-13]},
                {'number':32,'pitch':17.987, 'nominal angle': 44.95859,'offset':[27,-27]},
                {'number':42,'pitch':18.013, 'nominal angle':-45.04146,'offset':[13,-27]},
                {'number':33,'pitch':30.0362,  'nominal angle':-45.06914,'offset':[-13,-13]},
                {'number':43,'pitch':29.9639,  'nominal angle': 44.93102,'offset':[-27,-13]},
                {'number':34,'pitch':15.009, 'nominal angle': 45.03455,'offset':[-13,-27]},
                {'number':44,'pitch':14.991, 'nominal angle':-44.96549,'offset':[-27,-27]}]

def measure_moire(tiff_in,bounds,thetaran, threshold=.8,n=100,pix2um=24.,peaktol=5):
    all_theta=np.linspace(thetaran[0],thetaran[1],n)
    imarr=im2ndarray(tiff_in)[bounds[0]:bounds[1],bounds[2]:bounds[3]]
    marr=np.ma.masked_less(imarr,np.max(imarr)*threshold)

    fig,ax=plt.subplots()
    ax.imshow(marr)
    fig.show()

    def rot_matrix(marr,angle):
        rotarr=rotate(marr,angle)
        return np.max(np.gradient(np.sum(rotarr,axis=0)))

    bestangle=0.0
    maxgrad=0.0

    for angle in all_theta:
        thismax=rot_matrix(marr,angle)
        if thismax >=maxgrad:
            maxgrad=thismax
            bestangle=angle

    print bestangle, maxgrad
    bestarr=np.ma.masked_less(rotate(marr,bestangle),np.max(imarr)*threshold)
    bsumarr=np.sum(bestarr,axis=0)
    bgradarr=np.gradient(bsumarr) #normalize these
    bgradarr=bgradarr/np.max(bgradarr)
    bsumarr=bsumarr/np.max(bsumarr)


    #now measure period...
    from scipy.signal import argrelextrema
    lmax=list(argrelextrema(bsumarr, np.greater)[0])
    #lsmooth=[]
    #for i in range(0,len(lmax)-1):
    #    if lmax[i+1]-lmax[i] <= peaktol:
    #        lsmooth.append(np.mean([lmax[i+1],lmax[i]]))
    #    else:
    #        lsmooth.append(lmax[i])
    #print lsmooth
    ldiff=[lmax[i+1]-lmax[i] for i in range(0,len(lmax)-1)]
    #apply peaktol limits
    group_hist,group_xvec=[],[]
    #first get the groups of True values -- this is different now, already have peaks
    group=False
    for i,es in enumerate(lmax[:-1]):
        #if irange:
            #print i+irange[0],es,group,group_hist,group_xvec
        #print group, es, lmax[i+1]
        if not group: #start the group for the first non-zero value
            lhist=[bsumarr[es]]
            lxvec=[es]
            group=True
        if group and lmax[i+1]-es < peaktol: #if aready in group, stay in group unless following tol pixels are all zeros
            lxvec.append(es)
            lhist.append(bsumarr[es])
        elif lmax[i+1]-es >= peaktol:
            group=False
        if not group:
            try:
                if lxvec not in group_xvec:
                    group_xvec.append(lxvec)
                    group_hist.append(lhist)
            except NameError:
                continue

    xpeak=[]
    #print group_xvec,group_hist
    for xg,yg in zip(group_xvec,group_hist):
        xpeak.append(np.average(xg, weights=yg))

    xpeakdiff=[xpeak[i+1]-xpeak[i] for i in range(0,len(xpeak)-1)]
    print np.mean(xpeakdiff)*pix2um
    #bgradarr=bgradarr/np.max(bgradarr)
    #bsumarr=bsumarr/np.max(bsumarr)


    fig,ax=plt.subplots()
    ax.imshow(bestarr,origin='lower left', cmap='gray',alpha=0.6)
    ax2=ax.twinx()
    ax2.plot(range(0,len(bsumarr)),bsumarr, linewidth=2, label='sum')
    ax2.plot(range(0,len(bsumarr)),bgradarr, linewidth=2, label='gradient of sum')
    ax2.set_ylim([-2,2])
    ax2.legend()
    fig.show()
    return xpeak,xpeakdiff



def calc_diff(mp,mo,observed):
    diffp=mp-observed[0]
    diffo=mo-observed[1]
    return diffp,diffo

def match_moire(observed, start_params, param_ranges, param_steps, plot=True,reflect=True):
    '''Fit for a single given moire period'''
    #start_params:[front_period,front_orientation,rear_period,rear_orientation, d1,d2,d_det, relative twist]
    #param_ranges:[(),(),(),(),(),(),()] tuples of ranges for each parameter
    #param_steps: [] list of step numbers for each parameter
    #observed: [moire_period,moire_orientation]
    ### build parameter list for pool
    all_params=[]
    for sp,pr,ps in zip(start_params,param_ranges,param_steps):
        all_params.append(list(np.linspace(sp+pr[0],sp+pr[1],ps)))
    #print np.shape(all_params),all_params[:5]
    flat_list=list(itertools.product(*all_params))
    #print 'number of possible parameter combinations: ',len(flat_list)
    #distribute to workers, start the clock
    start=time.time()
    pool=mpi.Pool()
    results,outlist=[],[]
    best_params={}
    for t in flat_list:
        #with mpi.Pool() as pool:
        results.append(pool.apply_async(pd.calc_moire_period_mpi, args=(t,)).get())
    pool.close()
    #print len(results)
    for r1 in results:
        mp,mo,inp=r1#.sort()#pool.map(r1,range(2))
        dp,do=calc_diff(mp,mo,observed)
        ll=[dp,do,mp,mo,inp]#[:-3]] #if we are only concered about the distance between things... modify this later if not
        #ll.extend(f for f in fl)
        outlist.append(ll)
    #print 'Moire period calculation took %.2f seconds' % (time.time() - start)

    #sort
    outarr=np.array(outlist)
    outarr=outarr[outarr[:,-1].argsort()] #sort by last column with the string key
    #print np.shape(outarr)
    #np.sort(outarr,0)
    print 'best params to fit period:', outarr[0]
    #best_params["period"]=outarr[0]
    if plot:
        pd.plot_moire_fancy([outarr[0][2]],[outarr[0][3]],reflect=reflect)
    #np.sort(outarr,1)
    print 'best params to fit orientation:', outarr[0]
    if plot:
        pd.plot_moire_fancy([outarr[0][2]],[outarr[0][3]],reflect=reflect)
    #np.sort(outarr,2)
    #best_params["orientation"]=outarr[0]
    return outarr.tolist()#, best_params

def get_score(inp):
    score=abs(float(inp[0]))*abs(float(inp[1]))
    return score

def moire_hist(win=11, exp_params=0.0, sigmas=0.1, plot=True,n=500):
    '''calculate moire patterns for given distributions of parameters, plot hist of period and orientation'''
    #grid params: period, orientation
    #exp params: in this case, the twist
    #sigmas: must match number of exp_params [later]
    front_params=[[w['pitch'],w['nominal angle']] for w in windows if w['number'] == win][0]
    rear_params=[[w['pitch'],w['nominal angle']] for w in windowsr if w['number'] == win][0]
    print front_params, rear_params
    #generate normal distribution with given sigma
    gparams,periods,angles=[],[],[]
    for i in range(0,n):
        tp=np.random.normal(loc=exp_params,scale=sigmas,size=None)
        #print tp
        gparams.append([front_params[0],front_params[1],rear_params[0],rear_params[1],1,0,0,tp])

    #calculate moire patterns for each ... instance
    for g in gparams:
        pp,oo,_=calc_moire_period_mpi(g)
        periods.append(pp)
        angles.append(oo)

    print np.mean(periods), np.std(periods)
    print np.mean(angles),np.std(angles)
    #plot
    fig,ax=plt.subplots(nrows=1,ncols=2)
    ax[0].hist(periods, normed=True)
    ax[1].hist(angles, normed=True)
    fig.show()

def match_segment(infile, plot=True,reflect=True):
    '''Match moire patterns in the entire segment. For now, optimize for change in distances only'''
    obs,start_params,pr,ps=[],[],[],[]
    with open(infile) as f: #read from csv
        lines=f.readlines()
        #f.close()
    for l in lines:
        if l.startswith("observed"):
            ll=l[:-1].split(' ')[2:]
            obs.append([float(aa) for aa in ll])
        elif l.startswith('start_params'):
            ll=l[:-1].split(' ')[2:]
            start_params.append([float(aa) for aa in ll])
        elif l.startswith('param_ranges'):
            ll=l[:-1].split(' ')[2:]
            #have to make into tuples
            tpls=[]
            for aa in ll:
                tpls.append((float(aa[1:aa.find(',')]),float(aa[aa.find(',')+1:-1])))
            pr.append(tpls)
        elif l.startswith('param_steps'):
            ll=l[:-1].split(' ')[2:]
            ps.append([int(aa) for aa in ll])
    allout=[]
    #print obs
    #print start_params
    #print pr
    #print ps
    for o,sp,r,s in zip(obs,start_params,pr,ps): #once for each window basically
        allout.append(match_moire(obs[0],start_params[0],pr[0],ps[0],plot=False)) #sorted with set of parameters an key, so should all be in the same order

    #allout=allout[0]
    #now score. multiply abs(period_diff) by abs(orientation_diff) for all windows. smallest values should give best match
    keys=[al[-1] for al in allout[0]] #should be the same keys for each window
    scores,mpvals,nmpvals=list(np.ones(len(allout[0]))),[],[]
    for win in allout:
        wscores=[get_score(w[:2]) for w in win] #scores of each set of parameters for that window
        scores=[scores[i]*wscores[i] for i in range(0,len(wscores))] #score over all windows
        #print len(scores)
    #mp=[w[2:4] for w in win for win in allout]
    #mpvals.append(mp)

    #return scores,mpvals
    #re-arrange mpvals
    #for i in range(0,len(allout)):
    #nmpvals.append([mpvals[i] for i in range(0,len(allout))])

    scores, keys=zip(*sorted(zip(scores, keys)))

    #find the lowest score
    print scores[0],keys[0]
    print scores[-1],keys[-1]

    keylist=keys[0].split(',')
    keylist=[float(k) for k in keylist]

    for i in range(0,len(obs)):
        #print start_params[i][:-4]+keylist
        mpvals.append(pd.calc_moire_period_mpi(start_params[i][:-4]+keylist)[:2])
    #print mpvals
    print 'Configuration which gives closest match to observed patterns: ', keylist,mpvals #mpvals doesn't take twist into account yet
    mps=[m[0] for m in mpvals]
    mos=[m[1] for m in mpvals]

    if plot:
        #mps=[m[0] for m in mpvals]
        #mos=[m[1] for m in mpvals]
        pd.plot_moire_fancy(mps,mos)


    return keylist,mps,mos

def seg_function(p,fparr,yarr):
    total_err=[]
    for fp,y in zip(fparr,yarr):
        Mx=1000*p[0]*np.sin(np.deg2rad(fp[1]+0.5*p[3]))/((p[0]+p[1])*fp[0])-1000*np.sin(np.deg2rad(fp[3]-0.5*p[3]))/fp[2]
        My=-1000*p[0]*np.cos(np.deg2rad(fp[1]+0.5*p[3]))/((p[1]+p[0])*fp[0])+1000*np.cos(np.deg2rad(fp[3]-0.5*p[3]))/fp[2]
        func_period = 1./np.sqrt(Mx**2+My**2)
        errfunc_period=func_period- y[0]
        #print(func_period,y[0])
        func_angle=np.rad2deg(np.arctan2(My,Mx))
        errfunc_angle=func_angle-y[1]
        #print(func_angle,y[1])
        total_err.append(np.abs(errfunc_period)+np.abs(errfunc_angle))
        #print errfunc_period,errfunc_angle,total_err
    return np.mean(total_err)


def match_segment_chisq(infile, plot=True,reflect=True,all_free=False):
    '''Match moire patterns in the entire segment. For now, optimize for change in distances only'''
    obs,start_params,pr,ps=[],[],[],[]
    with open(infile) as f: #read from csv
        lines=f.readlines()
        #f.close()
    for l in lines:
        if l.startswith("observed"):
            ll=l[:-1].split(' ')[2:]
            obs.append([float(aa) for aa in ll])
        elif l.startswith('start_params'):
            ll=l[:-1].split(' ')[2:]
            start_params.append([float(aa) for aa in ll])
        elif l.startswith('param_ranges'):
            ll=l[:-1].split(' ')[2:]
            #have to make into tuples
            #tpls=[]
            upper,lower=[],[]
            for aa in ll:
                #tpls.append((float(aa[1:aa.find(',')]),float(aa[aa.find(',')+1:-1])))
                lower.append(float(aa[1:aa.find(',')]))
                upper.append(float(aa[aa.find(',')+1:-1]))
            pr=[lower,upper]#.append(tpls)
        elif l.startswith('param_steps'):
            ll=l[:-1].split(' ')[2:]
            ps.append([int(aa) for aa in ll])
    allout=[]

    #fix ranges to include start params
    bounds=[]
    for ss in start_params:
        bounds.append([np.array(pr[0])+np.array(ss),np.array(pr[1])+np.array(ss)])

    ###initial fit params####
    #results=[]
    #for sp,nom,bb in zip(start_params,obs,bounds):
    p0=[sp[4:] for sp in start_params]# Initial guess for the parameters
    fp=[sp[:4] for sp in start_params]#fixed parameters
    #nom=[o[0] for o in obs]
    bb1=bounds[0][0][4:]
    bb2=bounds[0][1][4:]
    bb=[bb1,bb2]
    #    #print p0
    #    #print fp
    #    print nom
    #    nom_period=nom[0]
    #    nom_ang=nom[1]
    p1 = optimize.least_squares(seg_function, p0[0], bounds=bb, args=(fp,obs),xtol=1e-10)
    print p1
    lout=list(p1.x)
    mps,mos=[],[]
    for i in range(0,len(obs)):
        inp=fp[i]+lout
        print inp
        foo= pd.calc_moire_period_mpi(inp)
        mps.append(foo[0])
        mos.append(foo[1])

    print mps,mos
    if plot:
        pd.plot_moire_fancy(mps,mos)

    return p1

