
 #######################################
# plot_divergent.py
# Erica Lastufka 26/1/18
# Make the widget for getting effective area plots for the grids
#######################################

#######################################
# Usage:

#
######################################

from numpy import arange, sin, pi
from matplotlib.figure import Figure
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pickle
import os
import numpy as np
from scipy.ndimage import rotate
import itertools
global windows
global windowsr

#include offsets
windows=[{'number':11,'pitch':89.6752,  'nominal angle':-44.79357,'offset':[-27,27]},
                {'number':21,'pitch':90.326,  'nominal angle': 45.20793,'offset':[-27,13]},
                {'number':12,'pitch':22.4797,'nominal angle': 44.94825,'offset':[-13,27]},
                {'number':22,'pitch':22.5203,'nominal angle':-45.05184,'offset':[-13,13]},
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

def calc_eff_period1d(p,ang,duty=0.5):
    '''all units um or degrees'''
    #return np.abs((duty*2)*p*np.tan(np.deg2rad(ang)))
    return np.abs(p/np.sin(np.deg2rad(ang)))

def plot_eff_period1d(pf=90., pr=90.,fang=45.,rang=45.,size=1.2,ycen_off=20,phase=0.,duty=0.5, d1=50.,d2=55.,d_det=65.,plot=True):
    #phase = 0 => start with slit?
    pf_eff = calc_eff_period1d(pf,fang)
    pr_eff = calc_eff_period1d(pr,rang)

    ybot_off = ycen_off - ((size*10.)/2.) #ycen_off in mm
    npds = int((size*10000.)/(pf_eff*2.))#number of effective periods
    print npds
    yvec = np.linspace(ybot_off, ybot_off+(size*10.),npds) #1 point at the start of every slit
    first_grid_periods=[y + pf_eff*(d_det/d1) for y in yvec]
    first_grid_offsets=[(d_det/d1)*y for y in yvec]
    second_grid_periods=[y + pr_eff*(d_det/d2) for y in yvec]
    second_grid_offsets=[(d_det/d2)*y for y in yvec]

    #convert into plottable... do I have to exend the mesh? probably....
    #xvec=range(yvec[0],second_grid_offsets[-1],p/10.)
    xvec1=[o + p for o,p in zip(first_grid_offsets, first_grid_periods)]
    xvec2=[o + p for o,p in zip(second_grid_offsets, second_grid_periods)]
    first_grid_plot =[np.mod(n,2) for n in range(0,len(xvec1))]
    second_grid_plot=[np.mod(n,2) for n in range(0,len(xvec2))]

    #plot the patterns overlaid on each other for examination
    if plot:
        fig,ax=plt.subplots()
        ax.step(xvec1, first_grid_plot,'r')
        ax.step(xvec2, second_grid_plot,'g')
        #ax.plot(xvec, first_grid_offsets,'r--')
        #ax.plot(xvec, second_grid_offsets,'g--')
        fig.show()
    return first_grid_periods,second_grid_periods

def get_win_params(num, ptol=False,atol=False):
    ''' Get periods, angles, offsets etc from the window dictionary'''
    fwdict=[dd for dd in windows if dd['number'] == num]
    rwdict=[dd for dd in windowsr if dd['number'] == num]
    pf =fwdict[0]["pitch"]
    pr =rwdict[0]["pitch"]
    fang =fwdict[0]["nominal angle"]
    rang =rwdict[0]["nominal angle"]
    offset=fwdict[0]["offset"]
    return pf,pr,fang,rang,offset #rear offset should be the same

def plot_eff_period_seg(seg,phase=[0.,0.,0.,0.],size=1.2,duty=0.5, d1=50.,d2=55.,d_det=65.,plot=True):
    if seg == 'B':
        wins= [11,21,12,22]
    elif seg == 'C':
        wins= [31,41,32,42]
    else:
        wins= [33,43,34,44]

    paramlist= [get_win_params(w) for w in wins]

    if plot: # plot all 4 windows in segment together
        fig,ax=plt.subplots(2,2)

    for w,p in zip(wins,paramlist): #every set of parameters pf,pr,fang,rang,offset
        #phase = 0 => start with slit?
        pf,pr,fang,rang,offset = p[0],p[1],p[2],p[3],p[4]
        ycen_off = offset[1]
        #print w,ycen_off
        pf_eff = calc_eff_period1d(pf,fang) #um
        pr_eff = calc_eff_period1d(pr,rang) #um
        #print pf_eff,pr_eff
        ybot_off = (ycen_off - ((size*10.)/2.))*1000 #ycen_off in um
        #print ybot_off
        npds = int((size*10000.)/(pf_eff*2.))#number of effective periods
        #print npds
        yvec = np.linspace(ybot_off, ybot_off+(size*10000.),npds) #1 point at the start of every slit
        fg_s=d_det/d1
        rg_s=d_det/d2 #geometric factor
        fg_period=pf_eff*fg_s
        rg_period=pr_eff*rg_s
        first_grid_periods=[y + fg_period for y in yvec] #have to re-adjust to absolute period vs. effective period eventually
        first_grid_offsets=[fg_s*y for y in yvec]
        second_grid_periods=[y + rg_period for y in yvec]
        second_grid_offsets=[rg_s*y for y in yvec]
        #front_eff_period=first_grid_periods[0]-ybot_off
        #rear_eff_period=second_grid_periods[0]-ybot_off
        print "Projected periods, subcoll ", str(w)
        print "front: ",np.abs(fg_period*np.sin(np.deg2rad(fang)))," rear: ",np.abs(rg_period*np.sin(np.deg2rad(rang)))
        foo,bar=calc_moire_period(fg_period,rg_period,fang,rang)
        xvec1=[(o + p)/10000. for o,p in zip(first_grid_offsets, first_grid_periods)] #cm
        xvec2=[(o + p)/10000. for o,p in zip(second_grid_offsets, second_grid_periods)] #cm
        #print xvec1[0],xvec2[0]
        first_grid_plot =[np.mod(n,2) for n in range(0,len(xvec1))]
        second_grid_plot=[np.mod(n,2) for n in range(0,len(xvec2))]

        #plot the patterns overlaid on each other for examination
        if plot:
            axx=int(str(w)[0])-1
            axy=int(str(w)[1])-1
            #ax[axx,axy].step(xvec1, first_grid_plot,'r')
            #ax[axx,axy].step(xvec2, second_grid_plot,'g')
            ax[axx,axy].step(first_grid_plot,xvec1,'r')
            ax[axx,axy].step(second_grid_plot,xvec2,'g')
            ax[axx,axy].set_xlabel('Window '+str(w))
            ax[axx,axy].set_ylabel('Distance from y-center (cm)')
            #ax.plot(xvec, first_grid_offsets,'r--')
            #ax.plot(xvec, second_grid_offsets,'g--')

    if plot:
        fig.show()
    return first_grid_periods,second_grid_periods

def plot_moire_fancy(moire_period,moire_orientation,titles=False,size=1.5,realsize=1.0,reflect=True,xaxoff=False,yaxoff=False):
    nplots = len(moire_period)
    mim,xoff,yoff=[],[],[]
    nsimple=[[0]]
    m,bext,b=[],[],[]
    for p,o in zip(moire_period,moire_orientation):
        #determine how big the array needs to be to fit when rotated and cropped
        asize = np.sqrt(2.)*100.*size#*#(np.abs(np.cos(np.deg2rad(o))))#+np.abs(np.tan(np.deg2rad(o)))) #mm/10
        #print asize
        #off.append(100.*size*np.abs(np.cos(np.deg2rad(o)))#0)#int(size*100*abs(np.tan(np.deg2rad(o)))))#offset
        #draw x-axis vector such that you have 2pi every period
        npts=int(p*10.)
        npds=asize/(p*10.)
        if npds > 25:
            nsimple.append([1])
            #nplots=nplots-1
            print('Asking to plot ' + str(int(npds)) + ' periods! #This will not end well...skipping', p,o)
            mm,bb=plot_moire_simple(p,o,reflect=reflect)
            m.append(mm)
            b.append(bb[:])
            #if o < 90:
            #    ex=[bb[-1]+cc for cc in bb[1:]]
            #else:
            #    ex=[bb[1]-cc for cc in bb[1:]]
            #bb.extend(ex)
            bext=bb
            #bext.append(bb) #extend for better plotting
            continue
        print npts, npds
        xp = list(np.linspace(0,2*np.pi,npts))
        for n in range(0,int(npds)):
            xp.extend(xp)
        xvec=xp[:int(asize)+1]
        #draw the cos pattern
        cpattern=[np.cos(x) for x in xvec]
        #replicate it into an image
        im=np.array([cpattern for y in range(0,int(asize)+1)])
        #print 'original angle',o
        #rotate
        if o <0: #imrotate doesn't like negative angles
            o=-o
        print 'rotation angle', o
        finim=rotate(im,o,reshape=True)
        if reflect:
            mim.append(np.fliplr(finim))
            ang=-o
        else:
            mim.append(finim)
            ang=o
        ##print 'rotated image size: ',np.shape(finim)
        #determine offsets for plotting
        ydiff=100.*size#np.abs(100.*size*np.tan(np.deg2rad(ang))*np.sin(np.deg2rad(ang)))
        xdiff=np.abs(100.*size*np.cos(np.deg2rad(ang))*np.cos(np.deg2rad(ang)))
        yoff.append(50.*size*np.abs(np.tan(np.deg2rad(ang))))#asize*np.abs(np.sin(np.deg2rad(ang))) - ydiff)
        xoff.append(asize*np.abs(np.cos(np.deg2rad(ang))) - xdiff) #asize*np.abs(np.cos(np.deg2rad(ang)))-50.*size)
        nsimple.append([0])
        #print xoff, yoff
    sumsimple=[i for i,j in enumerate(nsimple) if j[0] !=0]

    if nplots == 2:
        fig,ax= plt.subplots(nplots)
        if type(ax) != np.ndarray:
            ax.imshow(mim[0],cmap='gray')
            ax.axhline(yoff[0],xmin=xoff[0],xmax=int(size*100.)+xoff[0],color='b')
            ax.axhline(int(size*100.)+yoff[0],xmin=xoff[0],xmax=int(size*100.)+xoff[0],color='b')
            ax.axvline(xoff[0],ymin=yoff[0],ymax=int(size*100.)+yoff[0],color='b')
            ax.axvline(xoff[0],ymin=yoff[0],ymax=int(size*100.)+yoff[0],color='b')
            #ax.set_xlim([xoff[0],int(size*100.)+xoff[0]])
            #ax.set_ylim([yoff[0],int(size*100.)+yoff[0]])
        else:
            jm=0
            for i,a in enumerate(ax):
                if i+1 not in sumsimple:
                    a.imshow(mim[jm],cmap='gray')
                    a.axhline(yoff[jm],xmin=xoff[jm],xmax=int(size*100.)+xoff[jm],color='b')
                    a.axhline(int(size*100.)+yoff[jm],xmin=xoff[jm],xmax=int(size*100.)+xoff[jm],color='b')
                    a.axvline(xoff[jm],ymin=yoff[jm],ymax=int(size*100.)+yoff[jm],color='b')
                    a.axvline(xoff[jm],ymin=yoff[jm],ymax=int(size*100.)+yoff[jm],color='b')
                    #a.set_xlim([xoff[jm],size*100.+xoff[jm]])
                    #a.set_ylim([yoff[jm],size*100.+yoff[jm]])
                    jm+=1
                else:
                    idx=sumsimple.index(i+1)#int(abs(sumsimple[0]-(i+1)))
                    for k,j in enumerate(bext[idx]):
                        if k%2 ==0: col = 'k'
                        else: col = 'w'
                        a.fill_between(bext[idx],m[idx]*np.array(bext[idx])+j,m[idx]*np.array(bext[idx])+j+1,color=col)
                    #a.set_xlim([b[idx][0],b[idx][-1]])
                    #a.set_ylim([b[idx][0],b[idx][-1]])
                a.set_aspect('equal')
    else:

        fig,ax= plt.subplots(2,nplots/2,figsize=[2*nplots,8])
        jm=0
        pidx=0
        for i,r in enumerate(itertools.product(range(0,nplots/2),(range(0,2)))):
            if i+1 not in sumsimple:
                ax[r[1]][r[0]].imshow(mim[jm],cmap='gray',origin='lower left')
                #rect=patches.Rectangle((xoff[jm],yoff[jm]),int(size*100.),int(size*100.),linewidth=5,edgecolor='b',facecolor='none')
                #ax[r[1]][r[0]].add_patch(rect)
                #print(xoff[jm],int(size*100.)+xoff[jm],yoff[jm],int(size*100.)+yoff[jm])
                if not xaxoff and not yaxoff:
                    ax[r[1]][r[0]].set_ylim([yoff[jm],int(realsize*100.)+yoff[jm]])
                    ax[r[1]][r[0]].set_xlim([xoff[jm],int(realsize*100.)+xoff[jm]])
                else:
                    ax[r[1]][r[0]].set_ylim([yaxoff[pidx],int(realsize*100.)+yaxoff[pidx]])
                    ax[r[1]][r[0]].set_xlim([xaxoff[pidx],int(realsize*100.)+xaxoff[pidx]])

                jm+=1

            else:
                idx=sumsimple.index(i+1)
                print('len(bext):', len(bext))
                for k,j in enumerate(bext[idx:-1]):
                    #print len(bext[i-sumsimple]),len(bext[i-sumsimple]+j)
                    #print(idx, len(j))
                    #print m[i-sumsimple]
                    if k%2 ==0: col = 'k'
                    else: col = 'w'
                    ax[r[1]][r[0]].fill_between(bext,m[idx]*np.array(bext)+j,m[idx]*np.array(bext)+bext[k+1],color=col)
                    #rect=patches.Rectangle((b[idx][0],b[idx][0]),b[idx][-1],b[idx][-1],linewidth=5,edgecolor='b',facecolor='none')
                    #ax[r[1]][r[0]].add_patch(rect)
                if not xaxoff and not yaxoff:
                    ax[r[1]][r[0]].set_xlim([b[idx][0],b[idx][0]+int(realsize*100.)])
                    ax[r[1]][r[0]].set_ylim([b[idx][0],b[idx][0]+int(realsize*100.)])
                else:
                    ax[r[1]][r[0]].set_xlim([xaxoff[pidx],xaxoff[pidx]+int(realsize*100.)])
                    ax[r[1]][r[0]].set_ylim([yaxoff[pidx],yaxoff[pidx]+int(realsize*100.)])
                    print(r[1],r[0],xaxoff[pidx],yaxoff[pidx])

            ax[r[1]][r[0]].set_aspect('equal')

            if titles:
                ax[r[1]][r[0]].set_title(titles[i],fontsize='small')
            ax[r[1]][r[0]].axis('off')
            pidx+=1
            print('pidx', pidx)
    plt.tight_layout()
    fig.show()

def plot_moire_simple(moire_period,moire_orientation,size=1.5,reflect=True):
    '''just draw lines because there's too many periods. Only return plottables, plot them together with other windows in plot_moire_fancy'''
    mim,xoff,yoff=[],[],[]
    nlines=int((10.*size/(moire_period/2.)))+1
    b=np.linspace(0,size*200.,2*nlines)
    #npix =
    xp = list(b)#np.linspace(0,nlines))
    if reflect:
        slope= np.tan(np.deg2rad(-1*moire_orientation))#calculate the slope
    else:
        slope= np.tan(np.deg2rad(moire_orientation))#calculate the slope
    #print len(xp)
    return slope,xp #list(b) #x,slope and offsets


def calc_moire_period(fp,rp,fa,ra,eff=False,quiet=True):
    if eff:#first adjust the effective period back to the actual period
        fp=np.abs(fp*np.sin(np.deg2rad(fa)))
        rp=np.abs(rp*np.sin(np.deg2rad(ra)))
    #front and rear period in um, front and rear orientation in degrees
    fx=1000*np.sin(np.deg2rad(fa))/fp #mm
    fy=-1000*np.cos(np.deg2rad(fa))/fp #mm
    rx=1000*np.sin(np.deg2rad(ra))/rp #mm
    ry=-1000*np.cos(np.deg2rad(ra))/rp #mm
    kx=fx+rx
    ky=fy+ry
    mx=fx-rx
    my=fy-ry
    moire_period=1./np.sqrt(mx**2+my**2)
    moire_orientation=np.rad2deg(np.arctan2(my,mx))
    if not quiet:
        print "Moire Period: ", np.mean(moire_period) #this is in mm...
        print "Moire Orientation: ",np.mean(moire_orientation)
    return moire_period,moire_orientation

def calc_moire_period_mpi(inp):
    d1,d2,d_det=inp[4],inp[5],inp[6]
    fpo,fao,rpo,rao=inp[0],inp[1],inp[2],inp[3]
    twist=inp[7]
    s1=1.+(d2/d1) #d2-d1 replaced by 1.3 (7mm red part, 3mm each lower cover) since this is fixed
    #d2-d1 is replaced by d2 where d2 is acting as delta d, ranging from parameters given
    s2=1. #this is never used? i think it might need to be in this case...
    fp=fpo*s1
    fa=fao+twist
    rp=rpo
    ra=rao
    #fa=fa+twist
    #print fp,fa,rp,ra
    #print calc_moire_period(fp,rp,fa,ra,eff=False,quiet=True)

    #front and rear period in um, front and rear orientation in degrees
    fx=1000*np.sin(np.deg2rad(fa))/fp #mm
    fy=-1000*np.cos(np.deg2rad(fa))/fp #mm
    rx=1000*np.sin(np.deg2rad(ra))/rp #mm
    ry=-1000*np.cos(np.deg2rad(ra))/rp #mm
    kx=fx+rx
    ky=fy+ry
    mx=fx-rx
    my=fy-ry
    moire_period=1./np.sqrt(mx**2+my**2)
    moire_orientation=np.rad2deg(np.arctan2(my,mx))

    #make inp into string key
    inpstr=','.join([str(d1),str(d2),str(d_det),str(twist)])
    return moire_period,moire_orientation,inpstr

def calc_all_moire(seg=True,wins=False,nominal=True, plot=False):
    if seg:
        if seg == "B":
            wins = [11,21,12,22]
        elif seg == "C":
            wins = [31,41,32,42]
        else:
            wins = [33,43,34,44]
    for w in wins:
        fwdict=[dd for dd in windows if dd['number'] == w]
        rwdict=[dd for dd in windowsr if dd['number'] == w]
        pf =fwdict[0]["pitch"]
        pr =rwdict[0]["pitch"]
        fang =fwdict[0]["nominal angle"]
        rang =rwdict[0]["nominal angle"]
        fg_periods,rg_periods=plot_eff_period1d(pf=pf,pr=pr,fang=fang,rang=rang,plot=plot)
        print "Mean projected periods, subcoll ", str(w)
        print "front: ", np.mean(fg_periods)," rear: ",np.mean(rg_periods)
        mp,mf=calc_moire_period(np.mean(fg_periods),np.mean(rg_periods),fang,rang)

def test_visibility(seg,size=1.2):
    '''Test if in a 2d overlay of nominal periods and angles, a moire pattern is visible by just plotting'''
    if seg == "B":
        wins = [11,21,12,22]
    elif seg == "C":
        wins = [31,41,32,42]
    else:
        wins = [33,43,34,44]
    paramlist= [get_win_params(w) for w in wins]

    fig,ax=plt.subplots(2,2)

    for w,p in zip(wins,paramlist): #every set of parameters pf,pr,fang,rang,offset
        #phase = 0 => start with slit?
        pf,pr,fang,rang,offset = p[0],p[1],p[2],p[3],p[4]
        xcen_off, ycen_off =offset[0], offset[1]
        #print w,ycen_off
        pf_eff = calc_eff_period1d(pf,fang) #um
        pr_eff = calc_eff_period1d(pr,rang) #um
        #print pf_eff,pr_eff
        ybot_off = (ycen_off - ((size*10.)/2.)) #ycen_off in um
        xleft_off = (xcen_off - ((size*10.)/2.)) #xcen_off in um
        #print ybot_off
        npds = int((size*10000.)/(pf_eff*2.))#number of effective periods
        #print npds

        yvec = np.linspace(ybot_off, ybot_off+(size),npds) #1 point at the start of every slit
        xvec = np.linspace(xleft_off, xleft_off+(size),npds)
        #first_grid_periods=[y + pf_eff for y in yvec] #have to re-adjust to absolute period vs. effective period eventually
        #first_grid_offsets=[((d_det/d1))*y for y in yvec]
        #second_grid_periods=[y + pr_eff for y in yvec]
        #second_grid_offsets=[((d_det/d2))*y for y in yvec]
        #front_eff_period=first_grid_periods[0]-ybot_off
        #rear_eff_period=second_grid_periods[0]-ybot_off

        #xvec1=[(o + p)/10000. for o,p in zip(first_grid_offsets, first_grid_periods)] #cm
        #xvec2=[(o + p)/10000. for o,p in zip(second_grid_offsets, second_grid_periods)] #cm
        #print xvec1[0],xvec2[0]
        #first_grid_plot =[np.mod(n,2) for n in range(0,len(xvec1))]
        #second_grid_plot=[np.mod(n,2) for n in range(0,len(xvec2))]
        mf = np.tan(np.deg2rad(fang))#slope of line in front grid
        mr = np.tan(np.deg2rad(rang))
        print mf,mr
        #plot the patterns overlaid on each other for examination
        for n in range(0,npds): #plot the lines
            axx=int(str(w)[0])-1 #only works for seg B actually....
            axy=int(str(w)[1])-1
            #ax[axx,axy].step(xvec1, first_grid_plot,'r')
            #ax[axx,axy].step(xvec2, second_grid_plot,'g')
            bf=n*(pf_eff/10000.)
            br=n*(pr_eff/10000.)
            ax[axx,axy].plot(xvec,mf*xvec+bf,'r')
            ax[axx,axy].plot(xvec,mr*xvec+br,'g')

        ax[axx,axy].set_xlabel('Window '+str(w))
        ax[axx,axy].set_xlim([xleft_off,xleft_off+size])
        ax[axx,axy].set_ylim([xleft_off,xleft_off+size])
        ax[axx,axy].set_ylabel('Distance from y-center (cm)')

    fig.show()
    return xvec

def tex2params(infile,test=1):
    '''read the latex output file and get the periods and orientations for the different tests so that they can be easily plotted'''
    periods,orientations=[],[]
    with open(infile,'r') as f:
        lines=f.readlines()
    for line in lines:
        pp=line[:-3].split('&')
        idx=(test*2)-1
        periods.append(float(pp[idx]))
        ang=float(pp[idx+1])
        if ang<1.0:
            ang=90.+abs(float(pp[idx+1]))
        orientations.append(ang)
    return periods,orientations

