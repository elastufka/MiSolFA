"""
===================
all_possible_moire.py
Erica  Lastufka 9.11.17
===================

Plot distribution of all possible moire periods and orientations within 1 sigma of the mean

"""
import glob
import pickle
import os


def rv(value,dec=3):
    return np.round([value],decimals=dec)[0]

def all_possible_moire(windict,windownum,dec=3,n=50,percent=10):
    '''Takes dictionary of statistics and window number to calculate and plot histograms of period and orientation distributions'''
    fdict=[ww for ww in windict[:11] if ww['number'] == windownum][0]
    rdict=[ww for ww in windict[11:] if ww['number'] == windownum][0]
    FmeanP=rv(fdict['pmean'],dec)
    FdevP=rv(fdict['pvar'],dec)#np.linspace(-1*rv(fdict['pvar'],dec),rv(fdict['pvar'],dec),n) ##vectorize
    RmeanP=rv(rdict['pmean'],dec)
    RdevP=rv(rdict['pvar'],dec)#np.linspace(-1*rv(rdict['pvar'],dec),rv(rdict['pvar'],dec),n)##vectorize
    FmeanA=rv(fdict['amean'],dec)*np.pi/180.
    FdevA=rv(fdict['avar'],dec)*np.pi/180.#np.linspace(-1*rv(fdict['avar'],dec)*np.pi/180,rv(fdict['avar'],dec)*np.pi/180,n)
    RmeanA=-1.*rv(rdict['amean'],dec)*np.pi/180.
    RdevA=rv(rdict['avar'],dec)*np.pi/180.#np.linspace(-1*rv(rdict['avar'],dec)*np.pi/180,rv(rdict['avar'],dec)*np.pi/180,n)
    #FgaussP=np.random.normal(FmeanP,FdevP,n)
    #FgaussA=np.random.normal(FmeanA,FdevA,n)
    #RgaussP=np.random.normal(RmeanP,RdevP,n)
    #RgaussA=np.random.normal(RmeanA,RdevA,n)

    FgaussP=np.random.normal(59.855,.06,n)
    FgaussA=np.random.normal(44.86*np.pi/180.,.02*np.pi/180.,n)
    RgaussP=np.random.normal(60.145,.06,n)
    RgaussA=np.random.normal(45.14*np.pi/180.,.02*np.pi/180.,n)
    
    
    Fx_vec,Fy_vec=calc_Svec(FgaussP,FgaussA)
    Rx_vec,Ry_vec=calc_Svec(RgaussP,RgaussA)

    fig,(ax1,ax2)=plt.subplots(2)
    ax1.hist(FgaussP,n/2)
    ax2.hist(RgaussP,n/2)
    fig.show()
    
    fig,(ax1,ax2)=plt.subplots(2)
    ax1.hist(FgaussA*180./np.pi,n/2)
    ax2.hist(RgaussA*180./np.pi,n/2)
    fig.show()

    print np.mean(FgaussP),np.mean(RgaussP),np.std(FgaussP),np.std(RgaussP)
    print np.mean(FgaussA)*180./np.pi,np.mean(RgaussA)*180./np.pi,np.std(FgaussA)*180./np.pi,np.std(RgaussA)*180./np.pi
   
    period,orientation=[],[]
    gp,go=[],[]
    pl=8.8-(percent/100.)*8.8
    ph=8.8+(percent/100.)*8.8
    ol=90.-(percent/100.)*90.
    oh=90.+(percent/100.)*90.
    #print pl,ph,ol,oh
    
    for Fx,Fy in zip(Fx_vec,Fy_vec): #do this for all possible combinations of Fi, Rj...that's kind of a lot... 
       #for Fy in Fy_vec:
        for Rx,Ry in zip(Rx_vec,Ry_vec):
               #for Ry in Ry_vec:
            p,o=calc_period_and_orientation(Fx,Fy,Rx,Ry)
            period.append(p)
            orientation.append(o)
            if p > pl and p < ph:
                if np.abs(o) < oh and np.abs(o) > ol:#p > pl and p < ph and
                    gp.append(p)
                    go.append(o)

    period=np.array(period)
    orientation=np.array(orientation)

    #plot
    fig,(axp,axo)=plt.subplots(2)
    pbins=axp.hist(period,np.linspace(0,17.6,50))
    #axp.set_xlim([0,17.6])
    obins=axo.hist(orientation, np.linspace(np.min(orientation),np.max(orientation),n))
    axo.set_xlim([np.min(orientation),np.max(orientation)])
    axp.set_title('All possible (?) moir\'e periods and orientations within 1$\sigma$, subcoll. '+str(windownum))
    axp.set_xlabel('Period (mm)')
    axo.set_xlabel('Orientation (deg)')
    axp.set_yscale('log')
    axo.set_yscale('log')
    fig.show()

    print 'Mean period: ' , np.mean(period)
    print 'Mean orientation: ' , np.mean(orientation)
   
    #calculate probability that moire is within 10% of nominal period AND orientation
    #print len(gp), len(period)
    print 'Probability that values are within ', str(percent), '% of nominal values: ',(float(len(gp))/float(len(period)))*100 , ' %'

    
    return period,orientation

def calc_period_and_orientation(Fx,Fy,Rx,Ry):
    '''for a single value of Fx,Fy, Rx,Ry'''
    Mx=Fx-Rx
    My=Fy-Ry
    period=1./np.sqrt(Mx**2+My**2)
    orientation=np.arctan2(My,Mx)*180./np.pi
    return period,orientation

def calc_Svec(gaussP,gaussA): #(meanP,devP,meanA,devA):
    '''mean and dev in radians'''
    Sx=1000.*np.sin(gaussA)/(gaussP) 
    Sy=-1000.*np.cos(gaussA)/(gaussP) 
    return Sx,Sy    


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
    
