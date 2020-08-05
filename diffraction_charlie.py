#using Charlie Lindsey's method to test that it gives the same result as Bar-Ziv ...

#import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import scipy.constants as sc

def calc_lindsey(maxr=100):
    #define the variables
    p= np.array([0.015,0.03,0.045,0.09,0.18,0.225]) #period and pitch are the same thing, right? for a uniform grid?
    mlambda =np.arange(1.239*10**-7,1.239*10**-6,1.239*10**-8,dtype=float)  # 10-100 keV in mm
    print np.shape(mlambda)
    thetax=0.
    k=2*np.pi/mlambda
    kappa= 2*np.pi/p #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2))
    D=154.7 #mm
    x=np.arange(0,225,dtype=float)/(1000*maxr) #range from 0 to largest period
    uprimef,xstar=[],[]

    for pp in p:
        xstar.append(x/pp)

    for kap in kappa:
        uprime = []
        for i,val in enumerate(x):
            h=HsubN(D,k[0],thetax,kap,gsubn,maxr)
            J=JofX(h,kap,val,maxr)
            uprime.append( A*np.exp(1j*k[0]*D)*J*np.exp(1j*k[0]*val*np.sin(thetax)))
        uprimef.append(uprime)

    uprimef_conj=np.conj(uprimef)
    uprimef_amp = (uprimef*uprimef_conj)
    print np.max(uprimef_amp)
    #print np.shape(uprimef_amp),np.shape(mlambda)
    labels=[str(p[0]),str(p[1]),str(p[2]),str(p[3]),str(p[4]),str(p[5])]
    plot(uprimef_amp,xstar,[0,100],[0,1],'test',labels)
    return uprimef_amp
    #plot U'

def calc_lindsey_dimensionless(maxr=100):
    #define the variables
    p= np.array([0.015,0.018,0.0225,0.03,0.045,0.09]) #period and pitch are the same thing, right? for a uniform grid?
    mlambda =np.arange(1.239*10**-8,1.239*10**-7,1.239*10**-8,dtype=float)  # 10-100 keV in mm
    thetax=0.
    k=2*np.pi/mlambda
    kappa= 2*np.pi/p #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2))
    #D=[1.,1.01,1.02,1.03,1.04,1.05] #mm
    D=154.7 #mm
    #Dstar=mlambda*D/p**2 #for each p
    xstar=np.arange(0,2*maxr,dtype=float)/(maxr) #range from -1 to 1
    #x=np.arange(0,maxr,dtype=float)/(maxr) #range from 0 to 1
    uprimef=[]
    f,axes = plt.subplots(3,2)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)

    for m,pp in enumerate(p): #kappa: #for each pitch
        uprime = []
        for i,val in enumerate(xstar):
            h=HsubN(D,k[0],thetax,kappa[m],gsubn,maxr)
            J=JofXstar(h,kappa[m],val,maxr)
            uprime.append(A*np.exp(1j*k[0]*D)*J*np.exp(1j*2*np.pi*(pp/mlambda[0])*val*np.sin(thetax)))
        uprimef.append(uprime*np.conj(uprime))
        axes.flat[m].plot(xstar,uprimef[m],'.b-', label='p='+str(1000.*p[m])+'$\mu$m')
        axes.flat[m].set_xlim([0.75,1.75])
        axes.flat[m].set_ylim([0,2])
        axes.flat[m].legend(loc='upper right',fontsize='small')

    plt.suptitle('Diffraction over one grating period,100 keV')
    f.show()

    #uprimef_conj=np.conj(uprimef)
    #uprimef_amp = (uprimef*uprimef_conj)
    labels=[str(p[0]),str(p[1]),str(p[2]),str(p[3]),str(p[4]),str(p[5])]
    #plot(uprimef,xstar,[0,1],[0,1],'test',labels)
    return uprimef
    #plot U'

def calc_lindsey_energy_intensity_G1d(maxr=100):
    #define the variables
    p= np.array([0.015,0.03,0.045,0.09,0.18,0.225]) #period and pitch are the same thing, right? for a uniform grid?
    mlambda =np.arange(1.239*10**-8,1.239*10**-6,1.239*10**-9,dtype=float)  # 1-100 keV in mm
    thetax=0.
    k=2*np.pi/mlambda
    kappa= 2*np.pi/p #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2))
    D=154.7 #mm
    xstar=.5
    EkeV = (sc.h*sc.c/mlambda)/sc.e
    print np.max(EkeV),np.min(EkeV)
    #xstar=np.arange(0,2*maxr,dtype=float)/(maxr) #range from -1 to 1
    uprimef=[]
    f,axes = plt.subplots(3,2)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)

    for m,pp in enumerate(p): #kappa: #for each pitch
        uprime = []
        for i,val in enumerate(mlambda):
            h=HsubN(D,k[i],thetax,kappa[m],gsubn,maxr)
            J=JofXstar(h,kappa[m],xstar,maxr)
            uprime.append(A*np.exp(1j*k[i]*D)*J*np.exp(1j*2*np.pi*(pp/val)*xstar*np.sin(thetax)))
        uprimef.append(uprime*np.conj(uprime))
        axes.flat[m].plot(EkeV,uprimef[m],'.r-', label='p='+str(p[m])+'mm')
        axes.flat[m].set_xlim([0,50])
        axes.flat[m].set_ylim([0,2])
        axes.flat[m].legend(loc='upper right',fontsize='small')

    plt.suptitle('Lindsey algorithm for Ronchi grid diffraction')
    f.show()

    #uprimef_conj=np.conj(uprimef)
    #uprimef_amp = (uprimef*uprimef_conj)
    labels=[str(p[0]),str(p[1]),str(p[2]),str(p[3]),str(p[4]),str(p[5])]
    #plot(uprimef,xstar,[0,1],[0,1],'test',labels)
    #return uprimef
    #plot U'

def calc_lindsey_energy_intensity_I(maxr=100):
    #define the variables
    p= np.array([0.015,0.03,0.045,0.09,0.18,0.225]) #period and pitch are the same thing, right? for a uniform grid?
    mlambda =np.arange(1.239*10**-8,1.239*10**-6,1.239*10**-9,dtype=float)  # 1-100 keV in mm
    thetax=0.
    k=2*np.pi/mlambda
    kappa= 2*np.pi/p #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2))
    D=154.7 #mm
    xstar=1
    EkeV = (sc.h*sc.c/mlambda)/sc.e
    intensity=[]
    f,axes = plt.subplots(3,2)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)

    for m,pp in enumerate(p): #kappa: #for each pitch
        iint=[]
        for i,val in enumerate(mlambda):
            h=HsubN(D,k[i],thetax,kappa[m],gsubn,maxr)
            iint.append(Intensity(xstar,A,D,k[i],kappa[m],thetax,gsubn,maxr))
        #print np.shape(intensity[m]),np.shape(EkeV)
        intensity.append(iint)
        axes.flat[m].plot(EkeV,intensity[m],'.r-', label='p='+str(p[m])+'mm')
        axes.flat[m].set_xlim([0,50])
        #axes.flat[m].set_ylim([0,2])
        axes.flat[m].legend(loc='upper right',fontsize='small')

    plt.suptitle('Lindsey algorithm for Ronchi grid diffraction: intensity I vs. E (keV)')
    f.show()

    labels=[str(p[0]),str(p[1]),str(p[2]),str(p[3]),str(p[4]),str(p[5])]

def calc_lindsey_multiple_1period(maxr=100):
    #define the variables
    p= np.array([0.015,0.018,0.0225,0.03,0.045,0.09]) #period and pitch are the same thing, right? for a uniform grid?
    #
    #p= np.array([0.02,0.024,0.03,0.04,0.06,0.12]) #period and pitch are the same thing, right? for a uniform grid?
    mlambda =np.linspace(1.239*10**-8,1.239*10**-7,dtype=float,num=maxr)  # 10-100 keV in mm
    thetax=0.
    k=2*np.pi/mlambda
    kappa= 2*np.pi/p #array size 6
    A=1.0
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2))
    D=154.7 #mm
    #D=206.28 #mm
    xstar=np.arange(0,maxr,dtype=float)/(.5*maxr) #range from -1 to 1
    uprimef=[]
    f,axes = plt.subplots(3,2)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)

    for m,pp in enumerate(p): #kappa: #for each pitch
        uprime,fgrid = [],[]
        #fgrid.append(np.sum(gsubn*np.exp(1j*n*2*np.pi*[x]/pp)) for x in xstar)
        print np.shape(fgrid)
        for i,val in enumerate(xstar):
            fgrid.append(np.sum(gsubn[i]*np.exp(1j*(i-maxr/2)*2*np.pi*val/pp)))
            h=HsubN(D,k[i],thetax,kappa[m],gsubn,maxr)
            #hprime=Hprime(h,stuff)
            J=JofXstar(h,kappa[m],val,maxr)
            uprime.append(A*np.exp(1j*k[i]*D)*J*np.exp(1j*2*np.pi*(pp/mlambda[m])*val*np.sin(thetax)))
        uprimef.append(uprime)#*np.conj(uprime))
        #axes.flat[m].plot(xstar,uprimef[m]*np.array(fgrid),'.r-', label='p='+str(p[m])+'mm')
        axes.flat[m].plot(xstar,uprimef[m]*np.conj(uprimef[m]),'.r-', label='p='+str(p[m])+'mm')
        axes.flat[m].plot(xstar,fgrid*np.conj(fgrid),'.b-', label='p='+str(p[m])+'mm')
        axes.flat[m].set_xlim([0.5,1.5])
        axes.flat[m].set_ylim([0,2])
        axes.flat[m].legend(loc='upper right',fontsize='small')
        print pp,np.min(uprimef[m]*np.array(fgrid)),np.max(uprimef[m]*np.array(fgrid))

    plt.suptitle('Lindsey algorithm for Ronchi grid diffraction')
    f.show()

def calc_difference_sep(maxr=100):
    p1= np.array([0.015,0.018,0.0225,0.03,0.045,0.09]) #period and pitch are the same thing, right? for a uniform grid?
    p2= np.array([0.02,0.024,0.03,0.04,0.06,0.12]) #period and pitch are the same thing, right? for a uniform grid?
    mlambda =np.linspace(1.239*10**-8,1.239*10**-7,dtype=float,num=maxr)  # 10-100 keV in mm
    thetax=0.
    k=2*np.pi/mlambda
    kappa1= 2*np.pi/p1 #array size 6
    kappa2= 2*np.pi/p2 #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2))
    D1=154.7 #mm
    D2=206.28 #mm
    xstar=np.arange(0,maxr,dtype=float)/(.5*maxr) #range from -1 to 1
    uprimef1,uprimef2=[],[]
    f,axes = plt.subplots(3,2)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)

    for m,pp1,pp2 in enumerate(zip(p1,p2)): #kappa: #for each pitch
        uprime = []
        for i,val in enumerate(xstar):
            h1=HsubN(D,k[i],thetax,kappa1[m],gsubn,maxr)
            h2=HsubN(D,k[i],thetax,kappa2[m],gsubn,maxr)
            #hprime=Hprime(h,stuff)
            J1=JofXstar(h1,kappa1[m],val,maxr)
            J2=JofXstar(2,kappa2[m],val,maxr)
            uprime1.append(A*np.exp(1j*k[i]*D1)*J1*np.exp(1j*2*np.pi*(pp1/mlambda[0])*val*np.sin(thetax)))
            uprime2.append(A*np.exp(1j*k[i]*D2)*J2*np.exp(1j*2*np.pi*(pp2/mlambda[0])*val*np.sin(thetax)))
        uprimef1.append(uprime1*np.conj(uprime1))
        uprimef2.append(uprime2*np.conj(uprime2))
        axes.flat[m].plot(xstar,uprimef2[m]-uprimef1[m],'.r-', label='p='+str(p[m])+'mm')
        axes.flat[m].set_xlim([0.5,1.5])
        axes.flat[m].set_ylim([0,2])
        axes.flat[m].legend(loc='upper right',fontsize='small')

    plt.suptitle('Difference in diffraction between D=206.28mm and D=154.7')
    f.show()


def calc_lindsey_multiple_energy_intensity(maxr=100):
    #define the variables
    p= np.array([0.015,0.03,0.045,0.09,0.18,0.225]) #period and pitch are the same thing, right? for a uniform grid?
    mlambda =np.arange(1.239*10**-8,1.239*10**-6,1.239*10**-9,dtype=float)  # 1-100 keV in mm
    thetax=0.
    k=2*np.pi/mlambda
    kappa= 2*np.pi/p #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2))
    D=154.7 #mm
    xstar=1
    EkeV = (sc.h*sc.c/mlambda)/sc.e
    print np.max(EkeV),np.min(EkeV)
    #xstar=np.arange(0,2*maxr,dtype=float)/(maxr) #range from -1 to 1
    uprimef=[]
    f,axes = plt.subplots(3,2)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)

    for m,pp in enumerate(p): #kappa: #for each pitch
        uprime = []
        for i,val in enumerate(mlambda):
            h=HsubN(D,k[i],thetax,kappa[m],gsubn,maxr)
            J=JofXstar(h,kappa[m],xstar,maxr)
            uprime.append(A*np.exp(1j*k[i]*D)*J*np.exp(1j*2*np.pi*(pp/val)*xstar*np.sin(thetax)))
        uprimef.append(uprime*np.conj(uprime))
        axes.flat[m].plot(EkeV,uprimef[m],'.r-', label='p='+str(p[m])+'mm')
        axes.flat[m].set_xlim([0,50])
        axes.flat[m].set_ylim([0,2])
        axes.flat[m].legend(loc='upper right',fontsize='small')

    plt.suptitle('Lindsey algorithm for Ronchi grid diffraction')
    f.show()

    #uprimef_conj=np.conj(uprimef)
    #uprimef_amp = (uprimef*uprimef_conj)
    labels=[str(p[0]),str(p[1]),str(p[2]),str(p[3]),str(p[4]),str(p[5])]
    #plot(uprimef,xstar,[0,1],[0,1],'test',labels)
    #return uprimef
    #plot U'

def test_Ronchi_grid(maxr=100):
    #define the variables
    #p= np.array([0.015,0.03,0.045,0.09,0.18,0.225]) #period and pitch are the same thing, right? for a uniform grid?
    #mlambda =np.arange(1.239*10**-7,1.239*10**-6,1.239*10**-8,dtype=float)  # 10-100 keV in mm
    #print np.shape(mlambda)
    thetax=0.
    k=2*np.pi
    kappa= 2*np.pi #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2)) #g_n = h sinc (nh), h=0.5 for ronchi
    D=[1.,1.01,1.02,1.03,1.04] #mm
    x=np.arange(0,maxr,dtype=float)/(maxr) #range from 0 to 1
    #x=[-1.,-.5,0,.5,.999]
    uprimef=[]
    f,axes = plt.subplots(5, sharex=True, sharey=True)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)
    #uprime = []

    #val=-1.
    for j,r in enumerate(D):
    #r=1.
        uprime = []
        #print 'test', np.exp(2*np.pi*1j)
        for i,val in enumerate(x):
            h=HsubN(r,k,thetax,kappa,gsubn,maxr)
            J=JofX(h,kappa,val,maxr)
            uprime.append( A*np.exp(1j*k*r)*J*np.exp(1j*k*val*np.sin(thetax)))
        #print J,np.real(uprime)
        #print np.shape(uprime),np.shape(np.arange(-maxr/2,maxr/2))
        axes[j].plot(x,uprime*np.conj(uprime))
        axes[j].set_xlim([0,1])
        axes[j].set_ylim([0,1.4])

    plt.suptitle('Lindsey algorithm for Ronchi grid diffraction over one period')
    f.show()

    uprimef.append(uprime)

    uprimef_conj=np.conj(uprimef)
    uprimef_amp = (uprimef*uprimef_conj)
    #labels=[str(D[0]),str(D[1]),str(D[2]),str(D[3]),str(D[4]),str(D[5])]
    #plot(uprimef_amp,x,[0,1],[0,10],'test',labels)

    #g,ax2=plt.subplots(1)
    #ax2.plot(x,uprime*np.conj(uprime))
    #plt.show()

    return uprimef_amp
    #plot U'

def JofX(h,kappa,x,maxr):
    nn=np.arange(-maxr/2,maxr/2)
    J=np.zeros(maxr)
    for i,n in enumerate(nn):
        J[i]=h[i]*np.exp(1j*kappa*n*x)
        #print J[i] #should be real
    sumJ=np.sum(J) #do I need to pick an axis
    #print sumJ
    return sumJ

def JofXstar(h,kappa,xstar,maxr):
    nn=np.arange(-maxr/2,maxr/2)
    J=np.zeros(maxr)
    for i,n in enumerate(nn):
        J[i]=h[i]*np.exp(1j*2*np.pi*n*xstar) #don't think it matters about h_n to be h_n(D) or h_n(D*)
        #print J[i] #should be real
    sumJ=np.sum(J) #do I need to pick an axis
    #print sumJ
    return sumJ

def HsubN(D,k,thetax,kappa,gsubn,maxr):
    nn=np.arange(-maxr/2,maxr/2)
    h=np.zeros(maxr)
    for i,n in enumerate(nn):
        h[i]= gsubn[i]*np.exp((-1j*D/(2*k))*(k*np.sin(thetax)+kappa*n)**2)
    return h

def HsubPprime(D,k,thetax,kappa,gsubn,maxr):
    nn=np.arange(-maxr/2,maxr/2)
    hprime=np.zeros(maxr)
    h=HsubN(D,k,thetax,kappa,gsubn,maxr)
    for i,n in enumerate(nn):
        hprime[i]= h[i]*gsubn[i]*np.exp(-1j*D*(kappa**2*n**2)/(2*k))*np.exp(-1j*kappa*D*n*np.sin(thetax))
    return hprime

def bigK(kappa,x,hprime,maxr):
    pp=np.arange(-maxr/2,maxr/2)
    K=np.zeros(maxr)
    for i,p in enumerate(pp):
        K[i]=np.exp(1j*kappa*p*x)*hprime[i]
    return np.sum(K)

def HsubP2(D1,D2,k,thetax,kappa,gsubn,gsubn0,maxr): #do gsubn1 and gsubn2 have to be different?
    #HsubN_theta is HsubN calculated for the previous grid's parameters D,k,theatx,kappa,etc
    #GsubP-N is gsubn for the previous grid?
    n1=np.arange(-maxr/2,maxr/2)
    n2=np.arange(-maxr/2,maxr/2)
    h=[]
    hm=[]
    deltad=D1-D2

    for j,m in enumerate(n2):
        hm.append(gsubn[j]*np.exp(-1j*kappa*m*D2*np.sin(thetax)))
    for i,n in enumerate(n1):
        h.append(gsubn[i]*np.exp(-1j*kappa*n*D1*np.sin(thetax))*np.exp(-1j*kappa**2*n*2*deltad/(2*k))*np.sum(hm))
    #print h[0],h[99]
    #fig,ax=plt.subplots()
    #ax.plot(n2,np.real(hm),'r')
    #ax.plot(n1,np.real(h),'b')
    #ax.plot(n2,np.imag(hm),'r-')
    #ax.plot(n1,np.imag(h),'b-')
    #fig.show()
    return h*gsubn0

def test_UDoublePrime(maxr=100,D=[1.,1.01,1.02,1.03,1.04]):
    import scipy.integrate as integrate
    thetax=0.
    k=2*np.pi
    kappa= 2*np.pi #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2)) #g_n = h sinc (nh), h=0.5 for ronchi
    T=1. #period of grid...
    #D=[1.,1.01,1.02,1.03,1.04] #mm
    #D=[1.05,1.1,1.15,1.2,1.25]
    x=np.arange(0,maxr,dtype=float)/(maxr) #range from 0 to 1

    f,axes = plt.subplots(5, sharex=True, sharey=True)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)

    for j,r in enumerate(D):
        intensity,hprime,K= [],[],[]
        UdoublePrime=np.zeros(maxr)
        for i,val in enumerate(x):
            hprime.append(HsubPprime(r,k,thetax,kappa,gsubn,maxr))
            K.append(bigK(kappa,val,hprime[i],maxr))
            UdoublePrime[i]=A*np.exp(1j*k*D[j])*K[i]*np.exp(1j*k*val*np.sin(thetax))
            gridFT=np.sum(gsubn[i]*np.exp(1j*i*2*np.pi*val/T))
        #print np.shape(x),np.shape(UdoublePrime),type(x),type(UdoublePrime)
            def integrand(UdoublePrime):
                return UdoublePrime*np.conj(UdoublePrime)#*gridFT#*np.conj(UdoublePrime)/np.conj(gsubn)
            intensity.append(integrate.quad(integrand,val,val+.5))
        print np.shape(x),np.shape(intensity)
        axes[j].plot(x,intensity,label='D='+str(D[j]))
        axes[j].set_xlim([0,1])
        axes[j].legend()
        #axes[j].set_ylim([0,0.5])
        print np.min(intensity*np.conj(intensity)),np.max(intensity*np.conj(intensity))

    plt.suptitle('Lindsey algorithm for Ronchi grid diffraction over one period')
    f.show()

def test_BZIntensity(maxr=100,D=[1.,1.01,1.02,1.03,1.04]):
    '''I=integral over one period of U'U'*g'''
    import scipy.integrate as integrate
    thetax=0.
    k=2*np.pi
    kappa= 2*np.pi #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2)) #g_n = h sinc (nh), h=0.5 for ronchi
    T=1. #period of grid...
    #D=[1.,1.01,1.02,1.03,1.04] #mm
    #D=[1.05,1.1,1.15,1.2,1.25]
    x=np.arange(0,maxr,dtype=float)/(maxr) #range from 0 to 1

    f,axes = plt.subplots(5, sharex=True, sharey=True)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)

    for j,r in enumerate(D):
        intensity,hsubn,J= [],[],[]
        UPrime=np.zeros(maxr)
        gridFT=np.zeros(maxr)
        for i,val in enumerate(x):
            hsubn.append(HsubN(r,k,thetax,kappa,gsubn,maxr))
            J.append(bigK(kappa,val,hsubn[i],maxr))
            UPrime[i]=A*np.exp(1j*k*D[j])*J[i]*np.exp(1j*k*val*np.sin(thetax))
            gridFT[i]=np.sum(gsubn[i]*np.exp(1j*i*2*np.pi*val/T))
        #print np.shape(x),np.shape(UdoublePrime),type(x),type(UdoublePrime)
            def integrand(UPrime):
                return UPrime*np.conj(UPrime)*gsubn#*gridFT#*gridFT#*np.conj(UdoublePrime)/np.conj(gsubn)
            #print np.shape(UPrime),np.shape(gsubn),type(UPrime[0]),np.mean(UPrime*np.conj(UPrime)*gsubn)
            #intensity.append(integrate.quad(integrand,val+0,val+1)) #need to make a new 'x' vector as a result of the integration?
        print np.min(gridFT),np.mean(gridFT),np.max(gridFT)
        axes[j].plot(x,UPrime*np.conj(UPrime)*gridFT,label='D='+str(D[j]))
        axes[j].set_xlim([0,1])
        axes[j].legend()
        #axes[j].set_ylim([0,0.5])
        #print np.min(intensity*np.conj(intensity)),np.max(intensity*np.conj(intensity))

    plt.suptitle('Lindsey algorithm for Ronchi grid diffraction over one period')
    f.show()

def largeTheta(maxr=100):
    '''Test the large theta limit of the algorithm. It should reproduce what is calculated here, that is: the intensity in the detector plane is the product of the Fresnel diffraction patterns of the individual grids onto the detector plane, ie: U(x) = U1(x)*U2(x)*..*Un(x)'''
    thetax=1. #fix
    k=2*np.pi
    kappa= 2*np.pi #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2)) #g_n = h sinc (nh), h=0.5 for ronchi
    D1=100 #mm
    D2=0
    x=np.arange(0,maxr,dtype=float)/(maxr) #range from 0 to 1
    #x=[-1.,-.5,0,.5,.999]
    uprime,uprime1,uprime2=[],[],[]
    for i,val in enumerate(x):
        h=HsubN(D1,k,thetax,kappa,gsubn,maxr)
        J=JofX(h,kappa,val,maxr)
        u1= A*np.exp(1j*k*D1)*J*np.exp(1j*k*val*np.sin(thetax))
        h=HsubN(D2,k,thetax,kappa,gsubn,maxr)
        J=JofX(h,kappa,val,maxr)
        u2= A*np.exp(1j*k*D2)*J*np.exp(1j*k*val*np.sin(thetax))
        uprime1.append(u1)
        uprime2.append(u2)
        uprime.append(u1*u2)

    uprimef_conj=np.conj(uprime)
    uprimef_amp = (uprime*uprimef_conj)

    f,ax =plt.subplots()# plt.subplots(5, sharex=True, sharey=True)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)
    #uprime = []
    ax.plot(x,uprime1*np.conj(uprime1),'r')
    ax.plot(x,uprime2*np.conj(uprime2),'b')
    ax.plot(x,uprime*np.conj(uprime),'k-')
    ax.set_xlim([0,1])
    #ax.set_ylim([0,1.4])

    plt.suptitle('Expected result in the limit of large theta, D1=100mm, D2=10mm, theta=1')
    f.show()

    x=[0,0.25,0.5,.75,1.]
    theta=np.arange(100,100*np.pi)/100.
    f,ax =plt.subplots(5, sharex=True, sharey=True)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)
    for i,val in enumerate(x):
        uprime,uprime1,uprime2=[],[],[]
        for j,thetax in enumerate(theta):
            h=HsubN(D1,k,thetax,kappa,gsubn,maxr)
            J=JofX(h,kappa,val,maxr)
            u1= A*np.exp(1j*k*D1)*J*np.exp(1j*k*val*np.sin(thetax))
            h=HsubN(D2,k,thetax,kappa,gsubn,maxr)
            J=JofX(h,kappa,val,maxr)
            u2= A*np.exp(1j*k*D2)*J*np.exp(1j*k*val*np.sin(thetax))
            uprime1.append(u1)
            uprime2.append(u2)
            uprime.append(u1*u2)

        uprimef_conj=np.conj(uprime)
        uprimef_amp = (uprime*uprimef_conj)
        #ax[i].plot(theta,uprime1*np.conj(uprime1),'r')
        #ax[i].plot(theta,uprime2*np.conj(uprime2),'b')
        ax[i].plot(theta,uprime*np.conj(uprime),'k-',label='x= '+str(x[i]))
        ax[i].set_xlim([1,np.pi])
        ax[i].legend()
        #ax.set_ylim([0,1.4])

    plt.suptitle('Expected result in the limit of large theta, D1=100mm, D2=10mm')
    f.show()

    return uprime,uprime1,uprime2

def TwentyThree(maxr=100):
    '''Calculate the intensity through multiple grids using lindsey eq 24a'''
    #define x, kappa
    d=np.arange(0,100)/10. #0-100 mm away from last grid
    D1=154.+d
    #D1=[154.,164.,254.,1154.]
    D2=d
    #x=
    thetax=.1
    k=2*np.pi
    kappa= 2*np.pi #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2))
    gsubn0=gsubn
    intensity=[]

    for n in range(0,len(d)):
        h=HsubP2(D1[n],D2[n],k,thetax,kappa,gsubn,gsubn0,maxr)
        hstar=np.conj(h)
        intensity.append(np.sum(h*hstar)) #should be a function of d
        #print D1,D2,h[0],h[99],intensity[n]
    #print intensity[0:10]
    #plot

    fig,ax=plt.subplots()
    ax.plot(d,np.real(intensity))
    ax.set_ylim([0,1])
    #ax.set_yscale('log')
    fig.show()

    return intensity

def test_alg_largeTheta(maxr=100):
    '''See if results from largeTheta() can be reproduced by the algorithm'''
    thetax=0 #let's see if we get the same thing as last time...
    k=2*np.pi
    kappa= 2*np.pi #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2)) #g_n = h sinc (nh), h=0.5 for ronchi
    gsubn0=gsubn
    D1=100 #mm
    D2=0
    D=0 #test with this
    n=np.arange(0,maxr)
    m=n
    x=np.arange(0,maxr,dtype=float)/(maxr) #range from 0 to 1
    #x=[-1.,-.5,0,.5,.999]
    uprime,uprime1,uprime2=[],[],[]
    for i,val in enumerate(x):
        h=HsubP2(D1,D2,k,thetax,kappa,gsubn,gsubn0,maxr)
        bigK=np.sum(h*np.exp(1j*kappa*i*val*n))
        uprime.append(A*np.exp(1j*k*D)*bigK*np.exp(1j*k*val*np.sin(thetax)))
        #intensity[i]=(np.sum(h*hstar)) #should be a function of d

    f,ax =plt.subplots()# plt.subplots(5, sharex=True, sharey=True)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)
    #uprime = []
    ax.plot(x,uprime*np.conj(uprime))
    ax.set_xlim([0,1])
    #ax.set_ylim([0,1.4])

    plt.suptitle('Expected result in the limit of large theta, D1=100mm, D2=10mm, theta=1')
    f.show()


def lindseyfig2(maxr=100):
    '''reproduce Lindsey figure 2 to test algorithm. 1m between grids with 50 um slits separated by 100 um => period 150 um, wavelengths 0,1 and 10 A'''
    D2=0.001
    D1=1000.+D2 #1m = 1000mm
    thetax=np.arange(0,20*np.pi)/10 #arcsec
    mlambda = np.array([10**-20,10**-7,10**-6]) # in mm
    k=2*np.pi/mlambda
    p=.1
    kappa= 2*np.pi/p #array size 6
    A=1.
    gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2))
    intensity=np.empty([3,len(thetax)])
    for j,kk in enumerate(k):
        for n in range(0,len(thetax)):
            h=HsubP2(D1,D2,kk,thetax[n],kappa,gsubn,maxr)
            hstar=np.conj(h)
            intensity[j,n]=(np.sum(h*hstar)) #should be a function of d

    #plot

    fig,ax=plt.subplots()
    ax.plot(thetax,np.real(intensity[0,:]),'r')
    #ax.plot(thetax,np.real(intensity[1,:]),'b')
    #ax.plot(thetax,np.real(intensity)[2,:],'g')
    #ax.plot(thetax,np.imag(intensity[0,:]),'-')
    #ax.plot(thetax,np.imag(intensity[1,:]),'b-')
    #ax.plot(thetax,np.imag(intensity)[2,:],'g-')

    ax.set_ylim([0,1])
    #ax.set_yscale('log')
    fig.show()

    return intensity


def Intensity_of_x(maxr=100):
    intensity=[]
    for x in range(0,maxr):
        intensity.append(TwentyFourA(x))

    #plot
    return intensity

def Intensity(xstar,A,D,k,kappa,thetax,gsubn,maxr,g2=False):
    '''Calculate I= integral_x*^x*+a G_1d(x*)G_2(x*) dx* for a single value of x* and a single wavelength'''
    integrand,intensity=[],[]
    xxstar=np.arange(xstar,xstar+1) #maxr)/maxr #integrate over one pitch
    #gsubn=0.5*np.sinc(0.5*np.arange(-maxr/2,maxr/2))
    #if not g2:
    #    g2= gsubn*np.exp(1j*2*np.pi*xxstar)#fourier series of Ronchi grid. Sum of g2 over 1 period = 1,max amp = 0.5

    h=HsubN(D,k,thetax,kappa,gsubn,maxr)

    for x in xxstar:
        g2= np.sum(gsubn*np.exp(1j*2*np.pi*x))
        J=JofXstar(h,kappa,x,maxr)
        integrand.append(A*np.exp(1j*k*D)*J*np.exp(1j*k*x*np.sin(thetax))*g2)  #G1 affected by diffraction, G2 is just good ol' sum(g_n e^ikx)
    intensity=np.sum(integrand) #still has to be I(x*) ... does this have to be normalized? Also, shouldn't it be real?
    #print intensity
    return intensity

def plot(uprimef,mlambda,xran,yran,title,label):
    # Three subplots sharing both x/y axes
    f, (ax1, ax2, ax3,ax4,ax5,ax6) = plt.subplots(6, sharex=True, sharey=True)
    fig = plt.gcf()
    fig.set_size_inches(3.5, 6.5)
    if np.shape(np.shape(mlambda))[0] == 1:
        ax1.plot(mlambda, uprimef[0], color="y",label=label,linewidth='2')
        ax2.plot(mlambda, uprimef[1], color="g",label=label,linewidth='2')
        ax3.plot(mlambda, uprimef[2], color="m",label=label,linewidth='2')
        ax4.plot(mlambda, uprimef[3], color="c",label=label,linewidth='2')
        ax5.plot(mlambda, uprimef[4], color="k",label=label,linewidth='2')
        ax6.plot(mlambda, uprimef[5], color="r",label=label,linewidth='2')
    else:
        ax1.plot(mlambda[0], uprimef[0], color="y",label=label,linewidth='2')
        ax2.plot(mlambda[1], uprimef[1], color="g",label=label,linewidth='2')
        ax3.plot(mlambda[2], uprimef[2], color="m",label=label,linewidth='2')
        ax4.plot(mlambda[3], uprimef[3], color="c",label=label,linewidth='2')
        ax5.plot(mlambda[4], uprimef[4], color="k",label=label,linewidth='2')
        ax6.plot(mlambda[5], uprimef[5], color="r",label=label,linewidth='2')

    #plt.xlim(xran)
    #plt.ylim(yran)
    plt.show()

def calc_functions(maxr,n,p,Rstar,mlambda, function):
    '''Calculate the diffraction (Bar-Ziv eqs 8 and 21)'''
    x=np.arange(0,maxr,dtype=float)/(maxr) #this is still in meters right?
    xstar = x/p
    sum_Asubn = np.zeros((0))
    if function == 'G1d':
        sum_Bsubn = np.zeros((0))
        sum_Asin = np.zeros((0))
        sum_Acos = np.zeros((0))

        for xx in xstar: #sum first - x* will remain the variable. This is for fixed R*.
            sum_Asin = np.append(sum_Asin, np.sum(np.sin(2*np.pi*xx*(2*n+1))))
            sum_Acos = np.append(sum_Acos, np.sum((np.cos(np.pi*Rstar*(2*n+1)**2))/(2*n+1)))
            sum_Asubn = np.append(sum_Asubn, np.sum(np.sin(2*np.pi*xx*(2*n+1))*np.cos(np.pi*Rstar*(2*n+1)**2)/(2*n+1)))
            sum_Bsubn = np.append(sum_Bsubn, np.sum(np.sin(2*np.pi*xx*(2*n+1))*np.sin(np.pi*Rstar*(2*n+1)**2)/(2*n+1)))

        G1d = .25 + (2/np.pi)*sum_Asubn + (4/np.pi**2)*(sum_Asubn**2+sum_Bsubn**2) #from eq 8
        var=G1d

    if function == 'G1d(R)':
        sum_Bsubn = np.zeros((0))
        sum_Asin = np.zeros((0))
        sum_Acos = np.zeros((0))
        xstar =0.20
        for rr in Rstar:
            sum_Asin = np.append(sum_Asin, np.sum(np.sin(2*np.pi*xstar*(2*n+1))))
            sum_Acos = np.append(sum_Acos, np.sum((np.cos(np.pi*rr*(2*n+1)**2))/(2*n+1)))
            sum_Asubn = np.append(sum_Asubn, np.sum(np.sin(2*np.pi*xstar*(2*n+1))*np.cos(np.pi*rr*(2*n+1)**2)/(2*n+1)))
            sum_Bsubn = np.append(sum_Bsubn, np.sum(np.sin(2*np.pi*xstar*(2*n+1))*np.sin(np.pi*rr*(2*n+1)**2)/(2*n+1)))

        G1d = .25 + (2/np.pi)*sum_Asubn + (4/np.pi**2)*(sum_Asubn**2+sum_Bsubn**2) #from eq 8
        var=G1d


    if function == 'Id(x)' : #moire image. xstar is now xi
        for xx in xstar:
            sum_Asubn = np.append(sum_Asubn, np.sum(np.cos(2*np.pi*xx*(2*n+1))*np.cos(np.pi*Rstar*(2*n+1)**2)/(2*n+1)**2))
        Id  = .25 + (2/np.pi**2)*sum_Asubn
        var=Id
    if function == 'Id(R)': #contrast - eq. 35. Now I is a function of R, not xi
        xstar = 0.0
        for rr in Rstar:
           sum_Asubn = np.append(sum_Asubn,np.sum(np.cos(2*np.pi*xstar*(2*n+1))*np.cos(np.pi*rr*(2*n+1)**2)/(2*n+1)))
        #print np.shape(sum_Asubn)
        Id = .25 + (2/np.pi**2)*sum_Asubn
        var=Id

    if function == 'C': #contrast - eq. 35. Now I is a function of R, not xi
        #C(R*)=[Id'(R*)-Id'(R*=1/2)]/Id'(R*=1/2)]
        xstar = 0.0

        for rr in Rstar:
           sum_Asubn = np.append(sum_Asubn,np.sum(np.cos(2*np.pi*xstar*(2*n+1))*np.cos(np.pi*rr*(2*n+1)**2)/(2*n+1)))

        Id_Rhalf = np.zeros((len(Rstar)))+ .25 + (2/np.pi**2)*np.sum(np.cos(2*np.pi*xstar*(2*n+1))*np.cos(np.pi*.5*(2*n+1)**2)/(2*n+1))
        Id = .25 + (2/np.pi**2)*sum_Asubn
        #print np.shape(sum_Asubn),np.shape(Id),np.shape(Id_Rhalf)
        C=np.abs(Id-Id_Rhalf)/Id_Rhalf
        var=C#[C,Id,Id_Rhalf]
        #print Id_Rhalf[0],Id[0],Id[maxr-1], np.sum(np.cos(np.pi*(2*n+1)**2)/(2*n+1)**2)
    return var



