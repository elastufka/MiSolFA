
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import scipy.constants as sc
import diffraction as d


def plot_tests(xstar,var,xran,yran,title,label):

    if label != 'C':
        # Three subplots sharing both x/y axes
        f, (ax1, ax2, ax3,ax4,ax5) = plt.subplots(5, sharex=True, sharey=True)
        fig = plt.gcf()
        fig.set_size_inches(3.5, 6.5)

        #ax1.plot(xstar, sum_Asubn[0], color="g",label="A_n",linewidth='2')
        #ax1.plot(xstar, sum_Bsubn[0], color="r",label="B_n",linewidth='2')
        ax1.plot(xstar, var[0], color="y",label=label,linewidth='2')

        #ax2.plot(xstar, sum_Asubn[1], color="g",label="A_n",linewidth='2')
        #ax2.plot(xstar, sum_Bsubn[1], color="r",label="B_n",linewidth='2')
        ax2.plot(xstar, var[1], color="y",label=label,linewidth='2')

        #ax3.plot(xstar, sum_Asubn[2], color="g",label="A_n",linewidth='2')
        #ax3.plot(xstar, sum_Bsubn[2], color="r",label="B_n",linewidth='2')
        ax3.plot(xstar, var[2], color="y",label=label,linewidth='2')

        #ax4.plot(xstar, sum_Asubn[3], color="g",label="A_n",linewidth='2')
        #ax4.plot(xstar, sum_Bsubn[3], color="r",label="B_n",linewidth='2')
        ax4.plot(xstar, var[3], color="y",label=label,linewidth='2')

        #ax5.plot(xstar, sum_Asubn[4], color="g",label="A_n",linewidth='2')
        #ax5.plot(xstar, sum_Bsubn[4], color="r",label="B_n",linewidth='2')
        ax5.plot(xstar, var[4], color="y",label=label,linewidth='2')

    if label == 'C':
        f = plt.figure()
        ax1=f.add_subplot(111)

        ax1.plot(xstar, var, color='r', label=label) #here xstar is Rstar
        #ax1.plot(xstar, var[1], color='g', label='Id') #here xstar is Rstar
        #ax1.plot(xstar, var[2], color='b', label='Id_Rhalf') #here xstar is Rstar
    ax1.set_title(title)
    #f.subplots_adjust(hspace=0)
    #plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.xlabel('x*')
    #plt.ylabel('Normalized Intensity')
    ax1.set_ylim(yran[0],yran[1])
    #ax1.set_xlim(xran[0],xran[1])
    ax1.legend(loc='upper right',fontsize='medium')

    f.show()
    #plt.savefig('nequals1000.png', bbox_inches='tight')
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    X,Y=np.meshgrid(xstar,Rstar)
#    ax.plot_surface(X, Y,var)
#    fig.show()

def reproduce_plots():
    '''reproduce plots from Bar Ziv 1985'''
    maxr=1000
    n=np.arange(0,maxr, dtype=float) #n for the sum
    p=1#p= 15*10**-6 #pitch of finest grids = .015 mm
    R=.15 #150-200 mm
    mlambda = 1#6.206*10**-11#6.206^-12 median wavelength in 2-200 keV
    delta = [0.0,0.01,0.02,0.03,0.04] #see eq. 33 and 34
    F1=np.array([p**2,p**2,p**2,p**2,p**2])/mlambda #first Fourier image plane
    m=[2,2,2,2,2]
    Fm = m*F1
    Rstar = Fm + delta*F1 #R*mlambda/(p**2) #but I don't get square waves when R* = 1,2, ... like I should!
    print Rstar, F1,Fm
    data,data2=[],[]
    for rr in Rstar:
        data.append(d.calc_functions(maxr,n,p,rr,mlambda, 'G1d'))
        data2.append(d.calc_functions(maxr,n,p,rr,mlambda, 'Id(x)'))

    x=np.arange(0,maxr,dtype=float)/(maxr) #this is still in meters right?
    xstar = x/p
    #data=calc_functions(maxr,n,p,R,mlambda, 'G1d')
    #plot_tests(xstar,data) #reproduce fig 2A
    #plot_tests(xstar,data2)

    #plot the contrast now
    Rstar = np.arange(0,maxr,dtype=float)/maxr
    Rstar[maxr-1]=1

    data3= d.calc_functions(maxr,n,p,Rstar,mlambda, 'C')
    print np.shape(data3[2]),np.shape(Rstar)
    #yran=[np.min(data3),np.max(data3)]
    plot_tests(Rstar, data3,[0,1],[0,1], 'Contrast of Moire image', 'C')

def misolfa_G1d(maxr=1000):
    '''Now with parameters from MiSolFA grids. Plot I(R*), and I(x) for a single grid -> G1d(x*)'''
    #maxr=1000
    n=np.arange(0,maxr, dtype=float) #n for the sum
    xmm=n/maxr
    p=np.array([0.015,0.03,0.045,0.09,0.18,0.225]) #pitch of finest grids = .015 mm
    #xran=[0,8.8*10**-3] #moire period of 8.8mm
    R=154.7 #150-200 mm
    mlambda = 1.239*10**-7  # 10 keV in mm
    #mlambda=6.199*10**-8 #20 keV to mm
    #mlambda = 4.132*10**-8  # 30 keV in mm
    Rstar1 = mlambda*R/p**2 # + delta#Fm + delta*F1 #R*mlambda/(p**2) #but I don't get square waves when R* = 1,2, ... like I should!
    xstar,Rstar,data,data2=[],[],[],[]
    x=n/maxr
    
    for rr,pp in zip(Rstar1,p):
        data.append(d.calc_functions(maxr,n,pp,rr,mlambda, 'G1d')) #might need to modify the xrange on this
        xstar.append(x/pp)
    for pp in p:
        Rstar.append((np.zeros((maxr)) +mlambda*R)/pp**2)
        data2.append(d.calc_functions(maxr,n,pp,Rstar,mlambda, 'Id(R)'))
    #yran=[0,np.max(data)]
    #print np.shape(p),np.shape(data)
    #print 'energy range', EkeV[0],EkeV[maxr-1],np.min(Rstar),np.max(Rstar)
    
    f,axes = plt.subplots(3,2, sharex=False,sharey=True)
    title='lambda=620 A (20 keV), R=154.7mm'
    for i,ax in enumerate(f.axes):
        ax.plot(xstar[i], data[i], '.r-',label='p='+str(p[i])+'mm') #here xstar is Rstar
        ax.set_xlim(0,1)
        #ax.set_ylim(0,1)
        ax.legend(loc='upper right',fontsize='small')
    plt.suptitle(r"$G_{1d}$ intensity distribution calculated for pitches of MiSolFA grids, R=154.7mm, 10keV")
    f.text(0.5, 0.04, 'x*', ha='center')
    f.text(0.04, 0.5, 'Intensity', va='center', rotation='vertical')
    #plt.savefig('G1d_R154_10keV.png')
    f.show()

def misolfa_G1dR():
    '''Now with parameters from MiSolFA grids. Plot I(R*), and I(x) for a single grid -> G1d(x*)'''
    maxr=1000
    n=np.arange(0,maxr, dtype=float) #n for the sum
    xmm=n/maxr
    p=np.array([0.015,0.03,0.045,0.09,0.18,0.225]) #pitch of finest grids = .015 mm
    #xran=[0,8.8*10**-3] #moire period of 8.8mm
  #  R=154.7 #150-200 mm
    R=500*n/maxr
    mlambda = 1.239*10**-7  # 10 keV in mm
    #mlambda=6.199*10**-8 #20 keV to mm
    #mlambda = 4.132*10**-8  # 30 keV in mm
    #Rstar1 = mlambda*R/p**2 # + delta#Fm + delta*F1 #R*mlambda/(p**2) #but I don't get square waves when R* = 1,2, ... like I should!
    xstar,Rstar,data,data2=[],[],[],[]
    x=n/maxr
    
    for pp,i in zip(p,range(0,6)):
        Rstar.append(np.zeros((maxr)) +mlambda*R/pp**2)
        data.append(d.calc_functions(maxr,n,pp,Rstar[i],mlambda, 'G1d(R)')) #might need to modify the xrange on this
        #xstar.append(x/pp)
        print np.min(Rstar[i]),np.max(Rstar[i]),np.min(data[i]),np.max(data[i])
    f,axes = plt.subplots(3,2, sharex=False,sharey=True)
    title='lambda=620 A (20 keV), R=154.7mm'
    for i,ax in enumerate(f.axes):
        oldr=Rstar[i]*p[i]**2/mlambda
        ax.plot(oldr, data[i], '.r-',label='p='+str(p[i])+'mm') #here xstar is Rstar
        ax.set_xlim(0,10)
        #ax.set_ylim(0,1)
        ax.legend(loc='upper right',fontsize='small')
    plt.suptitle(r"$G_{1d}(R*)$ intensity distribution calculated for pitches of MiSolFA grids, x*=0.2, 10keV")
    f.text(0.5, 0.04, 'r (mm)', ha='center')
    f.text(0.04, 0.5, 'G1d', va='center', rotation='vertical')
    plt.savefig('G1d_x02_10keV.png')
    f.show()
    
def misolfa_I_Rstar():
    '''Now with parameters from MiSolFA grids. Plot I(R*), and I(x) for a single grid -> G1d(x*)'''
    maxr=1000
    n=np.arange(0,maxr, dtype=float) #n for the sum
    xmm=n/maxr
    p=np.array([0.015,0.03,0.045,0.09,0.18,0.225]) #pitch of finest grids = .015 mm
    #xran=[0,8.8*10**-3] #moire period of 8.8mm
    R=500*n/maxr #154.7 #150-200 mm
    #mlambda = 1.239*10**-7  # 10 keV in mm
    #mlambda=6.199*10**-8 #20 keV to mm
    mlambda = 4.132*10**-8  # 30 keV in mm
    #mlambda = 1.239*10**-7 * np.arange((maxr)) # 100 keV in mm
    #mlambda[0]=1.239*10**-7 #to avoid divide by zero
    #Rstar = mlambda*R/p**2 # + delta#Fm + delta*F1 #R*mlambda/(p**2) #but I don't get square waves when R* = 1,2, ... like I should!
    xstar,Rstar,data=[],[],[]
    x=n/maxr
    
    for pp,i in zip(p,range(0,6)):
        Rstar.append(np.zeros((maxr)) +mlambda*R/pp**2)
        data.append(d.calc_functions(maxr,n,pp,Rstar[i],mlambda, 'Id(R)'))

    print R[309], R[416],data[5][309],data[5][416]
    #print np.shape(data)
    f,axes = plt.subplots(3,2, sharex=True,sharey=True)
    title='lambda=620 A (20 keV), R=154.7mm'
    for i,ax in enumerate(f.axes):
        print np.shape(data[i])
        ax.plot(R, data[i], '.r-',label='p='+str(p[i])+'mm') #here xstar is Rstar
        ax.set_xlim(0,500)
        ax.set_ylim(0,1)
        ax.legend(loc='lower right',fontsize='medium')
    plt.suptitle(r"Moir$\'{e}$ Intensity, 30 keV. Maximum "+ str(data[5][309]) + " (154 mm) and " +str(data[5][416]) + " (208 mm)")
    f.text(0.5, 0.04, 'R (mm)', ha='center')
    f.text(0.04, 0.5, 'Intensity', va='center', rotation='vertical')
    plt.savefig('Id(R)_30keV.png')
    f.show()
        
def plot_pitch_contrast():
    '''docstring'''
    p=np.array([0.015,0.03,0.045,0.09,0.18,0.225]) #mm. Let's us mm for everything!
    data=[]
    maxr=1000
    n=np.arange(0,maxr, dtype=float) #n for the sum
    R=154.7 #150-200 mm
    mlambda = 6.206*10**-8#6.206^-11 m median wavelength in 2-200 keV
    Rstar = (np.zeros((6)) +mlambda*R)/p**2
    print p,Rstar
    xran=[np.min(Rstar),np.max(Rstar)]
    for pp in p:
        data.append(d.calc_functions(maxr,n,pp,Rstar,mlambda, 'Id(R)'))
    #data=data*10**10
    yran=[0,np.max(data)]
    print np.shape(p),np.shape(data)
    print yran
    #plot_tests(p,data,xran,yran,'Contrast of Moire image, R=154.7mm','Id')
    
    f = plt.figure()
    ax1=f.add_subplot(111)
    title='lambda=620 A (20 keV), R=154.7mm'
    ax1.plot(p, data[0], '.r-') #here xstar is Rstar
    ax1.set_title(r'$\lambda=620 \AA$ (20 keV), R=154.7mm')
    plt.xlabel('pitch (mm)')
    plt.ylabel(r"Moir$\'{e}$ Intensity")
    ax1.set_ylim(yran[0],yran[1])
    #ax1.set_xlim(xran[0],xran[1])
    ax1.legend(loc='upper right',fontsize='medium')

    f.show()

def plot_energy_intensity(p,R):
    '''intensity as a function of energy for the pitches of the Misolfa grids. length units are mm'''
    #p=np.array([0.015,0.03,0.045,0.09,0.18,0.225]) #mm. Let's us mm for everything!
    p=np.array([0.0200,0.0240,0.0300,0.0400,0.0600,0.1200])

    data=[]
    maxr=1000
    n=np.arange(0,maxr, dtype=float) #n for the sum
    R = 206 #R=154.7 #150-200 mm
    mlambda = 1.239*10**-8 * np.arange((maxr)) # 100 keV in mm
    mlambda[0]=1.239*10**-8 #to avoid divide by zero
    EkeV = (sc.h*sc.c)/(mlambda*sc.e) #E=hc/lambda, E(keV) = (E/E_electron)/1000 (we're already in mm so lose the 1000)
    xran=[np.min(EkeV),np.max(EkeV)]
    for pp in p:
        Rstar = (np.zeros((maxr)) +mlambda*R)/pp**2
        #EkeV = np.arange(1,1001,dtype=float)/10
        #Rstar = sc.h*sc.c*R/(EkeV*pp**2)
        data.append(calc_functions(maxr,n,pp,Rstar,mlambda, 'C'))
    print np.shape(data),np.shape(data[0]),np.mean(data[0]),np.max(data[0])
    yran=[0,np.max(data)]
    #print np.shape(p),np.shape(data)
    print 'energy range', EkeV[0],EkeV[maxr-1],np.min(Rstar),np.max(Rstar)

    #interpolate to make axis uniform in energy? or do that earlier?
    
    f,axes = plt.subplots(3,2, sharex=True,sharey=True)
    title='lambda=620 A (20 keV), R=154.7mm'
    for i,ax in enumerate(f.axes):
        ax.plot(EkeV, data[i], '.r-',label='p='+str(p[i])+'mm') #here xstar is Rstar
        ax.set_xlim(0,50)
        ax.set_ylim(0,1)
        ax.legend(loc='lower right',fontsize='medium')
    plt.suptitle(r"Moir$\'{e}$ Intensity calculated for pitches of MiSolFA grids, R=206mm")
    f.text(0.5, 0.04, 'E (keV)', ha='center')
    f.text(0.04, 0.5, 'Intensity', va='center', rotation='vertical')

    f.show()

#p=np.array([0.015,0.03,0.045,0.09,0.18,0.225])
#plot_energy_intensity(p,154.7)
#p=np.array([0.0200,0.0240,0.0300,0.0400,0.0600,0.1200])
#plot_energy_intensity(p,206)

