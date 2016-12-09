
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import scipy.constants as sc

def calc_functions(maxr,n,p,Rstar,mlambda, function):
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

    if function == 'Id' : #moire image. xstar is now xi
        for xx in xstar:
            sum_Asubn = np.append(sum_Asubn, np.sum(np.cos(2*np.pi*xx*(2*n+1))*np.cos(np.pi*Rstar*(2*n+1)**2)/(2*n+1)**2))
        Id  = .25 + (2/np.pi**2)*sum_Asubn
        var=Id
    if function == 'C': #contrast - eq. 35. Now I is a function of R, not xi
        #C(R*)=[Id'(R*)-Id'(R*=1/2)]/Id'(R*=1/2)]
        xstar = 0.0

        for rr in Rstar:
           sum_Asubn = np.append(sum_Asubn,np.sum(np.cos(2*np.pi*xstar*(2*n+1))*np.cos(np.pi*rr*(2*n+1)**2)/(2*n+1))) #since xstar=0,simplify:
           #sum_Asubn = np.append(sum_Asubn,np.sum(np.cos(np.pi*rr*(2*n+1)**2)/(2*n+1)**2))
           
        Id_Rhalf = np.zeros((len(Rstar)))+ .25 + (2/np.pi**2)*np.sum(np.cos(2*np.pi*xstar*(2*n+1))*np.cos(np.pi*.5*(2*n+1)**2)/(2*n+1))
        #since xi=0, simplify:
        #Id_Rhalf = np.zeros((maxr))+ .25 + (2/np.pi**2)*np.sum(np.cos(np.pi*.5*(2*n+1)**2)/(2*n+1)**2)
        #well at least we verified that the results are the same. 
        Id = .25 + (2/np.pi**2)*sum_Asubn
        #print np.shape(sum_Asubn),np.shape(Id),np.shape(Id_Rhalf)
        C=np.abs(Id-Id_Rhalf)/Id_Rhalf #what they plot is the absolute value, but then they fail to say that. Unless [] means absolute value in Israel?
        var=Id#[C,Id,Id_Rhalf]
        #print Id_Rhalf[0],Id[0],Id[maxr-1], np.sum(np.cos(np.pi*(2*n+1)**2)/(2*n+1)**2)
    return var

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

def make_plots():
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
        data.append(calc_functions(maxr,n,p,rr,mlambda, 'G1d'))
        data2.append(calc_functions(maxr,n,p,rr,mlambda, 'Id'))

    x=np.arange(0,maxr,dtype=float)/(maxr) #this is still in meters right?
    xstar = x/p
    #data=calc_functions(maxr,n,p,R,mlambda, 'G1d')
    #plot_tests(xstar,data) #reproduce fig 2A
    #plot_tests(xstar,data2)

    #plot the contrast now
    Rstar = np.arange(0,maxr,dtype=float)/maxr
    Rstar[maxr-1]=1

    data3= calc_functions(maxr,n,p,Rstar,mlambda, 'C')
    print np.shape(data3[2]),np.shape(Rstar)
    #yran=[np.min(data3),np.max(data3)]
    plot_tests(Rstar, data3,[0,1],[0,1], 'Contrast of Moire image', 'C')

def make_misolfa_plots():
    '''Now with parameters from MiSolFA grids'''
    maxr=1000
    n=np.arange(0,maxr, dtype=float) #n for the sum
    p= 15*10**-6 #pitch of finest grids = .015 mm
    xran=[0,8.8*10**-3] #moire period of 8.8mm
    R=.1547 #150-200 mm
    mlambda = 6.206*10**-11#6.206^-12 median wavelength in 2-200 keV
    delta = np.array([0.0,0.01,0.02,0.03,0.04]) #see eq. 33 and 34
    #F1=np.array([p**2,p**2,p**2,p**2,p**2])/mlambda 
    #m=[2,2,2,2,2]
    #Fm = m*F1
    Rstar = mlambda*R/p**2 + delta#Fm + delta*F1 #R*mlambda/(p**2) #but I don't get square waves when R* = 1,2, ... like I should!
    print Rstar
    data,data2=[],[]
    for rr in Rstar:
        data.append(calc_functions(maxr,n,p,rr,mlambda, 'G1d'))
        data2.append(calc_functions(maxr,n,p,rr,mlambda, 'Id'))

    x=np.arange(0,maxr,dtype=float)/(1.5*10**9) #this is still in meters right? so let's use nano
    xstar = x/p
    #print np.min(data),np.max(data),np.min(data2),np.max(data2),np.max(x),np.max(xstar),xran[1]
    yran=[np.min(data2),np.max(data2)]
    #data=calc_functions(maxr,n,p,R,mlambda, 'G1d')
    #plot_tests(xstar,data) #reproduce fig 2A
    plot_tests(xstar,data2, [0,p*10],yran,'Moire intensity for MiSolFA grids','Id')

    Rstar =Rstar[0]/2+ np.arange(0,maxr,dtype=float)/maxr
    print np.min(Rstar),np.max(Rstar)

    data3= calc_functions(maxr,n,p,Rstar,mlambda, 'C')
    yran=[np.min(data3[0]),np.max(data3[0])]
    print yran
    plot_tests(Rstar, data3,[0,Rstar[0]],yran, 'Contrast of Moire image', 'C')

def plot_pitch_contrast():
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
        data.append(calc_functions(maxr,n,pp,Rstar,mlambda, 'C'))
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

def plot_energy_intensity():
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
        data.append(calc_functions(maxr,n,pp,Rstar,mlambda, 'C'))
    print np.shape(data),np.shape(data[0]),np.mean(data[0]),np.max(data[0])
    yran=[0,np.max(data)]
    #print np.shape(p),np.shape(data)
    print 'energy range', EkeV[0],EkeV[maxr-1]
    #plot_tests(p,data,xran,yran,'Contrast of Moire image, R=154.7mm','Id')
    
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
