#matplotlib = reload(matplotlib)
#import matplotlib
#matplotlib = reload(matplotlib)

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

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
            sum_Asubn = np.append(sum_Asubn, np.sum(np.cos(2*np.pi*xx*(2*n+1))*np.cos(np.pi*Rstar*(2*n+1)**2)/(2*n+1)))
        Id  = .25 + (2/np.pi**2)*sum_Asubn
        var=Id
    if function == 'C': #contrast - eq. 35. Now I is a function of R, not xi
        #C(R*)=[Id'(R*)-Id'(R*=1/2)]/Id'(R*=1/2)]
        xstar = 0.0

        for rr in Rstar:
           sum_Asubn = np.append(sum_Asubn,np.sum(np.cos(2*np.pi*xstar*(2*n+1))*np.cos(np.pi*rr*(2*n+1)**2)/(2*n+1)))
        Id_Rhalf = np.zeros((maxr))+ .25 + (2/np.pi**2)*np.sum(np.cos(2*np.pi*xstar*(2*n+1))*np.cos(np.pi*.5*(2*n+1)**2)/(2*n+1))
        Id = .25 + (2/np.pi**2)*sum_Asubn
        #print np.shape(sum_Asubn),np.shape(Id),np.shape(Id_Rhalf)
        C=(Id-Id_Rhalf)/Id_Rhalf
        var=[C,Id,Id_Rhalf]
        #print np.shape(var)
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

        ax1.plot(xstar, var[0], color='r', label=label) #here xstar is Rstar
        ax1.plot(xstar, var[1], color='g', label='Id') #here xstar is Rstar
        ax1.plot(xstar, var[2], color='b', label='Id_Rhalf') #here xstar is Rstar
    ax1.set_title(title)
    #f.subplots_adjust(hspace=0)
    #plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.xlabel('distance x*')
    plt.ylabel('Normalized Intensity')
    #ax1.set_ylim(yran[0],yran[1])
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
    maxr=100
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
    data3= calc_functions(maxr,n,p,Rstar,mlambda, 'C')
    print np.shape(data3[2]),np.shape(Rstar)
    #yran=[np.min(data3),np.max(data3)]
    plot_tests(Rstar, data3,[0,1],[0,1], 'Contrast of Moire image', 'C')

def make_misolfa_plots():
    '''Now with parameters from MiSolFA grids'''
    maxr=100
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
    print np.min(data),np.max(data),np.min(data2),np.max(data2),np.max(x),np.max(xstar),xran[1]
    yran=[np.min(data2),np.max(data2)]
    #data=calc_functions(maxr,n,p,R,mlambda, 'G1d')
    #plot_tests(xstar,data) #reproduce fig 2A
    plot_tests(xstar,data2, xran,yran,'Moire intensity for MiSolFA grids','Id')
