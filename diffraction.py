
#import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import scipy.constants as sc

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



