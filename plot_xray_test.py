"""
===================
sim_det_fit.py
Erica  Lastufka 9.3.18
===================

Given nxm pixel detector, degrade the moire image (or 'bin') then fit sines or triangle functions to calculate the moire period and orientation

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm




def calc_hex(radius):
    '''calculate points on a hexagon about (0,0) of a given radius'''
    xvec=[0]
    yvec=[0]
    dtheta=60.
    for p in range(0,6):
        dtc=-p*dtheta
        xinc=radius*np.sin(np.radians(dtc))
        yinc=radius*np.cos(np.radians(dtc))
        xvec.append(xinc)
        yvec.append(yinc)
    return xvec,yvec

#xvec=[0,3.666666667,5.008759814,1.342093147,-3.666666667,-5.008759814,-1.342093147]
#yvec=[0,3.666666667,-1.342093147,-5.008759814,-3.666666667,1.342093147,5.008759814]
xvec,yvec=calc_hex(3)
print(xvec)
print(yvec)

fig,ax=plt.subplots()

r=2./2.
circle1 = plt.Circle((0, 0), r, color='b')
circle2 = plt.Circle((xvec[1], yvec[1]), r, color='b')
circle3 = plt.Circle((xvec[2], yvec[2]), r, color='b')
circle4 = plt.Circle((xvec[3], yvec[3]), r, color='b')
circle5 = plt.Circle((xvec[4], yvec[4]), r, color='b')
circle6 = plt.Circle((xvec[5], yvec[5]), r, color='b')
circle7 = plt.Circle((xvec[6], yvec[6]), r, color='b')
ax.add_artist(circle1)
ax.add_artist(circle2)
ax.add_artist(circle3)
ax.add_artist(circle4)
ax.add_artist(circle5)
ax.add_artist(circle6)
ax.add_artist(circle7)
ax.text(xvec[0],yvec[0],'p0',fontdict={'color':'white','size':16},horizontalalignment='center',verticalalignment='center')
ax.text(xvec[1],yvec[1],'p1',fontdict={'color':'white','size':16},horizontalalignment='center',verticalalignment='center')
ax.text(xvec[2],yvec[2],'p2',fontdict={'color':'white','size':16},horizontalalignment='center',verticalalignment='center')
ax.text(xvec[3],yvec[3],'p3',fontdict={'color':'white','size':16},horizontalalignment='center',verticalalignment='center')
ax.text(xvec[4],yvec[4],'p4',fontdict={'color':'white','size':16},horizontalalignment='center',verticalalignment='center')
ax.text(xvec[5],yvec[5],'p5',fontdict={'color':'white','size':16},horizontalalignment='center',verticalalignment='center')
ax.text(xvec[6],yvec[6],'p6',fontdict={'color':'white','size':16},horizontalalignment='center',verticalalignment='center')
ax.set_xlim([-5.5,5.5])
ax.set_ylim([-5.5,5.5])
ax.set_xlabel('x (mm)')
ax.set_ylabel('y (mm)')
fig.show()
