import matplotlib.mlab as mlab
import math


data11=pickle.load(open('win11_width_data5.0.p','rb'))
data42=pickle.load(open('win42_width_data5.0.p','rb'))
bins11=np.linspace(88.5,91.1,30)
bins42=np.linspace(17,19,30)

fig,ax=plt.subplots(1,2,sharey=True)
foo11=ax[0].hist(data11,bins11,facecolor='b')
print np.max(foo11[0]),np.max(foo11[2])
ax[0].plot(bins11,320.*mlab.normpdf(bins11,np.mean(data11),np.std(data11)))
foo42=ax[1].hist(data42,bins42,facecolor='b')
ax[1].plot(bins42,np.max(foo42[0])*mlab.normpdf(bins42,np.mean(data42),np.std(data42)))
ax[0].set_xlim([85,95])
ax[1].set_xlim([15,21])
ax[0].set_xlabel('Period $\mu$m')
ax[0].set_ylabel('Counts')
ax[0].set_yscale('log')
#ax[0].set_ylim([10,10000])

fig.show()
