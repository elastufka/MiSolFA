import numpy as np
import q
@q

def flare_XR_dist():
    '''Returns the typical distribution of counts for a flare in the energy range ... keV'''
    import pidly #set up pidly 
    idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
    
    eedd=np.linspace(3,300.5,595) #findgen(600-5)/2.+3.

    #thermal emission for MIXI proposal flare
    this_em=0.07 #(7*10**47)/(1*10**49)          #in units of 10d49
    this_te=1.72414 #20./11.6           #;in keV   1 keV = 11.6 M
    #to get thermal spectrum in photons/s/cm2/keV
    #thth =idl.function('f_vth',eedd,[this_em,this_te])
    thth =idl.f_vth(eedd,[this_em,this_te])

    #nonthermal spectrum in photons/s/cm2/keV

    this_norm=3         #norm (value of lower power law at 50 keV)
    this_below=1.7      #power law index below break
    this_break=13.5     #break energy in keV
    this_gamma=3.7      #power law index above break

    
    ntnt=idl.f_bpow(eedd,[this_norm,this_below,this_break,this_gamma])
    idl.close()
    dist= ntnt[1:]+thth
    
    return eedd[1:],dist,ntnt,thth

def plot_flare_XR_dist():
    dist=flare_XR_dist()
    eedd = dist[0]
    total=dist[1]
    thth=dist[2]
    ntnt=dist[3]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.loglog(eedd, thth[1:], color="r",label="Thermal emission")
    ax1.loglog(eedd, ntnt, color="b",label="Non-thermal emission")
    ax1.loglog(eedd, total, color="g",label="Total emission")
    plt.xlabel('Energy (keV)')
    plt.ylabel('Photons (s^-1 cm^-2 keV^-1)')
    ax1.set_ylim([10e-4,10e4])
    ax1.set_xlim([0,100])
    #plt.title("Effective area of grids")
    ax1.legend(loc='upper right',fontsize='medium')

    ax1.plot()

    fig.show()

if __name__ == '__main__':
    plot_flare_XR_dist()
