import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import get_window

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)

## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                  CALCULATION OF SPACE-TIME FOURIER SPECTRUM               ##
##                                                                           ##
## ============================================================================



### DISPERSION RELATIONS ###


def varicose(k,w,ro):

    return np.sqrt((g*k/ro+sigma*k**3/rho/ro**3)*np.tanh(k*w*ro/2/rc**2))

def varicose_N(k,w,ro,N):

    return N*np.sqrt((g*(k/N)/ro+sigma*(k/N)**3/rho/ro**3)*np.tanh((k/N)*w*ro/2/rc**2))

def sloshing(k,om0):

    #return om0+10+varicose(k,w,ro)
    #return np.sqrt((om0**2+g*k**2/(0.076)))
    #return om0+np.sqrt((g*k/rc+sigma*k**3/rho/ro**3)*np.tanh(k*w*rc/2/rc**2))
    return np.sqrt((om0**2+(g*k**2/(ro)+0.5*sigma*k**3/rho/ro**3)))

def sinuous(k,w,ri):

    om0 = np.pi/2*np.sqrt(9.81*np.sin(9*np.pi/180)/2/w)

    return np.sqrt(om0**2+np.abs(sigma/rho*(k/ri)**3))

g = np.sin(4.5*np.pi/180)*9.81

sigma = 0.05
rho = 1000

rc = 0.095
ro = 0.08
w = (rc-ro)
ri =ro-w

om0 = 7.14*2*np.pi#44.5


### LOAD SIGNAL ###

path = "/data/torus_mi/29_9_2022/7_9hz/7_9Hz_4min_750.npy"
path = "/data/torus_mi/10_5_2022/sweep/Basler_acA2040-120um__23597830__20220510_122929533.npy"
path = "/data/torus_wt/09_23/22_09/Basler_acA2040-120um__23597830__20230922_095128310.npy"


fps = 100

signal = np.load(path,allow_pickle=True)
signal = np.roll(signal,900)
signal = signal[:, 120:-300]

w = 15


#np.convolve(signal, np.ones(w), 'valid',axis=0) / w

signal = np.apply_along_axis(lambda m: np.convolve(m, np.ones(w), 'valid')/w, axis=1, arr=signal)

plt.plot(signal[1000])

#signal = np.pad(signal,[(250,250),(0,0)])

#signal = np.array([np.convolve(j, np.ones(200), "valid") / 200 for j in signal])

#signal = np.diff(signal,axis=1)

#plt.plot(signal[500])


wk = get_window("blackman", len(signal[0]))
wo = get_window("blackman", len(signal))
pn = 500
#signal = np.pad(si1gnal,[(pn,pn),(pn,pn)])
signal = signal*np.outer(wo,wk)#np.multiply(signal,w)
NS = len(signal)

#signal = np.abs(signal)
# signal = np.where(signal<0,signal,0)
### SPECTRUM ###

ff2 = np.fft.fft2(signal[:])
#ff2 = ff2[::-1]
del(signal)
#ff2 = np.fft.fftshift(ff2)

numFrames=len(ff2)

kmax = 180
omax = int((NS-1)/2)+6000

om = np.arange(0,omax)*fps/numFrames
k = (np.arange(0,100))

#plt.figure()
#plt.plot(om,np.sum(np.log(np.abs(ff2[:omax,:kmax])),axis=1))

### PLOTTING ###

fig, ax = plt.subplots(figsize=[12,9])


#four = ax.pcolormesh(k,om/2/np.pi,np.log(np.abs(ff2[:omax,:kmax])),rasterized=True,cmap="turbo",vmin=13,vmax=19)
four = ax.imshow(np.log(np.abs(ff2[:omax,:kmax]))[::-1,:],
                 cmap="turbo",extent=[k[0],kmax,om[0],om[-1]],
                 aspect="auto",vmin=8,vmax=17,interpolation="none")

#ax.plot(k,varicose(k,w,ro)/2/np.pi,"w",lw=2)

#delta = 5
#ax.plot(k,varicose(k,w,ro)+delta,"--w",lw=2)
#ax.plot(k,varicose(k,w,ro)-delta,"--w",lw=2)

# ax.plot(k,varicose_N(k,w,ro,2)/2/np.pi,"w",lw=2)
#ax.plot(k,varicose_N(k,w,ro,3),"w",lw=2)

cbar = plt.colorbar(four,format="%1.1f",ax=ax)
cbar.set_label(r"log$|\tilde{\eta}(k_\theta,\omega)|$",fontsize=30)
cbar.ax.tick_params(labelsize=30)

ax.set_xlabel(r"$k_\theta$",fontsize=35)
ax.set_ylabel("$f$ [Hz]",fontsize=35)

ax.tick_params(labelsize=30)

#plt.xlim([0,60])
#plt.ylim([0.0,45])

### VECTORS ###

# V = np.array([[11.33,2.9], [6.2,7.9], [24.7,7.9],[17.5,5.0]])
# origin = np.array([[0, 0, 0],[0, 0, 0]]) # origin point

# ax.quiver(*origin, V[:3,0], V[:3,1], color=['r','b','g'], angles='xy', scale_units='xy',
#           scale=1,width=0.004)

# v23 = [15.1,10.8] 
# v23c = [29,10.8] 

# #ax.quiver(*V[1], V[0,0], V[0,1], color=['r','b','g'], angles='xy', scale_units='xy',
# #          scale=1,linewidth=0.7)

# #ax.quiver(*V[1], -V[2,0], V[2,1], color=['r','b','g'], angles='xy', scale_units='xy',
# #          scale=1,linewidth=0.7)

# ax.quiver(*V[2], -V[0,0], V[0,1], color=['r'], angles='xy', scale_units='xy',
#           scale=1,width=0.004)

# #ax.quiver(*V[2], V[3,0], V[3,1], color=['r','b','g'], angles='xy', scale_units='xy',
# #          scale=1,linewidth=0.7)

# ax.quiver(*v23, V[2,0], V[2,1], color=['g'], angles='xy', scale_units='xy',
#           scale=1,width=0.004)

# ax.quiver(*v23, V[1,0], V[1,1], color=['b'], angles='xy', scale_units='xy',
#           scale=1,width=0.004)


# ax.axhline(7.9,linestyle='--',color='w')
# ax.axhline(10.8,linestyle='--',color='w')
# ax.axhline(15.8,linestyle='--',color='w')
# ax.axhline(18.7,linestyle='--',color='w')

# ax.text(66,8.85,r'$f_1$',fontsize=25,color='w')
# ax.text(66,11.7,r'$g_-^0$',fontsize=25,color='w')
# ax.text(66,16.7,r'$2f_1$',fontsize=25,color='w')
# ax.text(66,19.7,r'$h_-^0$',fontsize=25,color='w')

kt = np.linspace(0,150,1000)

#ax.plot(kt,sloshing(kt,om0)/2/np.pi,"k-.",linewidth=2.5)
#ax.plot(kt,varicose(kt,w,ro)/2/np.pi,"k--",linewidth=2.5)
#ax.plot(kt,varicose_N(kt,w,ro,2)/2/np.pi,"k--",linewidth=2.5)
#ax.plot(kt,sinuous(kt,w,ri)/2/np.pi,"w--",linewidth=2.5)

#plt.savefig("images/ff2_inv.pdf",bbox_inches="tight")
#plt.savefig("images/ff2_"+path.split("/")[-1].split(".npy")[0]+".pdf",bbox_inches="tight")

plt.show()
