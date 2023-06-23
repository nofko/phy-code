import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal

from os import listdir
from os.path import isfile, join
from scipy.io import loadmat

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)

## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                     PHASE SPECTRUM OF TEMPORAL SIGNAL                     ##
##                                                                           ##
## ============================================================================


### LOAD SIGNAL ###

path = "/data/torus_mi/10_5_2022/sweep/7_4hz_4min_A_600mV.mat"
#path = "/data/torus_mi/31_5_2022/amplitude/7_4hz_4min_A_650mV.mat"

T = 60*2
fs = 2000
nop = T*fs

S = loadmat(path)
S = S["signal"]
S = S[:,0]
S = S-np.mean(S)

# n1 = 100000
# n2 = 200000
# nop = n2-n1
# S = S[n1:n2]
# 

print("#############################################")
print("-- FILENAME : ", path.split("/")[-1])


### COSINE TEST ###

# t = np.linspace(0,T,2000*T)

# S = 0.5*np.cos(2*np.pi*10*t+np.pi/6)+np.cos(2*np.pi*5*t+np.pi/2)

### PHASE CALCULATION ###

# FF = np.fft.fftshift(np.fft.fft(S))

# threshold = max(abs(FF))/100

# FF[np.where(abs(FF)<threshold)] = 0

# PH = np.angle(FF)

# F = np.arange(-int(nop-1)/2,int(nop-1)/2+1)*fs/nop

### PLOTTING ###

# fig, ax = plt.subplots(figsize=[8,6])

# plt.plot(F,PH,"o")

# ax.set_xlim([0,30])
# ax.set_xlabel("$f$ [Hz]",fontsize=35)
# ax.set_ylabel(r"$\phi$ [rad]",fontsize=35)
# ax.tick_params(labelsize=25)



### BANDPASS ###

print("-- CALCULATING BANDPASS --")


sos = signal.butter(4,(7.33,7.53),'bandpass',fs=int(fs),output="sos")
Sf_forcing = signal.sosfilt(sos,S)

sos = signal.butter(4,(4.33,4.54),'bandpass',fs=int(fs),output="sos")
Sf_high = signal.sosfilt(sos,S)

sos = signal.butter(4,(2.93,3.13),'bandpass',fs=int(fs),output="sos")
Sf_low = signal.sosfilt(sos,S)

T = len(S)/fs


### HILBERT ###

print("-- PERFORMING HILBERT TRANSFORM --")

H1 = signal.hilbert(Sf_forcing)
H2 = signal.hilbert(Sf_high)
H3 = signal.hilbert(Sf_low)

phase_1 = np.unwrap(np.angle(H1))
phase_2 = np.unwrap(np.angle(H2))
phase_3 = np.unwrap(np.angle(H3))


### PLOTTING ###


print("-- DONE!")

fig, ax = plt.subplots(figsize=[12,9])

ax.set_xlabel(r"$t$ [s]",fontsize=35)
ax.set_ylabel(r"$\sin(\varphi_1+\varphi_2-\varphi_3)$",fontsize=35)
ax.tick_params(labelsize=25)


crop = 12300

t = np.linspace(0,T,len(S))
t = t[crop:]-t[crop]

ax.plot(t,np.sin(-phase_1+phase_2+phase_3)[crop:],lw=2,color="b")

plt.axhline(1,linestyle="-.",color="gray",lw=2)

#plt.axhline(0.5,linestyle="--",color="gray",lw=2)


plt.savefig("images/phase_lock_"+path.split("/")[-1].split(".mat")[0]+".pdf",bbox_inches="tight")

plt.show()
