import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from scipy.optimize import fsolve
from scipy.signal import welch
from scipy.io import loadmat


## ============================================================================
## Calculation of four-wave tricoherence in the torus from experimental data ##
##                                                                           ##
## ============================================================================




### LOAD SIGNAL ###

print("##############################################")
print("################ TRICOHERENCE ################")
print("")

path = "/data/torus_mi/5_5_2022/cascade/6_1hz_15min_A_1100mV.mat"
path = "/data/torus_mi/21_6_2022/7_7hz_60min_A_1000mV.mat"

s = loadmat(path)
s = s["signal"]
s = s[:,0]

fs = 2000

print("-- FILENAME : ", path.split("/")[-1])


### PWELCH ###

f, Pxx = welch(s, fs, window="blackman",noverlap = fs/2, nfft=4*fs, 
               axis=0, scaling="density", detrend=False, nperseg=4*fs)


fig, ax = plt.subplots(figsize=[8,6])
ax.set_title(path.split("/")[-1])

ax.loglog(f,Pxx,lw=1)

ax.set_xlabel("$f$ [Hz]",fontsize=35)
ax.set_ylabel(r"$S_\eta(\omega)$ [a.u.]",fontsize=35)
ax.tick_params(labelsize=25)

#ax.set_xlim(1.7,1e3)

###  FOURIER CALCULATION ###

print("-- CALCULATING FOURIER --")

NW = 101                                                      # NUMBER OF TEMPORAL WINDOWS

LW = len(s)//NW + 1                                           # WINDOW LENGTH

if LW % 2 ==0:
    LW -= 1

SW = np.array([s[i:i+LW] for i in range(0,LW*(NW-1),LW)])     # THE MATRIX OF THE WINDOWED SIGNAL

FF = np.zeros_like(SW,dtype=complex)                         # FOURIER MATRIX

for i in range(NW-1):

    FF[i] = np.fft.fftshift(np.fft.fft(SW[i]))

F = np.arange(-int(LW-1)/2,int(LW-1)/2+1)*fs/LW 
F0 = int(LW/2)

### TRICOHERENCE CALCULATION ###

print("-- CALCULATING TRICOHERENCE, PLEASE HOLD --")

NTI = 2000

nti = int(NTI/2)

T = np.zeros((NTI,NTI))                                     # TRICOHERENCE MATRIX

F3 = 20                                                     # f_3 IN HZ
D3 = int(F3*LW/fs)

for i in range(-nti,nti):

    for j in range(-nti,nti):
    
        T[i+nti,j+nti] = np.abs(np.mean(np.conjugate(FF[:,F0+i]*FF[:,F0+j])*FF[:,F0+D3]*FF[:,F0+i+j-D3]))


### PLOTTING ###

print("-- DONE!")

fig, ax = plt.subplots(figsize=[8,6])

image = ax.imshow(np.log(T[::-1,:]),cmap="turbo",extent=[F[F0-nti],F[F0+nti],F[F0-nti],F[F0+nti]])

ax.set_xlabel("$\omega_1$",fontsize=35)
ax.set_ylabel(r"$\omega_2$ [a.u.]",fontsize=35)
ax.tick_params(labelsize=25)

plt.show()
input("Press ENTER")
