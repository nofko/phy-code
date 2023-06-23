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
##                          GROWTH OF THE AMPLITUDES                         ##
##                                                                           ##
## ============================================================================


### LOAD SIGNAL ###

path = "/data/torus_mi/10_5_2022/sweep/7_4hz_4min_A_600mV.mat"
path = "/data/torus_mi/31_5_2022/amplitude/7_4hz_4min_A_750mV.mat"

T = 60*15
fs = 2000
nop = T*fs

norm = 14.5/6  ## mm/V

S = loadmat(path)
S = S["signal"]
S = S[:,0]*norm
S = S-np.mean(S)

print("#############################################")
print("-- FILENAME : ", path.split("/")[-1])


### COSINE TEST ###

#t = np.linspace(0,T,2000*T)

#S = 0.5*np.cos(2*np.pi*10*t+np.pi/6)+np.cos(2*np.pi*5*t+np.pi/2)
#S[:10000] = 0

### BANDPASS ###

print("-- CALCULATING BANDPASS --")

sos = signal.butter(3,(7.3,7.5),'bandpass',fs=int(fs),output="sos")
Sf_forcing = signal.sosfilt(sos,S)

Sf_forcing_pks, _ = signal.find_peaks(Sf_forcing)

sos = signal.butter(4,(4.2,4.4),'bandpass',fs=int(fs),output="sos")
Sf_high = signal.sosfilt(sos,S)

Sf_high_pks, _ = signal.find_peaks(Sf_high)

sos = signal.butter(4,(3.06,3.26),'bandpass',fs=int(fs),output="sos")
Sf_low = signal.sosfilt(sos,S)

Sf_low_pks, _ = signal.find_peaks(Sf_low)

T = len(S)/fs

### AMPLITUDES ###

t_m = 60 ## TIME DELAY FOR THE MEAN IN S
n_m = int(t_m*len(Sf_forcing_pks)/T/2)

a1 = np.mean(Sf_forcing[Sf_forcing_pks][n_m:])
a2 = np.mean(Sf_low[Sf_low_pks][n_m:])
a3 = np.mean(Sf_high[Sf_high_pks][n_m:])

print("")
print(f"-- AMPLITUDES AFTER APPROX {t_m:n} s: ")
print(f"--      a1 = {a1:f}")
print(f"--      a2 = {a2:f}")
print(f"--      a3 = {a3:f}")

### PLOTTING ###


print("-- DONE!")

fig, ax = plt.subplots(2,figsize=[12,9])

ax[0].set_ylabel(r"log$(a)$ [mm]",fontsize=35)
ax[0].tick_params(labelsize=25)

ax[1].set_xlabel(r"$t$ [s]",fontsize=35)
ax[1].set_ylabel(r"$a$ [mm]",fontsize=35)
ax[1].tick_params(labelsize=25)

t = np.linspace(0,T,len(S))
t_fit = np.linspace(0,60)


ax[1].plot(Sf_forcing_pks/fs, Sf_forcing[Sf_forcing_pks], marker="o",color="b",)

ax[1].plot(Sf_low_pks/fs, Sf_low[Sf_low_pks], marker="o",color="r",)

ax[1].plot(Sf_high_pks/fs, Sf_high[Sf_high_pks], marker="o",color="g",)

ax[1].set_xlim(-1,50)


ax[0].plot(Sf_forcing_pks/fs, np.log(Sf_forcing[Sf_forcing_pks]), marker="o",color="b",)

ax[0].plot(Sf_low_pks/fs, np.log(Sf_low[Sf_low_pks]), marker="o",color="r",)

ax[0].plot(Sf_high_pks/fs, np.log(Sf_high[Sf_high_pks]), marker="o",color="g",)

ax[0].set_xlim(-1,50)

coords = plt.ginput(2)

x1,y1 = coords[0]
x2,y2 = coords[1]

slope = (y2-y1)/(x2-x1)
cut = y1-slope*x1

print("Slope: ",slope)

ax[0].plot(t_fit,t_fit*slope+cut)

plt.show()


