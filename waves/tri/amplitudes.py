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

#path = "/data/torus_mi/10_5_2022/sweep/7_4hz_4min_A_600mV.mat"
#path = "/data/torus_mi/31_5_2022/amplitude/7_4hz_4min_A_800mV.mat"
#path = "/data/torus_mi/27_7_2022/7_4Hz_60min_A1000.mat"
path = "/data/torus_mi/11_5_2022/cascade/7_9hz_10min_A_1100mV.mat"

T = 60*3
fs = 2000
nop = T*fs

norm = 14.5/6  ## mm/V

S = loadmat(path)
#S = S["signal"]
S = S[next(reversed(S))]
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

sos = signal.butter(3,(7.65,7.8),'bandpass',fs=int(fs),output="sos")
Sf_forcing = signal.sosfilt(sos,S)

Sf_forcing_pks, _ = signal.find_peaks(Sf_forcing)

sos = signal.butter(4,(4.7,4.85),'bandpass',fs=int(fs),output="sos")
Sf_high = signal.sosfilt(sos,S)

Sf_high_pks, _ = signal.find_peaks(Sf_high)

sos = signal.butter(4,(2.9,3.05),'bandpass',fs=int(fs),output="sos")
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

plt.figure()
plt.plot(Sf_forcing)

fig, ax = plt.subplots(figsize=[12,9])

ax.set_xlabel(r"$t$ [s]",fontsize=35)
ax.set_ylabel(r"$a$ [mm]",fontsize=35)
ax.tick_params(labelsize=25)


#crop = 12300

t = np.linspace(0,T,len(S))
#t = t[crop:]-t[crop]

#ax.plot(Sf_forcing,lw=2,color="b")
#ax.plot(t,S,lw=2,color="r")
ax.plot(Sf_forcing_pks/fs, Sf_forcing[Sf_forcing_pks], marker="o",color="b",ms=2)

ax.plot(Sf_low_pks/fs, Sf_low[Sf_low_pks], marker="o",color="r",)

ax.plot(Sf_high_pks/fs, Sf_high[Sf_high_pks], marker="o",color="g",)

#ax.set_xlim(-1,30)

plt.savefig("images/amplitudes_"+path.split("/")[-1].split(".mat")[0]+".pdf",bbox_inches="tight")

plt.show()


