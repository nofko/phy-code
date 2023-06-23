import matplotlib.pyplot as plt
import numpy as np

from os import listdir
from os.path import isfile, join
from scipy.io import loadmat

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)

## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                  GROWTH OF PARTICULAR FOURIER COMPONENT                   ##
##                                                                           ##
## ============================================================================


### LOAD SIGNAL ###

path = "/data/torus_mi/27_4_2022/6Hz/6hz_2min_A_800mV.mat"

T = 2*60
fs = 2000
nop = T*fs

S = loadmat(path)
S = S["signal"]
S = S[:,0]

print("#############################################")
print("-- FILENAME : ", path.split("/")[-1])


### SPECTROGRAM ###

print("-- CALCULATING SPECTROGRAM --")

plt.figure(figsize=[8,6])

SP,F,T,IM = plt.specgram(S,Fs=2000,NFFT=3*fs)

plt.title(path.split("/")[-1])
    
plt.ylim([0,10])
plt.ylabel("$f$ [Hz]")
plt.xlabel("$t$ [s]")


### FOURIER COMPONENT ###

f = 5

plt.figure(figsize=[8,6])

plt.plot(T,SP[int(f*2),:],"o")

plt.ylabel("$A_f$ [a.u.]")
plt.xlabel("$t$ [s]")

plt.show()

np.save("data/fourier_f_"+str(f)+"Hz_"+path.split("/")[-1].split(".mat")[0],SP[int(f*2),:])
np.save("data/T_f_"+str(f)+"Hz_"+path.split("/")[-1].split(".mat")[0],T)
