import matplotlib.pyplot as plt
import numpy as np

from os import listdir
from os.path import isfile, join
from scipy.io import loadmat
import scipy.signal as signal


plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                         AMPLITUDE VS FREQUENCY                            ##
##                                                                           ##
## ============================================================================



### LOAD SIGNALS ###


path = "/data/pendulum/10_20_2022/"

f1 = '20hz_A1500.mat'
f2 = '20hz_A1200.mat'
f3 = '20hz_A800.mat'
f4 = '20hz_A400.mat'

s1 = loadmat(path+f1)
s2 = loadmat(path+f2)
s3 = loadmat(path+f3)
s4 = loadmat(path+f4)

s1 = s1[next(reversed(s1))].transpose()[0]
s2 = s2[next(reversed(s2))].transpose()[0]
s3 = s3[next(reversed(s3))].transpose()[0]
s4 = s4[next(reversed(s4))].transpose()[0]

s1 = np.convolve(s1, np.ones(10), "valid") / 10
s2 = np.convolve(s2, np.ones(10), "valid") / 10
s3 = np.convolve(s3, np.ones(20), "valid") / 20
s4 = np.convolve(s4, np.ones(20), "valid") / 20

s2 = s2[:20000]
s3 = s3[:20000]
s4 = s4[:20000]

fs = 2000

T = 240

plt.plot(s2)

fig, axs = plt.subplots(nrows=2, ncols=2,figsize=[12,9])

n = 1
axs[0,0].plot(s1[:-1],np.diff(s1),lw=0.5)
axs[0,1].plot(s2[:-1],np.diff(s2),lw=0.5)
axs[1,0].plot(s3[:-1],np.diff(s3),lw=0.5)
axs[1,1].plot(s4[:-n],np.diff(s4,n),lw=0.5)


plt.figure()

Fxx, Pxx = signal.welch(s1, fs, window="blackman",noverlap = fs/2, nfft=4*fs, 
               axis=0, scaling="density", detrend=False, nperseg=4*fs)


plt.loglog(Fxx,Pxx,lw=1,color="b")

plt.show()
