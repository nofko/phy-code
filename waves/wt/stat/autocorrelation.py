import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch
from scipy.io import loadmat
from scipy import signal

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)

## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                        CALCULATION OF AUTOCORRELATION                     ##
##                                                                           ##
## ============================================================================



### LOAD SIGNAL ###

path = "/data/torus_mi/10_5_2022/cascade/bruitCTP_6_8hz_15min_A1100.mat"
path = "/data/torus_mi/20_6_2022/7_4hz_60min_A_1000mV.mat"
path = "/data/torus_mi/21_6_2022/7_7hz_60min_A_1000mV.mat"

#path = "/data/torus_mi/27_7_2022/cascade_amp/7_7Hz_15min_A1500.mat"
#path = "/data/torus_mi/28_7_2022/7_5Hz_40min_A900.mat"

T = 60*30
fs = 2000
nop = T*fs

norm = 14.5/6*1e-3  ## mm/V

S = loadmat(path)
S = S["signal"]
S = S[:,0]*norm

print("#############################################")
print("-- FILENAME : ", path.split("/")[-1])

plt.plot(S,lw=0.2)

### PDF FUNCTIONS ###

def autocorr(x):

    A = np.correlate(x, x, mode='full')    
    A /= A[A.argmax()]  # Normalize autocorrelation

    return A


ac = autocorr(S[2400000:2600000])

### PLOTTING ###


fig, ax = plt.subplots(figsize=[12,9])

ax.plot(ac[int(len(ac)/2):])

ax.set_xlabel(r"$t$",fontsize=30)
ax.set_ylabel(r"$A.C$",fontsize=30)

ax.tick_params(labelsize=25)


plt.show()
