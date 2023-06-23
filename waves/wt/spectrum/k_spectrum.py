import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch
from scipy.io import loadmat
from scipy.signal import get_window


plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)

## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##         CALCULATION OF PWELCH FOR SPACETIME SIGNAL - WAVENUMBER           ##
##                                                                           ##
## ============================================================================


### LOAD SIGNAL ###

path = "/data/torus_mi/5_5_2022/bruitCTP_5_7hz_15min_A1100.npy"
path = "/data/torus_mi/29_9_2022/7_9hz/7_9Hz_4min_600.npy"
path = "/data/torus_wt/12_22/9_12/bruitCTP_6_9Hz_A1100.npy"
path = "/data/torus_wt/02_23/17_02/bruitCTP_2_5hz_7min_A1200.npy"

fps = 100

S = np.load(path)

# S = np.roll(S,-430)
# S = S[:, 250:]

# S = np.roll(S,-400)
# S = S[:, 500:]
# S = np.pad(S,[(250,250),(0,0)])

S = np.roll(S,-850)
S = S[:,200:]

wk = get_window("blackman", len(S[0]))
wo = get_window("blackman", len(S))
pn = 500
#signal = np.pad(signal,[(pn,pn),(pn,pn)])
signal = S*np.outer(wo,wk)

print("#############################################")
print("-- FILENAME : ", path.split("/")[-1])


### SPECTRUM ###

print("-- CALCULATING FOURIER --")

FF = np.fft.fft2(S)

print("-- CALCULATING SPECTRUM --")

S_k = np.sum(np.abs(FF)**2,axis=0)

kmax = 500
    
k = np.arange(1,kmax+1)


### PLOTTING ###

print("-- DONE!")

fig, ax = plt.subplots(figsize=[12,9])
    
ax.loglog(k,S_k[:kmax],lw=1,color="b")

x = np.linspace(10,80)

slope = -3.1

ax.loglog(x,x**(slope)*2e19,"--k")

ax.set_xlabel(r"$k_\theta$",fontsize=35)
ax.set_ylabel(r"$S_\eta(k_\theta)$ [a.u.]",fontsize=35)
ax.tick_params(labelsize=25)

plt.savefig("images/pwelch_k_"+path.split("/")[-1].split(".mat")[0]+".pdf",bbox_inches="tight",dpi=500)

plt.show()
