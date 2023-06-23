import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch
from scipy.io import loadmat

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)

## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##             CALCULATION OF PWELCH FOR A TEMPORAL SIGNAL                   ##
##                                                                           ##
## ============================================================================



### LOAD SIGNAL ###

path = "/data/torus_mi/28_7_2022/7_7Hz_10min"

T = 10*60
fs = 2000
nop = T*fs

norm = 14.5/6*1e-3  ## mm/V

S = loadmat(path)
S = S["signal"]
S = S[:,0]*norm


print("#############################################")
print("-- FILENAME : ", path.split("/")[-1])



### PLOTTING ###

plt.figure()
plt.plot(S,lw=1)

fig, ax = plt.subplots(figsize=[13,10])
    
ax.loglog(Fxx,Pxx,lw=1,color="b")

ax.set_xlabel("$f$ [Hz]",fontsize=35)
ax.set_ylabel(r"$S_\eta(\omega)$ [m$^2$s]",fontsize=35)
ax.tick_params(labelsize=25)

plt.show()


