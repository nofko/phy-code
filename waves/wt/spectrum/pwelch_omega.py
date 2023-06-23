import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch,argrelextrema
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

path = "/data/torus_mi/10_5_2022/cascade/bruitCTP_6_8hz_15min_A1100.mat"

path = "/data/torus_mi/20_6_2022/7_4hz_60min_A_1000mV.mat"
path = "/data/torus_mi/21_6_2022/7_7hz_60min_A_1000mV.mat"

path = "/data/torus_mi/27_7_2022/7_5Hz_3min_A400.mat"

T = 10*60
fs = 2000
nop = T*fs

norm = 14.5/6*1e-3  ## mm/V

S = loadmat(path)
print(S)
S = S["signal"]
S = S[:,0]*norm


print("#############################################")
print("-- FILENAME : ", path.split("/")[-1])


### PWELCH ###
    
Fxx, Pxx = welch(S, fs, window="blackman",noverlap = fs/2, nfft=4*fs, 
               axis=0, scaling="density", detrend=False, nperseg=4*fs)



### PLOTTING ###

plt.figure()
plt.plot(S,lw=1)

fig, ax = plt.subplots(figsize=[13,10])
    
ax.loglog(Fxx,Pxx,lw=1,color="b")

ax.set_xlabel("$f$ [Hz]",fontsize=35)
ax.set_ylabel(r"$S_\eta(\omega)$ [m$^2$s]",fontsize=35)
ax.tick_params(labelsize=25)

### FIT ###

x = np.linspace(25,190)

slope = -3.2

#ax.loglog(x,x**(slope)*0.9e-6,"--r")

ax.text(50,1e-11,r"$f^{"+str(slope)+"}$",fontsize=30,color="r")

ax.set_xlim([4,1e3])

ax.axvspan(6, 8, alpha=0.5, color='gray')
ax.text(6.3, 8e-13, "Forcing", rotation=90, va='center',fontsize=25)


print("-- SLOPE FIT : ", slope)

plt.savefig("images/pwelch_"+path.split("/")[-1].split(".mat")[0]+".pdf",bbox_inches="tight")

plt.show()


