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
##                             CALCULATION OF PDF                            ##
##                                                                           ##
## ============================================================================



### LOAD SIGNAL ###

path = "/data/torus_mi/10_5_2022/cascade/bruitCTP_6_8hz_15min_A1100.mat"
path = "/data/torus_mi/20_6_2022/7_4hz_60min_A_1000mV.mat"

#path = "/data/torus_mi/27_7_2022/cascade_amp/7_7Hz_15min_A1500.mat"
path = "/data/torus_wt/02_23/21_02/bruitCTP_6_9hz_15min_A800.mat"

T = 60*15
fs = 2000
nop = T*fs

norm = 14.5/6*1e-3  ## mm/V

S = loadmat(path)
S = S["signal"]
S = S[:,0]*norm
#S = S[120000:]
print("#############################################")
print("-- FILENAME : ", path.split("/")[-1])

plt.plot(S)
### PDF FUNCTIONS ###

def flatness(X):

    return (np.mean((X-np.mean(X))**4) )/(np.mean((X-np.mean(X))**2)**2 )

def skewness(X):

    return (np.mean((X-np.mean(X))**3) )/(np.mean((X-np.mean(X))**2)**(3/2) )


def pdf(X,nbin):

    maximum = np.max((X-np.mean(X))/np.std(X))
    minimum = np.min((X-np.mean(X))/np.std(X))

    
    bins = np.linspace(minimum,maximum,nbin)
    
    hist,b = np.histogram((X-np.mean(X))/np.std(X),bins)

    hist = hist / ( np.sum(hist)*np.mean(np.diff(b)))

    return hist,b


### PLOTTING ###

xg = np.linspace(-4,4,1000)
xp = np.linspace(-5,0,1000)

fig, ax = plt.subplots(figsize=[12,9])

hs, be = pdf(S, 200)

ax.semilogy(be[:-1], hs)

ax.semilogy(xg,np.exp(-(xg**2)/2)/np.sqrt(2*np.pi),'-.k')
#ax.semilogy(xp,np.exp(1.5*xp))

ax.set_xlabel(r"$\eta/\sigma$",fontsize=30)
ax.set_ylabel(r"$p(\eta/\sigma)$",fontsize=30)

ax.tick_params(labelsize=25)

### PLOTTING OF FILTERED###

# sos = signal.butter(4,60,'low',fs=int(fs),output="sos")
# Sf = signal.sosfilt(sos,S)

# hsf, bef = pdf(Sf[1000:], 200)
# #fig, ax = plt.subplots(figsize=[12,9])
# ax.semilogy(bef[:-1], hsf)

# #ax.semilogy(xg,np.exp(-(xg**2)/2)/np.sqrt(2*np.pi),'-.k')
# #ax.semilogy(xp,np.exp(1.5*xp))

# ax.set_xlabel(r"$\eta/\sigma$",fontsize=30)
# ax.set_ylabel(r"$p(\eta/\sigma)$",fontsize=30)

# ax.tick_params(labelsize=25)

#plt.savefig("images/pdf_small.pdf",bbox_inches="tight")

plt.show()
