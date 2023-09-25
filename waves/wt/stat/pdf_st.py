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


path = "/data/torus_wt/09_23/22_09/Basler_acA2040-120um__23597830__20230922_103544711.npy"


fps = 100

signal = np.load(path,allow_pickle=True)
signal = np.roll(signal,900)
signal = signal[:, :-270]

S = signal[:,500]
print(len(signal))
n = 1000
# plt.imshow(signal[8000:9000],aspect=len(signal)/len(signal[0]))


print("#############################################")
print("-- FILENAME : ", path.split("/")[-1])

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

hs, be = pdf(signal, 100)

# np.save("data/hs_5.npy",hs)
# np.save("data/bin_5.npy",be)


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
