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

num_files = 5

# Initialize empty dictionaries to store the loaded data
h_dict = {}
b_dict = {}

# Loop through the files and load the data
for i in range(1, num_files + 1):
    h_var_name = f"h{i}"
    b_var_name = f"b{i}"

    h_dict[h_var_name] = np.load(f"data/hs_{i}.npy", allow_pickle=True)
    b_dict[b_var_name] = np.load(f"data/bin_{i}.npy", allow_pickle=True)

    

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
xp = np.linspace(2.5,6,1000)

fig, ax = plt.subplots(figsize=[12,9])

#hs, be = pdf(signal, 100)

for i in range(5, num_files + 1):
    h_var_name = f"h{i}"
    b_var_name = f"b{i}"
    
    h_data = h_dict[h_var_name]
    b_data = b_dict[b_var_name]
    
    # Plot the data
    ax.semilogy(b_data[:-1], h_data,)

ax.semilogy(xg,np.exp(-(xg**2)/2)/np.sqrt(2*np.pi),'-.k')
ax.semilogy(xp,np.exp(-1.2*xp)*2)

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
