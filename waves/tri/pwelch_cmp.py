import matplotlib.pyplot as plt
import numpy as np

from os import listdir
from os.path import isfile, join
from scipy.io import loadmat
from scipy.signal import welch


plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                          COMPARING PWELCH SPECTRA                         ##
##                                                                           ##
## ============================================================================



### LOAD SIGNALS ###


def myFunc(e):

    return e.split("hz")[0]
    #return int(e.split("A_")[1].split("mV")[0])



#path = "/data/torus_mi/29_4_2022/sweep/"
path = "/data/torus_mi/10_5_2022/sweep/"
#path = "/data/torus_mi/31_5_2022/amplitude/"
#path = "/data/torus_mi/21_6_2022/"
path = "/data/torus_wt/01_23/17_1/"

files = [f for f in listdir(path) if isfile(join(path, f))]

files.sort(key=myFunc)

### SIGNAL CHOICE ###

print("###############################")
print("CHOOSE SIGNALS TO COMPARE")
print("")

for i in range(len(files)):
    print(str(i)+" - "+files[i])

print("")
print("###############################")
print("")

choices = input("")

choices = choices.split(",")

### PHYSICS ###

N = len(files)
fs = 2000
T = N/fs
dt = 1/2000


def plot_n(choices):


    fig, ax = plt.subplots(figsize=[12,9])

    for i in choices:
    
        signal = loadmat(path+files[int(i)])
        signal = signal[next(reversed(signal))]
    
        signal = signal[80000:,0]
        signal = signal-np.mean(signal)

        ff, pp = welch(signal, fs, window="blackman",noverlap = fs/2, nfft=4*fs, 
                       axis=0, scaling="density", detrend=False, nperseg=4*fs)

        ax.loglog(ff,pp,lw=1)
        
    ax.set_xlabel("$f$ [Hz]",fontsize=35)
    ax.set_ylabel(r"$S_\eta(\omega)$ [m$^2$s]",fontsize=35)
    ax.tick_params(labelsize=25)

    ax.set_xlim([3,1e3])
        
    return


plot_n(choices)

def plot_l(slope,c):
    
    x = np.linspace(20,90)

    plt.loglog(x,c*x**(slope))


plot_l(-3.5,1e-2)
    
plt.show()
    


