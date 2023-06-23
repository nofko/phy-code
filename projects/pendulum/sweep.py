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


path_f = "/data/pendulum/10_20_2022/sweep_forward_240s_5_40hz_A300.mat"
path_b = "/data/pendulum/10_20_2022/sweep_back_240s_5_40hz_A300.mat"


fs = 2000

T = 240


def get_amps(path):

    files = [f for f in listdir(path) if isfile(join(path, f))]

    N = len(files)

    amps = np.zeros(len(files))

    freq =  np.zeros(len(files))

    for i in range(N):
        
        S = loadmat(path+files[i])
        S = S[next(reversed(S))]*1000

        S = S-np.mean(S)

        sos = signal.butter(3,(1,35),'bandpass',fs=int(fs),output="sos")
        sff = signal.sosfilt(sos,S)

        sff = sff.reshape(len(S))
        S = S.reshape(len(S))
        ff, pp = signal.welch(S, fs, window="blackman",noverlap = fs/2, nfft=4*fs, 
                           axis=0, scaling="density", detrend=False, nperseg=4*fs)

        #plt.loglog(ff,pp,lw=1)

        #plt.plot(S[:-1],np.diff(S,1))

        sff_pks, _ = signal.find_peaks(sff,height=0)
        #plt.plot(sff_pks, sff[sff_pks], "x")

        amps[i] = np.mean(sff[sff_pks])
        freq[i] = float(files[i].split('_')[2])


    return freq,amps

def sweep(path,fmin,fmax,backwards):

    S = loadmat(path)

    S = S[next(reversed(S))]*1000

    S = S-np.mean(S)

    S = S[:fs*T]

    freq = np.linspace(fmin,fmax,len(S))

    if backwards:
        S = S[::-1]

    sos = signal.butter(3,(1,45),'bandpass',fs=int(fs),output="sos")
    sff = signal.sosfilt(sos,S)

    sff = sff.reshape(len(S))
    S = S.reshape(len(S))    

    sff_pks, _ = signal.find_peaks(sff,height=0,distance=100)

    plt.figure(figsize=[12,9])
    sp,f,t,fig = plt.specgram(sff,Fs=2000,NFFT=10*fs)

    plt.ylim([0,60])
    plt.title(path.split('/')[-1])
    plt.ylabel("$f$ [Hz]",fontsize=35)
    plt.xlabel("$t$ [s]",fontsize=35)

    # plt.plot(t,(fmax-fmin)*t/T+fmin)
    
    # plt.tick_params(labelsize=25)

    # plt.figure()
    # plt.plot(S)
    
    amps = sff[sff_pks]
    amps = np.convolve(amps, np.ones(20), "valid") / 20

    freq = freq[sff_pks]
    freq = freq[:-19]
    
    return freq,amps

### PHYSICS ###

ff,af = sweep(path_f,5,40,False)
fb,ab = sweep(path_b,5,40,True)

fig, ax = plt.subplots(figsize=[12,9])

ax.set_xlabel(r"$f$ [Hz]",fontsize=35)
ax.set_ylabel(r"$A$ [mV]",fontsize=35)
ax.tick_params(labelsize=25)

plt.plot(ff-1,af,color='b')
plt.plot(fb,ab,color='b',ls='--')


plt.show()
    


