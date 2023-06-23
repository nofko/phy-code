from os import listdir
from os.path import isfile, join
from scipy.io import loadmat
from scipy.signal import welch

#path = "/data/torus_mi/29_4_2022/sweep/"
path = "/data/torus_mi/10_5_2022/sweep/"
#path = "/data/torus_mi/31_5_2022/amplitude/"
path = "/data/torus_mi/11_5_2022/cascade/"
path = "/data/pendulum/19_01_2023/"

def myFunc(e):

    return e.split("hz")[0]
    #return int(e.split("A_")[1].split("mV")[0])

files = [f for f in listdir(path) if isfile(join(path, f))]

files.sort(key=myFunc)

N = len(files)
T = 60*10
fs = 2000
dt = 1/2000

n = 2


def plot_n(n):
    
    signal = loadmat(path+files[n])
    signal = signal[next(reversed(signal))]
    #signal = signal["signal"]
    
    signal = signal[:,0]
    signal = signal-np.mean(signal)
    #signal = np.convolve(signal, np.ones(100), "valid") /100

    plt.figure(figsize=[8,6])
    sp,f,t,fig = plt.specgram(signal,Fs=2000,NFFT=10*fs)

    plt.title(files[n])
    
    plt.ylim([0,15])
    plt.ylabel("$f$ [Hz]")
    plt.xlabel("$t$ [s]")

    # plt.figure()
    # plt.plot(signal,lw=0.5)

    # f, Pxx = welch(signal, fs, window="blackman",noverlap = fs/2, nfft=4*fs, 
    #                axis=0, scaling="density", detrend=False, nperseg=4*fs)

    # plt.figure()
    # plt.loglog(f,Pxx,lw=1)
    # plt.title(files[n])

    return
    


