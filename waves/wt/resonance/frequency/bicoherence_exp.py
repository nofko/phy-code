import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch
from scipy.io import loadmat
from scipy.fftpack import next_fast_len
from scipy.signal import get_window

# from polycoherence import _plot_signal, polycoherence, plot_polycoherence

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)

## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##          CALCULATION OF 3-WAVE BICOHERENCE IN FREQUENCY DOMAIN            ##
##                                                                           ##
## ============================================================================

red = "#ff0054"
green = "#0d954c"

### LOAD SIGNAL ###

print("#############################################")
print("################ BICOHERENCE ################")
print("")

path = "/data/torus_mi/5_5_2022/cascade/6_1hz_40min_A_1100mV.mat"
path = "/data/torus_mi/21_6_2022/7_7hz_60min_A_1000mV.mat"

#path = "/data/torus_mi/31_5_2022/amplitude/7_4hz_4min_A_580mV.mat"
path = "/data/torus_wt/02_23/17_02/bruitCTP_2_5hz_15min_A1000.mat"


s = loadmat(path)
s = s["signal"]
s = s[:,0]
s = s[84000:]

fs = 2000
nop = len(s)

s = s[:fs*2000]

print("-- FILENAME : ", path.split("/")[-1])

####### PHYSICS #######

g = np.sin(4.5*np.pi/180)*9.81

sigma = 0.055     # SURFACE TENSION
rho = 1000        # FLUID DENSITY

rc = 0.07         # CENTRAL RADIUS IN M
ro = 0.0785       # TORUS OUTER RADIUS
w = 2*(ro-rc)     # TORUS WIDTH
ri =ro-w          # TORUS INNER RADIUS

q = 1             # DISPERSION RELATION WIDTH

#om0 = 42.4  # R=7.9  # CUTOFF FREQUENCY OF SLOSHING BRANCH IN RAD/S
om0 = 44.5   # R=8.1

####### DISPERSION RELATIONS #######

def varicose(k,w,ro):

    return np.sqrt((g*k/ro+sigma*k**3/rho/ro**3)*np.tanh(k*w*ro/2/rc**2))

def sloshing(k,om0):
    
    return np.sqrt((om0**2+(g*k**2/(0.114)+0.9*sigma*k**3/rho/ro**3)))

    # return np.sqrt(om0**2+g*k**2/0.076)


### PWELCH ###

f, Pxx = welch(s, fs, window="blackman",noverlap = fs/2, nfft=4*fs, 
               axis=0, scaling="density", detrend=False, nperseg=4*fs)


fig, ax = plt.subplots(figsize=[8,6])
ax.set_title(path.split("/")[-1])

ax.loglog(f,Pxx,lw=1)

ax.set_xlabel("$f$ [Hz]",fontsize=35)
ax.set_ylabel(r"$S_\eta(\omega)$ [a.u.]",fontsize=35)
ax.tick_params(labelsize=25)

#ax.set_xlim(1.7,1e3)

###  FOURIER CALCULATION ###

print("-- CALCULATING FOURIER --")

NW = 101                                                      # NUMBER OF TEMPORAL WINDOWS

LW = len(s)//NW + 1                                           # WINDOW LENGTH

if LW % 2 ==0:
    LW -= 1

SW = np.array([s[i:i+LW] for i in range(0,LW*(NW-1),LW)])     # THE MATRIX OF THE WINDOWED SIGNAL

FF = np.zeros_like(SW,dtype=complex)                         # FOURIER MATRIX

wo = get_window("blackman", len(SW[0]))

for i in range(NW-1):

    FF[i] = np.fft.fftshift(np.fft.fft(SW[i]*wo))

F = np.arange(-int(LW-1)/2,int(LW-1)/2+1)*fs/LW 
F0 = int(LW/2)

### BICOHERENCE CALCULATION ###

print("-- CALCULATING BICOHERENCE, PLEASE HOLD --")

NBI = 1000

nbi = int(NBI/2)

B = np.zeros((NBI,NBI))                                     # BICOHERENCE MATRIX

for i in range(0,NBI):

    for j in range(0,NBI):
        
        N = np.sqrt(np.mean(np.abs(FF[:,F0+i]*FF[:,F0+j])**2)*np.mean(np.abs(FF[:,F0+i+j])**2))
        B[i,j] = np.abs(np.mean(np.conjugate(FF[:,F0+i]*FF[:,F0+j])*FF[:,F0+i+j]))/N


# kw = dict(nperseg=nop // 17, noverlap=nop // 20, nfft=next_fast_len(nop // 2))
# freq1, freq2, bicoh = polycoherence(s, fs, **kw)
# plot_polycoherence(freq1, freq2, bicoh)


####### RESONANCE VAR+VAR=VAR #######

N = 1000

kmax = 40

k1 = np.linspace(0,kmax,N)
k2 = np.linspace(-kmax,0,N)

kk1,kk2 = np.meshgrid(k1,k2)


fig, ax = plt.subplots(figsize=[8,6])

ax.set_title(r"$\omega_V^1+\omega_V^2=\omega_V^3$")

ax.set_xlabel('$k_1$',fontsize=35)
ax.set_ylabel('$k_2$',fontsize=35)
ax.tick_params(labelsize=25)

cs = ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)-sloshing(kk1+kk2,om0),[0],colors="r")
#cs = ax.contour(kk1,kk2,varicose(kk1,w,ro)+sloshing(kk2,om0)-varicose(kk1+kk2,w,ro),[0],colors="b")


data = cs.allsegs[0][0]

kr1 = data[:,0]
kr2 = data[:,1]

or1 = varicose(kr1,w,ro)
or2 = varicose(kr2,w,ro)

fr1 = or1/2/np.pi
fr2 = or2/2/np.pi

or3 = or1+or2 
kr3 = kr1+kr2


### PLOTTING ###

print("-- DONE!")

f1 = 7.5

fig, ax = plt.subplots(figsize=[12,9])

image = ax.imshow(B[::-1,:],cmap="magma_r",extent=[F[F0]/f1,F[F0+NBI]/f1,F[F0]/f1,F[F0+NBI]/f1],vmin=0.65,vmax=1.5)

f = np.linspace(0,7.9)


ax.plot(f/f1,(f1-f)/f1,'--k',lw=2)
ax.axvline(1,ls='--',color='k',lw=2)
ax.axhline(1,ls='--',color='k',lw=2)

ax.plot(fr1/f1,fr2/f1,'--',color=red,lw=2)

cbar = plt.colorbar(image,format="%1.1f",ax=ax)
cbar.set_label(r"$B(\nu_1,\nu_2)$",fontsize=30)
cbar.ax.tick_params(labelsize=25)

ax.set_xlim([0,1.5])
ax.set_ylim([0,1.5])

ax.set_xlabel(r"$\nu_1/f_1$",fontsize=35)
ax.set_ylabel(r"$\nu_2/f_1$",fontsize=35)
ax.tick_params(labelsize=25)

plt.savefig("images/bicoherence_"+path.split("/")[-1].split(".mat")[0]+".pdf",bbox_inches="tight")

plt.show()


np.save("data/bicoherence_"+path.split("/")[-1].split(".mat")[0],B)
