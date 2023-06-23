import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import get_window

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)



## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##          CALCULATION OF 3-WAVE BICOHERENCE IN WAVENUMBER DOMAIN           ##
##                                                                           ##
## ============================================================================



### LOAD SIGNAL ###

print("#############################################")
print("################ BICOHERENCE ################")
print("")

path = "/data/torus_mi/5_5_2022/bruitCTP_5_7hz_15min_A1100.npy"
path = "/data/torus_mi/29_9_2022/7_9hz/7_9Hz_4min_600.npy"
path = "/data/torus_wt/02_23/17_02/bruitCTP_2_5hz_5min_A1000.npy"
path = "/data/torus_wt/02_23/17_02/bruitCTP_2_5hz_7min_A1000.npy"

S = np.load(path)

#S = np.roll(S,-430)
#S = S[:, 250:]

#S = np.roll(S,-2730)
#S = S[:, 500:3000]

S = np.roll(S,-400)
S = S[:, 500:]
S = np.pad(S,[(250,250),(0,0)])

print("-- FILENAME : ", path.split("/")[-1])

wk = get_window("blackman", len(S[0]))


####### DISPERSION RELATIONS #######

def varicose(k,w,ro):

    return np.sqrt((g*k/ro+sigma*k**3/rho/ro**3)*np.tanh(k*w*ro/2/rc**2))


def sinuous(k,w,ri):

    om0 = np.pi/2*np.sqrt(9.81*np.sin(9*np.pi/180)/2/w)

    return np.sqrt(om0**2+np.abs(sigma/rho*(k/ri)**3))
    
def sloshing(k,om0):
    
    return np.sqrt(om0**2+g*k**2/0.076)


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


###  FOURIER CALCULATION ###

print("-- CALCULATING FOURIER --")

fps = 100

FF = np.zeros_like(S,dtype=complex)

for i in range(len(S)):

    FF[i] = np.fft.fftshift(np.fft.fft(S[i]*wk))


    
K0 = int(len(FF[0])/2)


print("-- DUMPING SIGNAL FROM RAM --")

del(S)

### BICOHERENCE CALCULATION ###

print("-- CALCULATING BICOHERENCE, PLEASE HOLD --")

NB = 250

B = np.zeros((NB,NB))

nb = int(NB/2)

k3 = 20
k4 = 30
k5 = 40

for k1 in range(-nb,nb):

    for k2 in range(-nb,nb):
        print((k1+nb)/NB)
        #N = np.sqrt(np.mean(np.abs(FF[:,K0+k1]*FF[:,K0+k2]*FF[:,K0+k3])**2)*np.mean(np.abs(FF[:,K0+k1+k2+k3])**2))
        
        B[k1+nb,k2+nb] = np.abs(np.mean(np.conjugate(FF[:,K0+k1]*FF[:,K0+k2])*FF[:,K0+k4]*FF[:,K0+k5]*FF[:,K0+k1+k2+k3-k4-k5]))



print("-- DUMPING FOURIER FROM RAM --")

del(FF)


### PLOTTING ###

print("-- DONE!")


k = np.linspace(-int(NB/2),int(NB/2),NB)

fig, ax = plt.subplots(figsize=[12,9])

image = ax.imshow(np.log(B[::-1,:]),cmap="plasma",extent=[k[0],k[-1],k[0],k[-1]],aspect=1,interpolation='nearest')

cbar = plt.colorbar(image,format="%1.1f",ax=ax)
cbar.set_label(r"$B(k_1,k_2)$",fontsize=30)
cbar.ax.tick_params(labelsize=25)

kmax = 150

N = 1000

k1 = np.linspace(-kmax,kmax,N)
k2 = np.linspace(-kmax,kmax,N)


kk1,kk2 = np.meshgrid(k1,k2)

#ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)+varicose(k3,w,ro)-varicose(kk1+kk2-k3,w,ro),[0],colors="k")

#ax.contour(kk1,kk2,varicose(kk1,w,ro)+sloshing(kk2,om0)-varicose(kk1-kk2,w,ro),[0],colors="k")


ax.set_xlabel("$k_1$",fontsize=35)
ax.set_ylabel("$k_2$",fontsize=35)
ax.tick_params(labelsize=25)




#plt.savefig("images/tricoherence_"+path.split("/")[-1].split(".mat")[0]+".pdf",bbox_inches="tight")
#np.save("data/tricoherence_"+path.split("/")[-1].split(".mat")[0],B)

plt.show()
