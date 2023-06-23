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


path = "data/tricoherence_7_7Hz_A900.npy"

B = np.load(path)

print("-- FILENAME : ", path.split("/")[-1])


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

### PLOTTING ###

print("-- DONE!")

k3 = 20

NB = 250

k = np.linspace(-int(NB/2),int(NB/2),NB)

fig, ax = plt.subplots(figsize=[12,9])

image = ax.imshow(np.log(B[:,::-1]),cmap="turbo",extent=[k[0],k[-1],k[0],k[-1]],aspect=1)

cbar = plt.colorbar(image,format="%1.1f",ax=ax)
cbar.set_label(r"$B(k_1,k_2)$",fontsize=30)
cbar.ax.tick_params(labelsize=25)

kmax = NB/2

N = 1000

k1 = np.linspace(-kmax,kmax,N)
k2 = np.linspace(-kmax,kmax,N)

delta = 4

kk1,kk2 = np.meshgrid(k1,k2)

ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)+varicose(k3,w,ro)-varicose(kk1+kk2+k3,w,ro),[0],colors="k")

ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)-varicose(k3,w,ro)-varicose(kk1+kk2-k3,w,ro),[0],colors="k")

ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)-varicose(k3,w,ro)-varicose(kk1+kk2-k3,w,ro),[0],colors="k")

ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)-varicose(k3,w,ro)-varicose(kk1+kk2-k3,w,ro),[0],colors="k")

#ax.contour(kk1,kk2,varicose(kk1,w,ro)+sloshing(kk2,om0)-varicose(kk1-kk2,w,ro),[0],colors="k")


ax.set_xlabel("$k_1$",fontsize=35)
ax.set_ylabel("$k_2$",fontsize=35)
ax.tick_params(labelsize=25)

ax.set_xlim([-kmax,kmax])
ax.set_ylim([-kmax,kmax])


plt.show()
