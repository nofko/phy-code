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
path = "/data/torus_wt/02_23/17_02/bruitCTP_2_5hz_7min_A1000.npy"

S = np.load(path)

#S = np.roll(S,-430)
#S = S[:, 250:]

# S = np.roll(S,-2730)
# S = S[:, 500:]

S = np.roll(S,-400)
S = S[:, 500:]

#S = np.pad(S,[(250,250),(0,0)])



print("-- FILENAME : ", path.split("/")[-1])

wk = get_window("blackman", len(S[0]))

# x = np.linspace(0,2*np.pi,len(S[0]))
# t = np.arange(len(S))/50

# xx,tt = np.meshgrid(x,t)

# f1 = 5
# f2 = 4
# kk1 = 5
# kk2 = 15

# f3 = f1+f2
# kk3 = kk1+kk2

# S =  np.sin(2*np.pi*f1*tt+kk1*xx)+np.sin(2*np.pi*f2*tt+kk2*xx)*np.sin(2*np.pi*f3*tt+kk3*xx)

####### DISPERSION RELATIONS #######

def varicose(k,w,ro):

    return np.sqrt((g*k/ro+sigma*k**3/rho/ro**3)*np.tanh(k*w*ro/2/rc**2))


def sinuous(k,w,ri):

    om0 = np.pi/2*np.sqrt(9.81*np.sin(9*np.pi/180)/2/w)

    return np.sqrt(om0**2+np.abs(sigma/rho*(k/ri)**3))
    
def sloshing(k,om0):

    return np.sqrt((om0**2+(g*k**2/(0.114)+0.9*sigma*np.abs(k)**3/rho/ro**3)))
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

    #FF[i] = np.fft.fftshift(np.fft.fft(S[i]*wk))

    FF[i] = np.fft.fft(S[i]*wk)

K0 = int(len(FF[0])/2)

plt.plot(np.log(np.abs(FF[100])))

print("-- DUMPING SIGNAL FROM RAM --")

del(S)

### BICOHERENCE CALCULATION ###

print("-- CALCULATING BICOHERENCE, PLEASE HOLD --")

NB = 200

B = np.zeros((NB,NB))

nb = int(NB/2)

# for k1 in range(-nb,nb):
#     print((k1+nb)/NB)
#     for k2 in range(-nb,nb):
        
#         N = np.sqrt(np.mean(np.abs(FF[:,K0+k1]*FF[:,K0+k2])**2)*np.mean(np.abs(FF[:,K0+k1+k2])**2))
        
#         B[k1+nb,k2+nb] = np.abs(np.mean(np.conjugate(FF[:,K0+k1]*FF[:,K0+k2])*FF[:,K0+k1+k2]))

for k1 in range(0,NB):
    print((k1)/NB)
    nn = int(len(FF[0])/2)
    B[k1] = np.abs(np.mean(np.conjugate(FF[:,k1]*FF[:,:nn-k1].transpose())*FF[:,k1:nn].transpose(),axis=1))[:NB]
    #print(B[k1])
    
# print("-- DUMPING FOURIER FROM RAM --")

#fig, ax = plt.subplots(figsize=[12,9])

#image = ax.imshow(np.log(np.abs(FF[:,:100])),cmap="jet",extent=[0,100,0,len(FF)],aspect=100/len(FF))

del(FF)


### PLOTTING ###

print("-- DONE!")


k = np.linspace(0,NB)

fig, ax = plt.subplots(figsize=[12,9])

image = ax.imshow(np.log(B[::-1,:]),cmap="turbo",extent=[k[0],k[-1],k[0],k[-1]],aspect=1)

cbar = plt.colorbar(image,format="%1.1f",ax=ax)
cbar.set_label(r"$B(k_1,k_2)$",fontsize=30)
cbar.ax.tick_params(labelsize=25)

kmax = 150

N = 1000

k1 = np.linspace(-kmax,kmax,N)
k2 = np.linspace(-kmax,kmax,N)

kk1,kk2 = np.meshgrid(k1,k2)

# ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)-sloshing(kk1+kk2,om0),[0],colors="k")

# ax.contour(kk1,kk2,varicose(kk1,w,ro)+sloshing(kk2,om0)-varicose(kk1+kk2,w,ro),[0],colors="k")


ax.set_xlabel("$k_1$",fontsize=35)
ax.set_ylabel("$k_2$",fontsize=35)
ax.tick_params(labelsize=25)




#plt.savefig("images/bicoherence_"+path.split("/")[-1].split(".mat")[0]+".pdf",bbox_inches="tight")
#np.save("data/bicoherence_"+path.split("/")[-1].split(".mat")[0],B)

plt.show()
