import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import fsolve

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)

## ============================================================================
## THEORETICAL ANALYSIS                                                      ##
##                                                                           ##
##             CALCULATION OF DISPERSION OVER 1D PERIODIC BOTTOM             ##
##                                                                           ##
## ============================================================================


### PHYSICS ###

g = 9.81   # M/S^2
sigma = 0.0
rho = 1000

L = 1    # M
om0 = np.sqrt(g/L)

d = 0.4*L

h = L   # M
h1 = 0.05*L


print("####PHYSICAL SIZES####")
print("PERIOD L : ",L*100,"cm")
print("WATER DEPTH : ",h*100,"cm")
print("STEP WIDTH: ",d*100,"cm")
print("WATER ABOVE STEP: ",h1*100,"cm")
print("")
### DISPERSION RELATION ###


def disp_gravity(KK,HH):

    return np.sqrt((g*KK+sigma*KK**3/rho)*np.tanh(KK*HH))


####### FUNCTION TO SOLVE ########

def zero_gravity(KK,HH,OM):

    return disp_gravity(KK,HH)-OM


####### SCALAR SOLUTION #######

def inv_gravity(OM,HH):
    
    sol = fsolve(zero_gravity,0.2,args=(HH,OM))
        
    return sol


### RHS ###

def rhs(KK,KK1):

    zeta = -0.5*np.log(KK/KK1)
    
    return np.cos(KK1*d)*np.cos(KK*(L-d))-np.cosh(2*zeta)*np.sin(KK1*d)*np.sin(KK*(L-d))


### SOLVING ###

N = 500

om = np.linspace(0.001,4,N)*om0

k = np.zeros(N)

r = np.zeros(N)

for i in range(N):

    k = inv_gravity(om[i],h)
    k1 = inv_gravity(om[i],h1)

    r[i] = rhs(k,k1)


# fig, ax = plt.subplots(figsize=[12,8])

# plt.plot(om,r)

# ax.set_xlabel("$\omega/\omega_0$",fontsize=35)
# ax.set_ylabel(r"$\cos(KL)$",fontsize=35)
# ax.tick_params(labelsize=25)
    
# ax.axhline(1,linestyle="--",color="gray")
# ax.axhline(-1,linestyle="--",color="gray")


fig, ax = plt.subplots(figsize=[12,8])

plt.plot(np.arccos(r)/np.pi,om/om0,c="tab:red")
plt.plot(-np.arccos(r)/np.pi,om/om0,c="tab:red")

ax.set_xlabel(r"$K$ [m$^{-1}$]",fontsize=35)
ax.set_ylabel("$f$ [Hz]",fontsize=35)
ax.tick_params(labelsize=25)

ax.axhline(0.82)
ax.axhline(1.455)

plt.show()


