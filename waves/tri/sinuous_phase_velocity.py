import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
from scipy.optimize import fmin


plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)

N =1000

####### PHYSICS #######

rc = 0.07

ro = 0.0786 #np.linspace(0.077,0.089,N)
w = 2*(ro-rc)
ri =ro-w

sigma = 0.055
rho = 1000

om0 = 44.5

g = np.sin(4.5*np.pi/180)*9.81


####### DISPERSION RELATIONS #######

def varicose(k,w,ro):

    return np.sqrt((g*k/ro+sigma*k**3/rho/ro**3)*np.tanh(k*w*ro/2/rc**2))

def sinuous(k,w,ri):

    om0 = np.pi/2*np.sqrt(9.81*np.sin(9*np.pi/180)/2/w)

    return np.sqrt(om0**2+np.abs(sigma/rho*(k/ri)**3))

def sin_v(k, w, ri):
    
    return sinuous(k, w, ri) / k

def sloshing(k,om0):
    
    return np.sqrt(om0**2+g*k**2/rc)

def difference(k,w,ro,om0):

    return sloshing(k,om0)-varicose(k,w,ro)

def distance(k2,k1,w,ro,om0):

    om1 = sloshing(k1,om0)
    om2 = varicose(k2,w,ro)

    return np.sqrt((k2-k1)**2+(om2-om1)**2)

####### PLOTTING #######


fig, ax = plt.subplots(figsize=[8,6])

ax.set_xlabel('$R$',fontsize=35)
ax.set_ylabel('$f_{min}$',fontsize=35)
ax.tick_params(labelsize=25)

k_min = np.zeros(N)

k = np.linspace(0.1,30,N)
dk = k[1]-k[0]
# for i in range(N):

#     m = fmin(distance, 7, args=(k[i],w, ro,om0))[0]

#     k_min[i] = m
    

plt.plot(k[:-1],np.diff(varicose(k,w,ro))/dk)
plt.plot(k[:-1],np.diff(sloshing(k,om0))/dk)

plt.plot(k,varicose(k,w,ro)/k)

plt.plot(k[:-1],np.diff(sloshing(k,om0))/dk-np.diff(varicose(k,w,ro))/dk)
plt.axhline(0)

print(sloshing(6.41,om0)/2/np.pi)
# for i in range(N):
    
#     m = fmin(sin_v, 7, args=(w[i], ri[i]))[0]

#     f[i] = sinuous(m,w[i],ri[i])/2/np.pi

# plt.plot(ro,f)


plt.show()
