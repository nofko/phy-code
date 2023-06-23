import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
from scipy.optimize import fsolve, bisect
from scipy.signal import argrelmin,argrelmax
import mat73
import time

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)



## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                        PLOTTING SPACE-TIME RESULTS                        ##
##                                                                           ##
## ============================================================================




path_ei = "/data/canal/Signaux_Soliton_Sin_Noise/results/CumSum_one_soliton_A10mm_h3cm_L4m_30s_Signal_Tot_ei.npy"

path_e = "/data/canal/Signaux_Soliton_Sin_Noise/results/CumSum_one_soliton_A10mm_h3cm_L4m_30s_Signal_Tot_e.npy"
path_tr = "/data/canal/Signaux_Soliton_Sin_Noise/results/CumSum_one_soliton_A10mm_h3cm_L4m_30s_Signal_Tot_tr.npy"


Eigs = np.load(path_ei)
e = np.load(path_e)
tr = np.load(path_tr)

Eigs = np.sort(Eigs)

h = 3 #water depth
lam = 3/2/h**3

fps = 10

N, M = np.shape(Eigs)

sindex = np.zeros((N,int(len(Eigs)/2)-1))

sref = np.zeros(N)
eref = np.zeros(N)

amps = np.zeros((N,int(len(Eigs)/2)-1))


for i in range(N):
    
    for j in range(1,int(M/2)-1):
        
        sindex[i,j-1] = (Eigs[i][2*j]-Eigs[i][2*j-1])/(Eigs[i][2*j]-Eigs[i][2*j-2])


for i in range(N):

     r = [ n for n,j in enumerate(sindex[i]) if j>0.99]

     if len(r)>0:

         sref[i] = r[-1]

     else:
         
         sref[i] = -1
     
for i in range(N):
    
    for j in range(1,int(M/2)):

        r = int(sref[i])+1+(1-j)

        eref[i] = Eigs[i][2*r]
        
        if r>0:
            
            amps[i,j-1] = 2*(Eigs[i][2*r]-Eigs[i][2*j-1])/lam

        else:

            amps[i,j-1] = (Eigs[i][2*j]-Eigs[i][2*j-1])/lam/2


t = np.arange(N)/fps

    
fig, ax = plt.subplots(figsize=[8,6])


ax.set_xlabel('$t$ [s]',fontsize=35)
ax.set_ylabel('$s$',fontsize=35)
ax.tick_params(labelsize=25)
            
ax.plot(t,sindex[:,0])
ax.plot(t,sindex[:,1])
ax.plot(t,sindex[:,2])

ax.axhline(0.99,ls="--")


# fig, ax = plt.subplots(figsize=[8,6])

# print(sref[50])
# print(amps[50])

# ax.set_xlabel('$t$ [s]',fontsize=35)
# ax.set_ylabel('$A$ [cm]',fontsize=35)
# ax.tick_params(labelsize=25)
            
# ax.plot(e,tr[50,:])
# ax.scatter(Eigs[50][::2],np.ones_like(Eigs[50][::2]))
# ax.scatter(Eigs[50][1::2],-np.ones_like(Eigs[50][1::2]))

# plt.axhline(1,0,1,ls="--",color="tab:red")
# plt.axhline(-1,0,1,ls="--",color="tab:red")
# plt.axvline(0,0,1,ls="--",color="tab:red")
# plt.axhline(0,0,1,ls="--",color="tab:red")

fig, ax = plt.subplots(figsize=[8,6])

ax.plot(t,amps[:,0])
ax.plot(t,amps[:,1])

ax.set_xlabel('$t$ [s]',fontsize=35)
ax.set_ylabel('$A$ [cm]',fontsize=35)
ax.tick_params(labelsize=25)

fig, ax = plt.subplots(figsize=[8,6])

ax.plot(t,eref*(sref+1))

print(sref+1)

plt.show()
