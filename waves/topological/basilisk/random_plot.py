import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import h5py

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
## THEORETICAL ANALYSIS                                                      ##
##                                                                           ##
##                  LOCALIZATION LENGTH OVER A RANDOM 1D BOTTOM              ##
##                                                                           ##
## ============================================================================


############################         SETUP         ############################

N = 2048

L = 4

x = np.linspace(0,L,N)
dx = L/N

f = np.linspace(0.5,1.8,50)

f7 = np.linspace(0.5,1.8,70)


############################         LENGTH        ############################


data_1 = h5py.File('/data/topo/exp/Xi_f_Sweep_Aff105_Rand25.mat','r')

data_4 = h5py.File('/data/topo/exp/Xi_f_Sweep_Aff30_Rand25.mat','r')

theo = scipy.io.loadmat('/data/topo/exp/zeta_th_L25_N9.mat')

f_t = np.linspace(0,2.1,2000)
g = np.load("/data/topo/exp/random.npy",allow_pickle=True)

xs_01 = np.load("/data/topo/data_A01_4096/loc_len.npy",allow_pickle=True)

xs_03 = np.load("/data/topo/random_1024/A_03/loc_len.npy",allow_pickle=True)
xs_1 = np.load("/data/topo/random_1024/A_1/loc_len.npy",allow_pickle=True)
xs_15 = np.load("/data/topo/random_1024/A_15/loc_len.npy",allow_pickle=True)
xs_4 = np.load("/data/topo/random_1024/A_4/loc_len.npy",allow_pickle=True)


fe_1 = data_1["freq"]
ze_1 = np.array(data_1["zeta"])

fe_4 = data_4["freq"]
ze_4 = np.array(data_4["zeta"])


############################         PLOT         ############################


fig, ax = plt.subplots(figsize=[12,9])        


#plt.plot(f7,-1/xs_01,"o",color="y",label="$\epsilon=0.001$")
#plt.plot(f,-1/xs_1,"--",color="b",label="$\epsilon=0.001$")
plt.plot(f,-1/xs_15,"--",color="g",label="$\epsilon=0.01$")
plt.plot(f,-1/xs_4,"--",color="r",label="$\epsilon=0.04$")

plt.scatter(fe_1[:-1],ze_1[:-1]/100,color='g')
plt.scatter(fe_4[:-1],ze_4[:-1]/100,color='r')

plt.plot(f_t,0.01*g,color="k")

plt.ylim([-0.1,7])

# plt.plot(x,y)
# plt.plot(x,np.exp(a*x+b))

# plt.semilogy()

plt.ylim([-0.1,7])
plt.xlim([0.3,1.85])

plt.xlabel("$f$ [Hz]",fontsize=25)
plt.ylabel("$x_0$ [m]",fontsize=25)

ax.tick_params(labelsize=20,direction="in")

mdic = {"z15": -1/xs_15, "z4": -1/xs_4, "f": f}

scipy.io.savemat("simu_z_R.mat", mdic)


plt.savefig("random.pdf")



plt.show()
