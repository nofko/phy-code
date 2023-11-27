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
f4 = np.linspace(0.5,1.8,40)


############################         LENGTH        ############################


#data_1 = h5py.File('data_exp/Xi_f_Sweep_Aff4_Perio25.mat','r')

#data_4 = h5py.File('data_exp/Xi_f_Sweep_Aff30_Perio25.mat','r')

#theo = scipy.io.loadmat('data_exp/zeta_th_L25_N9.mat')

#f_t = theo['om']/2/np.pi

#g = theo['GG']

# for name in theo:
#     print(name)


xs_1 = np.load("loc_len_A4_4096_n.npy",allow_pickle=True)
xs_4 = np.load("loc_len_A4_4096.npy",allow_pickle=True)


# fe_1 = data_1["freq"]
# ze_1 = np.array(data_1["zeta"])

# fe_4 = data_4["freq"]
# ze_4 = np.array(data_4["zeta"])

print(xs_1)
plt.plot(f4,-1/xs_1,"o",color="r")
plt.plot(f,-1/xs_4,"+",color="r")

# plt.plot(fe_1[:-1],ze_1[:-1]/100,'o',color="b")
# plt.plot(fe_4[:-1],ze_4[:-1]/100,'+',color="b")

#plt.plot(f_t.transpose(),0.25/np.array(g),color="k")

plt.ylim([-0.1,40])

# plt.plot(x,y)
# plt.plot(x,np.exp(a*x+b))

# plt.semilogy()

plt.xlabel("$f$ [Hz]")
plt.ylabel("$x_0$ [m]")

#plt.savefig("local.pdf")


plt.show()
