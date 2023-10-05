import matplotlib.pyplot as plt
import numpy as np


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

L = 6

x = np.linspace(0,L,N+1)
dx = L/N

f = np.arange(5,18)/10


############################         LENGTH        ############################


x0 = np.zeros((len(f)))

for i in range(len(f)):

    data = np.genfromtxt("data/out_p"+str(i))

    y = np.max(data,axis=0)

    a,b = np.polyfit(x[:1000],np.log(y[:1000]),1)

    print(a,b)

    x0[i] = a

    
np.save("loc_len.npy",x0)

plt.plot(f,-1/x0,"o")

# plt.plot(x,np.exp(a*x+b))

# plt.semilogy()

plt.xlabel("$f$ [Hz]")
plt.ylabel("$x_0$ [m]")
plt.savefig("local.pdf")
#plt.show()
