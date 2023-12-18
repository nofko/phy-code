import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

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

N = 1024

L = 4

n1 = 128
n2 = 750

n2 = 700 # random

x = np.linspace(0,(n2-n1)*L/N,n2-n1)
#x = np.linspace(0,L,N)
dx = x[1]-x[0]

f = np.linspace(0.5,1.8,50)

folder = "/data/topo/random_1024/A_4/"

print(folder)

############################         LENGTH        ############################


x0 = np.zeros((len(f)))

for i in range(len(f)):

    data = np.genfromtxt(folder+"out_p"+str(i))
    # data = data[:6000]
    
    y = np.max(data[2000:],axis=0)

    #a,b = np.polyfit(x,np.log(y[n1:n2]),1)

    r = curve_fit(lambda t,a,b: a*np.exp(b*t),  x,  y[n1:n2])

    print(r[0])
    print(i)

    x0[i] = r[0][1]

    
np.save(folder+"loc_len.npy",x0)

# x0 = np.load("loc_len.npy",allow_pickle=True)

# plt.plot(f[1:],-1/x0[1:],"o")

# plt.plot(x,y)
# plt.plot(x,np.exp(a*x+b))

# plt.semilogy()

# plt.xlabel("$f$ [Hz]")
# plt.ylabel("$x_0$ [m]")
# plt.savefig("local.pdf")


#plt.show()
