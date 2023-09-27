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

L = 12

x = np.linspace(0,L,N)
dx = L/N

f = np.arange(5,18)/10

############################         STEEP         ############################


def steep(data,n):

    return np.mean(np.sqrt(1/4*(np.sum(np.diff(data[n:,:680],axis=1)**2/dx,axis=1))))



n = 10

epsilon = np.zeros((len(f),n))

for i in range(len(f)):

    data = np.genfromtxt("data/out_f"+str(i))

    T = len(data)

    index = np.arange(0,T,int(T/n))
    
    for j in range(n):

        epsilon[i,j] = steep(data,index[j])


        



############################         PLOT          ############################



for i in range(10):
    
    plt.plot(epsilon[i,:]/f[i]*6,"o")


print(epsilon[:,-1])
plt.figure()
plt.plot(f,epsilon[:,-1],"o")


    
plt.show()
