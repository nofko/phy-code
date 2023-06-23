import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ==========================================================================##
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                         RUNGE-KUTTA PHASE PLOT                            ##
##                                                                           ##
## ==========================================================================##




############################          SETUP        ############################


th = np.load("rkp2.npy")


AS = np.linspace(0.1,5,100)
fs = np.linspace(0.1,3,100)

ap = (AS[-1]-AS[0])/(fs[-1]-fs[0])


############################         PLOTTING      ############################


fig, ax = plt.subplots(figsize=[12,9])        

u = th[:,::-1,0].transpose()

image = ax.imshow(u,
                  extent=[fs[0],fs[-1],AS[0],AS[-1]],
                  aspect=1/ap,
                  cmap="binary")



ax.set_xlabel("$f$",fontsize=30)
ax.set_ylabel(r"$A$ ",fontsize=30)

ax.tick_params(labelsize=20,direction="in")


plt.savefig("phase-diagram.pdf",bbox_inches='tight')


plt.show()
