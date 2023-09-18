import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ==========================================================================##
## NANTES CAMPAIGN                                                           ##
##                                                                           ##
##                                PLOTTING MSD                               ##
##                                                                           ##
## ==========================================================================##


############################          SETUP        ############################


big = np.genfromtxt("run2_big.csv",delimiter=',')

small = np.genfromtxt("run2_small.csv",delimiter=',')



############################          PLOT         ############################


fig, ax = plt.subplots(figsize=[12,9])        

ax.plot(big[:,0],big[:,1]*1e-6,"o", color="r")
ax.plot(small[:,0],small[:,1]*1e-6,"o",color= "b")


ax.set_xlabel(r"log$(t)$",fontsize=30)
ax.set_ylabel(r'$\langle \Delta r^2 \rangle$ [m$^2$]',fontsize=30)

ax.tick_params(labelsize=20,direction="in")


ax.set_xscale('log')
ax.set_yscale('log')

plt.show()
