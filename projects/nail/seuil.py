import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)




## ==========================================================================##
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                                 THRESHOLD                                 ##
##                                                                           ##
## ==========================================================================##




############################        LOAD FILE      ############################


path = "/data/nail/01_06/seuil.csv"


data = np.genfromtxt(path,delimiter=",")


f = data[1:,0]

A = data[1:,2]


plt.errorbar(f,A,yerr=0.05*A,fmt="o",color="b")

plt.xlabel(r"$1/f$ [s]", fontsize=25)
plt.ylabel(r"$B$ [mT]", fontsize=25)

plt.savefig("seuil.pdf",bbox_inches="tight")

plt.show()
