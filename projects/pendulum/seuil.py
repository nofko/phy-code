import matplotlib.pyplot as plt
import numpy as np

from os import listdir
from os.path import isfile, join
from scipy.io import loadmat
import scipy.signal as signal


plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                         AMPLITUDE VS FREQUENCY                            ##
##                                                                           ##
## ============================================================================



### LOAD SIGNALS ###


path = "/data/pendulum/10_25_2022/seuil.csv"

data = np.genfromtxt(path,delimiter=',')

f = data[1:,0]
a = data[1:,1]
a_l = data[1:,2]

era = data[1:,3]
erl  = data[1:,4]

fig, ax = plt.subplots(figsize=[12,9])

plt.errorbar(f,a,yerr=era,fmt="o",ms=8,ecolor="lightgray",label="$A$")
plt.errorbar(f,a_l,yerr=erl,fmt="o",ms=8,ecolor="lightgray",label="$A$")

#plt.fill_between(f, a, color='#539ecd')

ax.set_xlabel(r"$f$ [Hz]",fontsize=35)
ax.set_ylabel(r"$A$ [mV]",fontsize=35)
ax.tick_params(labelsize=25)

plt.show()
