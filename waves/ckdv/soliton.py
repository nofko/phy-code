import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
## ============================================================================



############################          SETUP        ############################


path = "/data/ckdv/19_06/14-Jun-2023-Result-tau10dms-a300mV"

fs = 1000

norm = 14.5/6*1e-3  ## mm/V

data = np.loadtxt(path)


fig, ax = plt.subplots(figsize=[12,9])        


solit = data[5001:10000,:] - np.mean(data[:4990,:],axis=0)

#ax.plot(solit[300,:])
ap = len(solit[:,0])/len(solit[0])
         
plt.imshow(solit,aspect=1/ap)


plt.figure()

plt.plot(solit[260])

plt.show()


