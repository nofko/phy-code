import matplotlib.pyplot as plt
import numpy as np
import pywt
from scipy.signal import welch

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)

## ============================================================================
## THEORETICAL ANALYSIS                                                      ##
##                                                                           ##
##                       TESTING OUT WAVELET ANALYSIS                        ##
##                                                                           ##
## ============================================================================



T = 1
fs = 2000

s = 129

N = T*fs

x = np.linspace(0,T,N)

y = 10*np.sin(x*5*2*np.pi)

coef, freqs=pywt.cwt(y,np.arange(1,s),'gaus1',T)


plt.pcolormesh(coef)
print(freqs)

plt.figure()
plt.plot(x,y)

### BICOHERENCE ###




plt.show()
