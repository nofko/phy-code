import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipj

mm = 0.9999
Lambda = np.sqrt(mm) * ellipk(mm)

Xmin = 0
Xmax = 2*Lambda
Xhalf = (Xmax-Xmin)/2.

print(Xmax)

plt.figure()
XX = np.linspace(Xmin, Xmax, 1001)
YYper = 2*ellipj((XX-Xhalf)/np.sqrt(mm),mm)[3] + np.pi

YYkink = 4*np.arctan(np.exp(XX-Xhalf))

plt.plot(XX, YYper , label ='am function')
plt.plot(XX, YYkink, label ='1-kink function')

plt.grid()
plt.xlabel('x')
plt.legend()
plt.title('Static am$\,(x,m)/\pi+1$, m = '+str(round(mm,4)))

plt.show()

