import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import get_window


N = 5000

L = 10
T = 10

c = 1

x = np.linspace(-L/2,L/2,N)
t = np.linspace(0,T,N)

k = np.arange(0,N)*2*np.pi/L
omega = c*k

k_wave = t/T*np.pi/L/2
o_wave = c*k_wave

xx, tt = np.meshgrid(x,t)

v = 0.4

u = 1/np.cosh((xx-v*tt+L/4)/0.3)**2


plt.plot(u[0])
plt.plot(u[-1])

plt.figure()

#wk = get_window("hann", N)
#wo = get_window("hann", N)

#u = u*np.outer(wo,wk)


ff = np.fft.fft2(u)

kmax = 200
fmax = 100

plt.imshow(np.log(np.abs(ff[-fmax:-1,:kmax])),aspect=kmax/fmax)


plt.colorbar()

plt.show()
