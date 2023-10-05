import matplotlib.pyplot as plt
import numpy as np



data = np.genfromtxt("out")

#data = data 

#data = data[1900:4000]

N = 2048

L = 12

x = np.linspace(0,L,N)
dx = L/N

steep = np.mean(np.sqrt(1/4*(np.sum(np.diff(data[2000:,:680],axis=1)**2/dx,axis=1))))

print(steep)

#plt.imshow(data)

plt.plot(x,np.max(data,axis=0))
#plt.imshow(data)

plt.plot(x,data[0])
plt.plot(x,data[50])
plt.plot(x,data[1000])


plt.xlabel("$x$ [m]")
plt.ylabel("$\eta$ [mm]")

plt.show()
