import matplotlib.pyplot as plt
import numpy as np



data = np.genfromtxt("out")

data = data * 1000

data = data[1900:4000]

N = 2048

L = 6

x = np.linspace(0,L,N)

plt.plot(x,np.max(data,axis=0))
#plt.imshow(data)

plt.plot(x,data[0])
plt.plot(x,data[50])
plt.plot(x,data[100])


plt.xlabel("$x$ [m]")
plt.ylabel("$\eta$ [mm]")

plt.show()
