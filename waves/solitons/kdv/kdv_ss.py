N = 1024
T = 10


t = np.linspace(0,T,3000)

dt = t[1]-t[0]

x = np.linspace(-10,10,N);

s = np.zeros((N,N))

dx = x[1]-x[0]

dk = 2*np.pi/(N*dx)

k = np.arange(int(N/2)+1)*dk

k = np.concatenate((k,-np.flip(k[1:-1]))) 

#u = 0.2*np.exp(-x**2)
u = 1/np.cosh(x/np.sqrt(2))**2+8/np.cosh(2*(x+2))**2

v = np.fft.fft(u)


for i in range(N):

    s[i] = u

    v = v*np.exp(1j*k**3*dt)
    v -= 3j*k*dt*np.fft.fft(np.real(np.fft.ifft(v))**2)
    
    u = np.real(np.fft.ifft(v))
    

#plt.figure()

for i in s:

   plt.clf()
   plt.plot(i)
   plt.pause(0.00001)
   plt.draw()
