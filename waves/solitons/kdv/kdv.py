N = 1024
T = 1
L = 500

t = np.linspace(0,T,N)
x = np.linspace(-L,L,N);

dt = t[1]-t[0]
dx = x[1]-x[0]
dk = 2*np.pi/(N*dx)



k = np.arange(int(N/2)+1)*dk
k = np.concatenate((k,-np.flip(k[1:-1]))) 

u = np.exp(-x**2)
u = 1/np.cosh(x/np.sqrt(2))**2#+8/np.cosh(2*(x+2))**2
v = np.fft.fft(u)


s = np.zeros((N,N))

ik3 = 1j*k**3

for i in range(N):

    s[i] = u
    
    g = -3*1j*k*dt
    
    E = np.exp(dt*ik3)

    a = g*np.fft.fft(np.real(np.fft.ifft(v))**2)
    b = g*np.fft.fft(np.real(np.fft.ifft(v+E*a/2))**2)
    c = g*np.fft.fft(np.real(np.fft.ifft(v+E*b/2))**2)
    d = g*np.fft.fft(np.real(np.fft.ifft(v+E*c))**2)

    v = v+(a+2*(b+c)+d)/6
    u = np.real(np.fft.ifft(v))

#plt.figure()

for i in s:

   plt.clf()
   plt.plot(i)
   plt.pause(0.00001)
   plt.draw(
