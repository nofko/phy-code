
N=1000

x = np.linspace(-3,3,N)
y = np.linspace(-3,3,N)

xx,yy = np.meshgrid(x,y)

l = 1

psi = 4*np.arctanh(l*np.cos(np.sqrt(l**2+1)*xx)/(np.sqrt(1+l**2)*np.cosh(l*yy)))

plt.pcolormesh(xx,yy,psi)
