


def omega(k,h):

    return np.sqrt(np.tanh(k*h)*(g*k+T*k**3/rho))

def omega_prime(k,h):

    return 0.5/omega(k,h)*((g+3*T*k**2/rho)*np.tanh(k*h)+(g*k+T*k**3/rho)/np.cosh(k*h)**2*h)


def mu(k,h):

    om = omega(k,h)
    sigma = np.tanh(k*h)
    
    result = -1/8/om**3*((g**2-6*g*T*k**2-3*T**2*k**4)*sigma**2-2*(g+T*k**2)*(g+3*T*k**2)*k*h*(sigma-sigma**3)+(g+T*k**2)**2*k**2*h**2*(1+2*sigma**2-3*sigma**4))

    return result


def nu(k,h):

    om = omega(k,h)
    sigma = np.tanh(k*h)

    omp = omega_prime(k,h)

    result = g**2*(9*sigma**4-10*sigma**2+9)+g*T*k**2*(15*sigma**4-44*sigma**2+30)+T**2*k**4*(6*sigma**4-25*sigma**2+21)
    result *= k/(g*sigma**2+T*k**2*(sigma**2-3))
    result += 2*(g+T*k**2)/(omp**2-g*h)*(6*(g+T*k**2)*sigma-2*(g+3*T*k**2)*sigma**3+3*(g+T*k**2)*kh)*(1-sigma**2)**2
    result *= -k**2/4/om/sigma

    return result

    

g = 9.81
rho = 1000
gamma = 0.75
T = gamma/rho

k = np.linspace(0.1,300,4000)
h = np.linspace(0.0001,0.2,4000)

sigma = np.tanh(k*h)


kh = k*h

xx,yy = np.meshgrid(kh,k/np.sqrt(g/T))

kk, hh = np.meshgrid(k,h)

#plt.imshow(np.log(np.abs(np.tanh(xx)**2 + yy**2 * (np.tanh(xx) ** 2 - 3))))

#plt.imshow(np.log(np.abs(mu(kk,hh))))
fig, ax = plt.subplots(figsize=[12,9])

cs1 = ax.contour(kk*hh,kk/np.sqrt(g/T),mu(kk,hh)*nu(kk,hh),[0],colors="r")
