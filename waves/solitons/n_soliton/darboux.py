
## Grid setup

N = 5000

x = np.linspace(-10,10,N)

t = np.linspace(-1,1,N)

xx, tt = np.meshgrid(x,t)


def psi(lam):

    omega = np.sqrt(lam)

    eta = omega*(xx-tt*(4*lam))

    return np.array([[np.exp(eta),np.exp(-eta)],[omega*np.exp(eta),-omega*np.exp(-eta)]])


u_0 = 0.0

nsol = 15

lam = [3.6,3.5,3,2.5,2,1]

lam = np.random.rand(nsol)*7+1

mus = np.exp(np.random.rand(nsol)*4)

lam = sorted(lam,reverse=True)

if nsol%2 == 0:
    mu = -1
else:
    mu = 1

psi0 = psi(lam[0])
psi1 = psi(lam[1])

sigma=(mu*psi0[1,1]+psi0[1,0])/(psi0[0,0]+mu*psi0[0,1])

u = 2*lam[0]-u_0-2*sigma**2

sigmas=np.zeros((nsol,N,N))

sigmas[0] = sigma

for i in range(1,nsol):

    if nsol%2 == 0:
        mu = (-1)**(i+1)*mus[i]
    else:
        mu = (-1)**(i)*mus[i]

    D = np.array([[-sigma,np.ones_like(sigma)],[sigma**2+lam[i]-lam[i-1],-sigma]])

    psi_t = psi(lam[i])

    Dn = D
    
    for j in range(0,i-1):

        n=i-2-j
        print(n)
        Dt = np.array([[-sigmas[n],np.ones_like(sigmas[n])],[sigmas[n]**2+lam[i]-lam[n],-sigmas[n]]])

        Dn = np.einsum('ijmn,jkmn->ikmn', Dn, Dt)

    psi0 = np.einsum('ijmn,jkmn->ikmn', Dn, psi_t)
    
    sigma=(psi0[1,0]+mu*psi0[1,1])/(psi0[0,0]+mu*psi0[0,1])

    sigmas[i] = sigma
    
    u = 2*lam[i]-u-2*sigma**2



#om1 = np.sqrt(lam[0])
#om2 = np.sqrt(lam[1])

#et1 = om1*(xx-tt*(4*lam[0]))
#et2 = om2*(xx-tt*(4*lam[1]))

#u = -2*(lam[1]-lam[0])*((om1/np.sinh(et1))**2+(om2/np.cosh(et2))**2)/(om1/np.tanh(et1)-om2*np.tanh(et2))**2

plt.imshow(np.where(np.abs(u)<30,u,0),cmap="inferno",vmin="0",vmax="14",extent=[-1,1,-10,10],aspect=0.1)

#plt.pcolormesh(xx,tt,np.where(np.abs(u)<30,u,0),cmap="inferno",vmin="0",vmax="14",shading="flat")
plt.xlabel("x")
plt.ylabel("t")
plt.savefig("15_solitons.png",bbox_inches="tight")
