from scipy.linalg import toeplitz
import numpy.linalg as la
## Setup

N = 200*4
L = 100
x = np.linspace(-L/2,L/2,N)
dx = L/N

Nf = int(N/4)

k0 = 2*np.pi/L

## Soliton

A = 1
x0 = 50


mu=1; B=-1/np.sqrt(1+16*mu/3);
u=np.sqrt(-4*B*mu/(B+np.cosh(2*np.sqrt(mu)*x)));
G1=-mu-u**2+u**4; G2=-mu-3*u**2+5*u**4;

c1 = np.zeros(2*Nf+1,dtype="complex")
c2 = np.zeros(2*Nf+1,dtype="complex")


for i in range(-Nf,Nf+1):

    c1[i+Nf] = dx*sum(G1*np.exp(-k0*i*x*1j))/L
    c2[i+Nf] = dx*sum(G2*np.exp(-k0*i*x*1j))/L

D = -(k0*np.diag(np.arange(-Nf,Nf+1)))**2;

row11 = np.concatenate((c1[Nf:],np.zeros(Nf)))
c1t = c1[:Nf+1]
row12 = np.concatenate((c1t[::-1],np.zeros(Nf)))

row21 = np.concatenate((c2[Nf:],np.zeros(Nf)))
c2t = c2[:Nf+1]
row22 = np.concatenate((c2t[::-1],np.zeros(Nf)))

C1 = toeplitz(row11,row12)

C2 = toeplitz(row21,row22)

C0 = np.zeros((2*Nf+1,2*Nf+1))

H = np.block([[C0,D+C1],[D+C2,C0]])

wc, vc = la.eig(H)

plt.plot(wc.real,wc.imag,"o")
