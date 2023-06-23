import matplotlib.pyplot as plt
import numpy as np
import itertools
from scipy.io import savemat

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ==========================================================================##
## N-SOLITON GENERATION                                                      ##
##                                                                           ##
##                                HIROTA METHOD                              ##
##                                                                           ##
## ==========================================================================##


############################         SETUP         ############################


N = 1000
Nt= 1000

L = 400
T = 5

x = np.linspace(-L/2,L/2,N)
dx = L/N

t = np.linspace(-T/2,T/2,Nt)

xx,tt = np.meshgrid(x,t)

############################        PHYSICS       ############################


M = 1   ## soliton number

h = 5.5

lam = 3/2/h**3

c0 = np.sqrt(9.81*h/100)*100

eta0 = np.array([0.4])
#eta0 = np.random.rand(M)*7+1

fname = "h55cm_A04cm.mat"

phi = np.array([0])
#phi = np.random.rand(M)*2-1

k = np.sqrt(2*eta0*lam)

A = np.zeros((M,M))




############################       GENERATION    ############################


for i in range(M):

    for j in range(M):

        if k[i]==k[j]:
            A[i,j] = 0

        else:
            A[i,j] = np.log(((k[i]-k[j])/(k[i]+k[j]))**2)

eta = np.zeros((M,Nt,N))

for i in range(M):
    eta[i]=k[i]*(xx-k[i]**2*tt*h**2*c0/6-c0*tt)+phi[i]

combinations = [np.asarray(i) for i in itertools.product(range(0,2),repeat=M)]

F = np.zeros((Nt,N))

for nu in combinations:

    e = np.zeros((Nt,N))

    for j in range(len(nu)):
        
        e += nu[j]*eta[j]+0.5*sum([A[l][j]*nu[l] for l in range(j+1,len(nu))])*nu[j]

    F+=np.exp(e)

FL = np.log(F)

Fp = np.gradient(FL,axis=1)/dx
Fpp = np.gradient(Fp,axis=1)/dx

u = 2*Fpp/lam


#plt.plot(A0/np.cosh((x-x0)/delta)**2)
plt.imshow(u,cmap="inferno",extent=[-T/2,T/2,-L/2,L/2],aspect=0.01)
plt.xlabel("x")
plt.ylabel("t")

plt.figure()
plt.plot(u[0])

mdic = {"signal": u}
savemat("single/"+fname, mdic)


plt.show()
