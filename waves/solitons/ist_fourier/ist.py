from numpy import linalg as la
from scipy.linalg import toeplitz

## Setup

N = 1000
L = 100
x = np.linspace(0,L,N)
dx = L/N


## Physics
h = 5
A = 4
x0 = 50

delta = np.sqrt(4*h**3/A/3)
lam = 3/2/h**3
c0 = np.sqrt(9.81*h*100)
lam = 3/2/h**3*c0**2

cst=0

u=A/np.cosh((x-30)/delta)**2+4.6/np.cosh((x-70)/delta)**2



mat = scipy.io.loadmat("data/19-1-2022/CumSum_two_solitons_A05cm_A1cm_h5cm_L4m_Pos150_Mes1m_10s_plage_P.mat")
sig = mat["tableau"]

N0 = len(sig)

u0 = sig[:,1]*1/1.7
u0 = u0-np.mean(u0)

sos = signal.butter(4,5,'low',fs=2000,output="sos")
u = signal.sosfilt(sos,u0)
#sos2 = signal.butter(4,0.3,'hp',fs=2000,output="sos")
#u = signal.sosfilt(sos2,u)

u = u[1000:15000]
u = u[::4]
print(len(u))


N = len(u)
#N = 10000
T = 10
t = sig[:,0]
t = t[1000:15000:4]
#t = np.linspace(0,T,N)
dx = t[1]-t[0]


D2 = np.diag(np.ones(N-1),-1)+np.diag(np.ones(N-1),1)-2*np.diag(np.ones(N))
D2 /= dx**2

M = D2+np.diag(lam*u)

w, v = la.eig(M[1:-1,1:-1])

ws = np.sort(-w)

print("A: ",2/lam*(ws[0]))

print("A-h: ", 2/lam*(ws[0]-lam*cst))

fig, ax1 = plt.subplots(figsize=[12,9])



left, bottom, width, height = [0.59, 0.57, 0.3, 0.3]
ax2 = plt.axes([left, bottom, width, height])

ax1.plot(t,u,color="red",lw=2)

ax1.set_xlabel("$x$ [cm]",fontsize=45)
ax1.set_ylabel("$u$ [cm]",fontsize=45)
ax1.tick_params(labelsize=30)


ax2.plot(ws[:6],ws[:6]*0,"o",color="b")
ax2.axvline(0,ls="--",color="gray")
ax2.set_xlabel("$\lambda$ [cm$^{-2}$]",fontsize=35)

ax2.set_yticklabels([])
ax2.tick_params(labelsize=25)



## Collocation

# Nf = int(N/4)


# k0 = 2*np.pi/L;

# b = np.zeros(2*Nf+1,dtype="complex")

# for i in range(-Nf,Nf+1):

#     b[i+Nf] = dx*sum(u*np.exp(-k0*i*x*1j))/L


# D = -(k0*np.diag(np.arange(-Nf,Nf+1)))**2;


# row1 = np.concatenate((b[Nf:],np.zeros(Nf)))

# bt = b[:Nf+1]
# row2 = np.concatenate((bt[::-1],np.zeros(Nf)))

# B = toeplitz(row1,row2)

# H = D+lam*B

# wc, vc = la.eig(H)
