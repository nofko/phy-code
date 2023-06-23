#Key Lengths

pxcm = 19/1525


# 1_3_f10 data

Ro = 7.87
Ri = 5.96
W = Ro-Ri

# n solitons data
Ro =1291/2*pxcm
Ri =935/2*pxcm
W =Ro-Ri

# spectrum 7
pxcm = 19/1506
Ro =8.2
Ri =5.63
W =Ro-Ri


# sine
#pxcm = 19/1604
#Ro =8.25
#Ri =5.53
#W =Ro-Ri

# video 14 of small torus
#Ro = 5.61
#W = 3.73
#pxm = 16/1528

Rc = 7

L = 2*Ro*np.pi

## Setup

N = 1000
Ns = 3801

x = np.linspace(0,L,N)
dx = L/Ns

## Physics

mu = 3

Bo = 0.85**2*Rc**4/W**2/Ro**4
deltaBo = 1/23.6 - Bo

lam = mu*Rc**4/3/W**3/Ro**4/abs(deltaBo)

#u = -np.load("data/1_3_f10.npy")*pxcm    #taken from solitons spectrum video 9
#u = -np.load("data/1_solit_2.npy")*pxcm
u1 = -np.load("data/8_soliton.npy")*pxcm
u = -np.load("data/7_spec_550.npy")*pxcm
#u = np.load("data/14_small_torus.npy")*pxcm

u = u-np.mean(u)

## Potential and energies

E_min = -lam*(np.max(u)-np.mean(u))

E_max =  E_min + 6*np.abs(E_min) # It can also be set to the Nyquist frequency np.pi/dx

E = np.linspace(E_min,E_max,N)


## Useful plotting function

def linlog(y):

    N = len(y)
    out = np.zeros(N)

    for i in range(N):
        if np.abs(y[i])<=1:
            out[i] = y[i]

        else:
            out[i] = np.sign(y[i])*(1+np.log(np.abs(y[i])))

    return out


## Calculating M

M = np.zeros((N,2,2))
S = np.zeros(N)


## The M matrix for a given energy E

def B(q):

    Mi = np.identity(2)    
    
    for i in range(len(q)):
    
        if q[i]<0:
            sq =np.sqrt(-q[i])
            
            cos = np.cosh(dx*sq)
            sh = np.sinh(dx*sq)
            
            exp12 = sh/sq
            exp21 = sh*sq
            
        else:
            sq =np.sqrt(q[i])
            
            cos = np.cos(dx*sq)
            sin = np.sin(dx*sq)

            exp12 = sin/sq
            exp21 = -sin*sq

            
        exp = np.array([[cos,exp12],[exp21,cos]])

        Mi = np.dot(exp,Mi)
        
    return Mi


## Calculating the matrix M for the range of values between Emin and Emax

for j in range(len(E)):
    print(j/len(E))
    q = lam*u+E[j]
    
    M[j] = B(q)
    

trM = np.trace(M,axis1=1,axis2=2)

trace = linlog(trM/2)
M11 = linlog(M[:,0,0])
M12 = linlog(M[:,0,1])


plt.figure(2)
plt.title("Trace")
plt.plot(E,trace)
plt.plot(E,S)
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("Trace")

