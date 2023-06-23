
import time
## Setup

#N = 3000
#L = 400
#x = np.linspace(0,L,N)
#dx = L/N


## Physics
h = 3
#A = 0.8
#x0 = 50

#delta = np.sqrt(4*h**3/A/3)
#lam = 3/2/h**3

#u=A/np.cosh((x-100)/delta)**2
#u+=0.5/np.cosh((x-50)/delta)**2

#u = np.load("2_sol.npy")

#h = 3
tst = mat73.loadmat("data/canal/signal_soliton_cam/CumSum_one_soliton_A5mm_h3cm_20s_Signal_Tot.mat")
sig = tst["SignalT"]
sig = np.array(sig["Tot"])

u0 = sig[80]
x = tst["XX"]
N = len(u)
L = x[-1]
dx = x[1]-x[0]

FS = N/L 
sos = signal.butter(4,5,'low',fs=int(FS),output="sos")
u = signal.sosfilt(sos,u0)


u = u[3800:9200]

lam = 3/2/h**3

## Potential and energies

E_min = -lam*(np.max(u)-np.mean(u))

E_max =  E_min + 4*np.abs(E_min) # It can also be set to the Nyquist frequency np.pi/dx

E = np.linspace(E_min,E_max,2000)


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

M = np.zeros((len(E),2,2))
S = np.zeros(N)


## The M matrix for a given energy E

start = time.time()

def B(q):

    Mi = np.identity(2)    
    
    for i in range(len(q)):

        
        
        if q[i]<0:
            sq =np.sqrt(-q[i])
            
            cos = np.cosh(dx*sq)
            sh = np.sinh(dx*sq)
            
            exp12 = sh/sq
            exp21 = sh*sq
            
        elif q[i]>0:
            sq =np.sqrt(q[i])
            
            cos = np.cos(dx*sq)
            sin = np.sin(dx*sq)

            exp12 = sin/sq
            exp21 = -sin*sq

        else:

            cos = 1
            exp12 = 1
            exp21 = 1
            
        exp = np.array([[cos,exp12],[exp21,cos]])

        #Mi1 = Mi
        Mi = np.dot(exp,Mi)
        
        #s += np.abs(np.sign(Mi[0,0])-np.sign(Mi1[0,0]))/2

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

end = time.time()


plt.title("Trace")
plt.plot(E,trace)
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("Trace")

print("The time of execution of above program is :", end-start)
