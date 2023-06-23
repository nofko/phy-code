from scipy.optimize import fsolve

N = 1000

L = 20

x = np.linspace(0,L,N)

dx = x[1]-x[0]

epsilon = 3

u = np.array([1.8,0,-1.8*np.sqrt(2-1.8),0])

sol = np.zeros((N,4))

def rhs(PP):

    x1_dot = PP[1]
    x2_dot = PP[2]
    x3_dot = PP[3]
    x4_dot = -PP[0]+0.75*PP[0]**2+epsilon*PP[2]
    
    return np.array([x1_dot,x2_dot,x3_dot,x4_dot])


def fun(y,y_old):

    return y_old+dx*rhs(y)-y
    

for i in range(N):

    sol[i] = u

    u = fsolve(fun,u,args=(u))


plt.figure(1)
plt.plot(sol[:200,0])
