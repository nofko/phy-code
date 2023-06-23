
N = 1000

L = 3

x = np.linspace(0,L,N)

dx = x[1]-x[0]

epsilon = 3

u0 = 1.8

u = np.array([u0,0,-u0*np.sqrt(2-u0),0])

sol = np.zeros((N,4))

def rhs(PP):

    x1_dot = PP[1]
    x2_dot = PP[2]
    x3_dot = PP[3]
    x4_dot = -PP[0]+0.75*PP[0]**2+epsilon*PP[2]
    
    return np.array([x1_dot,x2_dot,x3_dot,x4_dot])


for i in range(N):

    sol[i] = u

    a = dx*rhs(u)
    b = dx*rhs(u+a/2)
    c = dx*rhs(u+b/2)
    d = dx*rhs(u+c)

    u = u+(a+2*(b+c)+d)/6

plt.figure(1)
plt.plot(x,sol[:,0])


