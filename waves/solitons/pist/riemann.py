

E = [-0.02398377,-0.023983698,4.69805653e-06]

E = [-0.01595537, -0.01595289, -0.00555603, -0.00540625,  0.]






N = 20 # Degrees of freedom i.e. number of modes

phi = (np.random.rand(N)*2*np.pi).tolist()

E = np.sort(np.random.rand(N*2+1)).tolist()

M = 5   # Limits for theta function sum

p = 2000    # Number of Fourier components

L=100

k = np.arange(-int(p/2),int(p/2)+1)*2*np.pi/L

x = np.linspace(-L/2,L/2,p)


h = 5
A0 = 4
x0 = 50


delta = np.sqrt(4*h**3/A0/3)
lam = 3/2/h**3
u_o=A0/np.cosh((x)/delta)**2


def fun(x,E,i,j):

    e=(E[j+1]-E[j])/2*np.cos(x)+(E[j+1]+E[j])/2

    product = 1

    temp = E.copy()
    temp.pop(j)
    temp.pop(j)

    for k in temp:
        product *= (e-k)
    
    return 2*e**i/np.sqrt(np.abs(product))


G = np.zeros((2*N,N))
Gs = np.zeros((2*N,N))


for i in range(N):   ## Degree of freedom
    
    for j in range(0,2*N):   ## Boundaries i.e. gaps - 2N of these
        print(j)
        
        if j%2==0:
            sign = (-1)**(j/2)
        else:
            sign = (-1)**((j+1)/2)
            
        g1, err = quad(fun,0,np.pi,args=(E,i,j))

        #g2, err = quad(fun,0,np.pi,args=(E,i,j+1))
        
        G[j,i] = g1
        Gs[j,i] = sign*g1


K = np.zeros((N,N))
Kp = np.zeros((N,N))

Gt = Gs.transpose()

K = Gt[:,1::2]

K = K.transpose()

for i in range(N):
    
    for j in range(N):

        Kp[i,j] = sum(Gs[:2*i+1:2,j])


B = 2*np.pi*np.matmul(Kp,np.linalg.inv(K))



def I(m):

    s = 0

    for n in range(1,len(m[0])+1):
        s+=n*m[:,n-1]
    
    return s


m = np.array([list(i) for i in itertools.product(range(-M,M+1),repeat=N)])

l = I(m)



a = np.zeros(p)
b = np.zeros(p)

ns = np.arange(-1000,1000)

for i in range(p):

    an = 0
    bn = 0
    
    for j in range(len(l)):
        
        if l[j]==ns[i]:

            an+=np.exp(0.5*np.dot(m[j,:],np.dot(B,m[j,:].transpose())))*np.cos(np.dot(phi,m[j,:]))
            bn+=np.exp(0.5*np.dot(m[j,:],np.dot(B,m[j,:].transpose())))*np.sin(np.dot(phi,m[j,:]))

    a[i]=an
    b[i]=-bn

    

theta = np.zeros(p)

for i in range(p):
    theta+=a[i]*np.cos(k[i]*x)+b[i]*np.sin(k[i]*x)

FL = np.log(theta)

dx = L/p

Fp = np.gradient(FL)/dx
Fpp = np.gradient(Fp)/dx


u = 2*Fpp/lam


plt.plot(x,u)
#plt.plot(x,u_o-np.mean(u_o))


sindex = []

for i in range(1,int(len(E)/2)):
    sindex.append((E[2*i]-E[2*i-1])/(E[2*i]-E[2*i-2]))
