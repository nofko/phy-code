import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ==========================================================================##
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                           MATHIEU STABILITY CURVE                         ##
##                                                                           ##
## ==========================================================================##



############################          SETUP        ############################



N = 10    # NUMBER OF DETERIMANT ROWS

ll2 = np.arange(-N,N)
rr2 = np.arange(-N+1,N+1)

ll4 = np.concatenate((np.flip(np.arange(1,N)),np.arange(1,N+1)))
rr4 = np.concatenate((np.flip(np.arange(1,N+1)),np.arange(1,N)))

kappa = 0.0

S = 1

############################       FUNCTIONS       ############################


def E1_p1(b):
    
    # FIRST FOURIER COEFFICIENT #
    # b = b-1

    f1 = 4*(np.sqrt(1-b)+2*(b-1)*np.sqrt(1+b))*np.pi
    f2 = -8*(b-1)*np.sqrt(1+b)*np.arctan(np.sqrt((1-b)/(1+b)))
    f3 = -8*(b-1)*np.sqrt(1+b)*np.arctan(np.sqrt((1+b)/(1-b)))
    
    return (f1+f2+f3)/(2*np.sqrt(1-b)*b*np.pi)

def E2_p1(b):

    # SECOND FOURIER COEFFICIENT #

    b2 = b**2
    
    return 2*(-2+b2+2*np.sqrt(1-b2))/b2


def DD(alpha,beta):

    """
    DETERMINANT FOR 2PI PERIODIC FUNCTIONS FOR GIVEN VALUES OF ALPHA AND BETA 
    """
    
    M = len(alpha)

    DT = np.zeros((M,M),dtype="complex")
    
    for i in range(M):

        for j in range(M):

            A = 1/alpha[i,j]**2
            B = beta[i,j]
            C = -1/alpha[i,j]*0.4

            e1 = E1_p1(B)
            e2 = E2_p1(B)

            delta_l = C*e2/(A-ll2**2+1j*kappa*ll2)/1j
            delta_r = -C*e2/(A-rr2**2+1j*kappa*rr2)/1j
            
            gamma_l = C*e1/(A-ll2**2+1j*kappa*ll2)
            gamma_r = C*e1/(A-rr2**2+1j*kappa*rr2)
    
            D = np.diag(delta_r[1:],-2)+np.diag(gamma_r,-1)+np.diag(gamma_l,1)+np.identity(2*(N)+1)+np.diag(delta_l[:-1],2)
            
            DT[i,j] = np.linalg.det(D)

    return np.abs(DT)



def KK(alpha,beta):
    
    """
    DETERMINANT FOR 4PI PERIODIC FUNCTIONS FOR GIVEN VALUES OF ALPHA AND BETA 
    """

    M = len(alpha)

    DT = np.zeros((M,M),dtype="complex")
    
    for i in range(M):

        for j in range(M):

            A = 1/alpha[i,j]**2
            B = beta[i,j]
            C = -1/alpha[i,j]*0.4
            
            e1 = E1_p1(B)
            e2 = E2_p1(B)
            
            delta_l = -C*e2/(A-(2*ll4-1)**2/4+1j*kappa*(2*ll4-1)/2)/1j
            delta_r = C*e2/(A-(2*rr4-1)**2/4+1j*kappa*(2*rr4-1)/2)/1j
    
            gamma_l = C*e1/(A-(2*ll4-1)**2/4+1j*kappa*(2*ll4-1)/2)
            gamma_r = C*e1/(A-(2*rr4-1)**2/4+1j*kappa*(2*rr4-1)/2)
    
            D = np.diag(delta_r[:-1],2)+np.diag(gamma_r,1)+np.diag(gamma_l,-1)+np.identity(2*(N))+np.diag(delta_l[1:],-2)
            DT[i,j] = np.linalg.det(D)

    return np.abs(DT)


############################       PLOTTING       ############################



a = np.linspace(0.0,4,500)
b = np.linspace(0.0,1,500)


aa,bb = np.meshgrid(a,b)

TT = np.log(np.abs(KK(aa,bb))) + np.log(np.abs(DD(aa,bb)))

fig, ax = plt.subplots(figsize=[8,6])

#ax.contour(aa,bb,DD(aa,bb),[0],colors="b")

ar = abs((a[-1]-a[0])/(b[-1]-b[0]))

image = ax.imshow(TT[::-1,:],extent=[a[0],a[-1],b[0],b[-1]],aspect=ar,cmap="turbo")

c = 0.4*2

def E1(b):
    
    # FIRST FOURIER COEFFICIENT #

    return -2*b/(1-b**2)**(1.5) 

def E2(b):

    # SECOND FOURIER COEFFICIENT #

    b2 = b**2
    
    return 8/b2+4/b2*(3*b2-2)/(1-b2)**(1.5) 


def T1(x,y):

    e1 = E1_p1(y)

    gamma = c*np.sqrt(x)
    
    return (x-0.25)**2+0.25*kappa**2-0.25*gamma**2*e1**2


def T2(x,y):

    e1 = E1_p1(y)
    e2 = E2_p1(y)

    gamma = c*np.sqrt(x)
    
    return x*(4*(x-1)**2+4*kappa**2-gamma**2*e2**2)-2*gamma**2*e1**2*(x-1)

ax.contour(aa,bb,T1(1/aa**2,bb),[0],colors="r",linewidths=0.4)

ax.contour(aa,bb,T2(1/aa**2,bb),[0],colors="b",linewidths=0.4)

plt.colorbar(mappable=image)


plt.show()
