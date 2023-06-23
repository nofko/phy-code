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



N = 10          # NUMBER OF DETERIMANT ROWS

kappa = 0.0    # DISSIPATION

ll2 = np.arange(-N,N)
rr2 = np.arange(-N+1,N+1)

ll4 = np.concatenate((np.flip(np.arange(1,N)),np.arange(1,N+1)))
rr4 = np.concatenate((np.flip(np.arange(1,N+1)),np.arange(1,N)))


############################       FUNCTIONS       ############################


def DD(alpha,beta):

    """
    DETERMINANT FOR 2PI PERIODIC FUNCTIONS FOR GIVEN VALUES OF ALPHA AND BETA 
    """
    
    M = len(alpha)

    DT = np.zeros((M,M))
    
    for i in range(M):

        for j in range(M):
    
            gamma_l = beta[i,j]/(alpha[i,j]-ll2**2)/2
            gamma_r = beta[i,j]/(alpha[i,j]-rr2**2)/2
    
            D = np.diag(gamma_r,-1)+np.diag(gamma_l,1)+np.identity(2*(N)+1)

            DT[i,j] = np.linalg.det(D)

    return DT



def KK(alpha,beta):
    
    """
    DETERMINANT FOR 4PI PERIODIC FUNCTIONS FOR GIVEN VALUES OF ALPHA AND BETA 
    """

    M = len(alpha)

    DT = np.zeros((M,M),dtype=complex)
    
    for i in range(M):

        for j in range(M):
    
            gamma_l = beta[i,j]/(alpha[i,j]-(2*ll4-1)**2/4)/2
            gamma_r = beta[i,j]/(alpha[i,j]-(2*rr4-1)**2/4)/2
    
            D = np.diag(gamma_r,-1)+np.diag(gamma_l,1)+np.identity(2*(N))

            DT[i,j] = np.linalg.det(D)

    return DT



############################       PLOTTING       ############################


a = np.linspace(0.1,6,500)
b = np.linspace(0.1,8,500)


aa,bb = np.meshgrid(a,b)

fig, ax = plt.subplots(figsize=[8,6])

ax.contour(aa,bb,DD(aa**2/4,bb),[0],colors="b")

ax.contour(aa,bb,KK(aa**2/4,bb),[0],colors="r")


#TT = np.log(np.abs(KK(aa,bb)))+np.log(np.abs(DD(aa,bb)))

ar = abs((a[-1]-a[0])/(b[-1]-b[0]))

#image = ax.imshow(TT[::-1,:],extent=[a[0],a[-1],b[0],b[-1]],aspect=ar,cmap="turbo")


#plt.colorbar(mappable=image)


plt.show()
