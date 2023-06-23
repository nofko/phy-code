import numpy as np
import matplotlib.pyplot as plt
import copy

# ===========================================================================
# Boundary conditions :
# Dirichlet at ±L/2 i.e. index 0 and -2
# Dirichlet at x=0  i.e. index N//2  (N is supposed even)
# ---------------------------------------------------------------------------
#       y(x = +L/2) = y(x = -L/2) =  CstDir
#       y(x = 0) = A
# Unknowns :
# ---------------------------------------------------------------------------
# - values at nodes ---> N+1 points
# - value of celerity c that must be related to the amplitude value at x=0
#
# ===========================================================================

######################### build_matrices ###########################
def build_matrices(N, dxs):
    """ Build the three matrices with proper boundary conditions
    ZZ = full unknowns = [u_0, ... , u_{N}, c]
    !!! LAST ELEMENT is c !!!
    """
    # SPATIAL OPERATORS hence N+1 points x_i = 1st diagonal block

    M1 = np.zeros((N+2, N+2))
    M3 = np.zeros_like(M1)
    M5 = np.zeros_like(M1)

    # LAST COLUMN is NON-ZERO !
    one_m3 = np.diag(np.ones(N-1), -3)
    one_m2 = np.diag(np.ones(N),   -2)
    one_m1 = np.diag(np.ones(N+1), -1)

    one_p3 = np.diag(np.ones(N-1), +3)
    one_p2 = np.diag(np.ones(N),   +2)
    one_p1 = np.diag(np.ones(N+1), +1)

    # LAST COLUMN is now ZERO !
    one_m3[:, -1] = 0
    one_m2[:, -1] = 0
    one_m1[:, -1] = 0

    one_p3[:, -1] = 0
    one_p2[:, -1] = 0
    one_p1[:, -1] = 0

    # GENERIC MATRICES
    M1 = (-1/2) * one_m1 + (+1/2) * one_p1 
    M3 = (-1/2) * one_m2 + one_m1 - one_p1 + (+1/2) * one_p2 
    M5 = (-1/2) * one_m3 + (+2) * one_m2 + (-5/2) * one_m1 + \
        (+5/2) * one_p1 + (-2) * one_p2 + (+1/2) * one_p3

    # LAST LINES (concerning velocity c)
    # ==================================
    M1[-1, :] = 0
    M3[-1, :] = 0
    M5[-1, :] = 0
    # BOUNDARY POINTS x_0 and x_N
    # ===========================
    M1[0,:] = 0
    M1[-2,:] = 0
    #--------------------
    M3[0,:] = 0
    M3[-2,:] = 0
    #--------------------
    M5[0,:] = 0
    M5[-2,:] = 0

    # NEIGHBORING BOUNDARY POINTS x_1 and x_{N-1}
    # ===========================================
    # M1 remains the same -> NO ALTERATION
    # ------------------------------------
    M3[1, 0] = -(+3/2) 
    M3[1, 1] = -(-5)   
    M3[1, 2] = -(+6)   
    M3[1, 3] = -(-3)   
    M3[1, 4] = -(+1/2) 
    #
    M3[-3,-6] = +(+3/2)         
    M3[-3,-5] = +(-5)            
    M3[-3,-4] = +(+6)           
    M3[-3,-3] = +(-3)              
    M3[-3,-2] = +(+1/2)                  
    # -------------------
    M5[1, 0] = -(+5/2 )
    M5[1, 1] = -(-28/2)
    M5[1, 2] = -(+65/2)
    M5[1, 3] = -(-80/2)
    M5[1, 4] = -(+55/2)
    M5[1, 5] = -(-20/2)
    M5[1, 6] = -(+3/2 )
    M5[-3,-8] =  +(+5/2 )
    M5[-3,-7] =  +(-28/2)
    M5[-3,-6] =  +(+65/2)
    M5[-3,-5] =  +(-80/2)
    M5[-3,-4] =  +(+55/2)
    M5[-3,-3] =  +(-20/2)
    M5[-3,-2] =  +(+3/2 )
    # NEIGHBORING BOUNDARY POINTS x_2 and x_{N-2}
    # -------------------------------------------
    M5[2, 0] =  -(+3/2) 
    M5[2, 1] =  -(-16/2)
    M5[2, 2] =  -(+35/2)
    M5[2, 3] =  -(-40/2)
    M5[2, 4] =  -(+25/2)
    M5[2, 5] =  -(-8/2) 
    M5[2, 6] =  -(+1/2) 
    M5[-4,-8] =  +(+3/2) 
    M5[-4,-7] =  +(-16/2)
    M5[-4,-6] =  +(+35/2)
    M5[-4,-5] =  +(-40/2)
    M5[-4,-4] =  +(+25/2)
    M5[-4,-3] =  +(-8/2) 
    M5[-4,-2] =  +(+1/2) 

    # DIMENSIONAL DERIVATIVES
    M1 /= dxs
    M3 /= dxs**3
    M5 /= dxs**5
    
    # # DIRICHLET CONDITIONS
    # M1[0,0] = 1
    # M1[-2,-2] = 1
    # M3[0,0] = 1
    # M5[0,0] = 1
    return M1, M3, M5
######################## end of build_matrices ##############################


################################# G(Z) ######################################
def G(Z):
    """ KdV + BC
    """
    GG = np.zeros_like(Z)
    # KdV equation regardless of lines 0, -2, -1:
    GG = -Z[-1]*(M1@Z) + adv_coeff * Z*(M1@Z) + \
         s*(M3@Z) + epsilon2*(M5@Z)
    # Dirichlet BC:
    GG[0]  = Z[0]  - CstDir
    GG[-2] = Z[-2] - CstDir
    # central BC:
    GG[-1]  = Z[N//2] - Au
    return GG
############################### end of G(Z) #################################


################################## dG(Z) ####################################
def dG(Z):
    """ computes the MATRIX full_result of operator dG(Z_n) \cdot zeta
    """
    full_result = np.zeros((len(Z), len(Z)))
    # spatial part
    block1  = -Z[-1] * M1  
    block1 += adv_coeff * (np.diag(M1@Z) + np.diag(Z)@M1) 
    block1 += s * M3 + epsilon2 * M5
    full_result[:-1, :-1] = block1[:-1, :-1]
    # velocity part (LAST COLUMN)
    full_result[:-1, -1]  = -(M1@Z)[:-1]
    #
    # ENFORCING DIRICHLET BC so that dG is INVERTIBLE
    # @ line N+2
    full_result[-1,:] = 0
    full_result[-1,int(N/2)] = 1
    # @ line 0
    full_result[0, :-1]  = 0
    full_result[0, 0]  = 1
    # @ line N+1
    full_result[-2, :-1] = 0
    full_result[-2,-2] = 1
    return full_result
############################### end of dG(Z) ################################


def check_MATRICES(M1, M3, M5):
    # **********   CHECK MATRICES    ***********
    print("2 * M1 = ")
    print(2*M1)
    print("2 * M3 = ")
    print(2*M3)
    print("2 * M5 = ")
    print(2*M5)
    input("PAUSE ! Press ENTER")
    return None



#****************************************************************************
#****************************************************************************
#
#*****************************   MAIN PROGRAM   *****************************
#
#****************************************************************************
#****************************************************************************

#############----   GLOBAL VARIABLES   ----##############
#

# length of the domain    
L = 160
# number of sub-intervals (!!! N MUST BE EVEN !!!)
N = 2000
if N%2 == 1:
    print("WARNING : Adding one supplementary segment")
    N+=1

# spatial discretization
dx = L/N

tau = 0.3

XX  = np.linspace(-L/2,L/2,N+1)
XXs = copy.deepcopy(XX) #*np.sqrt(6/(1-3*tau))

dxs = XXs[1]-XXs[0]

# PARAMETERS of KdV equation: Eq. (9)
# -c u' + adv_coeff * u u' + s u''' + eps^2 u''''' = 0
# =================================================================
adv_coeff = +1   # convective term  (-6 for standard elevation KdV)
s = 1.          # 3rd derivative
epsilon2 = -30   # 5th derivative

# AMPLITUDE at x=0 and Dirichlet BC
# =======================================
Au = -0.4
CstDir = -0.1  #1e-5 * np.sign(Au)

iShow = True

M1, M3, M5 = build_matrices(N, dxs)
#check_MATRICES(M1, M3, M5)

############################### INITIAL GUESS ##########################
Z0 = np.zeros(N+2)

Yguess = Au / np.cosh(XXs*np.sqrt(abs(Au/12)))**2
Z0[:-1] = Yguess

# Initial velocity guess : 
Z0[-1] = -0.05

tol = 1.e-6
itermax = 20


Z = copy.deepcopy(Z0)
dZ = np.zeros_like(Z)

Z[0]  = Yguess[0]
Z[-2] = Z[0]

iteration = 0
print('Error value before entering Newton loop:')
err = np.linalg.norm(G(Z))
print(f' Iteration = {iteration:4d}\
    Error = {err:e}    Celerity = {Z[-1]:e}')
print()

print(" --> NEWTON loop <-- ")
###############################  NEWTON  ###############################
while (err > tol) and (iteration < itermax):
    RHS = G(Z)
    # DIRICHLET BC  
    RHS[0]   = 0.
    RHS[-2]  = 0.
    RHS[int(N/2)] = 0.
    # print('{0:e} {1:e} {2:e}'.format(RHS[0], RHS[-2], RHS[N//2]))
    dZ = np.linalg.solve(dG(Z), RHS)
    # print('{0:e} {1:e} {2:e}'.format(dZ[0], dZ[-2], dZ[N//2]))
    Z -= dZ
    # REDUNDANT Dirichlet BC enforcement
    Z[0]   = CstDir
    Z[-2]  = CstDir
    Z[N//2] = Au
    #
    err = np.linalg.norm(G(Z))
    iteration += 1
    print(f' Iteration = {iteration:4d}\
    Error = {err:e}    Celerity = {Z[-1]:e}')
    

print()
print("Values at index 1 and -3:")
print(Z[1], "   ", Z[-3])
print()

if iShow:
    #plt.close('all')
    plt.figure(1)
    plt.plot(XX,Z0[:-1], '--', label="Initial guess")
    plt.plot(XX,Z[:-1]       , label="Newton solution")
    plt.legend(loc="lower left",fontsize=14)
    plt.plot(XX[0], Z[0], 'sk', XX[-1], Z[-2], 'sk', XX[N//2], Z[N//2], '*r')
    plt.grid()
    plt.show(block=False)
    #input("Press ENTER")


