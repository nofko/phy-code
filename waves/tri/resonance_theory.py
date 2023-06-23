import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import fsolve



## ===========================================================================
## Calculation of resonant frequencies for the first three branches of the  ##
## torus:  the varicose, sinuous, sloshing.                                 ##
## The frequency of the forcing is imposed onto the resonance conditions    ##
## ===========================================================================


####### PHYSICS #######

g = np.sin(4.5*np.pi/180)*9.81

sigma = 0.055     # SURFACE TENSION
rho = 1000        # FLUID DENSITY

rc = 0.07         # CENTRAL RADIUS IN M
ro = 0.0785       # TORUS OUTER RADIUS
w = 2*(ro-rc)     # TORUS WIDTH
ri =ro-w          # TORUS INNER RADIUS

f = np.linspace(7.9,16,100)           # FORCING FREQUENCY IN HZ
q = 1             # DISPERSION RELATION WIDTH

om0 = 44.5  # R=7.9  # CUTOFF FREQUENCY OF SLOSHING BRANCH IN RAD/S
#om0 = 34   # R=8.1


####### DISPERSION RELATIONS #######

def varicose(k,w,ro):

    return np.sqrt((g*k/ro+sigma*k**3/rho/ro**3)*np.tanh(k*w*ro/2/rc**2))

def sinuous(k,w,ri):

    om0 = np.pi/2*np.sqrt(9.81*np.sin(9*np.pi/180)/2/w)

    return np.sqrt(om0**2+np.abs(sigma/rho*(k/ri)**3))

def sloshing(k,om0):
    
    return np.sqrt(om0**2+g*k**2/rc)


####### FUNCTIONS TO SOLVE ########

def func(k,om,w,ro):

    return varicose(k,w,ro)-om

def fsin(k,om,w,ri):

    return sinuous(k,w,ri)-om

def fslosh(k,om,w,om0):

    return sloshing(k,om0)-om

####### SCALAR SOLUTIONS #######

def var_sc(om,w,ro):
    

    s = fsolve(func,10,args=(om,w,ro))
        
    return s

def sin_sc(om,w,ri):

    s = fsolve(fsin,10,args=(om,w,ri))
        
    return s

def slo_sc(om,w,om0):

    s = fsolve(fslosh,10,args=(om,w,om0))
        
    return s

####### VECTOR SOLUTIONS #######

def var_v(om,w,ro):

    sol = np.zeros(len(om))

    for i in range(len(om)):

        s = fsolve(func,10,args=(om[i],w,ro))
            
        sol[i] = s
        
    return sol

def sin_v(om,w,ri):

    sol = np.zeros(len(om))

    for i in range(len(om)):

        s = fsolve(fsin,20,args=(om[i],w,ri))

        sol[i] = s
        
    return sol

def slo_v(om,w,om0):

    sol = np.zeros(len(om))

    for i in range(len(om)):
            
        s = fsolve(fslosh,20,args=(om[i],w,om0))
        
        sol[i] = s
        
    return sol

####### MATRIX SOLUTIONS #######

def var_m(om,w,ro):

    sol = np.zeros((len(om),len(om[0])))

    for i in range(len(om)):

        for j in range(len(om[0])):
            
            s = fsolve(func,10,args=(om[i][j],w,ro))
            
            sol[i][j] = s
        
    return sol

def sin_m(om,w,ri):

    sol = np.zeros((len(om),len(om[0])))

    for i in range(len(om)):

        for j in range(len(om[0])):
            
            s = fsolve(fsin,20,args=(om[i][j],w,ri))

            sol[i][j] = s
        
    return sol


def slo_m(om,w,om0):

    sol = np.zeros((len(om),len(om[0])))

    for i in range(len(om)):

        for j in range(len(om[0])):
   
            s = fsolve(fslosh,20,args=(om[i][j],w,om0))

            sol[i][j] = s
        
    return sol


####### PLOTTING #######

N = 100

om1 = np.linspace(0,14,N)*2*np.pi
om2 = np.linspace(0,14,N)*2*np.pi

k = np.linspace(0,30,N)

xx,yy = np.meshgrid(om1,om2)

fig, ax = plt.subplots(figsize=[8,6])

ax.set_xlabel('$f_1$',fontsize=35)
ax.set_ylabel('$f_2$',fontsize=35)
ax.tick_params(labelsize=25)



### FRACTAL - NOT UNDERSTOOD ###

#cs = ax.contour(xx/2/np.pi,yy/2/np.pi,-var_m(xx,w,ro)-sin_m(yy,w,ri)-var_m(xx+yy,w,ro),[0],colors="r")

### QUASI RESONANCE ###

#cs = ax.contour(xx/2/np.pi,yy/2/np.pi,var_m(xx,w,ro)+sin_m(yy,w,ri)-var_sc(2*np.pi*f,w,ro)+q,[0],colors="r")

#cs = ax.contour(xx/2/np.pi,yy/2/np.pi,var_m(xx,w,ro)+sin_m(yy,w,ri)-var_sc(2*np.pi*f,w,ro)-q,[0],colors="b")

#cs = ax.contour(xx/2/np.pi,yy/2/np.pi,xx+yy-2*np.pi*f,[0],colors="g")

### SLOSHING = SINUOUS + VARICOSE ###

# cs = ax.contour(xx/2/np.pi,yy/2/np.pi,var_m(xx,w,ro)-sin_m(yy,w,ri)-slo_sc(2*np.pi*f,w,om0),[0],colors="b")

# cs = ax.contour(xx/2/np.pi,yy/2/np.pi,xx+yy-2*np.pi*f,[0],colors="r")

#cs = ax.contour(xx/2/np.pi,yy/2/np.pi,var_m(xx,w,ro)-sin_m(yy,w,ri)-slo_sc(2*np.pi*f,w,om0),[0],colors="b")

### VARICOSE = SINUOUS + VARICOSE ###

#cs = ax.contour(xx/2/np.pi,yy/2/np.pi,var_m(xx,w,ro)+sin_m(yy,w,ri)-var_sc(2*np.pi*f,w,ro),[0],colors="b")

#cs = ax.contour(xx/2/np.pi,yy/2/np.pi,xx+yy-2*np.pi*f,[0],colors="r")

#cs = ax.contour(xx/2/np.pi,yy/2/np.pi,var_m(xx,w,ro)-sin_m(yy,w,ri)-slo_sc(2*np.pi*f,w,om0),[0],colors="b")


### SLOSHING = VARICOSE + VARICOSE ###

cs = ax.contour(xx/2/np.pi,yy/2/np.pi,var_m(xx,w,ro)-var_m(yy,w,ro)-slo_m(xx+yy,w,om0),[0],colors="b")

#cs = ax.contour(xx/2/np.pi,yy/2/np.pi,xx+yy-2*np.pi*f,[0],colors="r")

#cs = ax.contour(xx/2/np.pi,yy/2/np.pi,var_m(xx,w,ro)-sin_m(yy,w,ri)-slo_sc(2*np.pi*f,w,om0),[0],colors="b")


### VAR+VAR=VAR ###


# om3 = 80
# om4 = 60

# cs = ax.contour(xx/2/np.pi,yy/2/np.pi,var_m(xx,w,ro)+var_m(yy,w,ro)+var_sc(om3,w,ro)-var_sc(om4,w,ro)-var_m(xx+yy+om3-om4,w,ro),[0],colors="b")


### DISPERSION RELATION CHECK ###


#plt.plot(om,var_v(om,w,ro))

#plt.plot(om,sin_v(om,w,ri))

#plt.plot(om,slo_v(om,w,om0))

#plt.plot(sinuous(k,w,ri),k)

#plt.plot(varicose(k,w,ro),k)

#plt.plot(sloshing(k,om0),k,"--")



####

#plt.pcolormesh(xx/2/np.pi,yy/2/np.pi,np.log(abs(-var_m(xx,w,ro)+var_m(yy,w,ro)+var_m(xx+yy,w,ro))),rasterized=True)

plt.show()
