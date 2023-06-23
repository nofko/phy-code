import matplotlib.pyplot as plt
import numpy as np


from scipy.optimize import fsolve

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


N = 100
CP_lighter_blue = "#94d2bd"

colors = ["#390099","#9e0059","#ff0054","#ff5400","#f5cb5c","#ffbd00"]
light_gray = "b5b5b5"
dark_gray = "#000723"

####### PHYSICS #######

g = np.sin(4.5*np.pi/180)*9.81

sigma = 0.055
rho = 1000

rc = 0.07
ro = 0.0786
w = 2*(ro-rc)
ri =ro-w



om0 = 44.5  # R=7.9
#om0 = 34.5   # R=8.1


####### DISPERSION RELATIONS #######

def varicose(k,w,ro):

    return np.sqrt((g*k/ro+sigma*k**3/rho/ro**3)*np.tanh(k*w*ro/2/rc**2))

def sinuous(k,w,ri):

    om0 = np.pi/2*np.sqrt(9.81*np.sin(9*np.pi/180)/2/w)

    return np.sqrt(om0**2+np.abs(sigma/rho*(k/ri)**3))

def sloshing(k,om0):
    
    return np.sqrt(om0**2+g*k**2/rc)

####### GROUP VELOCITY #######

def group_slo(k,om0):

    return g*k/sloshing(k,om0)/rc

def group_var(k,w,ro):

    a = (g/ro+3*sigma*k**2/rho/ro**3)*np.tanh(k*w*ro/2/rc**2)

    b = (g*k/ro+sigma*k**3/rho/ro**3)/np.cosh(k*w*ro/2/rc**2)**2*w/2*ro/rc**2
    
    return 0.5/varicose(k,w,ro)*(a+b)

####### FUNCTIONS TO SOLVE ########

def func(k,om,w,ro):

    return varicose(k,w,ro)-om

def fsin(k,om,w,ri):

    return sinuous(k,w,ri)-om

def fslosh(k,om,om0):

    return sloshing(k,om0)-om

####### SCALAR SOLUTIONS #######

def var_sc(om,w,ro):
    

    s = fsolve(func,10,args=(om,w,ro))
        
    return s

def sin_sc(om,w,ri):

    s = fsolve(fsin,10,args=(om,w,ri))
        
    return s

def slo_sc(om,om0):

    s = fsolve(fslosh,10,args=(om,om0))
        
    return s

####### RESONANCE CONDITION VAR+SIN  #######


def fres_var(om,f0,w,ri,ro,om0,k0):

    return var_sc(om,w,ro)-sin_sc(2*np.pi*f0-om,w,ri)-k0

def fres_sin(om,f0,w,ri,ro,om0,k0):

    return var_sc(2*np.pi*f0-om,w,ro)-sin_sc(om,w,ri)-k0


####### RESONANCE CONDITION VAR+VAR #######

def fres_var_vv(om,f0,w,ri,ro,om0,k0):

    return var_sc(om,w,ro)-var_sc(2*np.pi*f0-om,w,ro)-k0

def fres_sin_vv(om,f0,w,ri,ro,om0,k0):
    
    return var_sc(2*np.pi*f0-om,w,ro)-var_sc(om,w,ro)-k0


####### SOLVING RESONANCE CONDITION #######

#f = np.linspace(5.45,7.2,N)

f = np.linspace(6.8,15,N)

om_r_var = np.zeros(N)
om_r_sin = np.zeros(N)

om_r_var_vv = np.zeros(N)
om_r_sin_vv = np.zeros(N)

k_r_v1 = np.zeros(N)
k_r_v2 = np.zeros(N)
k_f = np.zeros(N)

gdif_1 = np.zeros(N)
gdif_2 = np.zeros(N)

for i in range(N):

    k0_slo = slo_sc(2*np.pi*f[i],om0)
    k0_var = var_sc(2*np.pi*f[i],w,ro)
    k_f[i] = k0_slo
    #s_v = fsolve(fres_var,2.5*np.pi*f[i],args=(f[i],w,ri,ro,om0,k0_slo))

    s_s_vv = fsolve(fres_var_vv,2.5*np.pi*f[i],args=(f[i],w,ri,ro,om0,k0_slo))

    #s_s = fsolve(fres_sin,2.5*np.pi*f[i],args=(f[i],w,ri,ro,om0,k0_slo))

    #om_r_var[i] = 2*np.pi*f[i]-s_s
    
    #om_r_sin[i] = s_s

    om_r_var_vv[i] = 2*np.pi*f[i]-s_s_vv
    
    om_r_sin_vv[i] = s_s_vv

    k_r_v1[i] = var_sc(s_s_vv,w,ro)
    k_r_v2[i] = var_sc(2*np.pi*f[i]-s_s_vv,w,ro)

    gdif_1[i] = abs(group_slo(k0_slo,om0)-group_var(k_r_v1[i],w,ro))
    gdif_2[i] = abs(group_slo(k0_slo,om0)-group_var(k_r_v2[i],w,ro))

####### PLOTTING #######


fig, ax = plt.subplots(figsize=[12,9])

# plt.plot(f,om_r_var/2/np.pi,lw=2,color="r")
# plt.plot(f,om_r_sin/2/np.pi,lw=2,color="b")

plt.plot(f,om_r_var_vv/2/np.pi,lw=2,color=colors[0],ls="--")
plt.plot(f,om_r_sin_vv/2/np.pi,lw=2,color=colors[2],ls="--")

plt.plot(f,(om_r_sin_vv+om_r_var_vv)/2/np.pi,lw=2,color=colors[-2])


ax.set_xlabel('$f_0$ [Hz]',fontsize=35)
ax.set_ylabel('$f_{1,2}$ [Hz]',fontsize=35)
ax.minorticks_on()
ax.tick_params(which='minor', direction='in')
ax.tick_params(labelsize=25,colors=dark_gray,direction="in")


# plt.figure()
# plt.plot(f,-om_r_var_vv/2/np.pi+om_r_sin/2/np.pi,lw=2,color="b",ls="--")
# plt.plot(f,om_r_sin_vv/2/np.pi-om_r_var/2/np.pi,lw=2,color="r",ls="--")


#####

# fig, ax = plt.subplots(figsize=[12,9])

# plt.plot(f,gdif_1)
# plt.plot(f,gdif_2)

# plt.plot(f,np.sqrt((k_r_v1-k_f)**2+(om_r_var_vv-2*np.pi*f)**2),"r")
# plt.plot(f,k_r_v2,"b")


plt.show()
