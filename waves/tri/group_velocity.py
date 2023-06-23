import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
from scipy.optimize import fsolve

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)



## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                          GROUP VELOCITY MATCHING                          ##
##                                                                           ##
## ============================================================================




####### PHYSICS #######

N =1000

rc = 0.07

ro = 0.0785 #np.linspace(0.077,0.089,N)
w = 2*(ro-rc)
ri =ro-w

sigma = 0.055
rho = 1000

om0 = 44.5

g = np.sin(4.5*np.pi/180)*9.81

k = np.linspace(0.1,50,N)

####### DISPERSION RELATIONS #######

def varicose(k,w,ro):

    return np.sqrt((g*k/ro+sigma*k**3/rho/ro**3)*np.tanh(k*w*ro/2/rc**2))

# def sloshing(k,om0):
    
#     return np.sqrt(om0**2+g*k**2/0.076)


def sloshing(k,om0):
    
    return np.sqrt((om0**2+(g*k**2/ro+0.5*sigma*k**3/rho/ro**3)))


def group_slo(k,om0):

    return 0.5/sloshing(k,om0)*(2*g*k/ro+1.5*sigma*k**2/rho/ro**3)

def group_var(k,w,ro):

    a = (g/ro+3*sigma*k**2/rho/ro**3)*np.tanh(k*w*ro/2/rc**2)

    b = (g*k/ro+sigma*k**3/rho/ro**3)/np.cosh(k*w*ro/2/rc**2)**2*w/2*ro/rc**2
    
    return 0.5/varicose(k,w,ro)*(a+b)

def difference(k,om0,w,ro):
    
    return group_slo(k,om0)-group_var(k,w,ro)

def difference_phase(k,om0,w,ro):
    
    return group_slo(k,om0)-varicose(k,w,ro)/k

def difference_branch(k,om0,w,ro):
    
    return sloshing(k,om0)-varicose(k,w,ro)

####### PLOTTING #######

res = fsolve(difference,6,args=(om0,w,ro))

res_phase = fsolve(difference_phase,6,args=(om0,w,ro))

res_branch = fsolve(difference_branch,100,args=(om0,w,ro))

f_sigma = sloshing(res_phase[0],om0)/2/np.pi

print("k of group velocity match: ", res[0])
print("f of group velocity match (sloshing): ", sloshing(res[0],om0)/2/np.pi)
print()
print("k of group-phase match: ", res_phase[0])
print("f of group-phase velocity match: ", f_sigma)
print("fgc of group-phase velocity match: ", varicose(3.28,w,ro)/2/np.pi)
print()
print("k of branch match: ", res_branch[0])
print("f of branch match: ", varicose(105,w,ro)/2/np.pi)

fig, ax = plt.subplots(figsize=[12,8])

ax.set_xlabel('$k$',fontsize=35)
ax.set_ylabel(r'$c_g$',fontsize=35)
ax.tick_params(labelsize=25)

plt.plot(k,group_slo(k,om0),label="$c_g^{\Sigma}$",color="r")
plt.plot(k,group_var(k,w,ro),label="$c_g^{V}$",color="g")
plt.plot(k,varicose(k,w,ro)/k,label="$c_p^{V}$",color="b")

plt.axhline(0,linestyle="-.",color="gray")
plt.axvline(res_phase[0],linestyle="--",color="k")
plt.text(6.5,0.5,"$f_\Sigma=$"+str(round(f_sigma,3)))

leg = ax.legend(loc="upper left",fontsize=20,frameon=False,handletextpad=0.2,labelspacing=0.5)
leg.get_frame().set_edgecolor('k')

#plt.savefig("images/velocity_matching.pdf",bbox_inches="tight")

plt.figure()

k = np.linspace(0.1,150,N)

#plt.plot(k,sloshing(k,om0))
#plt.plot(k,varicose(k,w,ro))
plt.plot(k[:-1],np.diff(group_var(k,w,ro)),label="$c_g^{\Sigma}$",color="r")

print(varicose(178.7,w,ro)/2/np.pi)

plt.show()
