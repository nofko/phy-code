import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
from scipy.optimize import fsolve, bisect
from scipy.signal import argrelmin,argrelmax
import mat73
import time

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)



## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                         PIST OF SPACE-TIME SIGNAL                         ##
##                                                                           ##
## ============================================================================


def linlog(y):

    N = len(y)
    out = np.zeros(N)

    for i in range(N):
        if np.abs(y[i])<=1:
            out[i] = y[i]

        else:
            out[i] = np.sign(y[i])*(1+np.log(np.abs(y[i])))

    return out



def trace_finder_plus(E,u,dx):
    Mt = B(lam*u+E,dx)

    tr = np.trace(Mt)/2
    
    if np.abs(tr)<=1:
        return tr-1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))-1
    
def trace_finder_minus(E,u,dx):
    Mt = B(lam*u+E,dx)

    tr = np.trace(Mt)/2
    
    if np.abs(tr)<=1:
        return tr+1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))+1


def B(q,dx):

    Mi = np.identity(2)    
    
    for i in range(len(q)):

        sq =np.sqrt(abs(q[i]))
        
        if q[i]<0:
            
            cos = np.cosh(dx*sq)
            sh = np.sinh(dx*sq)
            
            exp12 = sh/sq
            exp21 = sh*sq
            
        else:

            cos = np.cos(dx*sq)
            sin = np.sin(dx*sq)

            exp12 = sin/sq
            exp21 = -sin*sq
            
        exp = np.array([[cos,exp12],[exp21,cos]])

        #Mi1 = Mi
        Mi = np.dot(exp,Mi)
        
        #s += np.abs(np.sign(Mi[0,0])-np.sign(Mi1[0,0]))/2

    return Mi


def frame(u,E,dx,ns):

    N = len(u)

    M = np.zeros((N,2,2))
    
    for j in range(N):
        q = lam*u+E[j]
    
        M[j] = B(q,dx)


    trM = np.trace(M,axis1=1,axis2=2)

    trL = linlog(trM/2)


    Eigs_plus = []
    Eigs_minus = []

    mins = argrelmin(trL)[0]
    maxs = np.concatenate(([int(mins[0]/2)],argrelmax(trL)[0]))


    Emaxs = E[maxs]
    Emins = E[mins]

    
    Ebr = np.zeros(int(len(Emaxs)+len(Emins)))
    
    for i in range(len(Emaxs)):

        Ebr[2*i] = Emaxs[i]

    for i in range(len(Emins)):
    
        Ebr[2*i+1] = Emins[i]

    for i in range(len(Ebr)-1):

        if (np.sign(trace_finder_plus(Ebr[i],u,dx))!=np.sign(trace_finder_plus(Ebr[i+1],u,dx))):

            Eigs_plus.append(bisect(trace_finder_plus,Ebr[i],Ebr[i+1],xtol=1e-13,args=(u,dx)))

        if (np.sign(trace_finder_minus(Ebr[i],u,dx))!=np.sign(trace_finder_minus(Ebr[i+1],u,dx))):
            
            Eigs_minus.append(bisect(trace_finder_minus,Ebr[i],Ebr[i+1],xtol=1e-13,args=(u,dx)))

    Eigs = np.concatenate((Eigs_plus,Eigs_minus))
    Eigs = np.sort(Eigs)

    sindex = []

    for i in range(1,int(len(Eigs)/2)):
        sindex.append((Eigs[2*i]-Eigs[2*i-1])/(Eigs[2*i]-Eigs[2*i-2]))

    return trL,Eigs[:ns]



def movie(name):

    tst = mat73.loadmat(folder_data+name)
    signal = tst["SignalT"]
    signal = signal["Tot"]    
    signal = np.array(signal)
    signal = signal[:,::3]

    x = tst["XX"]
    N = len(signal[0])
    Nf = len(signal)
    Ns = 12
    L = x[-1]

    

    dx = x[1]-x[0]
    
    
    E_min = -lam*(np.max(signal)-np.mean(signal[0]))
    E_max =  E_min + 6*np.abs(E_min)
    E = np.linspace(E_min,E_max,N)

    trace = np.zeros((Nf,N)) 
    ei = np.zeros((Nf,Ns))
    
    for i in range(Nf):
        print(i/Nf)
        start = time.time()
        
        t,s = frame(signal[i],E,dx,Ns)
        
        trace[i] = t
        ei[i] = s
        
        end = time.time()
        print("Time:", end-start)

    n2 = name.split(".")
    np.save(folder_res+n2[0]+"_tr.npy",trace)
    np.save(folder_res+n2[0]+"_ei.npy",ei)
    np.save(folder_res+n2[0]+"_e.npy",E)

    return

h = 3 #water depth
lam = 3/2/h**3

folder_data = "/data/canal/Signaux_Soliton_Sin_Noise/"
folder_res = "/data/canal/Signaux_Soliton_Sin_Noise/results/"

name = "CumSum_one_soliton_A10mm_h3cm_L4m_30s_Signal_Tot.mat"

movie(name)
