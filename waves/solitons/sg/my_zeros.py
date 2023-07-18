import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import bisect
from scipy.signal import argrelmin,argrelmax

import real_trace as sg


## ============================================================================
## SINE-GORDON DIRECT SCATTERING                                             ##
##                                                                           ##
##                        SINE-GORDON REAL LINE ZEROS                        ##
##                                                                           ##
## ============================================================================




def trace_finder_plus(E):
    
    Mt = sg.B(u,ux,ut,E+0j,dx,dt)

    tr = np.trace(Mt.real)/2

    if np.abs(tr)<=1:
        return tr-1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))-1
    
def trace_finder_minus(E):
    
    Mt = sg.B(u,ux,ut,E+0j,dx,dt)

    tr = np.trace(Mt.real)/2
    
    if np.abs(tr)<=1:
        return tr+1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))+1


def find_zeros(E,trace):

    Eigs_plus = []
    Eigs_minus = []

    mins = argrelmin(trace)[0]
    maxs = np.concatenate(([int(mins[0]/2)],argrelmax(trace)[0]))


    Emaxs = E[maxs]
    Emins = E[mins]

    Ebr = np.zeros(int(len(Emaxs)+len(Emins)))

    for i in range(len(Emaxs)):

        Ebr[2*i] = Emaxs[i].real

    for i in range(len(Emins)):
    
        Ebr[2*i+1] = Emins[i].real


    for i in range(len(Ebr)-1):

        if (np.sign(trace_finder_plus(Ebr[i]))!=np.sign(trace_finder_plus(Ebr[i+1]))):
    
            Eigs_plus.append(bisect(trace_finder_plus,Ebr[i],Ebr[i+1],xtol=1e-16))

        if (np.sign(trace_finder_minus(Ebr[i]))!=np.sign(trace_finder_minus(Ebr[i+1]))):
        
            Eigs_minus.append(bisect(trace_finder_minus,Ebr[i],Ebr[i+1],xtol=1e-16))


    return Eigs_plus,Eigs_minus
