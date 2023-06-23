import matplotlib.pyplot as plt
import numpy as np
import sys

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)



## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                               BROWNIAN MOTION                             ##
##                                                                           ##
## ============================================================================



N = 5000

M = 50

x = np.zeros(N)

y = np.zeros(N)


def trajectory():

    # for i in range(N-1):

    #     th = np.random.randn(1)*2*np.pi
    #     dx = np.cos(th)
    #     dy = np.sin(th)

    #     x[i+1] = x[i]+dx
    #     y[i+1] = y[i]+dy

    x = np.cumsum(np.random.randn(N))
    y = np.cumsum(np.random.randn(N))
        
    return x,y

def compute_MSD(path):
   totalsize=len(path)
   msd=[]
   for i in range(totalsize-1):
       j=i+1
       msd.append(np.sum((path[0:-j]-path[j::])**2)/float(totalsize-j))

   msd=np.array(msd)
   
   return msd


def max_rolling1(a, window,axis =1):
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        rolling = np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
        return np.max(rolling,axis=axis)

def compute_MME(path):

    totalsize = len(path)
    mme = []

    r = np.sqrt((path[:,0]-path[0,0])**2+(path[:,1]-path[0,1])**2)

    for i in range(1660):
        
        j = i+1
        
        rr = r[:-j]

        mi = np.sum(max_rolling1(rr[:-j],j)**2)
        
        mme.append(mi/float(totalsize-j))

    mme = np.array(mme)
   
    return mme


def pdf(X,nbin):

    maximum = np.max((X-np.mean(X))/np.std(X))
    minimum = np.min((X-np.mean(X))/np.std(X))

    
    bins = np.linspace(minimum,maximum,nbin)
    
    hist,b = np.histogram((X-np.mean(X))/np.std(X),bins)

    hist = hist / ( np.sum(hist)*np.mean(np.diff(b)))

    return hist,b

MS = np.zeros(N-1)
ME = np.zeros(1660)

for j in range(M):

    x,y = trajectory()

    ms = compute_MSD(np.array(list(zip(x,y))))

    me = compute_MME(np.array(list(zip(x,y))))

    MS += ms
    ME += me


plt.plot(MS/M)

plt.figure()

plt.plot(ME/M)

#plt.plot(np.sqrt(np.cumsum(x**2+y**2))/N)

plt.show()
