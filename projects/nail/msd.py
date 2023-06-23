import matplotlib.pyplot as plt
import numpy as np
#import Ccolors

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ==========================================================================##
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                            MEAN SQUARE DISTANCE                           ##
##                                                                           ##
## ==========================================================================##




############################        LOAD FILE      ############################


path = "/data/nail/05_06/2_Hz_A9Vpp.npy"


coord = np.load(path,allow_pickle=True)


############################        FUNCTIONS      ############################


def compute_MSD(path):

    totalsize = len(path)
    msd = []

    for i in range(totalsize-1):
        
        j = i+1
        msd.append(np.sum((path[0:-j]-path[j::])**2)/float(totalsize-j))

    msd = np.array(msd)
   
    return msd


def compute_MME(path):

    totalsize = len(path)
    msd = []

    for i in range(totalsize-1):
        
        j = i+1
        msd.append(np.sum((path[0:-j]-path[j::])**2)/float(totalsize-j))

    msd = np.array(msd)
   
    return msd

############################        PLOTTING      ############################

pxm = 0.15/1169

fps = 100

dt = 1/fps

coord = coord*pxm*1000

ms = compute_MSD(coord)*dt


t = np.arange(len(coord))/fps

ts = np.arange(len(coord)-1)/fps
 

plt.loglog(ts,ms)



#plt.plot(ts,ms)
plt.axvline(1/2)

r = 40
mt = r**2*np.pi/3
 
#plt.axhline(mt)
plt.loglog(t,t**0.75*0.2)


x = coord[:,0]
y = coord[:,1]

w = 10
xm = np.convolve(x, np.ones(w), 'valid') / w
ym = np.convolve(y, np.ones(w), 'valid') / w

cm = np.array(list(zip(xm,ym)))
mm = compute_MSD(cm)

plt.figure()

distance = np.sqrt(x**2+y**2)

mme = np.maximum.accumulate(distance)

# n = 2000
plt.plot(x)
plt.plot(y)

# plt.plot(mm)


plt.show()
