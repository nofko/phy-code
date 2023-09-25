import matplotlib.pyplot as plt
import numpy as np
import trackpy as tp
import pims
import sys

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ==========================================================================##
## NANTES CAMPAIGN                                                           ##
##                                                                           ##
##                                BALL TRACKING                              ##
##                                                                           ##
## ==========================================================================##


############################          SETUP        ############################

@pims.pipeline
def gray(image):
    return image[:, :, 0]  # Take just the green channel


## SMALL BALL SETUP ##

diam   = 9
thd    = 30
mmass  = 1500
mxsize = 2.2


## LARGE BALL SETUP ##

diam   = 11
thd    = 30
mmass  = 3500
mxsize = 15


fps  = 1/60
pxm  = 8e6/359

frames = gray(pims.open('/data/nantes/DCIM/run4_resized/*.JPG'))


arg = sys.argv[1]


############################         LOCATE       ############################


if arg == "s":

    n = 500
    
    f = tp.locate(frames[n], diam,
                  threshold = thd,
                  minmass = mmass,
                  maxsize = mxsize,
                  invert = False,
                  #processes='auto'
                  )



    tp.annotate(f, frames[n]);


    fig, ax = plt.subplots()
    ax.hist(f['size'], bins=20)
    plt.title("Size")

    fig, ax = plt.subplots()
    ax.hist(f['mass'], bins=20)
    plt.title("Mass")

    plt.show()



############################         BATCH        ############################


if arg == "b":


    f = tp.batch(frames, diam,
                  threshold = thd,
                  minmass = mmass,
                  maxsize = mxsize,
                  invert = False,
                  processes='auto'
                  )

    t = tp.link(f, 5, memory=3)


    t = tp.filter_stubs(t, 5)

    
    d = tp.compute_drift(t)


    t_nd = tp.subtract_drift(t.copy(), d)
    

    # im = tp.imsd(t_nd, pxm, fps)

    em_nd = tp.emsd(t_nd, 8e6/359, 1/60,detail=True)
    em = tp.emsd(t, 8e6/359, 1/60,detail=True)
    
    fig, ax = plt.subplots()
    ax.plot(em.index, em.msd, 'o')
    ax.plot(em_nd.index, em_nd.msd, 'o')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
           xlabel='lag time $t$')


    #pw = tp.utils.fit_powerlaw(em.msd)

    em.to_csv('run4_big.csv', sep=',', header=True)
    em_nd.to_csv('run4_big_nd.csv', sep=',', header=True)

    t_nd.to_csv('run4_big_track.csv', sep=',', header=True)

    #print(em.index)
    plt.show()
    #print(pw)

