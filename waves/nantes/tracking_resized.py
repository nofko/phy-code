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

# diam   = 11
# thd    = 30
# mmass  = 3500
# mxsize = 15


fps  = 1/60
pxm  = 8e6/359

frames = gray(pims.open('/data/nantes/DCIM/run2_resized/*.JPG'))


arg = sys.argv[1]


############################         LOCATE       ############################


if arg == "s":

    
    f = tp.locate(frames[0], diam,
                  threshold = thd,
                  minmass = mmass,
                  maxsize = mxsize,
                  invert = False,
                  #processes='auto'
                  )



    tp.annotate(f, frames[0]);


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


    t1 = tp.filter_stubs(t, 5)



    d = tp.compute_drift(t)


    tm = tp.subtract_drift(t.copy(), d)
    im = tp.imsd(tm, pxm, fps)

    em = tp.emsd(tm, 8e6/359, 1/60)



    fig, ax = plt.subplots()
    ax.plot(em.index, em, 'o')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
           xlabel='lag time $t$')


    pw = tp.utils.fit_powerlaw(em)

    em.to_csv('run2_small.csv', sep=',', header=True)

    print(em.index)
    
    print(pw)

