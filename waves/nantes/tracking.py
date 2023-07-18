import matplotlib.pyplot as plt
import numpy as np
import trackpy as tp
import pims

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ==========================================================================##
## NANTES CAMPAIGN                                                           ##
##                                                                           ##
##                                BALL TRACKING                              ##
##                                                                           ##
## ==========================================================================##




@pims.pipeline
def gray(image):
    return image[:, :, 0]  # Take just the green channel

frames = gray(pims.open('/data/nantes/DCIM/test/*.JPG'))

f = tp.batch(frames, 11,
              threshold = 30,
              minmass = 100,
              maxsize = 50,
              invert=False,
              percentile=70,
              topn = 200)



t = tp.link(f, 100, memory=2)


t1 = tp.filter_stubs(t, 3)

1
plt.figure()
tp.plot_traj(t);



d = tp.compute_drift(t)


tm = tp.subtract_drift(t.copy(), d)
im = tp.imsd(tm, 1, 10)


fig, ax = plt.subplots()
ax.plot(im.index, im, 'k-', alpha=0.1)  # black lines, semitransparent
ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
       xlabel='lag time $t$')
ax.set_xscale('log')
ax.set_yscale('log')



# em = tp.emsd(tm, 1, 10) 


#tp.utils.fit_powerlaw(em)

plt.show()
