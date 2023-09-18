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

frames = gray(pims.open('/data/nantes/DCIM/run2_300/*.JPG'))


f = tp.batch(frames, 105,
              threshold = 30,
              minmass = 50,
              maxsize = 140,
              invert=False,
              #percentile=70,
              #topn = 30
              )

tp.annotate(f, frames[0]);



fig, ax = plt.subplots()
ax.hist(f['size'], bins=20)

# Optionally, label the axes.
ax.set(xlabel='mass', ylabel='count');



t = tp.link(f, 50, memory=3)


t1 = tp.filter_stubs(t, 3)

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


em = tp.emsd(tm, 1, 24)


plt.figure()
plt.ylabel(r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]')
plt.xlabel('lag time $t$');
#tp.utils.fit_powerlaw(em) 


# em = tp.emsd(tm, 1, 10) 


#tp.utils.fit_powerlaw(em)

plt.show()
