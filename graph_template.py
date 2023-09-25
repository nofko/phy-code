import matplotlib.pyplot as plt
import numpy as np



plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ==========================================================================##
## DATA                                                                      ##
##                                                                           ##
##                              PLOTTING TEMPLATE                            ##
##                                                                           ##
## ==========================================================================##



############################          SETUP        ############################

N = 10

x = np.linspace(1,10,N)
y = 1/x+(np.random.rand(N)-0.5)/10

xl = np.linspace(0.75,10.5,1000)

xi = np.linspace(1,10,N)
yi = np.sqrt(xi+(np.random.rand(N)-0.5)/5)
yerr = np.random.rand(N)/2


############################        PLOTTING       ############################


fig, ax = plt.subplots(figsize=[12,9])        


ax.plot(x,y,"o",color="r",ms=10)

ax.plot(xl,1/xl,color="b",lw=1.5)


ax.set_xlabel(r"$x$ [m]",fontsize=30)
ax.set_ylabel(r"$y$ [m$^{-1}$]",fontsize=30)

ax.tick_params(labelsize=20,direction="in")



############################         INSET         ############################


left, bottom, width, height = [0.53, 0.5, 0.35, 0.35]
ax_i = fig.add_axes([left, bottom, width, height])

ax_i.errorbar(x,yi,yerr=yerr,fmt="o",ms=7,ecolor="gray",color="b")

ax_i.set_ylabel("$B_x(z)$ [T]",fontsize=22)
ax_i.set_xlabel(r"$z$ [cm]",fontsize=22)


ax_i.tick_params(labelsize=18,direction="in")
ax_i.tick_params(which='minor', direction='in',labelsize=18)



plt.savefig("example.pdf",bbox_inches="tight")

#plt.show()
