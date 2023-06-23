import numpy as np
import matplotlib.pyplot as plt

import matplotlib.animation as animation


plt.rcParams['savefig.bbox'] = 'tight' 

x = np.linspace(-10,10,1000)

t= np.linspace(-3.5,3.5,150)

xx, tt = np.meshgrid(x,t)

u = 12*(3+4*np.cosh(2*xx-8*tt)+np.cosh(4*xx-64*tt))/(3*np.cosh(xx-28*tt)+np.cosh(3*xx-36*tt))**2

k1 = 1.7
k2 = 1

th1 = 2*k1*xx-8*k1**3*tt
th2 = 2*k2*xx-8*k2**3*tt

et1 = np.exp(th1)
et2 = np.exp(th2)

u = 8*(k1**2*et1+k2**2*et2+2*(k1-k2)**2*et1*et2+(k1-k2)**2/(k1+k2)**2*(k1**2*et1*et2**2+k2**2*et2*et1**2))/(1+et1+et2+(k1-k2)**2/(k1+k2)**2*et1*et2)**2




fig, ax = plt.subplots(figsize=[5.3,4])
line, = ax.plot([], [], lw=2)

fig.patch.set_facecolor('white')

ax.grid(color="silver")

ax.set_xlim(-10, 10)
ax.set_ylim(0, 6)

ax.set_xlabel('$x$',fontsize=30)
ax.set_ylabel("$\phi$",fontsize=30)

plt.tight_layout()
def update_line(i):
    line.set_data(x,u[i])
    fig.patch.set_facecolor('white')

    return line,


line_ani = animation.FuncAnimation(fig, update_line, frames=range(len(u)),
                                   interval=1, blit=True)

#mywriter = animation.FFMpegWriter(fps=60)
#line_ani.save('mymovie.mpeg',writer=mywriter)

line_ani.save('myAnimation.gif', writer='imagemagick', fps=15,)


fig2, ax2 = plt.subplots(figsize=[10,8])

four = ax2.pcolormesh(xx,tt,u,shading="auto")

ax2.set_xlabel('$x$',fontsize=30)
ax2.set_ylabel("$t$",fontsize=30)

ax2.set_ylim(-2.8,2.8)

cbar = plt.colorbar(four,ax=ax2)
cbar.set_label(r"$\phi$",fontsize=30)
#cbar.ax.tick_params(labelsize=30)

plt.savefig("2_soliton_color.png",bbox_inches="tight",dpi=500)
