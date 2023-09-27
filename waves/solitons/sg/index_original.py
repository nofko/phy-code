import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz

x0 = -0.1
y0 = 1./np.sqrt(3.)

def funcFF(x,y):
	return ( ((x-x0) +1j* (y-y0)) )**2

def indice_contour(FF, XX, YY, ix1, ix2, iy1, iy2):
	""" 1/(2*j*np.pi) * Integrale de contour de F'/F
	FF est le resultat de F appliquee au meshgrid contenu dans XX et YY
	Le contour est un rectangle sur des points de discretisation
	d'indices ix1, ix2, iy1, iy2
	On calcule F' par differences finies centree
	sens de parcours :
	1: (ix1, iy1) --> (ix2, iy1)
	2: (ix2, iy1) --> (ix2, iy2)
	3: (ix2, iy2) --> (ix1, iy2)
	4: (ix1, iy2) --> (ix1, iy1)
	FF[i,j] = F(xi, yj)
	"""
	# integration horizontale
	dx = XX[0,1]-XX[0,0]
	DFF = (np.roll(FF, 1, axis=0)-np.roll(FF, -1, axis=0))/2/dx/FF
	VV = DFF[ix1 : ix2+1, iy1]
	I1 = +dx*(sum(VV)-(VV[0]+VV[-1])/2.) 
	VV = DFF[ix1 : ix2+1, iy2]		 
	I3 = -dx*(sum(VV)-(VV[0]+VV[-1])/2.) 
	
	# integration verticale
	dy = YY[1,0]-YY[0,0]
	DFF = (np.roll(FF, 1, axis=1)-np.roll(FF, -1, axis=1))/2/dy/FF
	VV = DFF[ix2, iy1 : iy2+1]
	I2 = +dy*(sum(VV)-(VV[0]+VV[-1])/2.) 
	VV = DFF[ix1, iy1 : iy2+1]		 
	I4 = -dy*(sum(VV)-(VV[0]+VV[-1])/2.) 
	#
	integrale = I1 + I2 + I3 + I4
	resultat = 1./(2.*1j*np.pi) * integrale
	return resultat


# domaine total
Xmax = 4
Xmin = -Xmax
Ymax = Xmax
Ymin = Xmin

# discretisations spatiales
Ndisc = 400   # en x
Mdisc = Ndisc # en y

Xdisc = np.linspace(Xmin, Xmax, Ndisc+1)
Ydisc = np.linspace(Ymin, Ymax, Mdisc+1)

XX, YY = np.meshgrid(Xdisc, Ydisc)

# valeurs complexes du champ et ses parties reelles et imaginaires 
MM = funcFF(XX, YY)
AA = np.real(MM)
BB = np.imag(MM)

X0 = [x0]
Y0 = [y0]


# indices des coins du contour
ix1 = 130
ix2 = 300
iy1 = 120
iy2 = 390

# coordonnes correspondantes des coins du contour
x1 = Xdisc[ix1]
x2 = Xdisc[ix2]
y1 = Ydisc[iy1]
y2 = Ydisc[iy2]

print(x1, x2, y1, y2)

Xcontour = [x1, x2, x2, x1, x1]
Ycontour = [y1, y1, y2, y2, y1]

# Min et Max pour la colorbar
Cmax = 2
Cmin = -Cmax

# ******************************************************
# CALCUL DE L'INDICE
#

zindice = indice_contour(MM, XX, YY, ix1, ix2, iy1, iy2)
print("\nIndice complexe = \n",zindice)
print("\n*** L'indice vaut : \n\n", int(round(np.real(zindice), 0)), "\n")

#
#
# ******************************************************

# -----------------------
# GRAPHES D'ILLUSTRATIONS
# -----------------------

#################################################################

plt.figure()
plt.pcolormesh(XX, YY, AA, cmap='bwr', shading = "auto")
plt.clim(Cmin, Cmax)
plt.colorbar()
plt.plot(X0, Y0, 'or')
plt.plot(Xcontour, Ycontour, '--')
plt.contour(XX, YY, AA)
plt.title("AA")

#################################################################

plt.figure()
plt.pcolormesh(XX, YY, BB, cmap='bwr', shading = "auto")
plt.clim(Cmin, Cmax)
plt.colorbar()
plt.plot(X0, Y0, 'or')
plt.plot(Xcontour, Ycontour, '--')
plt.contour(XX, YY, BB)
plt.title("BB")
plt.show()

#################################################################
