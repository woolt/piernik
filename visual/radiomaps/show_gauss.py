from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from gauss_beam import gauss_beam

nbeam=30
sigma=2

def show_gauss(nbeam,sigma):
# function drawing the 2D Gauss function representing the profile of radiotelescope beam

	Z=gauss_beam(nbeam,sigma)
	g_shape=np.shape(Z)[0]
	fig = plt.figure()
	plot_style = fig.gca(projection='3d')
	X = np.arange(0, g_shape)
	Y = np.arange(0, g_shape)
	X, Y = np.meshgrid(X, Y)
	surf = plot_style.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)

show_gauss(nbeam,sigma)
plt.show()
