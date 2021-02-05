from math import *
import numpy as np

def gauss_beam(nbeam,sigma):
# The function retorns a 2D Gauss function representing the angular sensitivity of the telescope beam
# nbeam - beam resolution
# sigma - (angular) width of the beam 
	g2d=np.zeros((nbeam+1,nbeam+1))
	x=np.zeros((nbeam+1,nbeam+1))
	y=np.zeros((nbeam+1,nbeam+1))

	l=np.linspace(0,nbeam,num=nbeam+1)

	for ii in range(nbeam+1):
		x[ii,:]=l-0.5*(nbeam)
		y[:,ii]=l-0.5*(nbeam)

	g2d[:,:]=(1/(sigma*sqrt(2*pi)))*np.exp(-((x/sigma)**2+(y/sigma)**2)/(2*sigma**2))
	maks=np.max(g2d)
	g2d=g2d/maks
	return g2d
