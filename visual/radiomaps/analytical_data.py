from math import *
import numpy as np
import settings as stg

def data_raw(nx,ny,nz):
# The function to create a simple model of an axisymmetric, magnetized disk for testing of the synchrotron mapping package
#nx, ny, nz - number of grid zones for x, y and z directions of the rectangular volume to be visualized
#x_ion - ionization degree of thermal gas needed for computing Faraday rotation  (stg)

   # inicjalizacja tablic
   shapeB =          (nx,ny,nz)
   shapeE = (stg.ncre,nx,ny,nz)
   Bx = np.zeros(shapeB)
   By = np.zeros(shapeB)
   Bz = np.zeros(shapeB)
   Ecrp = np.zeros(shapeB)
   Ecre = np.zeros(shapeE)
   rhoi = np.zeros(shapeB)
   i0 = 0.5*(nx-1)
   j0 = 0.5*(ny-1)
   k0 = 0.5*(nz-1)

   r0 = nx/4.0
   #srodek dysku
   for ii in range (nx):
      for jj in range (ny):
         r = sqrt((ii-i0)**2 + (jj-j0)**2)
         rfact = exp(-r/r0)
         for kk in range (nz):
            zfact = exp(-0.05*(kk-k0)**2)
#
            rhoi[ ii,jj,kk] = 1.0*zfact*stg.x_ion*rfact
            Ecrp[ii,jj,kk] = 1.0*zfact*rfact
            Ecre[:,ii,jj,kk] = 1.0*zfact*rfact
#
            Bx[ii,jj,kk] = -(jj-j0)/r*zfact+1e-5*rfact
            By[ii,jj,kk] =  (ii-i0)/r*zfact+1e-5*rfact
            Bz[ii,jj,kk] = 0.0
   return Bx, By, Bz, rhoi, Ecrp, Ecre
