from numpy import log, log10, zeros, asfarray, sqrt
from settings import MHz, p_min_fix, p_max_fix, const_synch, maxcren

def prepare_coeff(ncre):
# var_names now should be initialized with values from problem.par@hdf file
# Reconstruct momenta
   edges = []
   p_fix = []
   edges[0:ncre+1] = range(0,ncre+1, 1)
   p_fix[0:ncre+1] = zeros(ncre+1)

# WARNING !!! cutoff momenta are not precisely described here !!! WARNING
   log_width   = (log10(p_max_fix/p_min_fix))/(ncre-2.0)
   p_fix_ratio = 10.0 ** log_width

   for i in range(0,ncre-1):
      p_fix[i+1] = p_min_fix * p_fix_ratio**float(i)
   p_fix[0]    = ( sqrt(p_fix[1] * p_fix[2]) ) / p_fix_ratio               #let us rid of zero values
   p_fix[ncre] = ( sqrt(p_fix[ncre-2]* p_fix[ncre-1]) ) * p_fix_ratio
   p_fix = asfarray(p_fix)
   return p_fix

def nu_to_ind(nu, B_perp, ncre, p_fix):
   #print nu, MHz, B_perp
   p_nu = 1956.0 * sqrt(nu / (16.1*MHz)) / sqrt( B_perp +1.e-24 ) # [B_perp] = mGs - assumed at input, [nu] = 16.1 MHz, but already present
   #print p_nu, p_min_fix, p_max_fix, ncre
   nu_to_ind = int((log10(p_nu/p_min_fix)/log10(p_max_fix/p_min_fix)) * float(ncre - 2.0 )+1.0)

   if nu_to_ind > ncre: nu_to_ind = ncre
   if nu_to_ind < 0: nu_to_ind = 0
   if ( not p_nu < p_fix[nu_to_ind] and p_nu > p_fix[nu_to_ind-1]): # WARNING temporary fix
      nu_to_ind = min(nu_to_ind +1,ncre)
   return nu_to_ind, p_nu

def crenpp(nu_s, ncre, bperp, ecr, p_fix):
   p_ind, p_nu = nu_to_ind(nu_s, bperp, ncre, p_fix)
   n_ind = min(p_ind - 1, maxcren-1)
   #cren_i = 10**float(n_ind-11) #Ecr[n_ind,i3]
   cren_i = ecr[n_ind]
   #print 'numer binu: ', n_ind, i3
   p_i    = p_fix[p_ind]
   if (p_nu < p_fix[p_ind] and p_nu > p_fix[p_ind-1]):
      delta_p = p_i - p_fix[p_ind-1]
   else:
      delta_p = p_max_fix - p_i  # WARNING temporary fix
   #elfq = const_synch * cren_i / p_i # powinno byc podzielone przez delta_p = (p_i+1/2 - p_i-1/2)
   elfq = const_synch * cren_i / delta_p
   return elfq
   