!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"

!>
!! \brief Timestep particle module
!!
!! This module contains all particle related routines to determine dt limitations.
!<

module particle_timestep
! pulled by NBODY

   use types, only: value

   implicit none

   private
   public :: timestep_nbody, dt_nbody, pacc_max

   real        :: dt_nbody           !< timestep depends on particles
   type(value) :: pacc_max

contains

   subroutine timestep_nbody(dt)

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: big, one, pMIN, two, zero
      use dataio_pub,     only: msg, printinfo
      use func,           only: operator(.notequals.)
      use global,         only: dt_max
      use grid_cont,      only: grid_container
      use mpisetup,       only: piernik_MPI_Allreduce, master
      use particle_utils, only: max_pacc_3d
      use particle_pub,   only: lf_c, ignore_dt_fluid
#ifdef DUST_PARTICLES
      use constants,      only: ndims, xdim, zdim
      use particle_utils, only: max_pvel_1d
#endif /* DUST_PARTICLES */

      implicit none

      real,            intent(inout) :: dt
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      real                           :: eta, eps, factor, dt_hydro
#ifdef DUST_PARTICLES
      integer(kind=4)                :: cdim
      real, dimension(ndims)         :: max_v
#endif /* DUST_PARTICLES */

#ifdef VERBOSE
      call printinfo('[particle_timestep:timestep_nbody] Commencing timestep_nbody')
#endif /* VERBOSE */

      eps      = 1.0e-1
      factor   = one
      dt_nbody = big

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         eta = minval(cg%dl)   !scale timestep with cell size

         call max_pacc_3d(cg, pacc_max)

         if (pacc_max%val .notequals. zero) then
            dt_nbody = sqrt(two*eta*eps/pacc_max%val)

#ifdef DUST_PARTICLES
            call max_pvel_1d(cg, max_v)
            if (any(max_v*dt_nbody > cg%dl)) then
               factor = big
               do cdim = xdim, zdim
                  if ((max_v(cdim) .notequals. zero)) then
                     factor = min(cg%dl(cdim)/max_v(cdim), factor)
                  endif
               enddo
            endif
#endif /* DUST_PARTICLES */

            dt_nbody  = lf_c * factor * dt_nbody
         endif

         cgl => cgl%nxt
      enddo

      dt_nbody = min(dt_nbody, dt_max)
      call piernik_MPI_Allreduce(dt_nbody, pMIN)
      pacc_max%assoc = dt_nbody
      dt_hydro = dt

      if (ignore_dt_fluid) then
         dt = dt_nbody      !IGNORE HYDRO TIMESTEP
      else
         !> \todo verify this condition
         if (dt_nbody .notequals. 0.0) dt = min(dt, dt_nbody)
      endif

      write(msg,'(a,3g12.5)') '[particle_timestep:timestep_nbody] dt for hydro, nbody and both: ', dt_hydro, dt_nbody, dt
      if (master) call printinfo(msg)

#ifdef VERBOSE
      call printinfo('[particle_timestep:timestep_nbody] Finish timestep_nbody')
#endif /* VERBOSE */

   end subroutine timestep_nbody

end module particle_timestep
