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
!! \brief Module of routines that correspond to resistivity
!!
!! In this module following namelist of parameters is specified:
!! \copydetails resistivity::init_resistivity
!<
module resistivity
! pulled by RESISTIVE
   use constants, only: dsetnamelen
   use types,     only: value

   implicit none

   private
   public  :: init_resistivity, timestep_resist, cleanup_resistivity, etamax, diffuseb, cu2max, deimin, eta1_active, diffuse_mag

   real                                  :: cfl_resist                     !< CFL factor for resistivity effect
   real                                  :: eta_0                          !< uniform resistivity
   real                                  :: eta_1                          !< anomalous resistivity
   real                                  :: j_crit                         !< critical value of current density
   real                                  :: jcrit2                         !< squared critical value of current density
   real                                  :: deint_max                      !< COMMENT ME
   real                                  :: eta_weight                     !< weight for smoothing eta; no smoothing if negative
   real                                  :: d_eta_factor
   type(value)                           :: etamax, cu2max, deimin
   logical, save                         :: eta1_active = .true.           !< resistivity off-switcher while eta_1 == 0.0
   character(len=dsetnamelen), parameter :: eta_n = "eta", jcu_n = "jcu2", dei_n = "dei"

contains

   subroutine cleanup_resistivity
      implicit none
   end subroutine cleanup_resistivity

!>
!! \brief Routine to set parameters values from namelist RESISTIVITY
!!
!! \n \n
!! @b RESISTIVITY
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>cfl_resist</td><td>0.4  </td><td>real value   </td><td>\copydoc resistivity::cfl_resist</td></tr>
!! <tr><td>eta_0     </td><td>0.0  </td><td>real value   </td><td>\copydoc resistivity::eta_0     </td></tr>
!! <tr><td>eta_1     </td><td>0.0  </td><td>real value   </td><td>\copydoc resistivity::eta_1     </td></tr>
!! <tr><td>eta_weight</td><td>4    </td><td>integer value</td><td>\copydoc resistivity::eta_weight</td></tr>
!! <tr><td>j_crit    </td><td>1.0e6</td><td>real value   </td><td>\copydoc resistivity::j_crit    </td></tr>
!! <tr><td>deint_max </td><td>0.01 </td><td>real value   </td><td>\copydoc resistivity::deint_max </td></tr>
!! </table>
!! The list is active while \b "RESISTIVE" is defined.
!! \n \n
!<
   subroutine init_resistivity

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use constants,        only: PIERNIK_INIT_GRID, GEO_XYZ, wcu_n, zero
      use dataio_pub,       only: die, code_progress, nh
      use domain,           only: dom
      use func,             only: operator(.notequals.)
      use mpisetup,         only: rbuff, master, slave, piernik_MPI_Bcast
      use named_array_list, only: qna
      use types,            only: value
#if !defined(IONIZED) && defined(ISO)
      use dataio_pub,       only: warn
#endif /* !IONIZED || ISO */

      implicit none

      type(cg_list_element), pointer :: cgl

      namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, eta_weight, j_crit, deint_max

      if (code_progress < PIERNIK_INIT_GRID) call die("[resistivity:init_resistivity] grid not initialized.")
      if (dom%geometry_type /= GEO_XYZ)      call die("[resistivity:init_resistivity] Unsupported geometry")

      cfl_resist = 0.4
      eta_0      = 0.0
      eta_1      = 0.0
      eta_weight = 4.0
      j_crit     = 1.0e6
      deint_max  = 0.01

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=RESISTIVITY)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=RESISTIVITY, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "RESISTIVITY")
         read(nh%cmdl_nml,nml=RESISTIVITY, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "RESISTIVITY", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=RESISTIVITY)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = cfl_resist
         rbuff(2) = eta_0
         rbuff(3) = eta_1
         rbuff(4) = j_crit
         rbuff(5) = deint_max
         rbuff(6) = eta_weight

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         cfl_resist = rbuff(1)
         eta_0      = rbuff(2)
         eta_1      = rbuff(3)
         j_crit     = rbuff(4)
         deint_max  = rbuff(5)
         eta_weight = rbuff(6)

      endif

      call all_cg%reg_var(wcu_n)
      call all_cg%reg_var(eta_n)
#if !defined(ISO) && defined(IONIZED)
      call all_cg%reg_var(jcu_n)
      call all_cg%reg_var(dei_n)
#endif /* !ISO && IONIZED */

      cgl => leaves%first
      do while (associated(cgl))
         cgl%cg%q(qna%ind(eta_n))%arr = eta_0
         cgl => cgl%nxt
      enddo
      etamax = value(eta_0, 0., [0., 0., 0.], [0, 0, 0], 0_4)

      eta1_active = (eta_1 .notequals. zero)

      if (eta1_active) then
         jcrit2 = j_crit**2
         d_eta_factor = 1./(2.*dom%eff_dim + eta_weight)
#if !defined(IONIZED) && defined(ISO)
         call warn("[resistivity:init_resistivity] eta_1 is set, but IONIZED gas is not included or ISO is set.")
         eta1_active = .false.
#endif /* !IONIZED || ISO */
      endif

   end subroutine init_resistivity

#if !defined(ISO) && defined(IONIZED)
   subroutine compute_resist

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: zero
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real, dimension(:,:,:), pointer :: eta, jc2

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call compute_current_sq(cg)

         if (eta1_active) then
            eta => cg%q(qna%ind(eta_n))%arr
            jc2 => cg%q(qna%ind(jcu_n))%arr

!           eta(:,:,:) = eta_0 + eta_1 * sqrt( max(0.0,jc2(:,:,:) - jcrit2 ))
!           the above may cause FPE because compiler may transform it to max(0.0, sqrt(jc2(:,:,:)- jcrit2))
            where (jc2(:,:,:) - jcrit2 > zero)
               eta(:,:,:) = eta_0 + eta_1 * sqrt(jc2(:,:,:) - jcrit2)
            elsewhere
               eta(:,:,:) = eta_0
            endwhere

            if (eta_weight < 0.) call smooth_eta(cg)
         endif

         cgl => cgl%nxt
      enddo

   end subroutine compute_resist

!> \brief square current computing in cell corner step by step
   subroutine compute_current_sq(cg)

      use constants,        only: LO, HI, oneq, xdim, ydim, zdim
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      type(grid_container),   pointer :: cg
      real, dimension(:,:,:), pointer :: jc2
      real, dimension(cg%lhn(xdim,LO):cg%lhn(xdim,HI), cg%lhn(ydim,LO):cg%lhn(ydim,HI), cg%lhn(zdim,LO):cg%lhn(zdim,HI)) :: db, rotb

      jc2 => cg%q(qna%ind(jcu_n))%arr

!--- current_z **2
      if (dom%has_dir(xdim)) then
         rotb(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:) = (cg%b(ydim,cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:) - cg%b(ydim,cg%lhn(xdim,LO):cg%lhn(xdim,HI)-1,:,:)) * cg%idl(xdim)
         rotb(cg%lhn(xdim,LO),:,:) = rotb(cg%lhn(xdim,LO)+1,:,:)
      else ; rotb = 0.0 ; endif
      if (dom%has_dir(ydim)) then
         db(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI),:) = (cg%b(xdim,:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI),:) - cg%b(xdim,:,cg%lhn(ydim,LO):cg%lhn(ydim,HI)-1,:)) * cg%idl(ydim)
         db(:,cg%lhn(ydim,LO),:) = db(:,cg%lhn(ydim,LO)+1,:)
         rotb = rotb - db
      endif

      if (dom%has_dir(zdim)) then
         jc2(:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)) =                                              oneq*(rotb(:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)) + rotb(:,:,cg%lhn(zdim,LO):cg%lhn(zdim,HI)-1))**2
         jc2(:,:,cg%lhn(zdim,LO)) = jc2(:,:,cg%lhn(zdim,LO)+1)
      else
         jc2 = rotb**2
      endif

!--- current_x **2
      if (dom%has_dir(ydim)) then
         rotb(:,cg%lhn(ydim,LO)+1:cg%lhn(xdim,HI),:) = (cg%b(zdim,:,cg%lhn(ydim,LO)+1:cg%lhn(xdim,HI),:) - cg%b(zdim,:,cg%lhn(ydim,LO):cg%lhn(xdim,HI)-1,:)) * cg%idl(ydim)
         rotb(:,cg%lhn(ydim,LO),:) = rotb(:,cg%lhn(ydim,LO)+1,:)
      else ; rotb = 0.0 ; endif
      if (dom%has_dir(zdim)) then
         db(:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)) = (cg%b(ydim,:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)) - cg%b(ydim,:,:,cg%lhn(zdim,LO):cg%lhn(zdim,HI)-1)) * cg%idl(zdim)
         db(:,:,cg%lhn(zdim,LO)) = db(:,:,cg%lhn(zdim,LO)+1)
         rotb = rotb - db
      endif

      if (dom%has_dir(xdim)) then
         jc2(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:) = jc2(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:) + oneq*(rotb(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:) + rotb(cg%lhn(xdim,LO):cg%lhn(xdim,HI)-1,:,:))**2
         jc2(cg%lhn(xdim,LO),:,:) = jc2(cg%lhn(xdim,LO)+1,:,:)
      else
         jc2 = jc2 + rotb**2
      endif

!--- current_y **2
      if (dom%has_dir(zdim)) then
         rotb(:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)) = (cg%b(xdim,:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)) - cg%b(xdim,:,:,cg%lhn(zdim,LO):cg%lhn(zdim,HI)-1)) * cg%idl(zdim)
         rotb(:,:,cg%lhn(zdim,LO)) = rotb(:,:,cg%lhn(zdim,LO)+1)
      else ; rotb = 0.0 ; endif
      if (dom%has_dir(xdim)) then
         db(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:) = (cg%b(zdim,cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI),:,:) - cg%b(zdim,cg%lhn(xdim,LO):cg%lhn(xdim,HI)-1,:,:)) * cg%idl(xdim)
         db(cg%lhn(xdim,LO),:,:) = db(cg%lhn(xdim,LO)+1,:,:)
         rotb = rotb - db
      endif

      if (dom%has_dir(ydim)) then
         jc2(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI),:) = jc2(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI),:) + oneq*(rotb(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI),:) + rotb(:,cg%lhn(ydim,LO):cg%lhn(ydim,HI)-1,:))**2
         jc2(:,cg%lhn(ydim,LO),:) = jc2(:,cg%lhn(ydim,LO)+1,:)
      else
         jc2 = jc2 + rotb**2
      endif

   end subroutine compute_current_sq

   subroutine smooth_eta(cg)

      use constants,        only: LO, HI, xdim, ydim, zdim, zero
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      type(grid_container),   pointer :: cg
      real, dimension(:,:,:), pointer :: eta
      real, dimension(cg%lhn(xdim,LO):cg%lhn(xdim,HI), cg%lhn(ydim,LO):cg%lhn(ydim,HI), cg%lhn(zdim,LO):cg%lhn(zdim,HI)) :: eh

      eta => cg%q(qna%ind(eta_n))%arr

      eh = zero
      if (dom%has_dir(xdim)) then
         eh(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI)-1,:,:) = eh(cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI)-1,:,:) + eta(cg%lhn(xdim,LO):cg%lhn(xdim,HI)-2,:,:) + eta(cg%lhn(xdim,LO)+2:cg%lhn(xdim,HI),:,:)
         eh(cg%lhn(xdim,LO),:,:) = eh(cg%lhn(xdim,LO)+1,:,:) ; eh(cg%lhn(xdim,HI),:,:) = eh(cg%lhn(xdim,HI)-1,:,:)
      endif
      if (dom%has_dir(ydim)) then
         eh(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI)-1,:) = eh(:,cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI)-1,:) + eta(:,cg%lhn(ydim,LO):cg%lhn(ydim,HI)-2,:) + eta(:,cg%lhn(ydim,LO)+2:cg%lhn(ydim,HI),:)
         eh(:,cg%lhn(ydim,LO),:) = eh(:,cg%lhn(ydim,LO)+1,:) ; eh(:,cg%lhn(ydim,HI),:) = eh(:,cg%lhn(ydim,HI)-1,:)
      endif
      if (dom%has_dir(zdim)) then
         eh(:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)-1) = eh(:,:,cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)-1) + eta(:,:,cg%lhn(zdim,LO):cg%lhn(zdim,HI)-2) + eta(:,:,cg%lhn(zdim,LO)+2:cg%lhn(zdim,HI))
         eh(:,:,cg%lhn(zdim,LO)) = eh(:,:,cg%lhn(zdim,LO)+1) ; eh(:,:,cg%lhn(zdim,HI)) = eh(:,:,cg%lhn(zdim,HI)-1)
      endif
      eh = (eh + eta_weight * eta) * d_eta_factor

      where (eta > eta_0) eta = eh

   end subroutine smooth_eta
#endif /* !ISO && IONIZED */

!-----------------------------------------------------------------------

   subroutine timestep_resist(dt)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: big, zero, pMIN, DIVB_CT
      use grid_cont,        only: grid_container
      use func,             only: operator(.notequals.)
      use global,           only: divB_0_method
      use mpisetup,         only: piernik_MPI_Allreduce, piernik_MPI_Bcast
#if !defined(ISO) && defined(IONIZED)
      use constants,        only: MAXL, MINL, small, xdim, ydim, zdim
      use fluidindex,       only: flind
      use func,             only: ekin, emag
      use named_array_list, only: qna, wna
#endif /* !ISO && IONIZED */

      implicit none

      real, intent(inout)               :: dt
      type(cg_list_element),  pointer   :: cgl
      type(grid_container),   pointer   :: cg
      real                              :: dt_eta, dt_eint
#if !defined(ISO) && defined(IONIZED)
      real, dimension(:,:,:),   pointer :: eta, jc2, dei
      real, dimension(:,:,:,:), pointer :: uu, bb
#endif /* !ISO && IONIZED */

      dt_eta = big ; dt_eint = big
#if !defined(ISO) && defined(IONIZED)
      if (divB_0_method == DIVB_CT) then
         call compute_resist
         if (eta1_active) then
            call leaves%get_extremum(qna%ind(eta_n), MAXL, etamax)
            call piernik_MPI_Bcast(etamax%val)
            call leaves%get_extremum(qna%ind(jcu_n), MAXL, cu2max)
            call piernik_MPI_Bcast(cu2max%val)
         endif
      endif
#endif /* !ISO && IONIZED */

      if (etamax%val .notequals. zero) then
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            dt_eta = min(dt_eta, cfl_resist * cg%dxmn2 / (2. * etamax%val))
#if !defined(ISO) && defined(IONIZED)
            if (divB_0_method == DIVB_CT) then
               eta => cg%q(qna%ind(eta_n))%span(cg%ijkse)
               jc2 => cg%q(qna%ind(jcu_n))%span(cg%ijkse)
               dei => cg%q(qna%ind(dei_n))%span(cg%ijkse)
               uu => cg%w(wna%fi)%span(cg%ijkse)
               bb => cg%w(wna%bi)%span(cg%ijkse)
               dei = (uu(flind%ion%ien,:,:,:) - ekin(uu(flind%ion%imx,:,:,:), uu(flind%ion%imy,:,:,:), uu(flind%ion%imz,:,:,:), uu(flind%ion%idn,:,:,:)) - &
                     emag(bb(xdim,:,:,:), bb(ydim,:,:,:), bb(zdim,:,:,:)))/ (eta(:,:,:) * jc2 + small)
               dt_eint = min(dt_eint, deint_max * abs(minval(dei)))
            endif
#endif /* !ISO && IONIZED */
            cgl => cgl%nxt
         enddo
      endif

      call piernik_MPI_Allreduce(dt_eta, pMIN)
      if (divB_0_method == DIVB_CT) then
#if !defined(ISO) && defined(IONIZED)
         call piernik_MPI_Allreduce(dt_eint, pMIN)
         call leaves%get_extremum(qna%ind(dei_n), MINL, deimin)
         deimin%assoc = dt_eint
         cu2max%assoc = min(dt_eta, dt_eint)
#endif /* !ISO && IONIZED */
         etamax%assoc = dt_eta
      endif

      dt = min(dt, dt_eta, dt_eint)

   end subroutine timestep_resist

!-----------------------------------------------------------------------------
!>
!! \brief
!! \todo overload me or use class(*) if you dare
!<
   subroutine vanleer_limiter(f,a,b)

      implicit none

      real, dimension(:), intent(in)    :: a !< second order correction of left- or right- moving waves flux on the left cell boundary
      real, dimension(:), intent(in)    :: b !< second order correction of left- or right- moving waves flux on the right cell boundary
      real, dimension(:), intent(inout) :: f !< second order flux correction for left- or right- moving waves
      ! locals
      real, dimension(size(a,1))        :: c !< a*b

      c = a*b                                                                    !> \todo OPTIMIZE ME
      where (c > 0.0)
         f = f+2.0*c/(a+b)
      endwhere

   end subroutine vanleer_limiter

   subroutine tvdd_1d(b1d,eta1d,idi,dt,wcu1d)

      use constants,     only: half
      implicit none

      real, dimension(:), pointer, intent(in)  :: eta1d, b1d
      real, dimension(:), pointer, intent(out) :: wcu1d
      real, intent(in)                         :: idi,dt

      real, dimension(size(b1d))               :: w, wp, wm, b1
      integer                                  :: n

      n = size(b1d)
      w(2:n)    = eta1d(2:n) * ( b1d(2:n) - b1d(1:n-1) )*idi ;  w(1)  = w(2)
      b1(1:n-1) = b1d(1:n-1) + half*(w(2:n) - w(1:n-1))*dt*idi; b1(n) = b1(n-1)

      w(2:n)    = eta1d(2:n) * ( b1(2:n) - b1(1:n-1) )*idi   ; w(1)  = w(2)
      wp(1:n-1) = half*(w(2:n) - w(1:n-1))                   ; wp(n) = wp(n-1)
      wm(2:n)   = wp(1:n-1)                                  ; wm(1) = wm(2)

      call vanleer_limiter(w,wm,wp)
      wcu1d     = w*dt

   end subroutine tvdd_1d

!-------------------------------------------------------------------------------
!
! 6 routines have been substituted by one with parameters:
!   diffuseby_x  --> diffuseb(ibdir = ydim, sdir = xdim, etadir = zdim, emf = 'emfz', n1 = ydim, n2 = zdim)
!   diffusebz_x  --> diffuseb(ibdir = zdim, sdir = xdim, etadir = ydim, emf = 'emfy', n1 = ydim, n2 = zdim)
!   diffusebz_y  --> diffuseb(ibdir = zdim, sdir = ydim, etadir = xdim, emf = 'emfx', n1 = zdim, n2 = xdim)
!   diffusebx_y  --> diffuseb(ibdir = xdim, sdir = ydim, etadir = zdim, emf = 'emfz', n1 = zdim, n2 = xdim)
!   diffusebx_z  --> diffuseb(ibdir = xdim, sdir = zdim, etadir = ydim, emf = 'emfy', n1 = xdim, n2 = ydim)
!   diffuseby_z  --> diffuseb(ibdir = ydim, sdir = zdim, etadir = xdim, emf = 'emfx', n1 = xdim, n2 = ydim)

   subroutine diffuseb(ibdir, sdir)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, ndims, half, I_ONE, wcu_n, idm, INT4, LO, HI
      use domain,           only: dom
      use global,           only: dt
      use grid_cont,        only: grid_container
      use magboundaries,    only: bnd_emf
      use named_array_list, only: qna, wna

      implicit none

      integer(kind=4),  intent(in)      :: ibdir, sdir
      integer                           :: i1, i2
      integer(kind=4)                   :: n1, n2, etadir, dir, emf, wcu_i, eta_i
      integer(kind=4), dimension(ndims) :: idml, idmh
      real, dimension(:),    pointer    :: b1d, eta1d, wcu1d
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg

      n1 = I_ONE + mod(sdir    ,   ndims)
      n2 = I_ONE + mod(sdir+I_ONE, ndims)
      etadir = sum([xdim,ydim,zdim]) - ibdir - sdir

#if !defined(ISO) && defined(IONIZED)
      call compute_resist
#endif /* !ISO && IONIZED */

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         wcu_i = qna%ind(wcu_n)
         eta_i = qna%ind(eta_n)

         idmh(:) = cg%lhn(:,HI) - idm(:,etadir)
         idml(:) = cg%lhn(:,LO) + idm(:,etadir)
         cg%q(eta_i)%arr(cg%lhn(xdim,LO):idmh(xdim),cg%lhn(ydim,LO):idmh(ydim),cg%lhn(zdim,LO):idmh(zdim)) = half*(cg%q(eta_i)%span(cg%lhn(:,LO),idmh) + cg%q(eta_i)%span(idml,cg%lhn(:,HI)))

         do i1 = cg%lhn(n1,LO), cg%lhn(n1,HI)
            do i2 = cg%lhn(n2,LO), cg%lhn(n2,HI)
               b1d   => cg%w(wna%bi)%get_sweep(sdir,ibdir,i1,i2)
               eta1d => cg%q(eta_i )%get_sweep(sdir,      i1,i2)
               wcu1d => cg%q(wcu_i )%get_sweep(sdir,      i1,i2)
               call tvdd_1d(b1d, eta1d, cg%idl(sdir), dt, wcu1d)
            enddo
         enddo

         cgl => cgl%nxt
      enddo

      cgl => leaves%first
      do while (associated(cgl))
         do dir = xdim, zdim
            emf = idm(etadir,dir) + 2_INT4
            if (dom%has_dir(dir)) call bnd_emf(wcu_i, emf, dir, cgl%cg)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine diffuseb

   subroutine diffuse_mag

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: idm, ndims, xdim, zdim, two, I_TWO, LO, HI, DIVB_CT
      use domain,           only: dom
      use global,           only: dt, divB_0_method
      use grid_cont,        only: grid_container
      use named_array_list, only: wna

      implicit none

      integer(kind=4)                         :: dir
      integer(kind=4), dimension(ndims,LO:HI) :: i0, im, ip
      real, dimension(:,:,:,:), pointer       :: b0, bm, bp
      real                                    :: df
      type(cg_list_element),    pointer       :: cgl
      type(grid_container),     pointer       :: cg

      if (divB_0_method == DIVB_CT) return

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do dir = xdim, zdim
            if (.not.dom%has_dir(dir)) cycle

            df = two * eta_0 * dt / cg%dl(dir)**2

            i0(:,LO) = cg%lhn(:,LO) + idm(:,dir)
            i0(:,HI) = cg%lhn(:,HI) - idm(:,dir)
            im = cg%lhn ; im(:,HI) = im(:,HI) - I_TWO * idm(:,dir)
            ip = cg%lhn ; ip(:,LO) = ip(:,LO) + I_TWO * idm(:,dir)
            b0 => cg%w(wna%bi)%span(i0)
            bm => cg%w(wna%bi)%span(im)
            bp => cg%w(wna%bi)%span(ip)
            b0 = b0 + (bp + bm - two*b0) * df
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine diffuse_mag

end module resistivity
