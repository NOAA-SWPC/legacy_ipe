!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!sep2012: efield.f was separated into each routin for SMS compatibility.
!--------------------------------------------  
      module module_svbksb
!--------------------------------------------------------------------- 
! description: calculates the electric potential for a given year,
!      day of year,UT, F10.7, B_z(K_p)
! - low/midlatitudes electric potential is from an empirical model from
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! - the transition zone is smoothed
! - output is the horizontal global electric field in magnetic coordinates direction
!  at every magnetic local time grid point expressed in degrees (0 deg-0MLT; 360deg 24 MLT)
!
! input 
!      integer :: iday,     ! day number of year
!                 iyear     ! year
!      real:: ut,       ! universal time 
!                 F10.7,    ! solar flux       (see ionosphere module)
!                 bz        ! component of IMF (see ionosphere module)
! output
!      real ::               &
!       ed1(0:nmlon,0:nmlat),    &  ! zonal electric field Ed1  [V/m] 
!       ed2(0:nmlon,0:nmlat)        ! meridional electric field Ed2/sin I_m  [V/m]  
!
! notes:
!
! - !to be done (commented out): input S_a F10.7/ Kp from WACCM and calculate B_z 
!    from these inputs
! - assume regular geomagnetic grid 
! - uses average year 365.24 days/year 30.6001 day/mo s. Weimer
! - get_tilt works only for iyear >= 1900
! - Weimer model 1996, Dan Weimer (not with the updates from B.Emery)
! - fixed parameters: B_z, B_y units nT  CHANGE THIS
!                     F10.7
! - we assume that the reference height is 300km for the emperical potential model
! - as a first approximation the electric field is constant in height
!   WATCH what is the upper boundary condition in WACCM
! - for all the calculation done here we set the reference height to the same 
!   value as in tiegcm (hr=130km)
! - 12/15/03 input value iseasav : replaced by day -> month and day of month
! - 12/15/03 S_aM calculated according to Scherliess draft paper and added
!   S_aM(corrected) = 90*(S_aM+1) to get variation in fig 1 Scherliess draft
!
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

c     use shr_kind_mod,  only: r8 => shr_kind_r8
c     use physconst,     only: pi
c     use abortutils,    only: endrun
c     use cam_logfile,   only: iulog
   
      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: svbksb

      contains
!-----------------------------------------------------------------------
! purpose: solves a*x = b
!
! method:     
! solves a*x = b for a vector x, where a is specified by the arrays
! u,w,v as returned by svdcmp. m and n
! are the logical dimensions of a, and will be equal for square matrices.
! mp and np are the physical dimensions of a. b(1:m) is the input right-hand 
! side. x(1:n) is the output solution vector. no input quantities are 
! destroyed, so the routine may be called sequentially with different b
!
! author:  a. maute dec 2002   
! (* copyright (c) 1985 numerical recipes software -- svbksb *!
! from numerical recipes 1986 pp. 57 or can be find on web-sites
!-------------------------------------------------------------------------      

      subroutine svbksb( u, w, v, m, n, mp, np, b, x )
!------------------------------------------------------------------------- 
!	... dummy arguments
!------------------------------------------------------------------------- 
      implicit none
      integer, parameter :: nmax = 1600
      integer, intent(in)   :: m
      integer, intent(in)   :: n
      integer, intent(in)   :: mp
      integer, intent(in)   :: np
      real , intent(in)  :: u(mp,np)
      real , intent(in)  :: w(np)
      real , intent(in)  :: v(np,np)
      real , intent(in)  :: b(mp)
      real , intent(out) :: x(np)

!------------------------------------------------------------------------- 
!	... local variables
!------------------------------------------------------------------------- 
      integer  :: i, j, jj
      real :: s
      real :: tmp(nmax)

      do j = 1,n
        s = 0. 
        if( w(j) /= 0. ) then
          do i = 1,m
            s = s + u(i,j)*b(i)
          end do
          s = s/w(j)
        endif
        tmp(j) = s
      end do

      do j = 1,n
        s = 0. 
        do jj = 1,n
          s = s + v(j,jj)*tmp(jj)
        end do
        x(j) = s
      end do

      end subroutine svbksb

!-----------------------------------------------------------------------
      end module module_svbksb
