!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
!      module sub_efield
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
   
!      implicit none

!nm20121003:module parameters are separated into efield.f!
!      end module sub_efield
!
!      
! Purpose: 
! Subroutines to calculate the electric potentials from the Weimer '96 model of
! the polar cap ionospheric electric potentials.
!
! Method:
!
! To use, first call subroutine ReadCoef once.
! Next, call SetModel with the specified input parameters.
! The function EpotVal(gLAT,gMLT) can then be used repeatively to get the
! electric potential at the desired location in geomagnetic coordinates.
! Subroutines to calculate the electric potentials from the Weimer '96 model of
! the polar cap ionospheric electric potentials.
!
!
! Author: A. Maute Dec 2003  
! This code is protected by copyright and is
! distributed for research or educational use only.
! Commerical use without written permission from Dan Weimer/MRC is prohibited.
!

!================================================================================================

!======================================================================



!================================================================================================



!================================================================================================



!================================================================================================


!================================================================================================

!================================================================================================

	FUNCTION MLT(MagLong)
!
!-----------------------------------------------------------------------
! given magnetic longitude in degrees, return Magnetic Local Time
! assuming that TRANS has been called with the date & time to calculate
! the rotation matrices.
!
! btf 11/06/03:
! Call sub adjust instead of referencing it as a function
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        USE module_ADJUST, ONLY: ADJUST
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real mlt
!
!-------------------------------Commons---------------------------------
!
        real cx, st, ct, am
	COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)

!
!------------------------------Arguments--------------------------------
!
	REAL MagLong
!
!---------------------------Local variables-----------------------------
!
	REAL angle, rotangle
!
!-----------------------------------------------------------------------
!
	RotAngle=CX(7)
!       MLT=ADJUST(Maglong+RotAngle+180.)/15.
        angle = Maglong+RotAngle+180.
        call adjust(angle)
        mlt = angle/15.
	RETURN
	END FUNCTION MLT

!================================================================================================

	FUNCTION MagLong(MLT)
!
!-----------------------------------------------------------------------
! return magnetic longitude in degrees, given Magnetic Local Time
! assuming that TRANS has been called with the date & time to calculate
! the rotation matrices.
!
! btf 11/06/03:
! Call sub adjust instead of referencing it as a function
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        USE module_ADJUST, ONLY: ADJUST
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real MagLong
!
!-------------------------------Commons---------------------------------
!
        real cx, st, ct, am
        COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
!
!------------------------------Arguments--------------------------------
!
	REAL MLT
!
!---------------------------Local variables-----------------------------
!
	REAL angle, rotangle
!
!-----------------------------------------------------------------------
!
	RotAngle=CX(7)
	angle=MLT*15.-RotAngle-180.
!       MagLong=ADJUST(angle)
        call adjust(angle)
        MagLong = angle
	RETURN
	END FUNCTION MagLong

!================================================================================================

	SUBROUTINE SunLoc(SunLat,SunLong)
!
!-----------------------------------------------------------------------
! Return latitude and longitude of sub-solar point.
! Assumes that TRANS has previously been called with the
! date & time to calculate the rotation matrices.
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-------------------------------Commons---------------------------------
!
        real cx, st, ct, am
        COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
!
!------------------------------Arguments--------------------------------
!
	Real SunLat,SunLong
!
!-----------------------------------------------------------------------
!
	SunLong=CX(9)
	SunLat=CX(8)
	RETURN
	END SUBROUTINE SunLoc

!================================================================================================

      SUBROUTINE GECMP (AMLA,RMLT,ET,EP)
!
!-----------------------------------------------------------------------
!          Get Electric field components for the Weimer electrostatic
!          potential model.  Before use, first load coefficients (CALL
!          READCOEF) and initialize model conditions (CALL SETMODEL).
!
!          INPUTS:
!            AMLA = Absolute value of magnetic latitude (deg)
!            RMLT = Magnetic local time (hours).
!          RETURNS:
!            ET = Etheta (magnetic equatorward*) E field component (V/m)
!            EP = Ephi   (magnetic eastward)     E field component (V/m)
!
!          * ET direction is along the magnetic meridian away from the
!            current hemisphere; i.e., when ET > 0, the direction is
!              southward when RMLA > 0
!              northward when RMLA < 0
!
!          NCAR addition (Jan 97).  R.Barnes
!-----------------------------------------------------------------------
!
c     use shr_kind_mod, only: r8 => shr_kind_r8
      USE module_EpotVal ,ONLY: EpotVal
      implicit none 
!
!-------------------------------Commons---------------------------------
!
!          CECMP contains constants initialized in READCOEF
      real alamn, alamx, alamr, stpd, stp2, cstp, sstp
      COMMON /CECMP/ ALAMN,ALAMX,ALAMR,STPD,STP2,CSTP,SSTP
!
!------------------------------Arguments--------------------------------
!
      real amla, rmlt, et, ep
!
!-----------------------------Parameters------------------------------
!
      real d2r, r2d
      PARAMETER ( D2R =  0.0174532925199432957692369076847 , 
     &           R2D = 57.2957795130823208767981548147)
!
!---------------------------Local variables-----------------------------
!
      real p1, p2
      real xmlt, xmlt1, kpol, dphi, amla1
!
!-------------------------External Functions----------------------------
!
!nm20121012      real epotval
!nm20121012      external epotval
!
!-----------------------------------------------------------------------
!
      ET = -99999.
      EP = -99999.
      IF (AMLA .LT. 0.) GO TO 100

!          Calculate -(latitude gradient) by stepping 10 km along the
!          meridian in each direction (flipping coordinates when going
!          over pole to keep lat <= 90).
      KPOL  = 0
      XMLT  = RMLT
   10 XMLT1 = XMLT
      AMLA1 = AMLA + STPD
      IF (AMLA1 .GT. 90.) THEN
	AMLA1 = 180. - AMLA1
	XMLT1 = XMLT1 + 12.
      ENDIF
      P1 = EPOTVAL (AMLA1    ,XMLT1)
      P2 = EPOTVAL (AMLA-STPD,XMLT )
      IF (KPOL .EQ. 1) GO TO 20
      ET = (P1 - P2) / STP2

!          Calculate -(lon gradient).  For most latitudes, step along a
!          great circle.  However, limit minimum latitude to the model
!          minimum (distorting the path onto a latitude line).  Also,
!          avoid a divide by zero at the pole avoid by using Art's trick
!          where Ephi(90,lon) = Etheta(90,lon+90)
      IF (AMLA .LT. ALAMX) THEN
	AMLA1 = MAX (ASIN(SIN(AMLA*D2R)*CSTP) , ALAMR)
	DPHI  = ASIN (SSTP/SIN(AMLA1))*R2D
	AMLA1 = AMLA1*R2D
	P1 = EPOTVAL (AMLA1,XMLT+DPHI)
	P2 = EPOTVAL (AMLA1,XMLT-DPHI)
      ELSE
	AMLA = 90.
	XMLT = XMLT + 6.
	KPOL = 1
	GO TO 10
      ENDIF
   20 EP = (P2 - P1) / STP2
      IF (KPOL .EQ. 1) EP = -EP

!          Below model minimum lat, the potential is value at min lat
      IF (AMLA .LT. ALAMN) THEN
	ET = 0.
	EP = EP * COS(ALAMR)/COS(AMLA*D2R)
      ENDIF

  100 RETURN
      END SUBROUTINE GECMP

!=====================================================================
      subroutine svdcmp( a, m, n, mp, np, w, v )
!------------------------------------------------------------------------- 
! purpose: singular value decomposition
!
! method:
! given a matrix a(1:m,1:n), with physical dimensions mp by np,
! this routine computes its singular value decomposition,
! the matrix u replaces a on output. the
! diagonal matrix of singular values w is output as a vector
! w(1:n). the matrix v (not the transpose v^t) is output as
! v(1:n,1:n).
!
! author: a. maute dec 2003      
! (* copyright (c) 1985 numerical recipes software -- svdcmp *!
! from numerical recipes 1986 pp. 60 or can be find on web-sites
!------------------------------------------------------------------------- 
      implicit none
      integer, parameter :: nmax = 1600
!------------------------------------------------------------------------- 
!	... dummy arguments
!------------------------------------------------------------------------- 
      integer, intent(in)     :: m
      integer, intent(in)     :: n
      integer, intent(in)     :: mp
      integer, intent(in)     :: np
      real, intent(inout) :: a(mp,np)
      real, intent(out)   :: v(np,np)
      real, intent(out)   :: w(np)

!------------------------------------------------------------------------- 
!	... local variables
!------------------------------------------------------------------------- 
      integer  :: i, its, j, k, l, nm
      real :: anorm
      real  :: c
      real  :: f
      real  :: g
      real  :: h
      real  :: s
      real  :: scale
      real  :: x, y, z
      real  :: rv1(nmax)
      logical  :: cnd1
      logical  :: cnd2

      g     = 0.0
      scale = 0.0
      anorm = 0.0

      do i = 1,n  !loop1
        l = i + 1
        rv1(i) = scale*g
        g     = 0.0
        s     = 0.0
        scale = 0.0
        if( i <= m ) then
          do k = i,m
            scale = scale + abs(a(k,i))
          end do
          if( scale /= 0.0 ) then
            do k = i,m
              a(k,i) = a(k,i)/scale
              s = s + a(k,i)*a(k,i)
            end do
            f = a(i,i)
            g = -sign(sqrt(s),f)
            h = f*g - s
            a(i,i) = f - g
            if( i /= n ) then
              do j = l,n
                s = 0.0
                do k = i,m
                  s = s + a(k,i)*a(k,j)
                end do
                f = s/h
                do k = i,m
                  a(k,j) = a(k,j) + f*a(k,i)
                end do
              end do
            end if
            do k = i,m
              a(k,i) = scale*a(k,i)
            end do
          endif
        endif
        w(i) = scale *g
        g     = 0.0
        s     = 0.0
        scale = 0.0
        if( i <= m .and. i /= n ) then
          do k = l,n
            scale = scale + abs(a(i,k))
          end do
          if( scale /= 0.0 ) then
            do k = l,n
              a(i,k) = a(i,k)/scale
              s      = s + a(i,k)*a(i,k)
            end do
            f = a(i,l)
            g = -sign(sqrt(s),f)
            h = f*g - s
            a(i,l) = f - g
            do k = l,n
              rv1(k) = a(i,k)/h
            end do
            if( i /= m ) then
              do j = l,m
                s = 0.0
                do k = l,n
                  s = s + a(j,k)*a(i,k)
                end do
                do k = l,n
                  a(j,k) = a(j,k) + s*rv1(k)
                end do
              end do
            end if
            do k = l,n
              a(i,k) = scale*a(i,k)
            end do
          end if
        end if
        anorm = max( anorm,(abs(w(i)) + abs(rv1(i))) )
      enddo !loop1

      do i = n,1,-1
        if( i < n ) then
          if( g /= 0.0 ) then
            do j = l,n
              v(j,i) = (a(i,j)/a(i,l))/g
            end do
            do j = l,n
              s = 0.0
              do k = l,n
                s = s + a(i,k)*v(k,j)
              end do
              do k = l,n
                v(k,j) = v(k,j) + s*v(k,i)
              end do
            end do
          end if
          do j = l,n
            v(i,j) = 0.0
            v(j,i) = 0.0
          end do
        end if
        v(i,i) = 1.0
        g = rv1(i)
        l = i
      end do

      do i = n,1,-1
        l = i + 1
        g = w(i)
        if( i < n ) then
          do j = l,n
            a(i,j) = 0.0
          end do
        end if
        if( g /= 0.0  ) then
          g = 1./g
          if( i /= n ) then
            do j = l,n
              s = 0.0
              do k = l,m
                s = s + a(k,i)*a(k,j)
              end do
              f = (s/a(i,i))*g
              do k = i,m
                a(k,j) = a(k,j) + f*a(k,i)
              end do
            end do
          end if
          do j = i,m
            a(j,i) = a(j,i)*g
          end do
        else
          do j = i,m
            a(j,i) = 0.0
          end do
        end if
        a(i,i) = a(i,i) + 1.0
      end do

      do k = n,1,-1
        do its = 1,30 !loop2
          do l = k,1,-1
            nm = l - 1
            cnd1 = abs( rv1(l) ) + anorm == anorm
            if( cnd1 ) then
              cnd2 = .false.
              exit
            end if
            cnd2 = abs( w(nm) ) + anorm == anorm
            if( cnd2 ) then
              cnd1 = .true.
              exit
            else if( l == 1 ) then
              cnd1 = .true.
              cnd2 = .true.
            end if
          end do

          if( cnd2 ) then
            c = 0.0
            s = 1.0
            do i = l,k
              f = s*rv1(i)
              if( (abs(f) + anorm) /= anorm ) then
                g = w(i)
                h = sqrt(f*f + g*g)
                w(i) = h
                h = 1.0/h
                c = (g*h)
                s = -(f*h)
                do j = 1,m
                  y = a(j,nm)
                  z = a(j,i)
                  a(j,nm) = (y*c) + (z*s)
                  a(j,i) = -(y*s) + (z*c)
                end do
              end if
            end do
          end if

          if( cnd1 ) then
            z = w(k)
            if( l == k ) then
              if( z < 0.0 ) then
                w(k) = -z
                do j = 1,n
                  v(j,k) = -v(j,k)
                end do
              end if
c             exit loop2
              go to 20
            end if
          end if

          x = w(l)
          nm = k - 1
          y = w(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0*h*y)
          g = sqrt( f*f + 1.0 )
          f = ((x - z)*(x + z) + h*((y/(f + sign(g,f))) - h))/x
          c = 1.0
          s = 1.0
          do j = l,nm
            i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g
            z = sqrt( f*f + h*h )
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c)+(g*s)
            g = -(x*s)+(g*c)
            h = y*s
            y = y*c
            do nm = 1,n
              x = v(nm,j)
              z = v(nm,i)
              v(nm,j) = (x*c)+(z*s)
              v(nm,i) = -(x*s)+(z*c)
            end do
            z = sqrt( f*f + h*h )
            w(j) = z
            if( z /= 0.0 ) then
              z = 1.0/z
              c = f*z
              s = h*z
            end if
            f = (c*g)+(s*y)
            x = -(s*g)+(c*y)
            do nm = 1,m
              y = a(nm,j)
              z = a(nm,i)
              a(nm,j) = (y*c)+(z*s)
              a(nm,i) = -(y*s)+(z*c)
            end do
          end do
          rv1(l) = 0.0
          rv1(k) = f
          w(k)   = x
        end do  !loop2
   20 continue
      end do
      
      end subroutine svdcmp

!-------------------------------------------------------------------------      

