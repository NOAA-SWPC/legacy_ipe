!dbg20120501: add v// to perp transport
!20110911: note: openmp was tried on jet but did not work: only thread 0 was used not the other thread...although other threads did exist...needs more investigation...
! DATE: 08 September, 2011
!********************************************
!***      Copyright 2011 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!--------------------------------------------  
      MODULE module_sub_PLASMA
      USE module_precision
      USE module_IPE_dimension,ONLY: ISPEC,ISPET,ISPEV,IPDIM,NLP,NMP,ISTOT
      USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_3d,VEXBup,plasma_3d_old, MaxFluxTube
#ifdef TESTING
      USE ModelDataInstances_Class
      USE module_MDI
#endif
      IMPLICIT NONE
      include "gptl.inc"

      PRIVATE
      PUBLIC :: plasma 

      CONTAINS
!---------------------------
      SUBROUTINE plasma ( utime, timestamp_for_IPE )
      USE module_input_parameters,ONLY:mpstop,ip_freq_output,start_time,stop_time,&
&     sw_neutral_heating_flip,sw_perp_transport,lpmin_perp_trans,lpmax_perp_trans,sw_para_transport,sw_debug,        &
&     sw_dbg_perp_trans,sw_exb_up,parallelBuild,mype, &
& HPEQ_flip, ut_start_perp_trans
      USE module_physical_constants,ONLY:rtd,zero
      USE module_FIELD_LINE_GRID_MKS,ONLY:JMIN_IN,plasma_grid_3d,plasma_grid_GL,plasma_grid_Z,JMAX_IS,hrate_mks3d
      USE module_PLASMA,ONLY:utime_save,plasma_1d
      USE module_perpendicular_transport,ONLY:perpendicular_transport
      IMPLICIT NONE
!------------------------
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
!--- local variables ---
      INTEGER (KIND=int_prec) :: mp
      INTEGER (KIND=int_prec) :: lp
      INTEGER (KIND=int_prec) :: i,j,midpoint, i1d,k,ret  
      INTEGER (KIND=int_prec) :: jth  
      integer :: status
      character(len=13), INTENT(IN) :: timestamp_for_IPE
      REAL(KIND=real_prec)          :: localPlasma_3d(1:ISTOT,1:MaxFluxTube,1:NLP,1:NMP) 



! save ut so that other subroutines can refer to it
      utime_save=utime
      print *, 'GHGM IN PLASMA, UTIME : ', utime

      ret = gptlstart ('apex_lon_loop') !24772.857
!SMS$PARALLEL(dh, lp, mp) BEGIN
      plasma_3d_old = plasma_3d
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-1")
      ret = gptlstart ('exchange_barrier')
!sms$insert      call ppp_barrier(status)
      ret = gptlstop ('exchange_barrier')
      ret = gptlstart ('EXCHANGE')
!SMS$EXCHANGE(plasma_3d_old)
      ret = gptlstop  ('EXCHANGE')
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-2")

      IF ( sw_neutral_heating_flip==1 )  hrate_mks3d(:,:,:,:)=zero
!        if ( sw_debug )  WRITE (0,"('sub-p: mp=',I4)")mp
      IF( sw_dbg_perp_trans .AND. utime==start_time .AND. parallelBuild )THEN
        print*,'sw_dbg_perp_trans=true does not work for parallel runs'
        print*,'Stopping module_sub_plasma'
        STOP
      ENDIF

      DO mp = 1,mpstop
        DO lp = 1,NLP
          DO i=JMIN_IN(lp),JMAX_IS(lp)
             i1d=i-JMIN_IN(lp)+1
             DO jth=1,ISTOT
                localPlasma_3d(jth,i,lp,mp) = plasma_3d(i,lp,mp,jth)
             END DO !jth
          END DO !i
        ENDDO
      ENDDO

PRINT*, SHAPE( plasma_3d )
PRINT*, SHAPE( localplasma_3d )
STOP
      IF ( HPEQ_flip==0.5 .AND. utime==ut_start_perp_trans ) THEN

        print *, 'utime=',utime,&
                 ' plasma perp transport is not called when HPEQ_flip=0.5'// &
                 ' start_time=0 because initial profiles do not exist!'

      ELSE IF ( utime>0 ) THEN


        IF ( sw_perp_transport>=1 ) THEN

          DO mp = 1,mpstop
            DO lp = lpmin_perp_trans, lpmax_perp_trans
  
              DO i=JMIN_IN(lp),JMAX_IS(lp)
                 i1d=i-JMIN_IN(lp)+1
                 DO jth=1,ISTOT
                    plasma_1d(jth,i1d) = localPlasma_3d(jth,i,lp,mp)
                 END DO !jth
              END DO !i

              ret = gptlstart ('perp_transport')
              CALL perpendicular_transport ( utime,mp,lp )
              ret = gptlstop ('perp_transport')

              DO i=JMIN_IN(lp),JMAX_IS(lp)
                 i1d=i-JMIN_IN(lp)+1
                 DO jth=1,ISTOT
                    localPlasma_3d(jth,i,lp,mp) = plasma_1d(jth,i1d)
                 END DO !jth
              END DO !i
  
              

            ENDDO
          ENDDO

        ENDIF
      ENDIF

#ifdef TESTING
       CALL mdi % Update( "module_sub_plasma.f90", &           ! Module name -here
                          "plasma", &                          ! Subroutine name -
                          "perpendicular_transport update", &  ! Unique name of
                          133, &                               ! Line number near
                          SIZE(localPlasma_3d), &              ! The total number
                          PACK(localPlasma_3d,.TRUE.) )        ! The array that
#endif

! update the boundary conditions if the top of the flux tube is open
!t        CALL update_flux_tube_boundary_condition ( )

      IF ( sw_para_transport==1 ) THEN 
      DO mp = 1,mpstop
        DO lp = 1,NLP

          DO i=JMIN_IN(lp),JMAX_IS(lp)
             i1d=i-JMIN_IN(lp)+1
             DO jth=1,ISTOT
                plasma_1d(jth,i1d) = localPlasma_3d(jth,i,lp,mp)
             END DO !jth
          END DO !i

          ret = gptlstart ('flux_tube_solver')
          CALL flux_tube_solver ( utime,mp,lp )
          ret = gptlstop ('flux_tube_solver')

        ENDDO
      ENDDO

      ELSE IF ( sw_para_transport==0 ) THEN 

      DO mp = 1,mpstop
        DO lp = 1,NLP

          DO i=JMIN_IN(lp),JMAX_IS(lp)
            i1d=i-JMIN_IN(lp)+1
             DO jth=1,ISTOT
                plasma_3d(i,lp,mp,jth) = localPlasma_3d(jth,i,lp,mp)
             END DO !jth
           END DO !i

        END DO !: DO lp = 1
      END DO !: DO mp = 
      ENDIF
#ifdef TESTING
       CALL mdi % Update( "module_sub_plasma.f90", &      ! Module name -here
                          "plasma", &                     ! Subroutine name -
                          "flux_tube_solver update", &    ! Unique name of
                          178, &               ! Line number near
                          SIZE(plasma_3d), &              ! The total number
                          PACK(plasma_3d,.TRUE.) )        ! The array that
#endif


!SMS$PARALLEL END
      ret = gptlstop ('apex_lon_loop')
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-4")

!dbg20120228: debug how2validate the transport
!dbg20120501 if(sw_dbg_perp_trans) call dbg_estimate_trans_error (utime)

! output plasma parameters to a file
      ret = gptlstart ('io_plasma_bin')
write(6,*)'BEFORE MOD check output plasma',utime,start_time,ip_freq_output
! ghgm output every time for now......
!      IF ( MOD( (utime-start_time),ip_freq_output)==0 ) THEN 
!
write(6,*)'before call to output plasma',utime,start_time,ip_freq_output
!dbg20110923segmentation fault??? memory allocation run time error???
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-5")
        CALL io_plasma_bin ( 1, utime, timestamp_for_IPE)            
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-6")

      ret = gptlstop ('io_plasma_bin')
      
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-7")

      END SUBROUTINE plasma
      END MODULE module_sub_PLASMA
