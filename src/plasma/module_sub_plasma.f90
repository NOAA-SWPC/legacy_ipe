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
 USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_3d,VEXBup,plasma_3d_old
 USE module_input_parameters,ONLY:mpstop,ip_freq_output,start_time,stop_time,&
                                  sw_neutral_heating_flip,sw_perp_transport, &
                                  lpmin_perp_trans,lpmax_perp_trans, &
                                  sw_para_transport,sw_debug,        &
                                  sw_dbg_perp_trans,sw_exb_up,parallelBuild,mype, &
                                   HPEQ_flip, ut_start_perp_trans, dumpFrequency
 USE module_physical_constants,ONLY:rtd,zero
 USE module_FIELD_LINE_GRID_MKS,ONLY:JMIN_IN,plasma_grid_3d,plasma_grid_GL, &
                                     plasma_grid_Z,JMAX_IS,hrate_mks3d
 USE module_PLASMA,ONLY:utime_save
 USE module_perpendicular_transport,ONLY:perpendicular_transport

 IMPLICIT NONE

 CONTAINS

! ----------------------------------------------------------------------------------- !

   SUBROUTINE plasma ( utime_local )


     IMPLICIT NONE

     INTEGER (KIND=int_prec), INTENT(IN) :: utime_local !universal time [sec]

     INTEGER (KIND=int_prec) :: mp, lp, i, j, k, jth
     INTEGER (KIND=int_prec) :: midpoint, i1d, ret  
     INTEGER                 :: pppStatus

!---------------

      utime_save=utime_local

!SMS$PARALLEL(dh, lp, mp) BEGIN
      plasma_3d_old = plasma_3d
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-1")
!sms$insert      call ppp_barrier(pppStatus)
!SMS$EXCHANGE(plasma_3d_old)
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-2")


      IF( sw_neutral_heating_flip==1 )THEN
        ! SMS will likely hate this....
        hrate_mks3d(:,:,:,:)=0.0_real_prec
      ENDIF

      IF( sw_dbg_perp_trans .and. parallelBuild )THEN
        PRINT*,'sw_dbg_perp_trans=true does not work for parallel runs'
        PRINT*,'Stopping module_sub_plasma'
        STOP
      ENDIF

      IF ( sw_perp_transport>=1 ) THEN
        DO mp = 1,mpstop
          DO lp = lpmin_perp_trans,lpmax_perp_trans

#ifdef DEBUG
            WRITE (0,"('sub-p: lp=',I4)")lp
#endif

              CALL perpendicular_transport ( utime_local,mp,lp )

          END DO 
        END DO  

      END IF 
          

      IF ( sw_para_transport==1 ) THEN 

        !DO mp = 1,mpstop
        !  DO lp = 1,NLP
        !    CALL flux_tube_solver ( utime_local,mp,lp )
            CALL flux_tube_solver ( utime_local )
        !  END DO 
        !END DO  

      END IF           

!SMS$PARALLEL END

!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-4")
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-6")
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-7")

      END SUBROUTINE plasma
      END MODULE module_sub_PLASMA
