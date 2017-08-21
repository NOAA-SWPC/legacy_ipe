!note:20120207: v36: used only activating the perp.transport gradually only during daytime...
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
!
PROGRAM  test_plasma

 USE module_precision
 USE module_input_parameters
 USE module_FIELD_LINE_GRID_MKS
 USE module_init_plasma_grid
 USE module_NEUTRAL_MKS 
 USE module_sub_PLASMA
!SMS$IGNORE BEGIN
 USE module_init_ELDYN
 USE module_sub_ELDYN
!SMS$IGNORE END
 USE module_open_output_files
 USE module_output
 USE module_close_files
 USE module_IPE_dimension
!SMS$IGNORE BEGIN
#ifdef TESTING
 USE module_MDI
 USE module_eldyn
#endif
!SMS$IGNORE END
 USE module_IO,ONLY: PRUNIT
 USE module_open_file,ONLY: open_file

   IMPLICIT NONE
!SMS$INSERT   include "mpif.h"
   INTEGER(KIND=int_prec)           :: utime_driver ! Universal Time [sec]
   INTEGER(KIND=int_prec)           :: iterate
   CHARACTER(8)                     :: iterChar
   INTEGER(KIND=int_prec),parameter :: luntmp=300
   INTEGER(KIND=int_prec)           :: istat,mp,ret,OLDcomm,key,status=0,color=0,mypeWorld
   INTEGER(KIND=int_prec)           :: resultlen ! length return from MPI routines
!SMS$INSERT character(len=mpi_max_processor_name) :: mynode ! node name


!SMS$IGNORE BEGIN
#ifdef TESTING
     CALL mdi % Build( )
#endif
!SMS$IGNORE END


! set up input parameters
!SMS$INSERT         MPI_COMM_IPE = MPI_COMM_WORLD
     CALL read_input_parameters ( )


! open Input/Output files
     CALL open_file ( 'FLIP_ERR', PRUNIT,'formatted','unknown') !open by all processors

!SMS$SERIAL BEGIN
     CALL open_output_files ( )
!SMS$SERIAL END

! set up plasma grids by reading file
     CALL init_plasma_grid ( )

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-1")

     IF ( sw_output_plasma_grid ) THEN
       PRINT *, 'sub-init_p: output plasma_grid'
       CALL output_plasma_grid ( )
     END IF

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-2")

! initialise the flux tubes from previous runs
     IF ( HPEQ_flip==0.0 ) THEN
       PRINT *,'before CALL io_plasma_bin finished! READ: start_time=', start_time,stop_time
       CALL io_plasma_bin ( 2, start_time )
       PRINT *,'after CALL io_plasma_bin finished! READ: start_time=', start_time,stop_time
     END IF


! initialization of electrodynamic module:
! read in E-field

     IF ( sw_perp_transport>=1 ) THEN
       CALL init_eldyn ( )
     ENDIF

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-4")

     iterate = 0

     DO utime_driver = start_time, stop_time, time_step
       iterate = iterate + 1

       PRINT*,'utime_driver=',utime_driver

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-5")


!SMS$IGNORE BEGIN
#ifdef TESTING

       CALL mdi % Update( "driver_ipe.f90", &
                          "eldyn", &
                          "ed1_90 before eldyn", &
                          0, &
                          SIZE(ed1_90), &
                          PACK(ed1_90,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "eldyn", &
                          "ed2_90 before eldyn", &
                          0, &
                          SIZE(ed2_90), &
                          PACK(ed2_90,.TRUE.) )   

#endif
!SMS$IGNORE END
       IF ( sw_perp_transport>=1.AND. MOD( (utime_driver-start_time),ip_freq_eldyn)==0 ) THEN
         CALL eldyn ( utime_driver )
       ENDIF

!SMS$IGNORE BEGIN
#ifdef TESTING
       CALL mdi % Update( "driver_ipe.f90", &
                          "eldyn", &
                          "ed1_90 after eldyn", &
                          0, &
                          SIZE(ed1_90), &
                          PACK(ed1_90,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "eldyn", &
                          "ed2_90 after eldyn", &
                          0, &
                          SIZE(ed2_90), &
                          PACK(ed2_90,.TRUE.) )   
#endif
!SMS$IGNORE end

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-6")

!SMS$IGNORE BEGIN
#ifdef TESTING
        
       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "ON_m3 before neutral", &
                          0, &
                          SIZE(ON_m3), &
                          PACK(ON_m3,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "HN_m3 before neutral", &
                          0, &
                          SIZE(HN_m3), &
                          PACK(HN_m3,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "N2N_m3 before neutral", &
                          0, &
                          SIZE(N2N_m3), &
                          PACK(N2N_m3,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "O2N_m3 before neutral", &
                          0, &
                          SIZE(O2N_m3), &
                          PACK(O2N_m3,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "N4S_m3 before neutral", &
                          0, &
                          SIZE(N4S_m3), &
                          PACK(N4S_m3,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "TN_k before neutral", &
                          0, &
                          SIZE(TN_k), &
                          PACK(TN_k,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "TINF_k before neutral", &
                          0, &
                          SIZE(TINF_k), &
                          PACK(TINF_k,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "Un_ms1 before neutral", &
                          0, &
                          SIZE(Un_ms1), &
                          PACK(Un_ms1,.TRUE.) )   

#endif
!SMS$IGNORE END
        ! update neutral 3D structure: 
        ! use MSIS/HWM to get the values in the flux tube grid
       IF ( MOD( (utime_driver-start_time),ip_freq_msis)==0 ) THEN 

!SMS$IGNORE BEGIN
#ifdef DEBUG
         PRINT *,'CALL MSIS',utime_driver,start_time, &
                  ip_freq_msis,(utime_driver-start_time), &
                  MOD( (utime_driver-start_time),ip_freq_msis)
#endif
!SMS$IGNORE END
         CALL neutral ( utime_driver )

       END IF

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-7")

!SMS$IGNORE BEGIN
#ifdef TESTING
        
       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "ON_m3 after neutral", &
                          0, &
                          SIZE(ON_m3), &
                          PACK(ON_m3,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "HN_m3 after neutral", &
                          0, &
                          SIZE(HN_m3), &
                          PACK(HN_m3,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "N2N_m3 after neutral", &
                          0, &
                          SIZE(N2N_m3), &
                          PACK(N2N_m3,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "O2N_m3 after neutral", &
                          0, &
                          SIZE(O2N_m3), &
                          PACK(O2N_m3,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "N4S_m3 after neutral", &
                          0, &
                          SIZE(N4S_m3), &
                          PACK(N4S_m3,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "TN_k after neutral", &
                          0, &
                          SIZE(TN_k), &
                          PACK(TN_k,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "TINF_k after neutral", &
                          0, &
                          SIZE(TINF_k), &
                          PACK(TINF_k,.TRUE.) )   

       CALL mdi % Update( "driver_ipe.f90", &
                          "neutral", &
                          "Un_ms1 after neutral", &
                          0, &
                          SIZE(Un_ms1), &
                          PACK(Un_ms1,.TRUE.) )   

#endif
!SMS$IGNORE END

!SMS$IGNORE BEGIN
#ifdef TESTING
       CALL mdi % Update( "driver_ipe.f90", &
                          "plasma", &
                          "plasma_3d begin plasma", &
                          0, &
                          SIZE(plasma_3d), &
                          PACK(plasma_3d,.TRUE.) )   
#endif
!SMS$IGNORE END

        CALL plasma ( utime_driver, 'dummytimestam' )

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-8")

!SMS$IGNORE BEGIN
#ifdef TESTING
       CALL mdi % Update( "driver_ipe.f90", &
                          "plasma", &
                          "plasma_3d end plasma", &
                          0, &
                          SIZE(plasma_3d), &
                          PACK(plasma_3d,.TRUE.) )   
#endif
!SMS$IGNORE END

       IF( MOD(utime_driver,ip_freq_output)==0)THEN
          WRITE( iterChar, '(I8.8)' )utime_driver
          CALL io_plasma_bin ( 1, utime_driver, 'iter_'//iterChar )
       ENDIF
       CALL output ( utime_driver )


!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-9")

!SMS$IGNORE BEGIN
#ifdef TESTING
       CALL mdi % Write_ModelDataInstances( "ipe" )
       CALL mdi % CalculateStorageCost(  )
#endif
!SMS$IGNORE END
     END DO


    ! Deallocate arrays
     CALL allocate_arrays ( 1 )

     ! close all open files
     CALL close_files ( )


!SMS$IGNORE BEGIN
#ifdef TESTING
     CALL mdi % Trash( )
#endif
!SMS$IGNORE END
     CALL stop


END PROGRAM  test_plasma
