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
      USE module_input_parameters,ONLY: read_input_parameters,start_time,stop_time,time_step,HPEQ_flip,ip_freq_msis,sw_output_plasma_grid,sw_debug,sw_perp_transport,parallelBuild,mype,ip_freq_eldyn,SMScomm,lps,lpe,mps,mpe
      USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_3d
      USE module_init_plasma_grid,ONLY: init_plasma_grid
      USE module_NEUTRAL_MKS,ONLY: neutral 
      USE module_sub_PLASMA,ONLY: plasma
!SMS$IGNORE BEGIN
!nm20121003      USE module_ELDYN,ONLY: init_eldyn, eldyn
!nm20121003:
      USE module_init_ELDYN,ONLY: init_eldyn
      USE module_sub_ELDYN,ONLY: eldyn
!SMS$IGNORE END
      USE module_open_output_files,ONLY: open_output_files
      USE module_output,ONLY: output
      USE module_close_files,ONLY: close_files
      USE module_IPE_dimension,ONLY: NMP,NLP
      USE module_IO,ONLY: PRUNIT
      USE module_open_file,ONLY: open_file
      IMPLICIT NONE
!SMS$INSERT include "mpif.h"
      include "gptl.inc"

      INTEGER(KIND=int_prec)           :: utime !universal time [sec]
      INTEGER(KIND=int_prec),parameter :: luntmp=300
      INTEGER(KIND=int_prec)           :: istat,mp,ret,OLDcomm,key,status=0,color=0,mypeWorld
      INTEGER(KIND=int_prec)           :: resultlen ! length return from MPI routines
!SMS$INSERT character(len=mpi_max_processor_name) :: mynode ! node name

      call gptlprocess_namelist ('GPTLnamelist', 77, ret) 
      ret = gptlinitialize ()
      ret = gptlstart ('Total')

! set up input parameters
      ret = gptlstart ('read_input')
      CALL read_input_parameters ( ) !calls create_decomp
      ret = gptlstop  ('read_input')

!      key = mype
!      if(0<mype.and.mype<40)  then 
!        if(mod(mype,2) /= 0) then
!          key = mype+38
!        endif
!      endif
!      if(39<mype.and.mype<79)  then 
!        if(mod(mype,2) == 0) then
!          key = mype-38
!        endif
!      endif
!!SMS$ignore begin
!      print*,'mype,key=',mype,key
!!SMS$ignore end
!      key=mype
!SMS$INSERT parallelBuild=.true.
!!SMS$INSERT CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,SMScomm,status)
!!SMS$SET_COMMUNICATOR(SMScomm)
!      IF (status /= 0) THEN
!        print*,'Error in driver_ipe gettin new MPI commumicator',status
!        stop
!      endif
!!SMS$INSERT call NNT_ME    (mype  )
!if(mype == 40) then
!   call set_affinity (0)
!endif
!SMS$INSERT call mpi_get_processor_name (mynode, resultlen, status)
!SMS$ignore begin
!SMS$INSERT print"('mynode,mype,lps,lpe,mps,mpe',2x,a10,5i5)",mynode,mype,lps,lpe,mps,mpe
!SMS$ignore end
!SMS$INSERT call print_affinity (mype)

! open Input/Output files
      ret = gptlstart ('open_output_files')
!--- unit=8
!filename ='FLIP_ERROR_FLAG_'//TRIM(string_tmp)//'.log'
      CALL open_file ( 'FLIP_ERR', PRUNIT,'formatted','unknown') !open by all processors
!SMS$SERIAL BEGIN
      CALL open_output_files ( )
!SMS$SERIAL END
      ret = gptlstop  ('open_output_files')

! set up plasma grids by reading file
      ret = gptlstart ('init_plasma_grid')
      CALL init_plasma_grid ( )
      ret = gptlstop  ('init_plasma_grid')

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-1")
IF ( sw_output_plasma_grid ) THEN
  ret = gptlstart ('output_plasma_grid')
  print *, 'sub-init_p: output plasma_grid'
  CALL output_plasma_grid ( )
  ret = gptlstop  ('output_plasma_grid')
END IF
!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-2")

! initialise the flux tubes from previous runs
      IF ( HPEQ_flip==0.0 ) THEN
        print *,'before CALL io_plasma_bin finished! READ: start_time=', start_time,stop_time
        ret = gptlstart ('io_plasma_bin')
        CALL io_plasma_bin ( 2, start_time )
        ret = gptlstop  ('io_plasma_bin')
        print *,'after CALL io_plasma_bin finished! READ: start_time=', start_time,stop_time

      END IF
!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-3")
!20120215: CALL io_plasma_bin_readinit ( start_time )
!20120215: print *,'CALL io_plasma_bin_readinit finished!'

! initialization of electrodynamic module:
! read in E-field
      ret = gptlstart ('init_eldyn')
      IF ( sw_perp_transport>=1 ) THEN
        CALL init_eldyn ( )
      ENDIF
      ret = gptlstop  ('init_eldyn')
!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-4")

      ret = gptlstart ('time_loop')
      time_loop: DO utime = start_time, stop_time, time_step
      print"('UTime=',2i7,f11.4)",utime,(MOD(utime,86400)),(MOD(utime,86400)/3600.)
!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-5")
! updates auroral precipitation

! interplate from plasma to neutral grid: Nei,Tei,Vi,NHEAT, auroral heating?

!nm20110907:moved here because empirical Efield is needed for both neutral &plasma
      ret = gptlstart ('eldyn')
      IF ( sw_perp_transport>=1.AND. MOD( (utime-start_time),ip_freq_eldyn)==0 ) THEN
        CALL eldyn ( utime )
      ENDIF
      ret = gptlstop  ('eldyn')
!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-6")

! update neutral 3D structure: use MSIS/HWM to get the values in the flux tube grid
        IF ( MOD( (utime-start_time),ip_freq_msis)==0 ) THEN 
 IF ( sw_debug )  print *,'call MSIS',utime,start_time,ip_freq_msis,(utime-start_time),MOD( (utime-start_time),ip_freq_msis)
          ret = gptlstart ('neutral')
          CALL neutral ( utime )
          ret = gptlstop  ('neutral')
        END IF
!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-7")

! interpolate from neutral to plasma grid:Tn,Un,[O,N2,O2],EHT(1,k), auroral heating?

! update plasma
        ret = gptlstart ('plasma')
        CALL plasma ( utime )
        ret = gptlstop  ('plasma')
!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-8")

! update self-consistent electrodynamics
!t        CALL eldyn ( utime )

! output UTIME to a FLIP_ERR
!        ret = gptlstart ('output_barrier')
!!sms$insert      call ppp_barrier(istat)
!        ret = gptlstop  ('output_barrier')
        ret = gptlstart ('output')
        CALL output(utime)
        ret = gptlstop  ('output')

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-9")

      END DO  time_loop !: DO utime = start_time, stop_time, time_step
      ret = gptlstop  ('time_loop')

! DEallocate arrays
      ret = gptlstart ('allocate_arrays1')
      CALL allocate_arrays ( 1 )
      ret = gptlstop  ('allocate_arrays1')

! close all open files
      ret = gptlstart ('close_files')
      CALL close_files ( )
      ret = gptlstop  ('close_files')


!dbg20120509: no longer need 
!20120207 I have to output the sw_perp_transport to a file for the next run...
!IF ( sw_rw_sw_perp_trans ) THEN
!open(unit=luntmp, file='fort.300',status='unknown',form='formatted',iostat=istat)
!DO mp=1,NMP
!write(unit=luntmp, fmt='(2i3)') mp,sw_perp_transport(mp)
!print *,'mp=',mp,' sw_p',sw_perp_transport(mp)
!END DO
!close(unit=luntmp)
!END IF !( sw_tmp_sw_perp_trans ) THEN

      ret = gptlstop  ('Total')
      call stop

END PROGRAM  test_plasma
