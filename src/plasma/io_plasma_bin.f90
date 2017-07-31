!CAUTION!!!!!!!the plasma i-o orders/file names have been changed !!!

!CAUTION!!!!!!!the plasma i-o orders/file names have been changed !!!
!date:Thu Sep 29 18:58:05 GMT 2011
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
SUBROUTINE io_plasma_bin ( switch, utime, timestamp_for_IPE )
USE module_precision
USE module_IO,ONLY: LUN_PLASMA1,LUN_PLASMA2,lun_min1,lun_min2,lun_ut,record_number_plasma,lun_max1

USE module_FIELD_LINE_GRID_MKS,ONLY: JMIN_IN,JMAX_IS,plasma_3d,JMIN_ING,JMAX_ISG,VEXBup &
&, Un_ms1,tn_k,on_m3,n2n_m3,o2n_m3
USE module_IPE_dimension,ONLY: NMP,NLP,NPTS2D,ISPEC,ISPEV,IPDIM,ISPET,ISTOT
USE module_input_parameters,ONLY:sw_debug,record_number_plasma_start,mype &
&,sw_record_number,stop_time,start_time,duration,mpstop, sw_output_wind, sw_use_wam_fields_for_restart
USE module_physical_constants,ONLY:zero
! ghgm - now need the open_file module....
USE module_open_file,ONLY: open_file

!USE module_myIPE_Init, only: model_start_time, model_stop_time

IMPLICIT NONE

INTEGER (KIND=int_prec ),INTENT(IN) :: switch !2:read; 1:write
INTEGER (KIND=int_prec ),INTENT(IN) :: utime !universal time [sec]
character(len=13), INTENT(IN) :: timestamp_for_IPE
INTEGER (KIND=int_prec ),PARAMETER  :: n_max=10000
INTEGER (KIND=int_prec )            :: stat_alloc
INTEGER (KIND=int_prec )            :: jth,mp,lp,npts
INTEGER (KIND=int_prec )            :: lun,in,is
INTEGER (KIND=int_prec )            :: n_read,n_read_min, utime_dum,record_number_plasma_dum
INTEGER (KIND=int_prec )            :: n_count
INTEGER (KIND=int_prec )            :: ipts !dbg20120501
REAL    (KIND=real_prec)            :: dumm(NPTS2D,NMP)
CHARACTER(len=80) :: restart_directory
CALL getenv("RESDIR", restart_directory)
print *, 'GHGM IO_PLASMA restart_directory ',TRIM(restart_directory)

IF ( switch<1.or.switch>2 ) THEN
  print *,'sub-io_plasma:!STOP! INVALID switch',switch
  STOP
END IF

if(sw_debug)  print *,'sub-io_plasma_bin: switch=',switch,' utime[sec]' ,utime 

!output time dependent plasma parameters

! array initialization
dumm=zero

IF ( switch==1 ) THEN !1:Output the 16 plasma* files


  record_number_plasma = record_number_plasma+1
!SMS$SERIAL(<plasma_3d,IN>:default=ignore) BEGIN

! ghgm - open the new plasma file unit......
call open_file('ipe_grid_plasma_params.'//timestamp_for_IPE,5999,'unformatted','unknown')

  j_loop1: DO jth=1,ISTOT !=(ISPEC+3+ISPEV)
    mp_loop1:do mp=1,mpstop
      lp_loop1:do lp=1,nlp
        IN=JMIN_IN(lp)
        IS=JMAX_IS(lp)
        npts = IS-IN+1 
        dumm(JMIN_ING(lp):JMAX_ISG(lp),mp) = plasma_3d(IN:IS,lp,mp,jth)
      end do lp_loop1!lp
    end do mp_loop1!mp

    LUN = LUN_PLASMA1(jth-1+lun_min1)
    if(sw_debug) print *,'jth=',jth,' LUN=',LUN
    WRITE (UNIT=lun) (dumm(:,mp),mp=1,mpstop)
    if(sw_debug) print *,'!dbg! output dummy finished'

! ghgm - also write the 16 plasma files to a single unit.....
write (unit=5999) (dumm(:,mp),mp=1,mpstop)

  END DO j_loop1!jth
! ghgm - and then close the file
close (unit=5999)

!SMS$SERIAL END
  LUN = LUN_PLASMA1(lun_max1)
!SMS$SERIAL(<VEXBup,IN>:default=ignore) BEGIN
  WRITE (UNIT=LUN) (VEXBup(:,mp),mp=1,mpstop)
  WRITE (UNIT=lun_ut,FMT=*) record_number_plasma, utime
!SMS$SERIAL END

!nm20141001: moved from neutral
      IF ( sw_output_wind ) THEN

!SMS$SERIAL(<tn_k,Un_ms1,on_m3,n2n_m3,o2n_m3,IN>:default=ignore) BEGIN
call open_file('ipe_grid_neutral_params.'//timestamp_for_IPE,5998,'unformatted','unknown')
    write (unit=5998) tn_k
    write (unit=5998) un_ms1
    write (unit=5998) on_m3 
    write (unit=5998) n2n_m3
    write (unit=5998) o2n_m3
close (unit=5998)
!SMS$SERIAL END

      END IF !( sw_output_wind ) THEN


  if(sw_debug) then
    print *,'LUN=',lun_ut,'!dbg! output UT  finished: utime=',utime,record_number_plasma
  endif

ELSE IF ( switch==2 ) THEN !2:RESTART: 

!SMS$SERIAL BEGIN
! ghgm - read in saved plasma_3d data.....
!restart_directory='/scratch3/NCEPDEV/swpc/noscrub/George.Millward/restart_directory/'
print *,' ghgm reading in plasma_3d'
call open_file(trim(restart_directory) // 'ipe_grid_plasma_params',5997,'unformatted','old')
j_loop3: DO jth=1,(ISPEC+3)
    read (unit=5997) dumm
    mp_loop3:do mp=1,NMP
      lp_loop3:do lp=1,NLP
        IN   = JMIN_IN(lp)
        IS   = JMAX_IS(lp)
        npts = IS-IN+1 
        plasma_3d(IN:IS,lp,mp,jth) = dumm(JMIN_ING(lp):JMAX_ISG(lp),mp) 
      end do lp_loop3
    end do mp_loop3
END DO j_loop3
close(5997)
print *,' ghgm done reading in plasma_3d'
!SMS$SERIAL END

!SMS$SERIAL BEGIN
! ghgm - read in saved WAM neutral parameters.....
!restart_directory='/scratch3/NCEPDEV/swpc/noscrub/George.Millward/restart_directory/'
print *,' ghgm reading in wam neutrals'
call open_file(trim(restart_directory) // 'ipe_grid_neutral_params',5996,'unformatted','old')
    read (unit=5996) tn_k
    read (unit=5996) un_ms1
    read (unit=5996) on_m3 
    read (unit=5996) n2n_m3
    read (unit=5996) o2n_m3
close(5996)
print *,' ghgm done reading in wam neutrals'
!SMS$SERIAL END


END IF !( switch==1 ) THEN


print *,'END sub-io_pl: sw=',switch,' uts=' ,utime 
END SUBROUTINE io_plasma_bin
