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
SUBROUTINE io_plasma_bin ( switch, utime )
USE module_precision
USE module_IO,ONLY: LUN_PLASMA1,LUN_PLASMA2,lun_min1,lun_min2,lun_ut,lun_ut2,record_number_plasma,lun_max1 &
, lun_ipe_grid_neut_params_ut &
, lun_ipe_grid_neut_wind &
, lun_ipe_grid_neut_temp &
, lun_ipe_grid_neut_O_den &
, lun_ipe_grid_neut_N2_den &
, lun_ipe_grid_neut_O2_den &
, LUN_WAM_RESTART0,LUN_WAM_RESTART1,lun_wam_tn,LUN_WAM_RESTART3,LUN_WAM_RESTART4,LUN_WAM_RESTART5
USE module_FIELD_LINE_GRID_MKS,ONLY: JMIN_IN,JMAX_IS,plasma_3d,JMIN_ING,JMAX_ISG,VEXBup &
&, VEXBe,VEXBth &
&, Un_ms1,tn_k,on_m3,n2n_m3,o2n_m3
USE module_IPE_dimension,ONLY: NMP,NLP,NPTS2D,ISPEC,ISPEV,IPDIM,ISPET,ISTOT
USE module_input_parameters,ONLY:sw_debug,record_number_plasma_start,mype &
&,sw_record_number,stop_time,start_time,duration,mpstop, sw_output_wind, sw_use_wam_fields_for_restart
USE module_physical_constants,ONLY:zero
IMPLICIT NONE

INTEGER (KIND=int_prec ),INTENT(IN) :: switch !2:read; 1:write
INTEGER (KIND=int_prec ),INTENT(IN) :: utime !universal time [sec]
INTEGER (KIND=int_prec ),PARAMETER  :: n_max=10000
INTEGER (KIND=int_prec )            :: stat_alloc
INTEGER (KIND=int_prec )            :: jth,mp,lp,npts
INTEGER (KIND=int_prec )            :: lun,in,is
INTEGER (KIND=int_prec )            :: n_read,n_read_min, utime_dum,record_number_plasma_dum
INTEGER (KIND=int_prec )            :: n_count
INTEGER (KIND=int_prec )            :: ipts !dbg20120501
REAL    (KIND=real_prec)            :: dumm(NPTS2D,NMP)


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
  END DO j_loop1!jth
!SMS$SERIAL END
  LUN = LUN_PLASMA1(lun_max1-2)
!SMS$SERIAL(<VEXBup,IN>:default=ignore) BEGIN
  WRITE (UNIT=LUN) (VEXBup(:,mp),mp=1,mpstop)
  WRITE (UNIT=lun_ut,FMT=*) record_number_plasma, utime
!SMS$SERIAL END
!nm20140218: output zonal drift
  LUN = LUN_PLASMA1(lun_max1-1)
!SMS$SERIAL(<VEXBe,IN>:default=ignore) BEGIN
  WRITE (UNIT=LUN) (VEXBe(:,mp),mp=1,mpstop)
!SMS$SERIAL END
!nm20140701: output meridional drift
  LUN = LUN_PLASMA1(lun_max1)
!SMS$SERIAL(<VEXBth,IN>:default=ignore) BEGIN
  WRITE (UNIT=LUN) (VEXBth(:,mp),mp=1,mpstop)
!SMS$SERIAL END

!nm20141001: moved from neutral
      IF ( sw_output_wind ) THEN

!nm20160711 debug wam field: 
!SMS$SERIAL(<tn_k,IN>:default=ignore) BEGIN
    write (UNIT=lun_ipe_grid_neut_temp) tn_k
    write (UNIT=lun_ipe_grid_neut_params_ut,FMT=*) utime
!SMS$SERIAL END
!SMS$SERIAL(<Un_ms1,IN>:default=ignore) BEGIN
           write (UNIT=lun_ipe_grid_neut_wind) Un_ms1
!SMS$SERIAL END
!SMS$SERIAL(<on_m3,IN>:default=ignore) BEGIN
           write (UNIT=lun_ipe_grid_neut_O_den) on_m3
!SMS$SERIAL END
!SMS$SERIAL(<n2n_m3,IN>:default=ignore) BEGIN
           write (UNIT=lun_ipe_grid_neut_N2_den) n2n_m3
!SMS$SERIAL END
!SMS$SERIAL(<o2n_m3,IN>:default=ignore) BEGIN
           write (UNIT=lun_ipe_grid_neut_O2_den) o2n_m3
!SMS$SERIAL END
      END IF !( sw_output_wind ) THEN
!
  if(sw_debug) then
    print *,'LUN=',lun_ut,'!dbg! output UT  finished: utime=',utime,record_number_plasma
  endif

ELSE IF ( switch==2 ) THEN !2:RESTART: Read from the 16 plasma* files

  print *,'sub-io_pl: start_uts=',utime
!SMS$PARALLEL(dh, lp, mp) BEGIN
  ! array initialization
  DO mp=1,NMP
    DO lp=1,NLP
      DO ipts=JMIN_IN(lp),JMAX_IS(lp)
        DO jth=1,ISTOT
          plasma_3d(ipts,lp,mp,jth) = zero
        END DO!j
      END DO!ipts
    END DO!lp
  END DO!mp
  !nm20120509: automatically keep running global run
  !0:you need to specify the record_number_plasma_start by your self...
  IF ( sw_record_number==0 ) THEN
!SMS$SERIAL BEGIN
    rd_loop0: DO n_read=1,n_max !=10000
      READ (UNIT=lun_ut2,FMT=*) record_number_plasma_dum, utime_dum
      print *,' record_# ', record_number_plasma_dum,' uts=', utime_dum
      IF (n_read==1) THEN
        n_read_min=record_number_plasma_dum
        print *,'n_read=',n_read,'n_read_min=', n_read_min
      END IF
      IF (record_number_plasma_dum==record_number_plasma_start) THEN
        IF (utime_dum==utime) THEN
          print *,'n_read=',n_read,' confirmed that ipe output exist at rec#=',record_number_plasma_dum,' at UT=', utime_dum
          CLOSE(UNIT=lun_ut2)
          EXIT  rd_loop0
        ELSE !IF (utime_dum/=utime) THEN
          print *,'n_read',n_read,'!STOP! INVALID start_time!',utime, record_number_plasma_dum, utime_dum
          STOP
        END IF
      ELSE IF (n_read==n_max.OR.record_number_plasma_dum>record_number_plasma_start) THEN
        print *,'n_read',n_read,'!STOP! INVALID record number!',record_number_plasma_dum, utime_dum
        STOP
      END IF
    END DO rd_loop0!: DO n_read=1,n_max          
!SMS$SERIAL END
  ELSE IF ( sw_record_number==1 ) THEN
  !1: automatically keep running global run: the code sets up the record_number from the very last record...
!SMS$SERIAL(<record_number_plasma_start,start_time,n_read_min,OUT>:default=ignore) BEGIN
    n_count=0
    rd_loop1: DO n_read=1,n_max !=10000
      READ (UNIT=lun_ut2,FMT=*,END=19) record_number_plasma_dum, utime_dum
      n_count=n_count+1
      IF (n_read==1) THEN
        n_read_min=record_number_plasma_dum
        print *,'n_read=',n_read,'n_read_min=', n_read_min
      END IF
    END DO rd_loop1!: DO n_read=1,n_max !=10000
19  CONTINUE
    record_number_plasma_start = record_number_plasma_dum
    start_time = utime_dum
    print *,'new record_number_plasma_start=',record_number_plasma_start
    print *,'new start_time=',start_time
!SMS$SERIAL END
    stop_time = start_time + duration
    print *,'new stop_time=',stop_time,' duration=', duration
  END IF !( sw_record_number==0 ) THEN        
!SMS$SERIAL BEGIN
  j_loop2: DO jth=1,(ISPEC+3)  !t  +ISPEV)
    LUN = LUN_PLASMA2(jth-1+lun_min2)
    if(sw_debug) print *,'jth=',jth,' LUN2=',LUN
    rd_loop: DO n_read=n_read_min, record_number_plasma_start
      if(jth==ISPEC+3) print *,'read dumm n_read=',n_read
      READ (UNIT=lun ) dumm
      if (sw_debug.and.jth==ISPEC+3) print *,'!dbg! read dummy finished jth',jth
    END DO rd_loop !: DO n_read=1,n_read_max
    mp_loop2:do mp=1,NMP
      lp_loop2:do lp=1,NLP
        IN   = JMIN_IN(lp)
        IS   = JMAX_IS(lp)
        npts = IS-IN+1 
        plasma_3d(IN:IS,lp,mp,jth) = dumm(JMIN_ING(lp):JMAX_ISG(lp),mp) 


if(jth==1.and.lp==57.and.mp==71)then
print*,'subIo:',plasma_3d(30:40,lp,mp,jth)
endif

      end do lp_loop2!lp
    end do mp_loop2!mp

    print *,'closing lun',LUN
    CLOSE(LUN)
  END DO j_loop2!jth
!SMS$SERIAL END

! Read in the WAM startup fields.... 
IF ( sw_use_wam_fields_for_restart ) THEN
!SMS$SERIAL BEGIN
print *,'****** GEORGE ******* READING WAM TN 1'
  read (UNIT=lun_wam_tn) tn_k
print *,'****** GEORGE ******* READING WAM TN 2'
  read (UNIT=lun_wam_tn) tn_k
print *,'***** THIS INNIT ',tn_k(JMIN_IN(10)+2,10,10)
print *,'***** THIS TOOOOO INNIT '
  close(lun_wam_tn)
!  read (UNIT=LUN_WAM_RESTART0,FMT=*) utime
!  close(LUN_WAM_RESTART0)
print *,'****** GEORGE ******* READING WAM WIND'
  read (UNIT=LUN_WAM_RESTART1) Un_ms1
print *,'***** THIS INNIT2 ',Un_ms1(JMIN_IN(10)+2,10,10,1),Un_ms1(JMIN_IN(10)+2,10,10,2),Un_ms1(JMIN_IN(10)+2,10,10,3)
  close(LUN_WAM_RESTART1)
print *,'****** GEORGE ******* READING WAM ON'
  read (UNIT=LUN_WAM_RESTART3) on_m3
print *,'***** THIS INNIT3 ',on_m3(JMIN_IN(10)+2,10,10)
  close(LUN_WAM_RESTART3)
print *,'****** GEORGE ******* READING WAM N2N'
  read (UNIT=LUN_WAM_RESTART4) n2n_m3
print *,'***** THIS INNIT4 ',n2n_m3(JMIN_IN(10)+2,10,10)
  close(LUN_WAM_RESTART4)
print *,'****** GEORGE ******* READING WAM O2N'
  read (UNIT=LUN_WAM_RESTART5) o2n_m3
print *,'***** THIS INNIT5 ',o2n_m3(JMIN_IN(10)+2,10,10)
  close(LUN_WAM_RESTART5)
print *,'****** GEORGE ******* DONE READING WAM FIELDS'
!SMS$SERIAL END
END IF

END IF !( switch==1 ) THEN

!SMS$PARALLEL END


print *,'END sub-io_pl: sw=',switch,' uts=' ,utime 
END SUBROUTINE io_plasma_bin
