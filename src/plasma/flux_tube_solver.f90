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
  SUBROUTINE flux_tube_solver ( utime )
    USE module_precision
    USE module_IPE_dimension,ONLY: ISPEC,ISPEV,IPDIM, NLP
    USE module_FIELD_LINE_GRID_MKS,ONLY: JMIN_IN,JMAX_IS,plasma_grid_3d,&
                                         plasma_grid_Z,plasma_grid_GL, &
                                         Pvalue,ISL,IBM,IGR,IQ,IGCOLAT,&
                                         IGLON,plasma_3d,ON_m3,HN_m3, &
                                         N2N_m3,O2N_m3,HE_m3,N4S_m3,TN_k,&
                                         TINF_k,un_ms1
    USE module_input_parameters,ONLY: mpstop,time_step,F107D,F107AV,DTMIN_flip, &
                                      sw_INNO,FPAS_flip,HPEQ_flip, &
                                      HEPRAT_flip,COLFAC_flip,sw_IHEPLS, &
                                      sw_INPLS,sw_debug,iout, start_time,&
                                      sw_wind_flip, sw_depleted_flip,&
                                      start_time_depleted,&
                                      sw_output_fort167,mpFort167,&
                                      lpFort167, sw_neutral_heating_flip, &
                                      ip_freq_output, parallelBuild,mype
    USE module_physical_constants,ONLY: pi
    USE module_IO,ONLY: PRUNIT,LUN_FLIP1,LUN_FLIP2,LUN_FLIP3,LUN_FLIP4
    USE module_unit_conversion,ONLY: M_TO_KM
    USE module_heating_rate,ONLY: get_neutral_heating_rate

    IMPLICIT NONE
!------------------------
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
! Local
      REAL (KIND=real_prec)   :: ltime !local time [hour]
      INTEGER (KIND=int_prec) :: mp    !longitude
      INTEGER (KIND=int_prec) :: lp    !latitude

!--- for CTIPINT
      INTEGER :: IN,IS !.. lcv + spatial grid indices
      INTEGER JMINX,JMAXX !.. lcv + spatial grid indices
      INTEGER CTIPDIM         !.. CTIPe array dimension, must equal to FLDIM
      DOUBLE PRECISION ::  PCO
      INTEGER ::  INNO            !.. switch to turn on FLIP NO calculation if <0
      DOUBLE PRECISION, DIMENSION(IPDIM) ::  ZX,GLX,SLX,BMX,GRX,OX,HX,N2X,O2X,HEX,N4SX,NNOX &
     &, UNX &     ! assumed a component parallel to a field line [m s-1]
     &, TNX &
       !.. TINFX has to be an array for grazing incidence column densities
     &, TINFX & !.. Exospheric temperature [k]
     &, SZA_dum
      DOUBLE PRECISION :: DTMIN !.. Minimum time step allowed (>=10 secs?)
      DOUBLE PRECISION DT

      REAL ::  F107D_dum,F107A_dum         !.. daily and 81 day average F10.7      !.. 
      DOUBLE PRECISION ::  FPAS,HEPRAT,COLFACX
      DOUBLE PRECISION HPEQ

      INTEGER ::  IHEPLS,INPLS  !.. switches He+ and N+ diffusive solutions on 
      DOUBLE PRECISION &
      !.. EHTX(3,J) = e heating rate, EHTX(1,J) = ion heating rate, EHTX(2,J) unused
     &  EHTX(3,IPDIM) &
      !.. TE_TI(3,J) = Te, TE_TIX(2,J) = Ti = TE_TIX(2,J)
     & ,TE_TIX(3,IPDIM) &
     & ,XIONNX(ISPEC,IPDIM),XIONVX(ISPEC,IPDIM) &
     & ,NHEAT(IPDIM) &  !.. Neutral heating rate [eV/cm^3/s]
     & ,hrate_cgs(22,IPDIM)   !.. heating rates [eV/cm^3/s]

      INTEGER EFLAG(11,11)    !.. error flags, check =0 on return from FLIP
      INTEGER :: PRUNIT_dum !.. Unit number to print results
      INTEGER (KIND=int_prec) :: midpoint
      INTEGER (KIND=int_prec) :: ipts,i,ret
      INTEGER (KIND=int_prec),PARAMETER :: ip_freq_output_fort=900      
      INTEGER (KIND=int_prec) :: jth !dbg20120501
!----------------------------------

        DO mp = 1,mpstop
          DO lp = 1,NLP

 
            IN = JMIN_IN(lp)
            IS = JMAX_IS(lp)

            ! make sure that JMINX equals to 1
            JMINX   = 1
            JMAXX   = IS - IN + 1
            CTIPDIM = IS - IN + 1   

            ZX(1:CTIPDIM)  = plasma_grid_Z(IN:IS,lp) * M_TO_KM !convert from m to km 
            PCO = Pvalue(lp)  !Pvalue is a single value
            SLX(1:CTIPDIM) = plasma_grid_3d(IN:IS,lp,mp,ISL)
            GLX(1:CTIPDIM) = pi/2. - plasma_grid_GL(IN:IS,lp)  ! magnetic latitude [radians]
            BMX(1:CTIPDIM) = plasma_grid_3d(IN:IS,lp,mp,IBM)   !Tesla
            GRX(1:CTIPDIM) = plasma_grid_3d(IN:IS,lp,mp,IGR)

            OX(1:CTIPDIM) = ON_m3(IN:IS,lp,mp) !(m-3)
            HX(1:CTIPDIM) = HN_m3(IN:IS,lp,mp)
            N2X(1:CTIPDIM) = N2N_m3(IN:IS,lp,mp)
            O2X(1:CTIPDIM) = O2N_m3(IN:IS,lp,mp)
            HEX(1:CTIPDIM) = HE_m3(IN:IS,lp,mp)
            N4SX(1:CTIPDIM) = N4S_m3(IN:IS,lp,mp)
            
            INNO = sw_INNO

            TNX(1:CTIPDIM) = TN_k(IN:IS,lp,mp)
            TINFX(1:CTIPDIM) = TINF_k(IN:IS,lp,mp)

            ! FLIP assumes positive SOUTHWARD along a field line
            IF ( sw_wind_flip == 1 ) THEN
              UNX(1:CTIPDIM)  = (-1.) * Un_ms1(IN:IS,lp,mp,3) 
            ELSE IF ( sw_wind_flip == 0 ) THEN
              UNX(1:CTIPDIM)  = 0.0
            END IF

            DT        = REAL(time_step)
            DTMIN     = DTMIN_flip
            F107D_dum = F107D
            F107A_dum = F107AV

            CALL Get_SZA ( utime,mp,lp, SZA_dum )

            FPAS      = FPAS_flip

            IF ( utime == start_time ) THEN
              HPEQ      = HPEQ_flip
            ELSE
              HPEQ      = 0.0

              IF ( sw_depleted_flip==1 .AND. utime == start_time_depleted ) THEN
                HPEQ      = - 0.1
              ENDIF
            ENDIF

            HEPRAT    = HEPRAT_flip
            COLFACX   = COLFAC_flip

            IHEPLS    = sw_IHEPLS
            INPLS     = sw_INPLS
  
            EFLAG(:,:)=0


!SMS$IGNORE BEGIN
#ifdef DEBUG
            print *,'sub-flux_tube_solver'
            print "('mp=',i6,' lp=',i6,' JMINX=',I6,' JMAXX=',I6)", mp,lp,jminx,jmaxx
            print "('CTIPDIM=',I6)", ctipdim
            print "('Z [km]     =',2F10.4)", ZX(jminx),ZX(jmaxx)
            print "('PCO        =',2F10.4)", PCO,Pvalue(lp)
            print "('SLX [m]    =',2E12.4)", SLX(jminx), SLX(jmaxx) 
            print "('GLX [deg]  =',2F10.4)",(GLX(jminx)*180./pi),(GLX(jmaxx)*180./pi)
            print "('BMX [Tesla]    =',2E12.4)", BMX(jminx), BMX(jmaxx)
            print "('GRX[m2 s-1]=',2E12.4)",GRX(jminx),GRX(jmaxx)
            !---neutral parameters
            print "('LOG10 OX [m-3]     =',2F10.4)",LOG10(OX(jminx)),LOG10(OX(jmaxx))
            print "('LOG10 HX [m-3]     =',2F10.4)",LOG10(HX(jminx)),LOG10(HX(jmaxx))
            print "('LOG10 N2X [m-3]    =',2F10.4)",LOG10(N2X(jminx)),LOG10(N2X(jmaxx))
            print "('LOG10 O2X [m-3]    =',2F10.4)",LOG10(O2X(jminx)),LOG10(O2X(jmaxx))
            print "('LOG10 HEX [m-3]    =',2F10.4)",LOG10(HEX(jminx)),LOG10(HEX(jmaxx))
            print "('LOG10 N4SX [m-3]   =',2F10.4)",LOG10(N4SX(jminx)),LOG10(N4SX(jmaxx))
            
            print "('INNO =',I6)",INNO
            IF ( INNO>=0 )  & !when CTIPe calculates NO
                 & print "('LOG10 NNOX [m-3]   =',2F10.4)",LOG10(NNOX(jminx)),LOG10(NNOX(jmaxx))
            print "('Tn [K]       =',2F10.4)",TNX(jminx),TNX(jmaxx)
            print "('TINF [K]     =',2F10.4)",TINFX(jminx),TINFX(jmaxx)
            print "('UNX [m s-1]  =',2F10.4)",UNX(jminx),UNX(jmaxx)
            
            print "('DT [sec]     =',F10.4)",DT
            print "('DTMIN [sec]  =',F10.4)",DTMIN
            print "('F107D_dum    =',F10.4)",F107D_dum
            print "('F107A_dum    =',F10.4)",F107A_dum
            print "('SZA [deg]    =',2F10.4)",SZA_dum(jminx)*180./pi,SZA_dum(jmaxx)*180./pi
            print "('FPAS         =',F10.4)",FPAS
            print "('HPEQ         =',F10.4)",HPEQ
            print "('HEPRAT       =',F10.4)",HEPRAT
            print "('COLFACX      =',F10.4)",COLFACX
            print "('IHEPLS       =',I6)",IHEPLS
            print "('INPLS        =',I6)",INPLS
#endif
!SMS$IGNORE END

            midpoint = IN + (IS-IN)/2
            IF ( lp>=1 .AND. lp<=6 )  midpoint = midpoint - 1
            ltime = REAL(utime)/3600.0 + (plasma_grid_3d(midpoint,lp,mp,IGLON)*180.0/pi)/15.0
            IF ( ltime > 24.0 )  ltime = MOD(ltime, 24.0)

            IF( sw_output_fort167.AND.mp==mpFort167.AND.lp==lpFort167 ) then

                  WRITE(UNIT=LUN_FLIP1,&
                        FMT="('mp=',i3,' lp=',i3,' U',i3,' North, UT=',2F10.3)")&
                         mp,lp,LUN_FLIP1,REAL(UTIME)/3600., ltime

                  WRITE(UNIT=LUN_FLIP2,&
                        FMT="('mp=',i3,' lp=',i3,' U',i3,' North, UT=',2F10.3)")&
                        mp,lp,LUN_FLIP2,REAL(UTIME)/3600., ltime

                  WRITE(UNIT=LUN_FLIP3,&
                        FMT="('mp=',i3,' lp=',i3,' U',i3,' South, UT=',2F10.3)")&
                        mp,lp,LUN_FLIP3,REAL(UTIME)/3600., ltime

                  WRITE(UNIT=LUN_FLIP4,&
                        FMT="('mp=',i3,' lp=',i3,' U',i3,' South, UT=',2F10.3)")&
                        mp,lp,LUN_FLIP4,REAL(UTIME)/3600., ltime
            
            END IF !( sw_output_fort167

!SMS$IGNORE BEGIN
#ifdef DEBUG
            print*,'sub-fl: UTs=',UTIME,' LThr=',ltime,' mp',mp,' lp',lp
            WRITE(UNIT=PRUNIT,FMT="('mp=',i6,' lp=',i6,' UT=',F10.2)") mp,lp,REAL(UTIME)/3600.,mype
#endif
!SMS$IGNORE END

            DO ipts=1,CTIPDIM
      
              DO jth=1,ISPEC
                XIONNX(jth,ipts) = plasma_3d(ipts-1+IN,lp,mp,jth)
              END DO !jth
      
              TE_TIX(3,ipts) = plasma_3d(ipts-1+IN,lp,mp,ISPEC+1)
              DO jth=1,2
                TE_TIX(jth,ipts) = plasma_3d(ipts-1+IN,lp,mp,jth+ISPEC+1)
              END DO !jth
      
              XIONVX(1,ipts) = plasma_3d(ipts-1+IN,lp,mp,1+ISPEC+3) 
              XIONVX(2,ipts) = plasma_3d(ipts-1+IN,lp,mp,2+ISPEC+3) 
              DO jth=3,ISPEC
                XIONVX(jth,ipts) = 0.0_real_prec 
              END DO !jth

              ! auroral electron heating-->EHTX(3)
              ! frictional heating for ions-->EHTX(1) 
              EHTX(1:3,ipts)=0.0_real_prec

              IF ( INNO<0 ) THEN
                NNOX(ipts)=0.0_real_prec
              ELSE
                print *,'CTIPe calculates NO'
              END IF
              NHEAT(ipts)          = 0.0_real_prec 
              hrate_cgs(1:22,ipts) = 0.0_real_prec 
     
           END DO 

           CALL CTIPINT( &
          &             JMINX, & !.. index of the first point on the field line
          &             JMAXX, & !.. index of the last point on the field line
          &           CTIPDIM, & !.. CTIPe array dimension, must equal to FLDIM
          &                ZX, & !.. array, altitude (km)
          &               PCO, & !.. p coordinate (L-shell)
          &               SLX, & !.. array, distance of point from northern hemisphere (meter)
          &               GLX, & !.. array, magnetic latitude (radians)
          &               BMX, & !.. array, magnetic field strength, (Tesla)
          &               GRX, & !.. array, gravity, m2 s-1
          &                OX, & !.. array, O density (m-3)
          &                HX, & !.. array, H density (m-3)
          &               N2X, & !.. array, N2 density (cm-3)
          &               O2X, & !.. array, O2 density (cm-3)
          &               HEX, & !.. array, He density (cm-3)
          &              N4SX, & !.. array, N(4S) density (cm-3)
          &              INNO, & !.. switch to turn on FLIP NO calculation if <0
          &              NNOX, & !.. array, NO density (cm-3)
          &               TNX, & !.. array, Neutral temperature (K)
          &             TINFX, & !.. array, Exospheric Neutral temperature (K)
          &               UNX, & !.. array, Neutral wind (m/s), positive northward,horizontal component in the magnetic meridian??? or component parallel to a field line???
          &                DT, & !.. CTIPe time step (secs)
          &             DTMIN, & !.. Minimum time step allowed (>=10 secs?)
          &         F107D_dum, & !.. Daily F10.7
          &         F107A_dum, & !.. 81 day average F10.7
          &           SZA_dum, & !.. Solar Zenith angle (radians)
          &              FPAS, & !.. Pitch angle scattering fraction
          &              HPEQ, & !.. Sets initial equatorial H+ density. See declaration below
          &            HEPRAT, & !.. Intial He+/H+ ratio (.01 to 1.0)
          &           COLFACX, & !.. O+ - O collision frequency Burnside factor (1.0 to 1.7)
          &            IHEPLS, & !.. switches He+ diffusive solution on if > 0
          &             INPLS, & !.. switches N+ diffusive solution on if > 0
          &              EHTX, & !.. IN/OUT 2D array, Electron & ion heating rate (eV cm-3 s-1)
          &            TE_TIX, & !.. IN/OUT: 2D array, Electron and ion temperatures (K) (see below)
          &     XIONNX,XIONVX, & !.. IN/OUT: 2D array, Storage for ion densities and velocities
          &             NHEAT, & !.. OUT: array, Neutral heating rate (eV/cm^3/s) 
          &             EFLAG, & !.. OUT: 2D array, Error Flags
          &                mp, &
          &                lp, &
          &         hrate_cgs  ) !.. heating rates [eV/cm^3/s] !nm20121020


          DO ipts=1,CTIPDIM
    
             DO jth=1,ISPEC
                plasma_3d(ipts+IN-1,lp,mp,jth) = XIONNX(jth,ipts)
             END DO 
    
             plasma_3d(ipts+IN-1,lp,mp,ISPEC+1) = TE_TIX(3,ipts)
    
             DO jth=1,2
                plasma_3d(ipts+IN-1,lp,mp,jth+ISPEC+1) = TE_TIX(jth,ipts)
             END DO !jth
    
             DO jth=1,ISPEV
               plasma_3d(ipts+IN-1,lp,mp,jth+ISPEC+3) = XIONVX(jth,ipts)
             END DO !jth
    
             IF ( sw_neutral_heating_flip==1 .AND. &
                &  MOD( (utime-start_time),ip_freq_output)==0) THEN
                if(parallelBuild) then
                   print*,'sw_neutral_heating_flip=1 does not work in parallel'
                   print*,'Stopping in Neut_heating'
                   stop
                endif
                   CALL get_neutral_heating_rate ( hrate_cgs , lp,mp )
    
             END IF
    
          END DO

          PRUNIT_dum = PRUNIT
          CALL WRITE_EFLAG(PRUNIT_dum, &  !.. Unit number to print results
         &                      EFLAG, &  !.. Error flag array
         &                         mp, &
         &                         lp,utime,ltime )

        ENDDO
     ENDDO    

  END SUBROUTINE flux_tube_solver
