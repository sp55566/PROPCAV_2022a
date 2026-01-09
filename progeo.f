      SUBROUTINE PROGEO
************************************************************************
*     PROGEO: PROpeller GEOmetry                                       *
*      --- Generate a propeller geometry                               *
*                                                                      *
*     Date     Comment or Revision                                     *
*     -------- -------------------                                     *
*     CM010798 Added if statement to call gblade2 when icon = 5        *
*     CM010798 Added if statement to call gwake2 when icon = 5         *
*     CM010898 Commented out gwake2 generation in above comment        *
*     CM051898 Square root of negative number was being taken in Do 30 *
*              loop.  Put an if statement to prevent in the future     *
*     JY112198 Added a new camber and a new thickness option:          *
*                     ICAM=99  User input camber distribution          *
*                     ITHK=99  User input thickness distribution       *
*      JY120498 Added a new parameter in *.adm                         *
*                     ISC=0   blade section have zero T.E. thickness   *
*                     ISC=1   blade section have finite T.E. thickness *
*                                                                      *
************************************************************************

      use m_INPGEO2
      use su_inp
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'INPGEO.INC'
      INCLUDE 'PUFCAVC.INC'
      COMMON /HUB5TMP/ RHULT1
      COMMON /MODBLD/ XMOD(80,50),YMOD(80,50),ZMOD(80,50)
      DIMENSION XID1(201),ETAD1(201),DCUBICX1(800)
      LOGICAL IDSTP
      REAL*4 AREATUN1,AREATUN2,AREAHUB1,AREAHUB2
C-----------------------------------------------------------------------
C     Parameters for generating the geometry
C-----------------------------------------------------------------------

      NH=NC/2
      NHP=NH+1
      NHM=NH-1
      NCP=NC+1
      MRP=MR+1
      NPANB=NC*MR
      DELK=TWOPI/NBLADE

C-----------------------------------------------------------------------
C     Create the blade geometry
C-----------------------------------------------------------------------
C
C.....Radial spacing: IRSPAC = 0:cosine; 1:constant; 2:half cosine

      IF(IRSPAC.EQ.0) THEN
         DTH=PI/MR
         DO 5 M=1,MRP
            RZ(M)=RHUB+0.5*(1.0-RHUB)*(1.0-COS(DTH*(M-1)))
5        CONTINUE
      ELSE IF(IRSPAC.EQ.1) THEN
         DELR=(1.0-RHUB)/MR
         DO 10 M=1,MRP
            RZ(M)=RHUB+DELR*(M-1)
10       CONTINUE
      ELSE IF(IRSPAC.EQ.2) THEN
         DTH=HALFPI/MR
         DO 20 M=1,MRP
            RZ(M)=RHUB+(1.0-RHUB)*SIN(DTH*(M-1))
20       CONTINUE
      END IF

C-----------------------------------------------------------------------
C     Change RZ spacing according to the panel input (IDRT = 1)
C                                         Yiran Su 2017 08 18
C-----------------------------------------------------------------------
      IF (IDRT.EQ.1) THEN
         OPEN(198,FILE='panel_inp.dat',STATUS='OLD')
         READ(198,*) NN1,MM1
         DO M = 1,MM1
            DO N = 1,NN1
               READ(198,*) XXX,YYY,ZZZ
               IF (N.EQ.1) THEN
                  RZ(M) = SQRT(YYY**2 + ZZZ**2)
               END IF
            END DO
         END DO
         CLOSE(198)
      END IF

C-----------------------------------------------------------------------
C     Adjustment made here to make sure radical term is positive
C                                                              CM051898
C-----------------------------------------------------------------------
      DO 30 M=1,MRP
         TERM1  =1.0 - RZ(M)
         IF(TERM1.LT.0.0)THEN
            TERM1 = ABS(TERM1)
         ENDIF
        RZSQ(M)=1.0-SQRT(TERM1)
30    CONTINUE

C.....Interpolate geometrical parameters in the radial direction

      CALL EVALDK(NX,MRP,XR,RZ,PITCH,PICUB)
      CALL EVALDK(NX,MRP,XR,RZ,RAKE,RKCUB)
      CALL EVALDK(NX,MRP,XR,RZ,SKEW,SKCUB)
      CALL EVALDK(NX,MRP,XRSQ,RZSQ,CHORD,CHCUB)
      CALL EVALDK(NX,MRP,XR,RZ,CAMBR,CICUB)
      CALL EVALDK(NX,MRP,XR,RZ,THICK,TICUB)

C.....Tip chord for equal area of parabolic segment. Zero tip thickness

C-----------------------------------------------------------------------
C    The following line will be commented out and replaced with the
C    equivalent line from PSF10.  I'm trying to get the circulations
C    to match                                                  CM021598
C    Also, the following IF statment will be uncommented out.  It
C    was commented out in the original version that I received.
C    Actually, we need to be a little more careful for the hydrofoil
C    case.  If we are running a hydrofoil, then use neal's formulation,
C    otherwise use the PSF10 equivalent formulation.           CM021998
C-----------------------------------------------------------------------

c      IF((ICON.EQ.5).OR.(ICON.EQ.6))THEN
c         CHORD(MRP)=AMAX1(CHORD(MRP),CHORD(MR)/3.0)
c      ELSE
c         CHORD(MRP)=AMAX1(CHORD(MRP),7.0*CHORD(MRP-1)/15.0)
c         THICK(MRP)=ZERO
c      ENDIF

      IF(ITUN .EQ. 0) THEN
         IF(IDUCT .EQ. 0 .OR. DUCTGAP .NE. 0.0) THEN
            IF(ICON.NE.5.AND.ICON.NE.6) THICK(MRP)=ZERO
         ENDIF
      ENDIF

C-----------------------------------------------------------------------
C     Spline for user input camber.                             JY112798
C-----------------------------------------------------------------------
C.....XCS = F/C, XTS = T/C;  XCS1 = F/C, XTS1 = T/D
      IF(ICAM.EQ.99) THEN
         DO 32 JJ=1,15
            CALL UGLYDK(NX,1,1,XR,XCS(1,JJ),0.0,0.0,CUBTMP)
            CALL EVALDK(NX,MRP,XR,RZ,XCS0,CUBTMP)
            DO 34 MM=1,MRP
               XCS1(JJ,MM)=XCS0(MM)
 34         CONTINUE
 32      CONTINUE
      END IF

C-----------------------------------------------------------------------
C     Spline for user input thickness.                         JY112798
C-----------------------------------------------------------------------
C.....XCS = F/C, XTS = T/C;  XCS1 = F/C, XTS1 = T/D
      IF(ITHK.EQ.99) THEN
         DO 36 JJ=1,16
            CALL UGLYDK(NX,1,1,XR,XTS(1,JJ),0.0,0.0,CUBTMP)
            CALL EVALDK(NX,MRP,XR,RZ,XTS0,CUBTMP)
            DO 38 MM=1,MRP
               XTS1(JJ,MM)=XTS0(MM)*CHORD(MM)
 38         CONTINUE
            IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8)
     *           XTS1(JJ,MRP)=ZERO
 36      CONTINUE
      END IF


C.....Generate the blade geometry

C-----------------------------------------------------------------------
C     Either call the gblade subroutine or gblade2 depending on whether
C     or not we are running the hydrofoil case                  CM010798
C-----------------------------------------------------------------------
C....Same for ICON=8 as well (JY110100)

      IF((ICON.EQ.5).OR.(ICON.EQ.6).OR.(ICON.EQ.8))THEN
         CALL GBLADE2
C-----------------------------------------------------------------------
C     Call GBLADE3 for special fan-like SP propeller.           JY090100
C-----------------------------------------------------------------------
      ELSE IF(ICON.EQ.7) THEN
         CALL GBLADE3
      ELSE
         CALL GBLADE
      ENDIF

      NPAND = 0

C      IF(IDUCT .NE. 0) CALL DUCTINPUT

C-----------------------------------------------------------------------
C     Prepare parameters for hub and wake geometries
C-----------------------------------------------------------------------
      XHBLE=-XB(NC/2+1,1)
      XHBTE=XB(1,1)
      XHBNU=XHBU+XHBLE
      XHBND=XHBD+XHBTE
      XHBFD=XHBT+XHBND

C -- Begin Tip HSLEE(10/11/99)
      CALL CHECK_NWZ

      if(ian .eq. 2) then
         call gwaketip
      Else
         IF(IDUCT .NE. 0) THEN
            CALL GWAKE_DUCT
         ELSE
C            IF(ICON .EQ. 5) THEN
C               CALL GWAKE_FOIL
C            ELSE
               IF(ITUN .EQ. 0) THEN
C                  WRITE(*,*) 'FWA based on the blade-pitch wake alignment'
                  CALL GWAKE
                  IF(IHUB .EQ. 6) CALL PODGWAKE
               ELSE
                  CALL GWAKE_TUN
               ENDIF
C            ENDIF
         ENDIF
      Endif
C -- End Tip (10/11/99)

C Yiran Su read wake geometry for IAN=5 10/05/2016
      INQUIRE(FILE='wakinp.dat',EXIST=IDSTP)
      IF ((IDSTP.EQ.(.TRUE.)).AND.(IAN.EQ.1))  THEN
        open(1104,file='wakinp.dat',status='OLD')
        do m = 1 , mr+1
          read(1104,*) nsw(m)
          nwpanel = nsw(m)-1
          do n = 1 , nsw(m)
            read(1104,'(3F12.6)') xw(n,m),yw(n,m),zw(n,m)
          enddo
C/s S.N.KIM | Even when wake is copied from outsource, fisrt spanwise stip needs to
C           | be attached to blade T.E. forcibly. Otherwise there could be very 
C           | infinitesimal gap, breaking down the Kutta condition critically.
          xw(1,m) = xb(1,m)
          yw(1,m) = yb(1,m)
          zw(1,m) = zb(1,m)
C/e S.N.KIM | Nov. 2018.
        enddo
        close(1104)

        do k = 1 , nblade
          do m = 1, mr+1
            do n = 1, nsw(m)
              xww(n,m,k) = xw(n,m)
              yww(n,m,k) = sqrt(yw(n,m)**2+zw(n,m)**2)*cos(atan2(zw(n,m),yw(n,m))
     *                    -real(k-1)*twopi/nblade)
              zww(n,m,k) = sqrt(yw(n,m)**2+zw(n,m)**2)*sin(atan2(zw(n,m),yw(n,m))
     *                    -real(k-1)*twopi/nblade)
            enddo
            if (k.eq.1) then
              xww(1,m,k) = xb(1,m)
              yww(1,m,k) = yb(1,m)
              zww(1,m,k) = zb(1,m)
            endif 
          enddo
        enddo

      ENDIF
C Yiran end 10/05/2016

      NPANW=0
      DO 40 M=1,MR
         NPANW=NPANW+MIN(NSW(M),NSW(M+1))
         IF(NSW(M).EQ.NSW(M+1)) THEN
            NPANW=NPANW-1
         END IF
40    CONTINUE

!---------------- Let's set a gap between duct and blade and wake -------------
      IF(IDUCT .NE. 0) THEN
         CALL DUCTGEO
         IF(IREPANEL.EQ.1) THEN
           CALL DUCTWAKGEO_FIXEDPANEL ! ductwake length (ndwk=nwakep=ncdw-20) is determined here and depends on NWAKEP when repaneling is included. S.N.KIM | Aug. 2018.
         ELSE
           CALL DUCTWAKGEO ! ductwake length (ndwk=ncdw-20) in FWA is determined here and depends on DWAKEL when NO repaneling is included. S.N.KIM | Aug. 2018.
         ENDIF
         CALL GWAKES_DUCT
         CALL DUCTPLT

         NPAND = NDUCT*MDUCT
         IF(IDOPT .EQ. 0) THEN

            NN1 = NDDAT/2 + 1

            DO N = 1, NN1
               N1 = NN1 - N + 1
               XID1(N) = XID(N1)
               ETAD1(N) = ETAD(N1)
            ENDDO

            DO N = 1 , NCP
               CALL UGLYDK(NN1,1,1,XID1,ETAD1,ZERO,ZERO,DCUBICX1)
               CALL EVALDKs(NN1,1,XID1,XB(N,MRP),RR1,DCUBICX1)
               THETA = DANGLE(ZB(N,MRP),YB(N,MRP))

               RR2 = RR1 - DUCTGAP

               YB(N,MRP) = RR2 * COS(THETA)
               ZB(N,MRP) = RR2 * SIN(THETA)
            ENDDO

            XW(1,MRP) = XB(1,MRP)
            YW(1,MRP) = YB(1,MRP)
            ZW(1,MRP) = ZB(1,MRP)

            RREND = SQRT(YB(1,MRP)**2+ZB(1,MRP)**2)

            RTMP = RREND
            DO N = 2, NSW(MRP)
               IF(XW(N,MRP) .LT. XD(1,1)) THEN
                  THETA =  DANGLE(ZW(N,MRP),YW(N,MRP))
                  CALL EVALDKs(NN1,1,XID1,XW(N,MRP),RR2,DCUBICX1)
                  IF(RREND .LE. RR2) THEN
                     YW(N,MRP) = RREND * COS(THETA)
                     ZW(N,MRP) = RREND * SIN(THETA)
                     RTMP = RREND
                  ELSE
                     YW(N,MRP) = RR2 * COS(THETA)
                     ZW(N,MRP) = RR2 * SIN(THETA)
                     RTMP = RR2
                  ENDIF
               ELSE
                  THETA =  DANGLE(ZW(N,MRP),YW(N,MRP))
                  YW(N,MRP) = RTMP * COS(THETA)
                  ZW(N,MRP) = RTMP * SIN(THETA)
               ENDIF
            ENDDO

         ENDIF
      ENDIF
!---------------------End of putting a gap between blade and duct and wake----------

!s--YE TIAN 04/13/2012 add following lines
      if (ian.eq.6) then
c        write(*,*) 'nwpanel=',nwpanel
c        write(*,*) 'nwk=',nwk  ! nwk = nwpz = nwz+1
        nwpanel = 60
c        !         nwpanel = nwk-1
c        write(*,*) 'nwpanel=',nwpanel
      endif
!e--YE TIAN 04/13/2012 add following lines

        CALL gwakes(NSUB, NWSUB1)

C-----------------------------------------------------------------------
C     Create the hub geometry
C
C     Added new option to IHUB (JY052401):
C        IHUB=0: hub off
C            =1: hub on (open ended far upstream, closed far downstream)
C            =2: hub on (closed far upstream, closed far downstream)
C            =3: hub on (open ended far upstream, open far downstream)
C            =4: hub on (closed far upstream, open far downstream)
C            =5: hub on (open ended far upstream, open far downstream)
C-----------------------------------------------------------------------

      IF(IHUB.EQ.0)THEN
         NHBX=0
      ELSE

         if(ian .eq. 2 .and. iscav .ne. 0) go to 1111

         IF(IHUB.EQ.1) THEN
            CALL GHUB1
         ELSE IF(IHUB.EQ.2) THEN
            CALL GHUB2
         ELSE IF(IHUB.EQ.3) THEN
            CALL GHUB3
         ELSE IF(IHUB.EQ.4) THEN
            CALL GHUB4
         ELSE IF(IHUB.EQ.5) THEN
c            IF(IAN.NE.2) THEN
               RHULT=RHULT1
               NWPANEL=NSW(1)-1
c            END IF
            CALL GHUB5
         ELSEIF(IHUB .EQ. 6) THEN
            CALL GHUB6
         ELSE IF(IHUB.EQ.7) THEN
            CALL GHUB7
         END IF

C.......Next IF statement added to prevent out-of-bound errors(JY072400)
         IF(NHBX.GT.NHPZ) THEN
            WRITE(*,'(''!!! WARNING !!! NOT ENOUGH PANELS IN HUB'')')
            WRITE(*,'('' NHBX = '',I5,'' NHBZ = '',I5)') NHBX,NHBZ
            WRITE(*,'('' Increase NHBZ --> NHBX in PARAM.INC! '')')
            STOP
         ELSE IF(MHBT.GT.MHBZ) THEN
            WRITE(*,'(''!!! WARNING !!! NOT ENOUGH PANELS IN HUB'')')
            WRITE(*,'('' MHBT = '',I5,'' MHBZ = '',I5)') MHBT,MHBZ
            WRITE(*,'('' Increase MHBZ --> MHBT in PARAM.INC! '')')
            STOP
         END IF

      ENDIF

 1000 CONTINUE

      DO M = 1, MRP
         DO N = 1, NCP
            HRZ(N,M) = SQRT( YB(N,M)**2+ZB(N,M)**2 )
            TERM1  =1.0 - HRZ(N,M)
            IF(TERM1.LT.0.0)THEN
               TERM1 = ABS(TERM1)
            ENDIF
            HRZSQ(N,M)=1.0-SQRT(TERM1)

C            CALL EVALDK(NX,1,XR,HRZ(N,M),HPITCH(N,M),PICUB)
C            CALL EVALDK(NX,1,XR,HRZ(N,M),HRAKE(N,M),RKCUB)
C            CALL EVALDK(NX,1,XR,HRZ(N,M),HSKEW(N,M),SKCUB)
C            CALL EVALDK(NX,1,XRSQ,HRZSQ(N,M),HCHORD(N,M),CHCUB)
C            CALL EVALDK(NX,1,XR,HRZ(N,M),HCAMBR(N,M),CICUB)
C            CALL EVALDK(NX,1,XR,HRZ(N,M),HTHICK(N,M),TICUB)
         ENDDO
      ENDDO

      IF(IHUB .EQ. 0) NHBX = 0

      NPANH=NHBX*MHBT
      NPANEL=NPANB+NPANH

 1111 continue

C***********************************************************************
C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
C***********************************************************************
c      IF(IAN .NE. 2) THEN
         NPANT = 0
         NPANC = 0
c      ENDIF
C***********************************************************************
C/e S.N.KIM | Aug. 2018.
C***********************************************************************

      NPANEL = NPANEL + NPAND

      IF(ISP .NE. 0) THEN
         CALL TRANS_GEO1(0)
         CALL TRANS_GEO2(0)
      ENDIF

C-----------------------------------------------------------------------
C     Prepare an output file for hidden line plotting
C-----------------------------------------------------------------------

      CALL PROPLT

      IF(ISC.EQ.1) CALL PROPLT_SR

C-----------------------------------------------------------------------
C     Plot the wake panels.                                    JY012299
C-----------------------------------------------------------------------

      CALL WAKEPLT

C --- if icon .eq. 9 : tunnel

      NPANTN1 = 0
      NPANTN2 = 0
      NPANTN = 0

      IF(ITUN .EQ. 1) THEN
         CALL TUNINP
         IF(ITUNGEO .EQ. 1) CALL TUNNELGEO1
         IF(ITUNGEO .EQ. 2) CALL TUNNELGEO2
         IF(ITUNGEO .EQ. 3) CALL TUNNELGEO3
         CALL TUNPLT
         NPANTN1 = NAX * MTUNEL
         NPANTN2 = NSIDE * MTUNEL
      ENDIF

C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)
C      NPANTN = NPANTN1 + 2 * NPANTN2
      NPANTN = NPANTN1 + NPANTN2
C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)

      NPANEL = NPANEL + NPANTN

C --- S.H.CHANG 04/20/2010 FOR CONTINUITY INSIDE OF A TUNNEL
      AREATUN1 = 0.0
      AREATUN2 = 0.0
      AREAHUB1 = 0.0
      AREAHUB2 = 0.0
      IF(ITUN.NE.0) THEN
        RTUN1 = TUNRAD
        RTUN2 = TUNRAD
        AREATUN1 = RTUN1**2
        AREATUN2 = RTUN2**2

        IF(IHUB.NE.0) THEN
          IF(XH(1,1) .LE. XTUN(1,1)) THEN
            AREAHUB1 = YH(1,1)**2+ZH(1,1)**2
          END IF
          IF(XH(NHBX,1) .GE. XTUN(NAX+1,1)) THEN
            AREAHUB2 = YH(NHBX,1)**2+ZH(NHBX,1)**2
          END IF
        END IF

        AREAR = (AREATUN1-AREAHUB1)/(AREATUN2-AREAHUB2)

        WRITE(*,*)
        WRITE(*,*) ' AREA RATIO FOR INFLOW =',AREAR
        WRITE(*,*)
      END IF
C --- S.H.CHANG 04/20/2010 FOR CONTINUITY INSIDE OF A TUNNEL


C-----------------------------------------------------------------------
C     Calculate abscissae for plotting
C     (Approximate positions of control points spanwise and chordwise)
C-----------------------------------------------------------------------
C
C.....Chordwise abscissa on the blade
C.....  SBP:  from l.e. (0.0) to t.e. (1.0)
      SBP(1)=HALF*SB(1)
      DO 50 I=2,NH
        SBP(I)=HALF*(SB(I)+SB(I-1))
   50 CONTINUE

C.....Spanwise abscissa on the blade

      DO M=1,MR
         RZP(M)=HALF*(RZ(M)+RZ(M+1))
         DO N = 1 , NC
            HRZP(N,M) =0.25*(HRZ(N,M)+HRZ(N,M+1)
     *                     + HRZ(N+1,M)+HRZ(N+1,M+1))
            HRZPSQ(N,M)=1.0-SQRT(ABS(1.0-HRZP(N,M)))
         ENDDO
      ENDDO

C-----------------------------------------------------------------------
C     Compute number of radii smaller than 0.95R (MRTIP)
C-----------------------------------------------------------------------
      DO 120 M=1,MR
         RZPSQ(M)=1.0-SQRT(ABS(1.0-RZP(M)))
         IF(RZP(M).LT.rmrtip) THEN
            MRTIP=M
         END IF
120   CONTINUE

C-----------------------------------------------------------------------
C     No smoothing for hydrofoil case.                          JY112998
C-----------------------------------------------------------------------
      IF((ICON.EQ.5).OR.(ICON.EQ.6) .OR.(ICON.EQ.8)) MRTIP=MR

C-----------------------------------------------------------------------
C     Write geometries to the output file
C-----------------------------------------------------------------------
      WRITE(2,900) (RZ(M),PITCH(M),RAKE(M),SKEW(M),CHORD(M),CAMBR(M),
     *              THICK(M),M=1,MRP)
      WRITE(2,910) NC,MR,NPANB, NHBX,MHBT,NHBDT,NPANH, NPANEL
      WRITE(2,920) NHBDK,NUWDK,NPANW, (M,NSW(M),M=1,MRP)
C-----------------------------------------------------------------------
C--------------------------- FORMAT ------------------------------------
C-----------------------------------------------------------------------
900   FORMAT(//5X,19(1H-),'INTERPOLATED BLADE GEOMETRY',15(1H-),/5X,
     *' R/RO      P/D      XS/D      SKEW    C/D     FO/C     TO/D'
     * /(3X,F8.4,F10.4,F9.4,F9.3,F8.4,2F9.4) )
910   FORMAT(//5X,'----------------NO. OF PANELS-------------------'
     */5X,'NC   =',I4,'   MR   =',I4,'   NPANB=',I4,
     */5X,'NHBX =',I4,'   MHBT =',I4,'   NHBDT  =',I4,'   NPANH=',I4,
     */5X,34X,'NPANEL=',I5 )
920   FORMAT(5X,'NHBDK=',I4,'   NUWDK=',I4,14X,'   NPANW=',I4,
     * /5X,'    M    NSW(M)'/(5X,I5,I8) )
C-----------------------------------------------------------------------

      RETURN
C))))))))))))))))))))) End of subroutine PROGEO ((((((((((((((((((((((((
      END
