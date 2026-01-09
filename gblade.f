      SUBROUTINE GBLADE
************************************************************************
*     GBLADE: Geometry of the propeller BLADE                          *
*      --- Generate a blade surface geometry                           *
*                                                                      *
*  Date     Comment or Revision                                        *
*  -------- -------------------                                        *
*  JY010798 Added a parabolic camber and thickness option for          *
*           blade geometry definition.  Here are the new               *
*           additions:                                                 *
*           ITHK = 10   Parabolic thickness form                       *
*           ICAM = 0    NACA A=0.8 mean line                           *
*           ICAM = 1    Parabolic camber line                          *
*  JY112198 Added a new camber and a new thickness option:             *
*                     ICAM=99  User input camber distribution          *
*                     ITHK=99  User input thickness distribution       *
*                                                                      *
************************************************************************

      use m_BLADEM
      use m_INPGEO2
      use su_inp
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'INPGEO.INC'
      DOUBLE PRECISION DDZ,DDY
      DIMENSION XIMG(2),YIMG(2),XYIMG(4)
      DIMENSION XIM(NBHZP),ETAM(NBHZP),XIOLD(NBPZ,MBPZ),
     *     ETAOLD(NBPZ,MBPZ),XIMOLD(NBHZP,MBPZ),ETAMOLD(NBHZP,MBPZ)
      DIMENSION XBRIN(MBPZ),YBRIN(MBPZ),ZBRIN(MBPZ)
      DIMENSION XBROUT(MBPZ),YBROUT(MBPZ),ZBROUT(MBPZ)
!     COMMON/BLADEM/XBM(NBHZP,MBPZ),YBM(NBHZP,MBPZ),ZBM(NBHZP,MBPZ)

      open (unit = 200, file = "bsurf.plt")
      open (unit = 201, file = "bcam.plt")

C-----------------------------------------------------------------------
C     Set up chordwise spacing
C-----------------------------------------------------------------------
      IF(ICSPAC.EQ.2.OR.ICSPAC.EQ.3) THEN
         RLET=RADLE(XTI(1),XCHD(1),ITHK)
      END IF
      CALL SPACE(NC,ICSPAC,RLET,SB)

C -- Create geometry for blade section when ITHK=99, ICAM =99, and IFORMAT=2

C-- If ITHK =99, ICAM=99, and IFORMAT=2
C   Directly interpolate blade section geometry using SPLINE

C-----------------------------------------------------------------------
C     Generate blade geometry
C-----------------------------------------------------------------------
      DO 50 M=1,MRP
C
C........Midchord centered local coordinate xi and eta
C........(nondimensionalized by R
C........CHORD=C/D,THK=T/D,CAMBR=F/C,YTC=Y/D,YCC=Y/C,SB=S/C,RAKE=XM/D

C-----------------------------------------------------------------------
C     ITHK=99, user input thickness distribution.               JY112798
C-----------------------------------------------------------------------

            IF(ITHK.EQ.99) THEN
               CALL THKINP(NH,XTS1(1,M),SB,YTC)
               GO TO 1000
            END IF

            IF(THICK(M).EQ.0.0) THEN
               DO 10 N=1,NH
                  YTC(N)=ZERO
 10            CONTINUE
            ELSE
               IF(ITHK.EQ.1) THEN
                  CALL NACA66(NH,THICK(M),RLE,SB,YTC)
               ELSE IF(ITHK.EQ.2) THEN
                  CALL RAE(NH,THICK(M),SB,YTC)
               ELSE IF(ITHK.EQ.3) THEN
                  CALL NACA65A(NH,THICK(M),RLE,SB,YTC)
               ELSE IF(ITHK.EQ.4) THEN
                  CALL NACA64A(NH,THICK(M),RLE,SB,YTC)
               ELSE IF(ITHK.EQ.5) THEN
                  CALL NACA00(NH,THICK(M),RLE,SB,YTC)
               ELSE IF(ITHK.EQ.6) THEN
                  CALL ELIPSE(NH,THICK(M),RLE,SB,YTC)
C--------------------------------------------------------------------------
C ITHK=8, NAVSEA TYPE_ONE (NOTE: NACA16 is NAVSEA TYPE_TWO) SNKIM 08232017
C--------------------------------------------------------------------------
               ELSE IF(ITHK.EQ.8) THEN
                  CALL NAVSEA1(NH,THICK(M),RLE,SB,YTC)
               ELSE IF(ITHK.EQ.9) THEN
                  CALL NACA16(NH,THICK(M),RLE,SB,YTC)
               ELSE IF(ITHK.EQ.10) THEN
                  CALL TKPB(NH,THICK(M),RLE,SB,YTC)
               END IF

            END IF

 1000       CONTINUE

C-----------------------------------------------------------------------
C        New if statement allows user to select camber distribution
C        using the variable ICAM.  The available cambers are:
C
C        ICAM = 0      A=0.8 Mean Line Camber
C        ICAM = 1      Parabolic Camber Distribution           JY010898
C        ICAM = 2      NACA 65 Section Meanline             SNKIM082317
C-----------------------------------------------------------------------
            IF(ICAM.EQ.0) THEN
               CALL A8ML(NH,CAMBR(M),SB,YCC,TCC)
            ELSE IF(ICAM.EQ.1) THEN
               CALL CBPB(NH,CAMBR(M),SB,YCC,TCC)
            ELSE IF(ICAM.EQ.2) THEN                        !SNKIM082317
               CALL NACA65(NH,CAMBR(M),SB,YCC,TCC)
C-----------------------------------------------------------------------
C       ICAM=99, user input thickness distribution.            JY112798
C-----------------------------------------------------------------------
            ELSE IF(ICAM.EQ.99) THEN
               CALL CAMINP(NH,XCS1(1,M),SB,YCC,TCC)
            END IF

C........local coordinate xi and eta
            XI(NHP)=-CHORD(M)
            ETA(NHP)=ZERO
            XIM(1)=XI(NHP)
            ETAM(1)=ZERO
            DO 20 N=1,NH

               IF(ITHK .EQ. 99 .AND. ICAM .EQ. 99
     %              .AND. IFORMAT .EQ. 2) THEN
                  DXT = 0.0
                  DYT = YTC(N)
               ELSE
                  DXT=YTC(N)*SIN(TCC(N))
                  DYT=YTC(N)*COS(TCC(N))
               ENDIF

               NBOT=NHP-N
               NTOP=NHP+N
               NMEAN=N+1
               XI(NTOP)=( (SB(N)-.5)*CHORD(M)-DXT )*TWO
               ETA(NTOP)=(YCC(N)*CHORD(M) +DYT)*TWO
               XI(NBOT)=( (SB(N)-.5)*CHORD(M)+DXT )*TWO
               ETA(NBOT)=(YCC(N)*CHORD(M) -DYT)*TWO
               XIM(NMEAN)=(SB(N)-.5)*CHORD(M)*TWO
               ETAM(NMEAN)=YCC(N)*CHORD(M)*TWO
 20         CONTINUE

         IF(ISC.EQ.1) THEN
C..........Store coordinates of original blades.
            DO N=1,NCP
               XIOLD(N,M)=XI(N)
               ETAOLD(N,M)=ETA(N)
               IF(N.LE.NHP) THEN
                  XIMOLD(N,M)=XIM(N)
                  ETAMOLD(N,M)=ETAM(N)
               END IF
            END DO
         END IF

C........Reference coordinate of blade element
         IF(RZ(M).EQ.0.)THEN
            RZ(M)=1
         ENDIF
         TANP=PITCH(M)/PI/RZ(M)

         COSP=1.0/SQRT(1+TANP*TANP)
         SINP=TANP*COSP
         DO 30 N=1,NCP
            DX= XI(N)*SINP-ETA(N)*COSP
            XB(N,M)=RAKE(M)*TWO +DX
            THETA=SKEW(M)*RAD+(XI(N)*COSP+ETA(N)*SINP)/RZ(M)
            IF(M.EQ.1) THEN
               THR(N)=THETA
            END IF
            IF(N.EQ.1) THEN
               THT(M)=THETA
            END IF
            YB(N,M)=RZ(M)*COS(THETA)
            ZB(N,M)=RZ(M)*SIN(THETA)
30       CONTINUE
         DO 40 N=1,NHP
            DX= XIM(N)*SINP-ETAM(N)*COSP
            XBM(N,M)=RAKE(M)*TWO +DX
            THETA=SKEW(M)*RAD+ (XIM(N)*COSP+ETAM(N)*SINP)/RZ(M)
            YBM(N,M)=RZ(M)*COS(THETA)
            ZBM(N,M)=RZ(M)*SIN(THETA)
40       CONTINUE

C.......Plot blade sectsions (JY071201)
         IF(M.EQ.1) THEN
            OPEN(55,FILE='bldsec.plt',STATUS='UNKNOWN')
            WRITE(55,*) 'VARIABLES="c/R","y/R"'
         END IF
         WRITE(55,*) 'ZONE T="M=',M,'"'
         DO N=1,NCP
            WRITE(55,*) XI(N)+CHORD(M),ETA(N)+RZ(M)
C...... Added XIV, ETAV for 2-D strip node coordinates calculation
C                    By  Hong Sun            11/04/05
            XIV(N,M) = XI(N)
            ETAV(N,M) = ETA(N)
         END DO
         IF(M.EQ.MRP) CLOSE(55)

50    CONTINUE


C-----Finished building the blade geometry
C-----If IDRT=1, use node points location on panel surface
C-----Yiran 6/3/14
      IF (IDRT.EQ.1) THEN
        CALL PANELINP

C Blade cut by the hub
        IF (IHUB.NE.6) THEN
          DO N=1,NCP
            XIMG(1)=XB(N,1)
            XIMG(2)=XB(N,2)
            YIMG(1)=SQRT(YB(N,1)**2+ZB(N,1)**2)
            YIMG(2)=SQRT(YB(N,2)**2+ZB(N,2)**2)
            IF (RHUB.GT.YIMG(2)) THEN
              WRITE(*,*) "Error!"
              WRITE(*,*) "RHUB not consistent with input panels!"
              STOP
            END IF
            TTMPS=(RHUB-YIMG(1))/(YIMG(2)-YIMG(1))
            IF ((ABS(TTMPS)).GT.(0.5)) THEN
              WRITE(*,*) "Warning:RHUB not consistent with input panels"
            END IF
            XB(N,1)=XIMG(1)*(1.0E0-TTMPS)+XIMG(2)*TTMPS
            YIMG(1)=ATAN2(ZB(N,1),YB(N,1))
            YIMG(2)=ATAN2(ZB(N,2),YB(N,2))
            TTMPS2=YIMG(1)*(1.0E0-TTMPS)+YIMG(2)*TTMPS
            YB(N,1)=RHUB*COS(TTMPS2)
            ZB(N,1)=RHUB*SIN(TTMPS2)
          END DO

          DO N=1,NHP
            XIMG(1)=XBM(N,1)
            XIMG(2)=XBM(N,2)
            YIMG(1)=SQRT(YBM(N,1)**2+ZBM(N,1)**2)
            YIMG(2)=SQRT(YBM(N,2)**2+ZBM(N,2)**2)
            IF (RHUB.GT.YIMG(2)) THEN
              WRITE(*,*) "Error!"
              WRITE(*,*) "RHUB not consistent with input panels!"
              STOP
            END IF
            TTMPS=(RHUB-YIMG(1))/(YIMG(2)-YIMG(1))
            IF ((ABS(TTMPS)).GT.(0.5)) THEN
              WRITE(*,*) "Warning:RHUB not consistent with input panels"
            END IF
            XBM(N,1)=XIMG(1)*(1.0E0-TTMPS)+XIMG(2)*TTMPS
            YIMG(1)=ATAN2(ZBM(N,1),YBM(N,1))
            YIMG(2)=ATAN2(ZBM(N,2),YBM(N,2))
            TTMPS2=YIMG(1)*(1.0E0-TTMPS)+YIMG(2)*TTMPS
            YBM(N,1)=RHUB*COS(TTMPS2)
            ZBM(N,1)=RHUB*SIN(TTMPS2)
          END DO
        ELSE  ! IHUB=6
          DO N=1,NCP
            XIMG(1)=XB(N,1)
            XIMG(2)=XB(N,2)
            YIMG(1)=SQRT(YB(N,1)**2+ZB(N,1)**2)
            YIMG(2)=SQRT(YB(N,2)**2+ZB(N,2)**2)
            CALL HUBTRIM(NHIN,XHIN,XYHINCUB,XIMG,YIMG,TTMPS,RTMPS)
            XB(N,1)=TTMPS
            TTMPS=(RTMPS-YIMG(1))/(YIMG(2)-YIMG(1))
            YIMG(1)=ATAN2(ZB(N,1),YB(N,1))
            YIMG(2)=ATAN2(ZB(N,2),YB(N,2))
            TTMPS2=YIMG(1)*(1.0E0-TTMPS)+YIMG(2)*TTMPS
            YB(N,1)=RTMPS*COS(TTMPS2)
            ZB(N,1)=RTMPS*SIN(TTMPS2)
          END DO
          DO N=1,NHP
            XIMG(1)=XBM(N,1)
            XIMG(2)=XBM(N,2)
            YIMG(1)=SQRT(YBM(N,1)**2+ZBM(N,1)**2)
            YIMG(2)=SQRT(YBM(N,2)**2+ZBM(N,2)**2)
            CALL HUBTRIM(NHIN,XHIN,XYHINCUB,XIMG,YIMG,TTMPS,RTMPS)
            XBM(N,1)=TTMPS
            TTMPS=(RTMPS-YIMG(1))/(YIMG(2)-YIMG(1))
            YIMG(1)=ATAN2(ZBM(N,1),YBM(N,1))
            YIMG(2)=ATAN2(ZBM(N,2),YBM(N,2))
            TTMPS2=YIMG(1)*(1.0E0-TTMPS)+YIMG(2)*TTMPS
            YBM(N,1)=RTMPS*COS(TTMPS2)
            ZBM(N,1)=RTMPS*SIN(TTMPS2)
          END DO
        END IF
      END IF
C-----Finished Yiran

      WRITE(200,*) NCP,MRP
      DO M=1,MRP
              WRITE(200,*)RZ(M)
        ENDDO
      DO N=1,NCP
              DO M=1,MRP
                     WRITE(200,'(3F10.6)')XB(N,M),YB(N,M),ZB(N,M)
              ENDDO
         ENDDO
      CLOSE(200)

      WRITE(201,*) NHP,MRP
      DO N=1,NHP
              DO M=1,MRP
                     WRITE(201,'(3F10.6)')XBM(N,M),YBM(N,M),ZBM(N,M)
              ENDDO
         ENDDO
      CLOSE(201)

      IF ((IHUB.EQ.6).AND.(IDRT.NE.1)) THEN
         write(*,*) ' Gblade : IHUB=6'
           DO N=1,NCP
             DO M=1,MRP
                XBRIN(M) = XB(N,M)
                YBRIN(M) = YB(N,M)
                ZBRIN(M) = ZB(N,M)
             END DO
             CALL PROPEXT(XBRIN,YBRIN,ZBRIN,MRP,IRSPAC,
     %                    XBROUT,YBROUT,ZBROUT,
     %                    XHIN,YHIN,XYHINCUB,NHIN,0)
             DO M =1,MRP
                XB(N,M) = XBROUT(M)
                YB(N,M) = YBROUT(M)
                ZB(N,M) = ZBROUT(M)
             END DO
          END DO

C-- Regenerate camber surface

           DO N=1,NHP
             DO M=1,MRP
                XBRIN(M) = XBM(N,M)
                YBRIN(M) = YBM(N,M)
                ZBRIN(M) = ZBM(N,M)
             END DO
             CALL PROPEXT(XBRIN,YBRIN,ZBRIN,MRP,IRSPAC,
     %                    XBROUT,YBROUT,ZBROUT,
     %                    XHIN,YHIN,XYHINCUB,NHIN,0)
             DO M =1,MRP
                XBM(N,M) = XBROUT(M)
                YBM(N,M) = YBROUT(M)
                ZBM(N,M) = ZBROUT(M)
             END DO
          END DO
      END IF

CRyan DRTINP
      IF (IDRT.NE.1) THEN
         OPEN(223,FILE='panel_opt.plt',STATUS='UNKNOWN')
         WRITE(223,*) NCP,MRP
         DO M=1,MRP
            DO N=1,NCP
               WRITE(223,"(3F24.20)") XB(N,M),YB(N,M),ZB(N,M)
            ENDDO
         ENDDO
         CLOSE(223)
      END IF
CRyan

C-----------------------------------------------------------------------
C     Create geometry for the separated region.                 JY081601
C-----------------------------------------------------------------------
      IF(ISC.EQ.1) CALL SRGEO(XIOLD,ETAOLD,XIMOLD,ETAMOLD)

C-----------------------------------------------------------------------
C     Compute the normal vectors of the camber surface
C-----------------------------------------------------------------------
      DO 90 N=1,NH
         DO 80 M=1,MR
            XG(1,1,1)=XBM(N,M)
            XG(1,1,2)=YBM(N,M)
            XG(1,1,3)=ZBM(N,M)
            XG(1,2,1)=XBM(N,M+1)
            XG(1,2,2)=YBM(N,M+1)
            XG(1,2,3)=ZBM(N,M+1)
            XG(1,3,1)=XBM(N+1,M+1)
            XG(1,3,2)=YBM(N+1,M+1)
            XG(1,3,3)=ZBM(N+1,M+1)
            XG(1,4,1)=XBM(N+1,M)
            XG(1,4,2)=YBM(N+1,M)
            XG(1,4,3)=ZBM(N+1,M)

            CALL GEOM3D(1,XG,CHRLEPS,IER)
            IF(IER.EQ.0) THEN
               WRITE(*,'(A)') ' UNACCEPTABLE PANELS ON CAMBER SURFACE'
               STOP
            END IF
            XCON(N,M)=VEL(1,1)
            YCON(N,M)=VEL(1,2)
            ZCON(N,M)=VEL(1,3)
80       CONTINUE
90    CONTINUE

      RETURN
C))))))))))))))))))))) End of subroutine GBLADE ((((((((((((((((((((((((
      END


      SUBROUTINE PANELINP
C-----Finished building the blade geometry
C-----If IDRT=1, use node points location on panel surface
C-----Yiran 6/3/14
      use m_BLADEM
      use m_INPGEO2
      use su_inp
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'INPGEO.INC'
      DOUBLE PRECISION XXX,YYY,ZZZ
c      REAL XXX,YYY,ZZZ
C/s S.N.KIM | Chord length calculation
      DOUBLE PRECISION CSUM, CDELL
      DIMENSION P1(3), P2(3)
      DIMENSION SC1(1201)
C/e S.N.KIM | Oct. 2018
      OPEN(197,FILE='panel_inp.dat',STATUS='OLD')
      READ(197,*) NN1,MM1
      DO M = 1,MM1
         DO N = 1,NN1
            READ(197,*) XXX, YYY, ZZZ
            XB(N,M) = REAL(XXX)
            YB(N,M) = REAL(YYY)
            ZB(N,M) = REAL(ZZZ)
            IF (N.EQ.NN1) THEN ! Trailing Edge Closure. Without this, Kutta Condition breaks down. Super important!
              XB(1,M) = XB(NN1,M)
              YB(1,M) = YB(NN1,M)
              ZB(1,M) = ZB(NN1,M)
            ENDIF
         END DO
      END DO
      CLOSE(197)

      DO M = 1, MM1
        DO N = 1, NN1
          IF (M.EQ.1) THR(N) = REAL(ATAN2(ZB(N,M),YB(N,M)))
          IF (N.EQ.1) THT(M) = REAL(ATAN2(ZB(N,M),YB(N,M)))
        ENDDO
      ENDDO

      DO M = 1,MRP
         DO N = 1,NHP
            XBM(N,M) = 0.5*(XB(NHP+1-N,M)+XB(NHP-1+N,M))
            YBM(N,M) = 0.5*(YB(NHP+1-N,M)+YB(NHP-1+N,M))
            ZBM(N,M) = 0.5*(ZB(NHP+1-N,M)+ZB(NHP-1+N,M))
         END DO
      END DO

C/s S.N.KIM | Chord length calculation. Note this calculation is not perfectly accurate
C           | since we go over mean camber line at each station, not the actual chord line. 
      DO M = 1, MRP
        CSUM = 0.0
        IF(M.EQ.1) SC1(1) = CSUM
        DO N = 1, NH
          P1(1) = XBM(N+1,M)
          P1(2) = YBM(N+1,M)
          P1(3) = ZBM(N+1,M)
          P2(1) = XBM(N,M)
          P2(2) = YBM(N,M)
          P2(3) = ZBM(N,M)
          CDELL = SQRT(DOT_PRODUCT(P1-P2,P1-P2))
          CSUM = CSUM + CDELL
          IF(M.EQ.1) THEN
            SC1(N+1) = CSUM
            IF(N.EQ.NH) THEN
              DO I = 1, NHP
                SC1(I) = SC1(I)/CSUM
              ENDDO
            ENDIF
          ENDIF
        ENDDO
        CHORD(M) = .5*CSUM
      ENDDO
C/e S.N.KIM | Oct. 2018.

      DO N = 1, NH
        SB(N) = SC1(N+1)
      ENDDO

      RETURN
      END


      SUBROUTINE HUBTRIM(NHIN,XHIN,XYHINCUB,XIMG,YIMG,TTMPS,RTMPS)
      IMPLICIT NONE
      INTEGER NHIN,I,J,K
      REAL XYHINCUB(*),XHIN(*),XIMG(2),YIMG(2),TTMPS,RTMPS
      REAL X1,X2,R1,R2,RR1,RR2,X0,R0,RR0

      X1=XIMG(1)*2.0E0-XIMG(2)
      X2=XIMG(2)
      R1=YIMG(1)*2.0E0-YIMG(2)
      R2=YIMG(2)
      CALL EVALDKs(NHIN,1,XHIN,X1,RR1,XYHINCUB)
      CALL EVALDKs(NHIN,1,XHIN,X2,RR2,XYHINCUB)

      IF (((R1-RR1)*(R2-RR2)).GT.(0.0E0)) THEN
        WRITE(*,*) "Error!"
        WRITE(*,*) "RHUB not consistent with input panels!"
        STOP
      END IF

      DO I=1,20
        X0=(X1+X2)/2.0E0
        R0=(R1+R2)/2.0E0
        CALL EVALDKs(NHIN,1,XHIN,X0,RR0,XYHINCUB)
        IF (RR0.GE.R0) THEN
          X1=X0
          R1=R0
        ELSE
          X2=X0
          R2=R0
        END IF
      END DO

      TTMPS=(X1+X2)/2.0E0
      RTMPS=(R1+R2)/2.0E0

      RETURN
      END
