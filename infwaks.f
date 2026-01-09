      SUBROUTINE INFWAKS
************************************************************************
*     INFWAKS: INFluence coefficients due to the WAKe Subpanels        *
*      --- Calculate the influence coefficients due to the             *
*          wake subpanels                                              *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)
C-----------------------------------------------------------------------
C     Prepare parameters
C-----------------------------------------------------------------------

      IDDUCT2 = NPANB + NPANH
      IDDUCT3 = NPANB + NPANH + NPAND
 
      IAAA1 = NPANB
      IF(IHUB .NE. 0) IAAA1 = IAAA1 + NPANH
      IF(IDUCT .NE. 0) IAAA1 = IAAA1 + NPAND

      DO 50 M=MR,1,-1
C-----------------------------------------------------------------------
C
C        Prepare panel geometries of     (N,M+1) 2 *-----* 3 (N+1,M+1)
C          the wake, clockwise viewing             |     |
C          from the suction side of the            |     |
C          blade                           (N,M) 1 *-----* 4 (N+1,M)
C
C-----------------------------------------------------------------------
        DO 20 N=1,NWSUB

          L=NWSUB*(M-1)+N
          XGW(L,1,1)=XWS(N,M)
          XGW(L,1,2)=YWS(N,M)
          XGW(L,1,3)=ZWS(N,M)
          XGW(L,2,1)=XWS(N,M+1)
          XGW(L,2,2)=YWS(N,M+1)
          XGW(L,2,3)=ZWS(N,M+1)
          XGW(L,3,1)=XWS(N+1,M+1)
          XGW(L,3,2)=YWS(N+1,M+1)
          XGW(L,3,3)=ZWS(N+1,M+1)
          XGW(L,4,1)=XWS(N+1,M)
          XGW(L,4,2)=YWS(N+1,M)
          XGW(L,4,3)=ZWS(N+1,M)
   20    CONTINUE
   50 CONTINUE

C.....NPWAKS is the total number of SUBpanels in the wake...............
      NPWAKS=NWSUB*MR
      CALL GEO3DW(NPWAKS,XGW,CHRLEWS,IER)
      IF(IER.EQ.0) THEN
        WRITE(*,*) 'UNACCEPTABLE PANEL IN INFWAKS'
        STOP
      END IF

C****** Geom info on the wake surface   by Hong Sun 04/06/2005  *******
      DO 60 J = 1, NPWAKS

        RCPW = SQRT(XCTW(J,2)**2+XCTW(J,3)**2)
        THPW = ATAN2(XCTW(J,3),XCTW(J,2))
        XCTPWs(J,1,1)=XCTW(J,1)
        XCTPWs(J,2,1)=XCTW(J,2)
        XCTPWs(J,3,1)=XCTW(J,3)

        DO K = 1, 3
         DO K1 = 1, 3
           DIRWsCP(J,K1,K,1) = DIRW(J,K1,K)     
         ENDDO
        ENDDO

      IF(IHULL .EQ. 1 .and. ICON .EQ. 5) THEN
         XCTPWs(J,1,2) = XCTW(J,1)
         XCTPWs(J,2,2) = 2. - XCTW(J,2)
         XCTPWs(J,3,2) = XCTW(J,3)
      ELSE
          IF(NBLADE.GT.1) THEN
           DO 55 KK=2,NBLADE
            XCTPWs(J,1,KK)=XCTW(J,1)
            XCTPWs(J,2,KK)=RCPW*COS(THPW+DELK*(KK-1))
            XCTPWs(J,3,KK)=RCPW*SIN(THPW+DELK*(KK-1))
 55        CONTINUE    
          ENDIF 
        ENDIF

 60   CONTINUE
C***********************************************************************

C/s S.N.KIM | XCTPWs for unsteady boundary layer correction.
C***********************************************************************
      IF (IAN.EQ.2.AND.NBLADE.GT.1) THEN
        CALL callgeow
        DO N = 1, NWPANEL + 1
          DO M = 1, MR + 1
            DO KK = 2, NBLADE
              IK = KK - 1
              THT1 = -IK * DELK  ! DELK = TWOPI/NBLADE
              XWW1(N,M,KK) = XWO(N,M,IK) 
              YWW1(N,M,KK) = YWO(N,M,IK)*COS(THT1) 
     &                     - ZWO(N,M,IK)*SIN(THT1)  
              ZWW1(N,M,KK) = YWO(N,M,IK)*SIN(THT1)
     &                     + ZWO(N,M,IK)*COS(THT1)
            ENDDO
          ENDDO
        ENDDO

        DO KK = 2, NBLADE  
          DO M = 1, MR + 1
            XWSO(1,M,KK) = XWW1(1,M,KK)
            YWSO(1,M,KK) = YWW1(1,M,KK)
            ZWSO(1,M,KK) = ZWW1(1,M,KK)
            DO N = 1, NSUB
              DWX = (XWW1(N+1,M,KK) - XWW1(N,M,KK))/FLOAT(NWSUB1)
              DWY = (YWW1(N+1,M,KK) - YWW1(N,M,KK))/FLOAT(NWSUB1)
              DWZ = (ZWW1(N+1,M,KK) - ZWW1(N,M,KK))/FLOAT(NWSUB1)
              DO L = 1, NWSUB1
                NIDX = (N-1)*NWSUB1 + L + 1
                XWSO(NIDX,M,KK) = XWSO(NIDX-1,M,KK) + DWX
                YWSO(NIDX,M,KK) = YWSO(NIDX-1,M,KK) + DWY
                ZWSO(NIDX,M,KK) = ZWSO(NIDX-1,M,KK) + DWZ
              ENDDO
            ENDDO
          ENDDO  

          DO M = MR, 1, -1
            DO N = 1, NWSUB
              L = NWSUB*(M-1) + N
              XGW(L,1,1) = XWSO(N,M,KK)
              XGW(L,1,2) = YWSO(N,M,KK)
              XGW(L,1,3) = ZWSO(N,M,KK)
              XGW(L,2,1) = XWSO(N,M+1,KK)
              XGW(L,2,2) = YWSO(N,M+1,KK)
              XGW(L,2,3) = ZWSO(N,M+1,KK)
              XGW(L,3,1) = XWSO(N+1,M+1,KK)
              XGW(L,3,2) = YWSO(N+1,M+1,KK)
              XGW(L,3,3) = ZWSO(N+1,M+1,KK)
              XGW(L,4,1) = XWSO(N+1,M,KK)
              XGW(L,4,2) = YWSO(N+1,M,KK)
              XGW(L,4,3) = ZWSO(N+1,M,KK)
            ENDDO
          ENDDO
  
          CALL GEO3DW(NPWAKS,XGW,CHRLEWS,IER)
          IF(IER.EQ.0) THEN
            WRITE(*,*) 'UNACCEPTABLE PANEL IN INFWAKS'
            STOP
          ENDIF

          DO J = 1, NPWAKS
            DO K = 1, 3
              DO K1 = 1, 3
                DIRWsCP(J,K1,K,KK) = DIRW(J,K1,K)
              ENDDO
            ENDDO
            XCTPWs(J,1,KK) = XCTW(J,1)
            XCTPWs(J,2,KK) = XCTW(J,2)
            XCTPWs(J,3,KK) = XCTW(J,3)
          ENDDO   

          IF(IVISC.NE.0) THEN
            DO M = 1, MR + 1
              DO N = 1, NWSUB
                XWVS3DO(N,M) = XWSO(N,M,KK)
                YWVS3DO(N,M) = YWSO(N,M,KK)
                ZWVS3DO(N,M) = ZWSO(N,M,KK)
              ENDDO
              DO N = NSUB + 1, NWPANEL + 1
                NIDX = N + NWSUB - NSUB
                XWVS3DO(NIDX,M) = XWW1(N,M,KK)
                YWVS3DO(NIDX,M) = YWW1(N,M,KK)
                ZWVS3DO(NIDX,M) = ZWW1(N,M,KK)
              ENDDO          
            ENDDO
            NWV = NWPANEL + NWSUB - NSUB 
            NPWAKEV = NWV*MR
            DO M = MR, 1, -1
              DO N = 1, NWV
                L = (MR - M)*NWV + N
                XGW(L,1,1)=xwvs3do(N,M)
                XGW(L,1,2)=ywvs3do(N,M)
                XGW(L,1,3)=zwvs3do(N,M)
                XGW(L,2,1)=xwvs3do(N,M+1)
                XGW(L,2,2)=ywvs3do(N,M+1)
                XGW(L,2,3)=zwvs3do(N,M+1)
                XGW(L,3,1)=xwvs3do(N+1,M+1)
                XGW(L,3,2)=ywvs3do(N+1,M+1)
                XGW(L,3,3)=zwvs3do(N+1,M+1)
                XGW(L,4,1)=xwvs3do(N+1,M)
                XGW(L,4,2)=ywvs3do(N+1,M)
                XGW(L,4,3)=zwvs3do(N+1,M)
              ENDDO 
            ENDDO
            CALL GEO3DW(NPWAKEV,XGW,CHRLEWS,IER)
            IF(IER.EQ.0) THEN
              WRITE(*,*) 'UNACCEPTABLE PANEL IN INFWAKS'
              STOP
            ENDIF
            DO J = 1, NPWAKEV
              XCTPWV(J,1,KK) = XCTW(J,1)
              XCTPWV(J,2,KK) = XCTW(J,2)
              XCTPWV(J,3,KK) = XCTW(J,3)
            ENDDO
          ENDIF

        ENDDO
      ENDIF
C***********************************************************************
C/e S.N.KIM | Aug. 2018.

C-----------------------------------------------------------------------
C        Calculate influence coefficients of the wake.
C        W(I,M):  induced potentials at I due to Mth helical strip of 
C                 dipoles
C-----------------------------------------------------------------------
      DO 160 M=MR,1,-1
        DO 90 I=1,NPANEL
          WSUBIF(I,M)=ZERO
   90   CONTINUE
        DO 150 N=1,NWSUB
          J=(M-1)*NWSUB+N
C.........Transfer data to the common block /GEOM/......................
          DO 100 K=1,4
            XV(K)=XVPW(J,K)
            YV(K)=YVPW(J,K)
            SIDE(K)=SIDW(J,K)
  100     CONTINUE
          DO 110 K=1,15
            S(K)=SSW(J,K)
  110     CONTINUE
C.........XM1,XM2,XM3,XM4 for Morino's formulation......................
          XM1(1)=XWS(N,M)
          XM1(2)=YWS(N,M)
          XM1(3)=ZWS(N,M)
          XM2(1)=XWS(N,M+1)
          XM2(2)=YWS(N,M+1)
          XM2(3)=ZWS(N,M+1)
          XM3(1)=XWS(N+1,M+1)
          XM3(2)=YWS(N+1,M+1)
          XM3(3)=ZWS(N+1,M+1)
          XM4(1)=XWS(N+1,M)             
          XM4(2)=YWS(N+1,M)             
          XM4(3)=ZWS(N+1,M) 
          IMR0=0    

          DO 140 I=1,NPANEL
C...........Transfer control points to local coordinate.................
             XLOC=ZERO
             YLOC=ZERO
             ZLOC=ZERO
             DO 120 K=1,3
                XLOC=XLOC+(XCTP(I,K,1)-XCTW(J,K))*DIRW(J,1,K)
                YLOC=YLOC+(XCTP(I,K,1)-XCTW(J,K))*DIRW(J,2,K)
                ZLOC=ZLOC+(XCTP(I,K,1)-XCTW(J,K))*DIRW(J,3,K)
 120         CONTINUE

C...........Compute the induced potentials due to the blades............
C...........modified 04-11-90...........................................

             IMR=IMR0
             CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,FSY,
     *            FDX,FDY,FDZ,0,IMR)
C...........Near field use Morino's formulation.........................
             IF(IMR.EQ.2) THEN
                DO 122 IXYZ=1,3
                   XMC(IXYZ)=XCTP(I,IXYZ,1)
 122            CONTINUE
                CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
             END IF

c............if not self-inf. functions, and fd > 6.28, then............
c............I will correct to be 0.....................................
CJY             IF(ABS(FD).GT.6.27) THEN
             IF(ABS(FD).GT.6.28) THEN
C                WRITE(*,*) 'INFWAKS-A',NTSTEP,I,M,N,IMR,FD
                FD=0.0
             END IF

             IF(ITUN .EQ. 1 .AND. I .GT. IAAA1) FD = -FD
             IF(IDUCT .EQ. 1 .AND. IDOPT .EQ. 1 .AND. 
     %            I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                FD = 0.0
             ENDIF
             WSUBIF(I,M)=WSUBIF(I,M)
     *            +(0.5-(FLOAT(N)-0.5)/FLOAT(NWSUB))*FD
 140         CONTINUE
 150   CONTINUE
 160  CONTINUE

cC------------------------------(S.H.CHANG 02/25/2010)------------------------- 
cC    OUTPUT GEOMETRIES IN THE WAKE SUBPANELS FOR HULLFPP 
c      WRITE(737,*) NPWAKS,NWSUB
c      DO L = 1,NPWAKS
c         DO M = 1,4
c            DO N = 1,3
c               WRITE(737,*) XG(L,M,N)
c            END DO
c         END DO
c      END DO 
cC------------------------------(S.H.CHANG 02/25/2010)------------------------- 

      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine INFWAKS>>>>>>>>>>>>>>>>>>>>>>>>
      END
