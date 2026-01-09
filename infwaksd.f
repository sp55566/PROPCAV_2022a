      SUBROUTINE INFWAKSD
************************************************************************
*     INFWAKSD: Generate Duct wake subpanel                            *
*               and INFluence coeff. due to the Duct WAKe Subpanels    *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)


      IDDUCT2 = NPANB + NPANH
      IDDUCT3 = NPANB + NPANH + NPAND
 
      IAAA1 = NPANB
      IF(IHUB .NE. 0) IAAA1 = IAAA1 + NPANH
      IF(IDUCT .NE. 0) IAAA1 = IAAA1 + NPAND

      DO 50 M=1, MDUCT
        DO 20 N=1,NWSUB

           L=NWSUB*(M-1)+N
           XGW(L,1,1)=XWSD(N,M)
           XGW(L,1,2)=YWSD(N,M)
           XGW(L,1,3)=ZWSD(N,M)
           XGW(L,2,1)=XWSD(N,M+1)
           XGW(L,2,2)=YWSD(N,M+1)
           XGW(L,2,3)=ZWSD(N,M+1)
           XGW(L,3,1)=XWSD(N+1,M+1)
           XGW(L,3,2)=YWSD(N+1,M+1)
           XGW(L,3,3)=ZWSD(N+1,M+1)
           XGW(L,4,1)=XWSD(N+1,M)
           XGW(L,4,2)=YWSD(N+1,M)
           XGW(L,4,3)=ZWSD(N+1,M)
 20     CONTINUE
 50   CONTINUE

      NPWAKSD=NWSUB*MDUCT

      CALL GEO3DW(NPWAKSD,XGW,CHRLEWS,IER)

      IF(IER.EQ.0) THEN
        WRITE(*,*) 'UNACCEPTABLE PANEL IN INFWAKS'
        STOP
      END IF
      

C******  Geom info on the wake surface   by Hong Sun 08/24/2005  ******
      DO 60 J = 1, NPWAKSD

        RCPWD = SQRT(XCTW(J,2)**2+XCTW(J,3)**2)
        THPWD = ATAN2(XCTW(J,3),XCTW(J,2))
        XCTPDWs(J,1,1)=XCTW(J,1)
        XCTPDWs(J,2,1)=XCTW(J,2)
        XCTPDWs(J,3,1)=XCTW(J,3)

        DO K = 1, 3
         DO K1 = 1, 3
           DIRDWsCP(J,K1,K,1) = DIRW(J,K1,K)     
         ENDDO
        ENDDO

      IF(IHULL .EQ. 1 .and. ICON .EQ. 5) THEN
         XCTPDWs(J,1,2) = XCTW(J,1)
         XCTPDWs(J,2,2) = 2. - XCTW(J,2)
         XCTPDWs(J,3,2) = XCTW(J,3)
      ELSE
          IF(NBLADE.GT.1) THEN
           DO 55 KK=2,NBLADE
            XCTPDWs(J,1,KK)=XCTW(J,1)
            XCTPDWs(J,2,KK)=RCPWD*COS(THPWD+DELK*(KK-1))
            XCTPDWs(J,3,KK)=RCPWD*SIN(THPWD+DELK*(KK-1))
 55        CONTINUE    
          ENDIF 
        ENDIF

 60   CONTINUE

C**********************************************************************


      DO 160 M=1, MDUCT
        DO 90 I=1,NPANEL
          WSUBIFD(I,M)=ZERO
   90   CONTINUE
        DO 150 N=1,NWSUB
          J=(M-1)*NWSUB+N
          DO 100 K=1,4
            XV(K)=XVPW(J,K)
            YV(K)=YVPW(J,K)
            SIDE(K)=SIDW(J,K)
  100     CONTINUE
          DO 110 K=1,15
            S(K)=SSW(J,K)
  110     CONTINUE
C.........XM1,XM2,XM3,XM4 for Morino's formulation......................
          XM1(1)=XWSD(N,M)
          XM1(2)=YWSD(N,M)
          XM1(3)=ZWSD(N,M)
          XM2(1)=XWSD(N,M+1)
          XM2(2)=YWSD(N,M+1)
          XM2(3)=ZWSD(N,M+1)
          XM3(1)=XWSD(N+1,M+1)
          XM3(2)=YWSD(N+1,M+1)
          XM3(3)=ZWSD(N+1,M+1)
          XM4(1)=XWSD(N+1,M)             
          XM4(2)=YWSD(N+1,M)             
          XM4(3)=ZWSD(N+1,M) 
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
     %            (I .LE. IDDUCT2 .OR. I .GT. IDDUCT3)) THEN
                FD = 0.0
             ENDIF

             WSUBIFD(I,M)=WSUBIFD(I,M)
     *            +(0.5-(FLOAT(N)-0.5)/FLOAT(NWSUB))*FD
 140         CONTINUE
 150   CONTINUE
 160  CONTINUE

      RETURN
      END
