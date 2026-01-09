      SUBROUTINE INFWAKS_IM
************************************************************************
*     INFWAKS: INFluence coefficients due to the WAKe Subpanels Images *
*      --- Calculate the influence coefficients due to the             *
*          wake subpanels                                              *
*                                                                      *
*     Date      Comments                                               *
*     --------- -----------------------------------------------------  *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)

C-----------------------------------------------------------------------
C     Prepare parameters
C-----------------------------------------------------------------------

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
 20      CONTINUE
 50   CONTINUE

C.....NPWAKS is the total number of SUBpanels in the wake...............
      NPWAKS=NWSUB*MR
      CALL GEO3DW(NPWAKS,XGW,CHRLEWS,IER)
      IF(IER.EQ.0) THEN
         WRITE(*,*) 'UNACCEPTABLE PANEL IN INFWAKS'
         STOP
      END IF

C-----------------------------------------------------------------------
C        Calculate influence coefficients of the wake.
C        W(I,M):  induced potentials at I due to Mth helical strip of 
C                 dipoles
C-----------------------------------------------------------------------
      DO 160 M=MR,1,-1
         DO 90 I=1,NPANEL
            WSUBIF_IM(I,M)=ZERO
 90      CONTINUE
         DO 150 N=1,NWSUB
            J=(M-1)*NWSUB+N
C.........Transfer data to the common block /GEOM/......................
            DO 100 K=1,4
               XV(K)=XVPW(J,K)
               YV(K)=YVPW(J,K)
               SIDE(K)=SIDW(J,K)
 100        CONTINUE
            DO 110 K=1,15
               S(K)=SSW(J,K)
 110        CONTINUE
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

               IF(ISUBM(I,IDXREV).EQ.1.AND.
     *              (ICW(M,IDXREV).EQ.1)) THEN

C...........Transfer control points to local coordinate.................
               XLOC=ZERO
               YLOC=ZERO
               ZLOC=ZERO
               DO 120 K=1,3
                  XLOC=XLOC+(XCTP_IM(I,K,1)-XCTW(J,K))*DIRW(J,1,K)
                  YLOC=YLOC+(XCTP_IM(I,K,1)-XCTW(J,K))*DIRW(J,2,K)
                  ZLOC=ZLOC+(XCTP_IM(I,K,1)-XCTW(J,K))*DIRW(J,3,K)
 120           CONTINUE
               
C...........Compute the induced potentials due to the blades............
C...........modified 04-11-90...........................................

               IMR=IMR0
               CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,FSY,
     *              FDX,FDY,FDZ,0,IMR)
               
C...........Near field use Morino's formulation.........................
               IF(IMR.EQ.2) THEN
                  DO 122 IXYZ=1,3
                     XMC(IXYZ)=XCTP_IM(I,IXYZ,1)
 122              CONTINUE
                  CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
               END IF
               
               IF(ABS(FD).GT.6.28) THEN
C                  WRITE(*,*) 'INFWAKS_IM-A',NTSTEP,I,M,N,IMR,FD
                  FD=0.0
               END IF
               
               WSUBIF_IM(I,M)=WSUBIF_IM(I,M)
     *              +(0.5-(FLOAT(N)-0.5)/FLOAT(NWSUB))*FD

               END IF

 140        CONTINUE
 150     CONTINUE
 160  CONTINUE
      
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine INFWAKS>>>>>>>>>>>>>>>>>>>>>>>>
      END
