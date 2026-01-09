      SUBROUTINE INFCOF_IM
************************************************************************
*     INFCOF: INFluence COeFficients due to blade IMages               *
*     Compute influence coefficients due to blade images and add them  *
*     to the original influence coefficients.                          *
*                                                                      *
*     Date      Comments                                               *
*     --------- -----------------------------------------------------  *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)

C-----------------------------------------------------------------------
C     Initialization
C-----------------------------------------------------------------------
C.....FILE 141 (IUMA1)  -- dipole inf. functions for each blade image
C.....FILE 142 (IUMB1)  -- source inf. functions for each blade image
      IUMA1=141
      IUMB1=142
      REWIND IUMA1
      REWIND IUMB1      

      DO 10 I=1,NPANZ
         TEMP1(I)=ZERO
 10   CONTINUE      

C-----------------------------------------------------------------------
C     Compute influence coefficients A(I,J) and B(I)
C       A(I,J): The induced potential at panel I due to a dipole with
C               unit strength at panel J
C       B(I):  The induced potential at panel I due to souces of all 
C               the panels
C-----------------------------------------------------------------------
C    The way this section works is that the vertices (corners) of each 
C    panel are put into the variables XV, and YV.  Then the control
C    point is put into XLOC,YLOC,ZLOC.  Rpan is then called to evaluate
C    source and dipole influence of the panel at that point.        cm
C-----------------------------------------------------------------------
      DO 160 M=MR,1,-1
         DO 150 N=1,NC
            J=INDEXB(N,M)

C.........Transfer data to the common block /GEOM/
            DO 100 K=1,4
               XV(K)=XVP(J,K)
               YV(K)=YVP(J,K)
               SIDE(K)=SID(J,K)
 100        CONTINUE
            DO 110 K=1,15
               S(K)=SS(J,K)
 110        CONTINUE

C.........XM1,XM2,XM3,XM4 for Morino's formulation......................
C ........ERRMAX is used to check if the wake panel is a triang. panel..
            XM1(1)=XB(N,M)
            XM1(2)=YB(N,M)
            XM1(3)=ZB(N,M)
            XM2(1)=XB(N,M+1)
            XM2(2)=YB(N,M+1)
            XM2(3)=ZB(N,M+1)
            XM3(1)=XB(N+1,M+1)
            XM3(2)=YB(N+1,M+1)
            XM3(3)=ZB(N+1,M+1)
            XM4(1)=XB(N+1,M)
            XM4(2)=YB(N+1,M)
            XM4(3)=ZB(N+1,M)   
            IMR0=0 
            ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))+ABS(XM2(3)-XM1(3))
            ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))+ABS(XM3(3)-XM2(3))
            ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))+ABS(XM4(3)-XM3(3))
            ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))+ABS(XM1(3)-XM4(3))
            ERRMIN=AMIN1(ER1,ER2,ER3,ER4) 

C.........Detect if triangular panels or not
            IF(ERRMIN.LE.1.0E-6) THEN
               IMR0=1
            END IF

            DO 140 KK=1,NBLADE

               IREC1=NTPOS(KK)

               DO 130 I=1,NPANEL

C.............Only calculate influence coefficients for fully...........
C.............submerged panels..........................................
                 IF(ISUBM(I,IDXREV).EQ.1.AND.ISUBM(J,IREC1).EQ.1) THEN

C.............Transfer control points to local coordinate...............
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO 120 K=1,3
                     XLOC=XLOC+(XCTP_IM(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                     YLOC=YLOC+(XCTP_IM(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                     ZLOC=ZLOC+(XCTP_IM(I,K,KK)-XCT(J,K))*DIR(J,3,K)
 120              CONTINUE

C.............Compute the induced potentials due to the blades..........
C.............modified 04-11-90.........................................

                  IMR=IMR0

                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *                 FDX,FDY,FDZ,0,IMR)
                  
C.............Near field use Morino's formulation.......................
                  IF(IMR.EQ.2) THEN
                     DO 122 IXYZ=1,3
                        XMC(IXYZ)=XCTP_IM(I,IXYZ,KK)
 122                 CONTINUE
                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                  END IF

                  IF(ABS(FD).GT.6.28) THEN
                     WRITE(66,*) 'INFCOF_IM',NTSTEP,KK,I,M,N,IMR,FD
                     FD=0.0
                  END IF

                  ELSE
                  FD=ZERO
                  FS=ZERO
                  END IF

C.................adding the effect of the image panels to the original
C.................influence coefficients.
                  A(I)=FD
                  TEMP1(I)=FS

 130           CONTINUE

C..............storing the modified influence coefficients.
               CALL WRITE1(IUMA1,A,NPANEL)
               CALL WRITE1(IUMB1,TEMP1,NPANEL)

 140        CONTINUE
 150     CONTINUE
 160  CONTINUE

C-----------------------------------------------------------------------
C     Compute influence coefficients due to the hub
C-----------------------------------------------------------------------
      IF(IHUB.NE.0) THEN

        DO 260 N=1,NHBX
          DO 250 M=1,MHBT
            J=INDEXH(N,M)

C...........Transfer data to the common block /GEOM/...................
            DO 200 K=1,4
              XV(K)=XVP(J,K)
              YV(K)=YVP(J,K)
              SIDE(K)=SID(J,K)
  200       CONTINUE
            DO 210 K=1,15
              S(K)=SS(J,K)
  210       CONTINUE

C...........XM1,XM2,XM3,XM4 for Morino's formulation
C...........ERRMAX is used to check if the wake panel is a triang. panel
            XM1(1)=XH(N,M+1)
            XM1(2)=YH(N,M+1)
            XM1(3)=ZH(N,M+1)
            XM2(1)=XH(N,M)
            XM2(2)=YH(N,M)
            XM2(3)=ZH(N,M)
            XM3(1)=XH(N+1,M)
            XM3(2)=YH(N+1,M)
            XM3(3)=ZH(N+1,M)
            XM4(1)=XH(N+1,M+1)             
            XM4(2)=YH(N+1,M+1)             
            XM4(3)=ZH(N+1,M+1)             
            IMR0=0 
            ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))
     *           +ABS(XM2(3)-XM1(3))
            ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))
     *           +ABS(XM3(3)-XM2(3))
            ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))
     *           +ABS(XM4(3)-XM3(3))   
            ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))
     *           +ABS(XM1(3)-XM4(3))
            ERRMIN=AMIN1(ER1,ER2,ER3,ER4) 

C..........Detect if triangular panels or not
            IF(ERRMIN.LE.1.0E-6) THEN
               IMR0=1
            END IF

            DO 240 KK=1,NBLADE
               
              IREC1=NTPOS(KK)

              DO 230 I=1,NPANEL

C.............Only calculate influence coefficients for fully...........
C.............submerged panels..........................................
                IF(ISUBM(I,IDXREV).EQ.1.AND.ISUBM(J,IREC1).EQ.1) THEN

C...............Transfer control points to the local coordinate.........
                XLOC=0.
                YLOC=0.
                ZLOC=0.
                DO 220 K=1,3
                  XLOC=XLOC+(XCTP_IM(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                  YLOC=YLOC+(XCTP_IM(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                  ZLOC=ZLOC+(XCTP_IM(I,K,KK)-XCT(J,K))*DIR(J,3,K)
  220           CONTINUE

C...............Compute the induced potentials due to the hub...........
                IMR=IMR0

                CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *                    FDX,FDY,FDZ,0,IMR)

C..............Near field use Morino's formulation
                IF(IMR.EQ.2) THEN
                   DO 221 IXYZ=1,3
                      XMC(IXYZ)=XCTP_IM(I,IXYZ,KK)
 221               CONTINUE
                   CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                END IF

                IF(ABS(FD).GT.6.28) THEN
                   WRITE(66,*) 'INFCOF_IM',NTSTEP,KK,I,M,N,IMR,FD
                   FD=0.0
                END IF

                ELSE
                FD=ZERO
                FS=ZERO
                END IF

C..............adding the effect of the image panels to the original
C..............influence coefficients.
                A(I)=FD
                TEMP1(I)=FS

  230         CONTINUE

C............storing the modified influence coefficients.
              CALL WRITE1(IUMA1,A,NPANEL)
              CALL WRITE1(IUMB1,TEMP1,NPANEL)

  240       CONTINUE
  250     CONTINUE
  260   CONTINUE
      END IF

      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine INFCOF_IM>>>>>>>>>>>>>>>>>>>>>>
      END



