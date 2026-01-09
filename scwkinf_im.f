      SUBROUTINE SCWKINF_IM
************************************************************************
*                                                                      *
*   Subroutine SCWKINF computes                                        *
*      >the potential influence at the foil control points due to the  *
*       wake image source panels                                       *
*      >the potential influence at the wake control points due to the  *
*       foil image source & dipole panels                              *
*      >the potential influence at the wake control points due to the  *
*       wake image source and dipole panels                            *
*                                                                      *
*     Date      Comments                                               *
*     --------- -----------------------------------------------------  *
*     JY011200  Subroutine created.  Please note that the influence    *
*               due to hub image panels are not included in the        *
*               calculation.                                           *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3) 

C.....NWMINFW is the number of wake panels at each radius that..........
C.....represent the shed vorticity from the fully wetted solution.......
      NWMIN=NWMINFW

C-----------------------------------------------------------------------
C    define the supercavity wake (using subdivided panels)
C-----------------------------------------------------------------------
      CALL GWAKESCS

C.....NTRA is the number of subpanels on each wake strip................
C.....NWSUB is the number of panels inside each macro panel.............
C.....NPWAKS is the total number of subpanels in the supercavity wake...

      NTRA=NWSUB
      NWSUB=NWSUB1
      NPWAKS=NTRA*MR

      DO 120 M=MR,1,-1

C-----------------------------------------------------------------------
C
C        Prepare panel geometries of     (N,M+1) 2 *-----* 3 (N+1,M+1)
C          the wake, clockwise viewing             |     |
C          from the suction side of the            |     |
C          blade                           (N,M) 1 *-----* 4 (N+1,M)
C
C-----------------------------------------------------------------------
         DO 110 N=1,NTRA
            L=NTRA*(MR-M)+N
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
 110     CONTINUE
 120  CONTINUE

      CALL GEO3DW(NPWAKS,XGW,CHRLEWS,IER)

      IF(IER.EQ.0) THEN
        WRITE(*,*) 'UNACCEPTABLE PANEL IN SCWKINF_IM'
        STOP
      END IF

C-----------------------------------------------------------------------
C  C(I,J) = the induced potential at wing panel I due to a unit strength
C            source at wake panel J
C  D(I,J) =the induced  potential at wake panel I due to a unit strength
C            dipole at wing panel J
C  E(I,J) = the induced potential at wake panel I due to a unit strength
C            source at wing panel J
C  F(I,J) = the induced potential at wake panel I due to a unit strength
C            source at wake panel J
C WK(I,M) = the induced potential at wake panel I due to a unit strength
C            dipole at wake panel J
C-----------------------------------------------------------------------

C.....set scratch file unit numbers for housing the above ic's..........
      ISCWC=210
      ISCWD=211
      ISCWE=212
      ISCWF=213

      REWIND ISCWC
      REWIND ISCWD
      REWIND ISCWE
      REWIND ISCWF
      DO KK=1,NBLADE
         IO=190+KK
         REWIND IO
      END DO

C-----------------------------------------------------------------------
C     compute the induced potential at the foil cp's, due to the
C     wake source panels, C(I,J)
C  C(I,J) = the induced potential at wing panel I due to a unit strength
C            source at wake panel J
C-----------------------------------------------------------------------
      IFL=1

      DO 201 MM=MR,1,-1
      DO 200 LL=1,NTRA
         J=(MR-MM)*NTRA+LL

         DO 150 K=1,4
            XV(K)=XVPW(J,K)
            YV(K)=YVPW(J,K)
            SIDE(K)=SIDW(J,K)
 150     CONTINUE
         DO 160 K=1,15
            S(K)=SSW(J,K)
 160     CONTINUE
         DO 190 KK=1,NBLADE

            IREC1=NTPOS(KK)
            NTMP=N1SUB+(MSW(MM,IREC1)-1)*NWSUB1

            DO 180 I=1,NPANEL

               IF(ISUBM(I,IDXREV).EQ.1.AND.
     *              (ICW(MM,IREC1).EQ.1.AND.LL.LE.NTMP)) THEN

C...........transfer control points to local coordinate.................
               XLOC=ZERO
               YLOC=ZERO
               ZLOC=ZERO
               DO 170 K=1,3
                  XLOC=XLOC+(XCTP_IM(I,K,KK)-XCTW(J,K))*DIRW(J,1,K)
                  YLOC=YLOC+(XCTP_IM(I,K,KK)-XCTW(J,K))*DIRW(J,2,K)
                  ZLOC=ZLOC+(XCTP_IM(I,K,KK)-XCTW(J,K))*DIRW(J,3,K)
 170           CONTINUE

C...........compute the induced velocities due to the blade.............
               CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,FSY,FDX,
     *              FDY,FDZ,0,IFL)

               ELSE
               FS=ZERO
               ENDIF

               TEMP1(I)=FS
               
 180        CONTINUE

            CALL WRITE1(ISCWC,TEMP1,NPANEL)

 190     CONTINUE
 200  CONTINUE
 201  CONTINUE

C-----------------------------------------------------------------------
C     compute the induced potentials at the wake control points
C  D(I,J) =the induced  potential at wake panel I due to a unit strength
C            dipole at wing panel J
C  E(I,J) = the induced potential at wake panel I due to a unit strength
C            source at wing panel J
C-----------------------------------------------------------------------

C-----1a.due to the blade source and dipole panels----------------------
      DO 300 MM=MR,1,-1
         DO 290 NN=1,NC
            J=INDEXB(NN,MM)
            
C.........XM1,XM2,XM3,XM4 for Morino's formulation
C.........ERRMAX is used to check if the wake panel is a triang. panel
            XM1(1)=XB(NN,MM)
            XM1(2)=YB(NN,MM)
            XM1(3)=ZB(NN,MM)
            XM2(1)=XB(NN,MM+1)
            XM2(2)=YB(NN,MM+1)
            XM2(3)=ZB(NN,MM+1)
            XM3(1)=XB(NN+1,MM+1)
            XM3(2)=YB(NN+1,MM+1)
            XM3(3)=ZB(NN+1,MM+1)
            XM4(1)=XB(NN+1,MM)             
            XM4(2)=YB(NN+1,MM)             
            XM4(3)=ZB(NN+1,MM)             
            
            IMR0=0
            
            ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))+ABS(XM2(3)-XM1(3))
            ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))+ABS(XM3(3)-XM2(3))
            ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))+ABS(XM4(3)-XM3(3))
            ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))+ABS(XM1(3)-XM4(3))
            ERRMIN=AMIN1(ER1,ER2,ER3,ER4) 

C.........Detect if triangular panels or not............................
            IF(ERRMIN.LE.1.0E-6) THEN
               WRITE(*,2901) MM,NN
 2901          Format(' ...... triangular panels (m,n) -- > ',2I4)
               IMR0=1
            END IF

C.........transfer data to the common block /GEOM/ .....................
            DO 220 K=1,4
               XV(K)=XVP(J,K)
               YV(K)=YVP(J,K)
               SIDE(K)=SID(J,K)
 220        CONTINUE
            DO 230 K=1,15
               S(K)=SS(J,K)
 230        CONTINUE
            
C.........loop over the all of the blades...............................
            DO 280 KK=1,NBLADE

               IREC1=NTPOS(KK)
               
C...........loop over the wake control points for each blade wake.......
               DO 271 M=MR,1,-1
               DO 270 L=1,NTRA
                  I=(MR-M)*NTRA+L

                  NTMP=N1SUB+(MSW(M,IDXREV)-1)*NWSUB1
                  IF((ICW(M,IDXREV).EQ.1.AND.L.LE.NTMP).AND.
     *                 ISUBM(J,IREC1).EQ.1) THEN

C.............transfer control points to local coordinate...............
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO 240 K=1,3
                     XLOC=XLOC+(XCPW_IM(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                     YLOC=YLOC+(XCPW_IM(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                     ZLOC=ZLOC+(XCPW_IM(I,K,KK)-XCT(J,K))*DIR(J,3,K)
 240              CONTINUE
                  
C.............compute the induced velocities due to the foil............
                  IMR=IMR0
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *                 FDX,FDY,FDZ,0,IMR)
 
C.............Near field use Morino's formulation.......................
                  IF(IMR.EQ.2) THEN
                     DO 250 IXYZ=1,3
                        XMC(IXYZ)=XCPW_IM(I,IXYZ,KK)
 250                 CONTINUE
                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                  END IF

                  ELSE
                  FD=ZERO
                  FS=ZERO
                  ENDIF
                  
                  TEMP4(I)=FD
                  TEMP5(I)=FS

 270           CONTINUE
 271           CONTINUE
               
               CALL WRITE1(ISCWD,TEMP4,NPWAKS)
               CALL WRITE1(ISCWE,TEMP5,NPWAKS)

 280        CONTINUE
 290     CONTINUE
 300  CONTINUE
      
C-----1b.due to the hub source and dipole panels-----------------------
      IF(IHUB.NE.0)THEN
        DO 370 NN=1,NHBX
          DO 360 MM=1,MHBT
            J=INDEXH(NN,MM)

C...........transfer data to the common block /GEOM/ ...................
            DO 310 K=1,4
              XV(K)=XVP(J,K)
              YV(K)=YVP(J,K)
              SIDE(K)=SID(J,K)
  310       CONTINUE
            DO 320 K=1,15
              S(K)=SS(J,K)
  320       CONTINUE

C...........loop over the all of the blades.............................
            DO 350 KK=1,NBLADE

               IREC1=NTPOS(KK)

C.............loop over the wake control points for each blade wake.....
               DO 340 M=MR,1,-1
               DO 341 L=1,NTRA
                  I=(MR-M)*NTRA+L

                NTMP=N1SUB+(MSW(M,IDXREV)-1)*NWSUB1
                IF((ICW(M,IDXREV).EQ.1.AND.L.LE.NTMP).AND.
     *               ISUBM(J,IREC1).EQ.1) THEN

C...............transfer control points to local coordinate.............
                XLOC=ZERO
                YLOC=ZERO
                ZLOC=ZERO
                DO 330 K=1,3
                  XLOC=XLOC+(XCPW_IM(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                  YLOC=YLOC+(XCPW_IM(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                  ZLOC=ZLOC+(XCPW_IM(I,K,KK)-XCT(J,K))*DIR(J,3,K)
  330           CONTINUE

C...............compute the induced velocities due to the foil..........
                IMR=1
                CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,FDX,
     *                    FDY,FDZ,0,IMR)

                ELSE
                FD=ZERO
                FS=ZERO
                ENDIF

                TEMP4(I)=FD
                TEMP5(I)=FS

 341          CONTINUE
 340          CONTINUE

              CALL WRITE1(ISCWD,TEMP4,NPWAKS)
              CALL WRITE1(ISCWE,TEMP5,NPWAKS)

  350       CONTINUE
  360     CONTINUE
  370   CONTINUE
      ENDIF

C-----2.due to the wake source and dipole panels in the supercavity wake
C  F(I,J) = the induced potential at wake panel I due to a unit strength
C            source at wake panel J
C-----------------------------------------------------------------------

      DO 500 MM=MR,1,-1
         DO 490 JJ=1,NSUB

            DO 390 KK=1,NBLADE
               DO 380 I=1,NPWAKS
                  STRGTH1(I,KK)=ZERO
 380           CONTINUE
 390        CONTINUE

            IF(JJ.EQ.1) THEN
               NN1=N1SUB
            ELSE
               NN1=NWSUB1
            END IF
            
            DO 470 II=1,NN1
               IF(JJ.EQ.1) THEN
                  J=(MR-MM)*NTRA+II
               ELSE
                  J=(MR-MM)*NTRA+N1SUB+(JJ-2)*NWSUB1+II
               END IF

C...........XM1,XM2,XM3,XM4 for Morino's formulation....................
C...........ERRMAX is used to check if the wake panel is a triang. panel
               ERRMAX=-999.
               IMR0=0 
               DO 400 KI=1,3
                  XM1(KI)=XGW(J,1,KI)
                  XM2(KI)=XGW(J,2,KI)
                  XM3(KI)=XGW(J,3,KI)
                  XM4(KI)=XGW(J,4,KI)
                  ERR=ABS(XM1(KI)-XM4(KI))
                  ERRMAX=AMAX1(ERRMAX,ERR) 
 400           CONTINUE

C...........Detect if triangular panels or not
               IF(ERRMAX.LE.1.0E-6) THEN
                  IMR0=1
               END IF
               DO 410 K=1,4
                  XV(K)=XVPW(J,K)
                  YV(K)=YVPW(J,K)
                  SIDE(K)=SIDW(J,K)
 410           CONTINUE
               DO 420 K=1,15
                  S(K)=SSW(J,K)
 420           CONTINUE
               DO 460 KK=1,NBLADE
                  
                  IREC1=NTPOS(KK)

                  DO 451 M=MR,1,-1
                  DO 450 L=1,NTRA
                     I=(MR-M)*NTRA+L

                     NTMP=N1SUB+(MSW(M,IDXREV)-1)*NWSUB1
                     IF((ICW(M,IDXREV).EQ.1.AND.L.LE.NTMP).AND.
     *                    JJ.LE.MSW(MM,IREC1)) THEN

C...............transfer control points to local coordinate.............
                     XLOC=ZERO
                     YLOC=ZERO
                     ZLOC=ZERO
                     DO 430 K=1,3
                        XLOC=XLOC+(XCPW_IM(I,K,KK)-XCTW(J,K))*
     *                       DIRW(J,1,K)
                        YLOC=YLOC+(XCPW_IM(I,K,KK)-XCTW(J,K))*
     *                       DIRW(J,2,K)
                        ZLOC=ZLOC+(XCPW_IM(I,K,KK)-XCTW(J,K))*
     *                       DIRW(J,3,K)
 430                 CONTINUE

C...............compute the induced velocities due to the foil..........
                     IMR=IMR0
                     CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,
     *                    FSY,FDX,FDY,FDZ,0,IMR)

C...............Near field use Morino's formulation
                     IF(IMR.EQ.2) THEN
                        DO 440 IXYZ=1,3
                           XMC(IXYZ)=XCPW_IM(I,IXYZ,KK)
 440                    CONTINUE
                        CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                     END IF

                     ELSE
                     FD=ZERO
                     FS=ZERO
                     ENDIF

                     STRGTH1(I,KK)=STRGTH1(I,KK)+FD                     
                     TEMP4(I)=FS
                     
 450              CONTINUE
 451              CONTINUE

                  CALL WRITE1(ISCWF,TEMP4,NPWAKS)

 460           CONTINUE
 470        CONTINUE

            DO 480 KK=1,NBLADE
               CALL WRITE1(190+KK,STRGTH1(1,KK),NPWAKS)
 480        CONTINUE

 490     CONTINUE
 500  CONTINUE

      RETURN
C<<<<<<<<<<<<<<<<<<<<<<end of Subroutine SCWKINF_IM>>>>>>>>>>>>>>>>>>>>>
      END







