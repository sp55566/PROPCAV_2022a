      SUBROUTINE SCWKINF
************************************************************************
*                                                                      *
*   Subroutine SCWKINF computes                                        *
*      >the potential influence at the foil control points due to the  *
*       wake source panels                                             *
*      >the potential influence at the wake control points due to the  *
*       foil source & dipole panels                                    *
*      >the potential influence at the wake control points due to the  *
*       wake source and dipole panels                                  *
*                                                                      *
*   09-26-91 Neal Fine                                                 *
*                                                                      *
*   Date of last Revision                   Revision                   *
*   ---------------------            ----------------------            *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3) 
!s-- YE TIAN ---- 08/19/2013 ----
!     DIMENSION XVPWBIG(NSCWZ,4),YVPWBIG(NSCWZ,4),XCTWBIG(NSCWZ,3)
!    *     ,    DIRWBIG(NSCWZ,3,3),SSWBIG(NSCWZ,15),SIDWBIG(NSCWZ,4)
!    *     ,    XGWBIG(NSCWZ,4,3),CHRLEWSBIG(NSCWZ)
      real,allocatable ::XVPWBIG(:,:),YVPWBIG(:,:),XCTWBIG(:,:)
     &     ,    DIRWBIG(:,:,:),SSWBIG(:,:),SIDWBIG(:,:)
     &     ,    XGWBIG(:,:,:),CHRLEWSBIG(:)

      if (.NOT. allocated(XVPWBIG)) then
        allocate(XVPWBIG(NSCWZ,4),YVPWBIG(NSCWZ,4),XCTWBIG(NSCWZ,3)
     &     ,    DIRWBIG(NSCWZ,3,3),SSWBIG(NSCWZ,15),SIDWBIG(NSCWZ,4)
     &     ,    XGWBIG(NSCWZ,4,3),CHRLEWSBIG(NSCWZ))
      end if
!e-- YE TIAN ---- 08/19/2013 ----


      IF(IDUCT .EQ. 0) THEN
         CALL TRANS_GEO1(-1)

C/s S.N.KIM | When FWA is adopted, wake does not need to be re-generated.
C           | if 'GWAKE' is called, wake panel numbers will be messed up.
C           | Same for tip loaded propellers (IAN=5).
         IF(IAN.NE.6.AND.IREADW.NE.1) THEN
           IF(ITUN .EQ. 0) THEN
              CALL GWAKE
              IF(IHUB .EQ. 6) CALL PODGWAKE
           ELSE
              CALL GWAKE_TUN
           ENDIF
         ENDIF
C/e S.N.KIM | Aug. 2018.

         IF(ISP .NE. 0) THEN
            DO M = 1, MRP
               DO N = 1 , NSW(M)
                  CALL ROTATE2(0,SPANGLE,XW(N,M),YW(N,M))
               ENDDO
            ENDDO            
            
            DO M = 1 , 2
               DO N = 1 , MHBT+1
                  CALL ROTATE2(0,SPANGLE,XWDK(N,M),YWDK(N,M))
               ENDDO
            ENDDO            
         ENDIF
         
         CALL TRANS_GEO1(0)
      ENDIF


C.....NWMINFW is the number of wake panels at each radius that..........
C.....represent the shed vorticity from the fully wetted solution.......

      NWMINFW=NWMIN

C.....NWMIN is the number of wake panels at each radius that represent..
C.....the shed voriticity...............................................

      NWMIN=999
      NWBIG=0
      DO 10 M=1,MR
        IF(NSW(M).EQ.NSW(M+1))THEN
          NWPAN(M)=NSW(M)-1
        ELSE
          NWPAN(M)=MIN(NSW(M),NSW(M+1))
          NWPAN(M)=NWPAN(M)-1
        ENDIF
        NWMIN=MIN(NWMIN,NSW(M))
        NWBIG=NWBIG+NWPAN(M)
   10 CONTINUE
      NWMIN=NWMIN-1

      DO 30 M=MR,1,-1
C-----------------------------------------------------------------------
C
C        Prepare panel geometries of     (N,M+1) 2 *-----* 3 (N+1,M+1)
C          the wake, clockwise viewing             |     |
C          from the suction side of the            |     |
C          blade                           (N,M) 1 *-----* 4 (N+1,M)
C
C-----------------------------------------------------------------------
        DO 20 N=1,NWPAN(M)
          L=INDEXW(N,M)
          XGWBIG(L,1,1)=XW(N,M)
          XGWBIG(L,1,2)=YW(N,M)
          XGWBIG(L,1,3)=ZW(N,M)
          XGWBIG(L,2,1)=XW(N,M+1)
          XGWBIG(L,2,2)=YW(N,M+1)
          XGWBIG(L,2,3)=ZW(N,M+1)
          XGWBIG(L,3,1)=XW(N+1,M+1)
          XGWBIG(L,3,2)=YW(N+1,M+1)
          XGWBIG(L,3,3)=ZW(N+1,M+1)
          XGWBIG(L,4,1)=XW(N+1,M)
          XGWBIG(L,4,2)=YW(N+1,M)
          XGWBIG(L,4,3)=ZW(N+1,M)
   20    CONTINUE
   30 CONTINUE

      CALL GEO3DW(NWBIG,XGWBIG,CHRLEWSBIG,IER)
      IF(IER.EQ.0) THEN
        WRITE(*,*) 'UNACCEPTABLE MACRO-PANEL IN SCWKINF'
        STOP
      END IF

C******  Geom info on the wake surface   by Hong Sun 01/20/2005  ******
      DO 35 J = 1, NWBIG

        RCPW = SQRT(XCTW(J,2)**2+XCTW(J,3)**2)
        THPW = ATAN2(XCTW(J,3),XCTW(J,2))
        XCTPW(J,1,1)=XCTW(J,1)
        XCTPW(J,2,1)=XCTW(J,2)
        XCTPW(J,3,1)=XCTW(J,3)

        DO K = 1, 3
         DO K1 = 1, 3
           DIRWCP(J,K1,K,1) = DIRW(J,K1,K)     
         ENDDO
        ENDDO

      IF(IHULL .EQ. 1 .and. ICON .EQ. 5) THEN
         XCTPW(J,1,2) = XCTW(J,1)
         XCTPW(J,2,2) = 2. - XCTW(J,2)
         XCTPW(J,3,2) = XCTW(J,3)
      ELSE
          IF(NBLADE.GT.1) THEN
           DO 25 KK=2,NBLADE
            XCTPW(J,1,KK)=XCTW(J,1)
            XCTPW(J,2,KK)=RCPW*COS(THPW+DELK*(KK-1))
            XCTPW(J,3,KK)=RCPW*SIN(THPW+DELK*(KK-1))
 25        CONTINUE    
          ENDIF 
        ENDIF

 35   CONTINUE


C-----------------------------------------------------------------------
C     store macro wake quantities in new variables
C-----------------------------------------------------------------------
      DO 80 M=MR,1,-1

C.......The 70 loop should be from 1 to NWPAN(M). (JY030400)
         DO 70 N=1,NWPAN(M)
          I=INDEXW(N,M)
          DO 50 J=1,4
            XVPWBIG(I,J)=XVPW(I,J)
            YVPWBIG(I,J)=YVPW(I,J)
            SIDWBIG(I,J)=SIDW(I,J)
            IF(J.LE.3)THEN
              XCTWBIG(I,J)=XCTW(I,J)
              DO 40 K=1,3
                DIRWBIG(I,J,K)=DIRW(I,J,K)
   40         CONTINUE
            ENDIF
   50     CONTINUE
          DO 60 J=1,15
            SSWBIG(I,J)=SSW(I,J)
   60     CONTINUE
   70   CONTINUE
   80 CONTINUE

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
  110    CONTINUE
  120 CONTINUE

      CALL GEO3DW(NPWAKS,XGW,CHRLEWS,IER)

      IF(IER.EQ.0) THEN
        WRITE(*,*) 'UNACCEPTABLE PANEL IN SCWKINF'
        STOP
      END IF



C******  Geom info on the wake surface   by Hong Sun 04/06/2005  ******
      DO 122 J = 1, NPWAKS

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
           DO 121 KK=2,NBLADE
            XCTPWs(J,1,KK)=XCTW(J,1)
            XCTPWs(J,2,KK)=RCPW*COS(THPW+DELK*(KK-1))
            XCTPWs(J,3,KK)=RCPW*SIN(THPW+DELK*(KK-1))
 121       CONTINUE    
          ENDIF 
        ENDIF

 122    CONTINUE


C-----put wake control points in array containing other blades----------

      DO 130 J=1,NPWAKS
         IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XCTW(J,1),XCTW(J,2))

         XCPW(J,1,1)=XCTW(J,1)
         XCPW(J,2,1)=XCTW(J,2)
         XCPW(J,3,1)=XCTW(J,3)
         RCP=SQRT(XCTW(J,2)*XCTW(J,2)+XCTW(J,3)*XCTW(J,3))
         THP=ATAN2(XCTW(J,3),XCTW(J,2))

C.......Calculate the control points on the other blades................
         DO 125 KK=2,NBLADE
            XCPW(J,1,KK)=XCTW(J,1)
            XCPW(J,2,KK)=RCP*COS(THP+DELK*(KK-1))
            XCPW(J,3,KK)=RCP*SIN(THP+DELK*(KK-1))
 125     CONTINUE

         IF(ISP .NE. 0) THEN
            DO KK = 1, NBLADE
               CALL ROTATE2(0,SPANGLE,XCPW(J,1,KK),XCPW(J,2,KK))
            ENDDO
            CALL ROTATE2(0,SPANGLE,XCTW(J,1),XCTW(J,2))
         ENDIF
         
 130  CONTINUE

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

      ISCWC=110
      ISCWD=111
      ISCWE=112
      ISCWF=113

C-----------------------------------------------------------------------
C     compute the induced potential at the wing cp's, due to the
C     wake source panels, C(I,J)
C  C(I,J) = the induced potential at wing (wing,hub,tunnel,duct,DTV..) panel I due to a unit strength
C            source at wake panel J
C-----------------------------------------------------------------------

      IFL=1
      DO 200 J=1,NPWAKS
        DO 150 K=1,4
          XV(K)=XVPW(J,K)
          YV(K)=YVPW(J,K)
          SIDE(K)=SIDW(J,K)
  150   CONTINUE
        DO 160 K=1,15
          S(K)=SSW(J,K)
  160   CONTINUE
        DO 190 KK=1,NBLADE
          DO 180 I=1,NPANEL
C...........transfer control points to local coordinate.................
            XLOC=ZERO
            YLOC=ZERO
            ZLOC=ZERO
            DO 170 K=1,3
              XLOC=XLOC+(XCTP(I,K,KK)-XCTW(J,K))*DIRW(J,1,K)
              YLOC=YLOC+(XCTP(I,K,KK)-XCTW(J,K))*DIRW(J,2,K)
              ZLOC=ZLOC+(XCTP(I,K,KK)-XCTW(J,K))*DIRW(J,3,K)
  170       CONTINUE

C...........compute the induced velocities due to the blade.............

            CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,FSY,FDX,FDY,
     *                FDZ,0,IFL)

            TEMP1(I)=FS

  180     CONTINUE

          CALL WRITE1(ISCWC,TEMP1,NPANEL)

  190   CONTINUE
  200 CONTINUE
 
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
 2901       Format(' ...... triangular panels (m,n) -- > ',2I4)
            IMR0=1
          END IF

C.........transfer data to the common block /GEOM/ .....................
          DO 220 K=1,4
            XV(K)=XVP(J,K)
            YV(K)=YVP(J,K)
            SIDE(K)=SID(J,K)
  220     CONTINUE
          DO 230 K=1,15
            S(K)=SS(J,K)
  230     CONTINUE
C.........loop over the all of the blades...............................
          DO 280 KK=1,NBLADE
C...........loop over the wake control points for each blade wake.......
            DO 270 I=1,NPWAKS
C.............transfer control points to local coordinate...............
              XLOC=ZERO
              YLOC=ZERO
              ZLOC=ZERO
              DO 240 K=1,3
                XLOC=XLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                YLOC=YLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                ZLOC=ZLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,3,K)
  240         CONTINUE
C.............compute the induced velocities due to the foil............
              IMR=IMR0
              CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,FDX,
     *                  FDY,FDZ,0,IMR)
C.............Near field use Morino's formulation.......................
              IF(IMR.EQ.2) THEN
                DO 250 IXYZ=1,3
                  XMC(IXYZ)=XCPW(I,IXYZ,KK)
  250           CONTINUE
                CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
              END IF
              TEMP4(I)=FD
              TEMP5(I)=FS
  270       CONTINUE
            CALL WRITE1(ISCWD,TEMP4,NPWAKS)
            CALL WRITE1(ISCWE,TEMP5,NPWAKS)
  280     CONTINUE
  290   CONTINUE
  300 CONTINUE

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
C.............loop over the wake control points for each blade wake.....
              DO 340 I=1,NPWAKS
C...............transfer control points to local coordinate.............
                XLOC=ZERO
                YLOC=ZERO
                ZLOC=ZERO
                DO 330 K=1,3
                  XLOC=XLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                  YLOC=YLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                  ZLOC=ZLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,3,K)
  330           CONTINUE
C...............compute the induced velocities due to the foil..........
                IMR=1
                CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,FDX,
     *                    FDY,FDZ,0,IMR)
                TEMP4(I)=FD
                TEMP5(I)=FS
  340         CONTINUE
              CALL WRITE1(ISCWD,TEMP4,NPWAKS)
              CALL WRITE1(ISCWE,TEMP5,NPWAKS)
  350       CONTINUE
  360     CONTINUE
  370   CONTINUE
      ENDIF

C-----1c.due to the duct source and dipole panels-----------------------

      IF(IDUCT .NE. 0)THEN
         DO MM=1,MDUCT
            DO NN=1,NDUCT
               J=INDEXD(NN,MM)
               
C...........transfer data to the common block /GEOM/ ...................
               DO K=1,4
                  XV(K)=XVP(J,K)
                  YV(K)=YVP(J,K)
                  SIDE(K)=SID(J,K)
               ENDDO
               DO K=1,15
                  S(K)=SS(J,K)
               ENDDO
               DO KK=1,NBLADE
                  DO I=1,NPWAKS
                     XLOC=ZERO
                     YLOC=ZERO
                     ZLOC=ZERO
                     DO K=1,3
                        XLOC=XLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                        YLOC=YLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                        ZLOC=ZLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,3,K)
                     ENDDO
                     IMR=1
                     CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *                    FDX,FDY,FDZ,0,IMR)
                     TEMP4(I)=FD
                     TEMP5(I)=FS
                  ENDDO
                  CALL WRITE1(ISCWD,TEMP4,NPWAKS)
                  CALL WRITE1(ISCWE,TEMP5,NPWAKS)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

C-----1d.due to the tunnel source and dipole panels-----------------------

      IF(ITUN.NE.0)THEN
         DO NN=1,NAXT
            DO MM=1,MTUNEL
               J=INDEXTN(NN,MM)
               
C...........transfer data to the common block /GEOM/ ...................
               DO K=1,4
                  XV(K)=XVP(J,K)
                  YV(K)=YVP(J,K)
                  SIDE(K)=SID(J,K)
               ENDDO
               DO K=1,15
                  S(K)=SS(J,K)
               ENDDO
               DO KK=1,NBLADE
                  DO I=1,NPWAKS
                     XLOC=ZERO
                     YLOC=ZERO
                     ZLOC=ZERO
                     DO K=1,3
                        XLOC=XLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                        YLOC=YLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                        ZLOC=ZLOC+(XCPW(I,K,KK)-XCT(J,K))*DIR(J,3,K)
                     ENDDO
                     IMR=1
                     CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *                    FDX,FDY,FDZ,0,IMR)
                     TEMP4(I)=FD
                     TEMP5(I)=FS
                  ENDDO
                  CALL WRITE1(ISCWD,TEMP4,NPWAKS)
                  CALL WRITE1(ISCWE,TEMP5,NPWAKS)
               ENDDO
            ENDDO
         ENDDO
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
  380       CONTINUE
  390     CONTINUE

          DO I=1,NPWAKS
             WKFACE(I,MM,JJ)=ZERO
          END DO

C-----------------------------------------------------------------------
C     The index needs to be change here because of changes in 
C     gwakescs.f.                                               JY081499
C-----------------------------------------------------------------------

C.... MODIFIED FOR THE FIRST 4 WAKE PANELS FOR VISCOUS RUN  by HONG SUN
          IF(IVISC.EQ.1) THEN
C          IF(IVISC.EQ.1.OR.IREADW.EQ.1) THEN
            IF(JJ.LE.4) THEN
               NN1=N1SUB
            ELSE
               NN1=NWSUB1
            END IF
          ELSE 
            IF(JJ.EQ.1) THEN
               NN1=N1SUB
            ELSE
               NN1=NWSUB1
            END IF 
          END IF   ! IVISC 
 

          DO 470 II=1,NN1
             IF(IVISC.EQ.1) THEN
C             IF(IVISC.EQ.1.OR.IREADW.EQ.1) THEN
               IF(JJ.LE.4) THEN
                  J=(MR-MM)*NTRA+(JJ-1)*N1SUB+II
               ELSE
                  J=(MR-MM)*NTRA+N1SUB*4+(JJ-5)*NWSUB1+II
               END IF
             ELSE
               IF(JJ.EQ.1) THEN
                  J=(MR-MM)*NTRA+II
               ELSE
                  J=(MR-MM)*NTRA+N1SUB+(JJ-2)*NWSUB1+II
               END IF
             END IF       !  IVISC  
 
C--------End of changes (JY081499)--------------------------------------

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
  400       CONTINUE
C...........Detect if triangular panels or not
            IF(ERRMAX.LE.1.0E-6) THEN
              IMR0=1
            END IF
            DO 410 K=1,4
              XV(K)=XVPW(J,K)
              YV(K)=YVPW(J,K)
              SIDE(K)=SIDW(J,K)
  410       CONTINUE
            DO 420 K=1,15
              S(K)=SSW(J,K)
  420       CONTINUE
            DO 460 KK=1,NBLADE
              DO 450 I=1,NPWAKS
C...............transfer control points to local coordinate.............
                XLOC=ZERO
                YLOC=ZERO
                ZLOC=ZERO
                DO 430 K=1,3
                  XLOC=XLOC+(XCPW(I,K,KK)-XCTW(J,K))*DIRW(J,1,K)
                  YLOC=YLOC+(XCPW(I,K,KK)-XCTW(J,K))*DIRW(J,2,K)
                  ZLOC=ZLOC+(XCPW(I,K,KK)-XCTW(J,K))*DIRW(J,3,K)
  430           CONTINUE
C...............compute the induced velocities due to the foil..........
                IMR=IMR0

                CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,FSY,FDX,
     *                    FDY,FDZ,0,IMR)
C...............Near field use Morino's formulation
                IF(IMR.EQ.2) THEN
                   DO 440 IXYZ=1,3
                     XMC(IXYZ)=XCPW(I,IXYZ,KK)
  440              CONTINUE
                   CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                END IF

                IF((J.EQ.I).AND.(KK.EQ.1))FD=-TWOPI

                IF((J.EQ.I).AND.(KK.EQ.1)) THEN
                   WKFACE(I,MM,JJ)=WKFACE(I,MM,JJ)+TWOPI
                ELSEIF (KK.EQ.1) THEN
                   WKFACE(I,MM,JJ)=WKFACE(I,MM,JJ)+FD
                END IF

                STRGTH1(I,KK)=STRGTH1(I,KK)+FD
                TEMP4(I)=FS
  450         CONTINUE
              CALL WRITE1(ISCWF,TEMP4,NPWAKS)
  460       CONTINUE
  470     CONTINUE
          DO 480 KK=1,NBLADE
            CALL WRITE1(90+KK,STRGTH1(1,KK),NPWAKS)
  480     CONTINUE
  490   CONTINUE
  500 CONTINUE

C
C.....influence of the wake beyond the supercavity wake on the wake cp's
c
      DO 600 MM=MR,1,-1
        DO 590 JJ=NSUB+1,NWMIN
          J=INDEXW(JJ,MM)
C.........XM1,XM2,XM3,XM4 for Morino's formulation......................
C.........ERRMAX is used to check if the wake panel is a triang. panel..
          ERRMAX=-999.
          IMR0=0 
          DO 510 KI=1,3
            XM1(KI)=XGWBIG(J,1,KI)
            XM2(KI)=XGWBIG(J,2,KI)
            XM3(KI)=XGWBIG(J,3,KI)
            XM4(KI)=XGWBIG(J,4,KI)
            ERR=ABS(XM1(KI)-XM4(KI))
            ERRMAX=AMAX1(ERRMAX,ERR) 
  510     CONTINUE
C.........Detect if triangular panels or not
          IF(ERRMAX.LE.1.0E-6) THEN
            IMR0=1
          END IF
          DO 520 K=1,4
            XV(K)=XVPWBIG(J,K)
            YV(K)=YVPWBIG(J,K)
            SIDE(K)=SIDWBIG(J,K)
  520     CONTINUE
          DO 530 K=1,15
            S(K)=SSWBIG(J,K)
  530     CONTINUE
          DO 570 KK=1,NBLADE
            DO 560 I=1,NPWAKS
C.............transfer control points to local coordinate...............
              XLOC=ZERO
              YLOC=ZERO
              ZLOC=ZERO
              DO 540 K=1,3
                XLOC=XLOC+(XCPW(I,K,KK)-XCTWBIG(J,K))*DIRWBIG(J,1,K)
                YLOC=YLOC+(XCPW(I,K,KK)-XCTWBIG(J,K))*DIRWBIG(J,2,K)
                ZLOC=ZLOC+(XCPW(I,K,KK)-XCTWBIG(J,K))*DIRWBIG(J,3,K)
  540         CONTINUE
C.............compute the induced velocities due to the foil............
              IMR=IMR0
              CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWSBIG(J),FS,FD,FSX,FSY,FDX,
     *                  FDY,FDZ,0,IMR)
C.............Near field use Morino's formulation
              IF(IMR.EQ.2) THEN
                 DO 550 IXYZ=1,3
                    XMC(IXYZ)=XCPW(I,IXYZ,KK)
  550            CONTINUE
                 CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
              END IF
              STRGTH1(I,KK)=FD
  560       CONTINUE

            CALL WRITE1(30+KK,STRGTH1(1,KK),NPWAKS)
  570     CONTINUE
  590   CONTINUE
  600 CONTINUE

C.....influence of the wake "tail", beyond the shed vortices............
      DO 700 MM=MR,1,-1
        DO 605 I=1,NPWAKS
           WUSINF1(I,MM)=ZERO
  605   CONTINUE
        DO 690 JJ=NWMIN+1,NWPAN(MM)
          J=INDEXW(JJ,MM)
C.........XM1,XM2,XM3,XM4 for Morino's formulation......................
C.........ERRMAX is used to check if the wake panel is a triang. panel..
          ERRMAX=-999.
          IMR0=0 
          DO 610 KI=1,3
            XM1(KI)=XGWBIG(J,1,KI)
            XM2(KI)=XGWBIG(J,2,KI)
            XM3(KI)=XGWBIG(J,3,KI)
            XM4(KI)=XGWBIG(J,4,KI)
            ERR=ABS(XM1(KI)-XM4(KI))
            ERRMAX=AMAX1(ERRMAX,ERR) 
  610     CONTINUE
C.........Detect if triangular panels or not
          IF(ERRMAX.LE.1.0E-6) THEN
            IMR0=1
          END IF
          DO 620 K=1,4
            XV(K)=XVPWBIG(J,K)
            YV(K)=YVPWBIG(J,K)
            SIDE(K)=SIDWBIG(J,K)
  620     CONTINUE
          DO 630 K=1,15
            S(K)=SSWBIG(J,K)
  630     CONTINUE
          DO 670 KK=1,NBLADE
            DO 660 I=1,NPWAKS
C.............transfer control points to local coordinate...............
              XLOC=ZERO
              YLOC=ZERO
              ZLOC=ZERO
              DO 640 K=1,3
                XLOC=XLOC+(XCPW(I,K,KK)-XCTWBIG(J,K))*DIRWBIG(J,1,K)
                YLOC=YLOC+(XCPW(I,K,KK)-XCTWBIG(J,K))*DIRWBIG(J,2,K)
                ZLOC=ZLOC+(XCPW(I,K,KK)-XCTWBIG(J,K))*DIRWBIG(J,3,K)
  640         CONTINUE
C.............compute the induced velocities due to the foil............
              IMR=IMR0
              CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWSBIG(J),FS,FD,FSX,FSY,FDX,
     *                  FDY,FDZ,0,IMR)
C.............Near field use Morino's formulation
              IF(IMR.EQ.2) THEN
                 DO 650 IXYZ=1,3
                    XMC(IXYZ)=XCPW(I,IXYZ,KK)
  650            CONTINUE
                 CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
              END IF
              WUSINF1(I,MM)=WUSINF1(I,MM)+FD
  660       CONTINUE
  670     CONTINUE
  690   CONTINUE
  700 CONTINUE


C
C.....influence of the duct wake beyond the supercavity wake on the wake cp's
c

      IF(IDUCT .NE. 0) THEN
         DO JJ=1, NDWK 
            DO MM=1, MDUCT
               J=INDEXWD(JJ,MM)
               ERRMAX=-999.
               IMR0=0 
               DO KI=1,3
                  XM1(KI)=XGWD(J,1,KI)
                  XM2(KI)=XGWD(J,2,KI)
                  XM3(KI)=XGWD(J,3,KI)
                  XM4(KI)=XGWD(J,4,KI)
                  ERR=ABS(XM1(KI)-XM4(KI))
                  ERRMAX=AMAX1(ERRMAX,ERR) 
               ENDDO

               IF(ERRMAX.LE.1.0E-6) THEN
                  IMR0=1
               END IF

               DO K=1,4
                  XV(K)=XVPWD(J,K)
                  YV(K)=YVPWD(J,K)
                  SIDE(K)=SIDWD(J,K)
               ENDDO
               DO K=1,15
                  S(K)=SSWD(J,K)
               ENDDO

               DO KK=1,NBLADE
                  DO I=1,NPWAKS
                     XLOC=ZERO
                     YLOC=ZERO
                     ZLOC=ZERO
                     DO K=1,3
                        XLOC=XLOC+(XCPW(I,K,KK)-XCTWD(J,K))*DIRWD(J,1,K)
                        YLOC=YLOC+(XCPW(I,K,KK)-XCTWD(J,K))*DIRWD(J,2,K)
                        ZLOC=ZLOC+(XCPW(I,K,KK)-XCTWD(J,K))*DIRWD(J,3,K)
                     ENDDO

                     IMR=IMR0
                     CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWSD(J),FS,FD,FSX,
     %                         FSY,FDX,FDY,FDZ,0,IMR)

                     IF(IMR.EQ.2) THEN
                        DO IXYZ=1,3
                           XMC(IXYZ)=XCPW(I,IXYZ,KK)
                        ENDDO
                        CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                     END IF
                     STRGTH1(I,KK)=FD
                  ENDDO
               ENDDO
               
               DO KK =1 , NBLADE
                  CALL WRITE1(520+KK,STRGTH1(1,KK),NPWAKS)
               ENDDO
               
            ENDDO
         ENDDO
      ENDIF

C---------For foil case, we don't want the ultimate wake disk-CM011598--
C....The wake disk should not be used if ISP=1 (JY112499)...............
C....Same for ICON=8 as well (JY110100)

      CALL CLEAR(CHDK1,NPWAKS)
      CALL CLEAR(CHDKDT1,NPWAKS)
      CALL CLEAR(WDK1,NPWAKS)
      
      IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8.AND.ISP.NE.1)
     *     CALL INFDISKSC 
      
C-----------------------------------------------------------------------
C     Total influence coefficient at 0.7R
C-----------------------------------------------------------------------
      do i = 1 ,  mr
         if(rzp(i) .ge. 0.7) then
            mpdk = i
            go to 1100
         endif
      enddo
 1100 continue

C/s S.N.KIM | For the option IAN=5 or 6, ultimate wake disk (WDK1) is already
C           | ignored in INFDISKSC.f, so loop 730 can remain for other
C           | options.
        DO 730 I=1,NPWAKS
          WUSINF1(I,MPDK)=WUSINF1(I,MPDK)-(0.8/RULT/DELK/TANBUW)*WDK1(I)
  730   CONTINUE
C/e S.N.KIM | Aug. 2019.

C.....restore wake quantities...........................................
      DO 750 M=MR,1,-1
         DO 740 N=1,NTRA
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
 740     CONTINUE
 750  CONTINUE
      CALL GEO3DW(NPWAKS,XGW,CHRLEWS,IER)
      IF(IER.EQ.0) THEN
         WRITE(*,*) 'UNACCEPTABLE PANEL IN SCWKINF'
         STOP
      END IF

      RETURN
      END







