C =========================================================
      SUBROUTINE SCWKINF2
C =========================================================
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)

      NWMINFW=NWPANEL

      DO M = 1 , MR
         NWPAN(M) = NWPANEL
      ENDDO

      NWMIN = NWPANEL

      DO KK = NBLADE,1,-1

         CALL GWAKESCS2(KK)

C.....NTRA is the number of subpanels on each wake strip................
C.....NWSUB is the number of panels inside each macro panel.............
C.....NPWAKS is the total number of subpanels in the supercavity wake...

         NTRA=NWSUB
         NWSUB=NWSUB1
         NPWAKS=NTRA*MR

         DO 120 M=MR,1,-1
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
 110        CONTINUE
 120     CONTINUE

         CALL GEO3DWO(NPWAKS,XGW,CHRLEWSO,KK,IER,NSCWZ,KZ)

         IF(IER.EQ.0) THEN
            WRITE(*,*) 'UNACCEPTABLE PANEL IN SCWKINF'
            STOP
         END IF

      ENDDO
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

      ISCWC=110
      ISCWD=111
      ISCWE=112
      ISCWF=113

C-----------------------------------------------------------------------
C     compute the induced potential at the foil cp's, due to the
C     wake source panels, C(I,J)
C  C(I,J) = the induced potential at wing panel I due to a unit strength
C            source at wake panel J
C-----------------------------------------------------------------------

      IMR0 = 1
      DO 200 J=1,NPWAKS
         DO 190 KK=1,NBLADE
            DO K=1,4
               XV(K)=XVPWO(J,K,KK)
               YV(K)=YVPWO(J,K,KK)
               SIDE(K)=SIDWO(J,K,KK)
            ENDDO
            DO K=1,15
               S(K)=SSWO(J,K,KK)
            ENDDO
            DO 180 I=1,NPANEL
               XLOC=ZERO
               YLOC=ZERO
               ZLOC=ZERO
               DO K=1,3
                  XLOC=XLOC+(XCTP(I,K,KK)-XCTWO(J,K,KK))
     %                 *DIRWO(J,1,K,KK)
                  YLOC=YLOC+(XCTP(I,K,KK)-XCTWO(J,K,KK))
     %                 *DIRWO(J,2,K,KK)
                  ZLOC=ZLOC+(XCTP(I,K,KK)-XCTWO(J,K,KK))
     %                 *DIRWO(J,3,K,KK)
               ENDDO
C...........compute the induced velocities due to the blade.............
               IMR = IMR0
               CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWSO(J,KK),FS,FD,
     %              FSX,FSY,FDX,FDY,FDZ,0,IMR)

C/s S.N.KIM: Let's turn off the effect from the first three columns in wake
C            (1 degree blade angel with N1SUB=3*NDLTAT in gwakescs2.f) on blade
C            in case the unsteady wake alignment scheme is coupled with cavity runs.
C            This is to avoid the square root singulrity of S.C. source when wake panels
C            get extremly closer to the blade T.E. It happens as the wake panels are free 
C            to move based on the local flow with UWA applied.  
               IF(MOD(J,NTRA).EQ.1.and.KK.EQ.1) FS=0.0
               IF(MOD(J,NTRA).EQ.2.and.KK.EQ.1) FS=0.0
               IF(MOD(J,NTRA).EQ.3.and.KK.EQ.1) FS=0.0
C               IF(MOD(J,NTRA).EQ.4.and.KK.EQ.1) FS=0.0
C/e S.N.KIM

               TEMP1(I)=FS
 180        CONTINUE
            CALL WRITE1(ISCWC,TEMP1,NPANEL)
 190     CONTINUE
 200  CONTINUE

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

          ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))
     %         +ABS(XM2(3)-XM1(3))
          ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))
     %         +ABS(XM3(3)-XM2(3))
          ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))
     %         +ABS(XM4(3)-XM3(3))
          ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))
     %         +ABS(XM1(3)-XM4(3))
          ERRMIN=AMIN1(ER1,ER2,ER3,ER4) 

          IF(ERRMIN.LE.1.0E-6) THEN
            WRITE(*,2901) MM,NN
 2901       Format(' ...... triangular panels (m,n) -- > ',2I4)
            IMR0=1
          END IF

          DO 220 K=1,4
            XV(K)=XVP(J,K)
            YV(K)=YVP(J,K)
            SIDE(K)=SID(J,K)
  220     CONTINUE
          DO 230 K=1,15
            S(K)=SS(J,K)
  230     CONTINUE

          DO 280 KK=1,NBLADE
             
             DO 270 I=1,NPWAKS
                
                XLOC=ZERO
                YLOC=ZERO
                ZLOC=ZERO
                DO 240 K=1,3
                   XLOC=XLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                   YLOC=YLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                   ZLOC=ZLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,3,K)
 240            CONTINUE

                IMR=IMR0
                CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,
     *               FSX,FSY,FDX,FDY,FDZ,0,IMR)

                IF(IMR.EQ.2) THEN
                   DO 250 IXYZ=1,3
                      XMC(IXYZ)=XCTWO(I,IXYZ,KK)
 250               CONTINUE
                   CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                END IF
                TEMP4(I)=FD
                TEMP5(I)=FS
 270         CONTINUE
             CALL WRITE1(ISCWD,TEMP4,NPWAKS)
             CALL WRITE1(ISCWE,TEMP5,NPWAKS)
 280      CONTINUE
 290   CONTINUE
 300  CONTINUE
      IF(IHUB.NE.0)THEN
         DO 370 NN=1,NHBX
            DO 360 MM=1,MHBT
               J=INDEXH(NN,MM)
               
               DO 310 K=1,4
                  XV(K)=XVP(J,K)
                  YV(K)=YVP(J,K)
                  SIDE(K)=SID(J,K)
 310           CONTINUE
               DO 320 K=1,15
                  S(K)=SS(J,K)
 320           CONTINUE

               DO 350 KK=1,NBLADE

                  DO 340 I=1,NPWAKS

                     XLOC=ZERO
                     YLOC=ZERO
                     ZLOC=ZERO
                     DO 330 K=1,3
                        XLOC=XLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                        YLOC=YLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                        ZLOC=ZLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,3,K)
330                 CONTINUE

                     IMR=1
                     CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,
     *                    FSY,FDX,FDY,FDZ,0,IMR)
                     TEMP4(I)=FD
                     TEMP5(I)=FS
 340              CONTINUE
                  CALL WRITE1(ISCWD,TEMP4,NPWAKS)
                  CALL WRITE1(ISCWE,TEMP5,NPWAKS)
 350           CONTINUE
 360        CONTINUE
 370     CONTINUE
      ENDIF
C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
C-----1c.due to the Tip hub source and dipole panels-----------------------

c      DO NN=1,NTHX
c         DO MM=1,MCVT
c            J=INDEXT(NN,MM)
c            
c            DO K=1,4
c               XV(K)=XVP(J,K)
c               YV(K)=YVP(J,K)
c               SIDE(K)=SID(J,K)
c            ENDDO
c            DO K=1,15
c               S(K)=SS(J,K)
c            ENDDO
c            DO KK=1,NBLADE
c               IF(KK .GE. 2) THEN
c                  DO I = 1 , NPWAKS
c                     TEMP4(I) = 0.0
c                     TEMP5(I) = 0.0
c                  ENDDO
c               ELSEIF(KK .EQ. 1) THEN
c                  DO I=1,NPWAKS
c                     XLOC=ZERO
c                     YLOC=ZERO
c                     ZLOC=ZERO
c                     DO K=1,3
c                        XLOC=XLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,1,K)
c                        YLOC=YLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,2,K)
c                        ZLOC=ZLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,3,K)
c                     ENDDO
c                     IMR=1
c                     CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),
c     %                    FS,FD,FSX,FSY,FDX,FDY,FDZ,0,IMR)
c                     TEMP4(I)=FD
c                     TEMP5(I)=FS
c                  ENDDO
c               ENDIF
c               CALL WRITE1(ISCWD,TEMP4,NPWAKS)
c               CALL WRITE1(ISCWE,TEMP5,NPWAKS)
c          ENDDO
c         ENDDO
c      ENDDO
c
cC-----1d.due to the Tip vortex source and dipole panels-----------------------
c
c        DO NN=1,NCVX
c          DO MM=1,MCVT
c            J=INDEXC(NN,MM)
c
c            DO K=1,4
c              XV(K)=XVP(J,K)
c              YV(K)=YVP(J,K)
c              SIDE(K)=SID(J,K)
c            ENDDO
c            DO K=1,15
c              S(K)=SS(J,K)
c            ENDDO
c          DO KK = 1 , NBLADE
c
c               IF(KK .GE. 2) THEN
c                  DO I = 1 , NPWAKS
c                     TEMP4(I) = 0.0
c                     TEMP5(I) = 0.0
c                  ENDDO
c               ELSEIF(KK .EQ. 1) THEN
c                  DO I=1,NPWAKS
c                     XLOC=ZERO
c                     YLOC=ZERO
c                     ZLOC=ZERO
c                     DO K=1,3
c                        XLOC=XLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,1,K)
c                        YLOC=YLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,2,K)
c                        ZLOC=ZLOC+(XCTWO(I,K,KK)-XCT(J,K))*DIR(J,3,K)
c                     ENDDO
c                     IMR=1
c                     CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),
c     %                    FS,FD,FSX,FSY,FDX,FDY,FDZ,0,IMR)
c                     TEMP4(I)=FD
c                     TEMP5(I)=FS
c                  ENDDO
c               ENDIF
c
c               CALL WRITE1(ISCWD,TEMP4,NPWAKS)
c               CALL WRITE1(ISCWE,TEMP5,NPWAKS)
c          ENDDO
c         ENDDO
c      ENDDO
C/e S.N.KIM | Aug. 2018.
      
C-----2.due to the wake source and dipole panels in the supercavity wake
C  F(I,J) = the induced potential at wake panel I due to a unit strength
C            source at wake panel J
C-----------------------------------------------------------------------

      DO 500 MM=MR,1,-1
         DO 490 JJ=1,NSUB
            CALL CLEAR(STRGTH1,NSCWZ*KZ)
            DO I=1,NPWAKS
               WKFACE(I,MM,JJ)=ZERO
            END DO
            
C......MODIFIED for the first 4 wake panels for viscous run  by Hong Sun 
        
            IF(IVISC.EQ.1) THEN
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
            END IF  !  IVISC

            DO 470 II=1,NN1
          
               IF(IVISC.EQ.1) THEN
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
               END IF   !  IVISC
               
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
               
               IF(ERRMAX.LE.1.0E-6) THEN
                  IMR0=1
               END IF
               
               DO 460 KK=1,NBLADE
                  
                  DO 410 K=1,4
                     XV(K)=XVPWO(J,K,KK)
                     YV(K)=YVPWO(J,K,KK)
                     SIDE(K)=SIDWO(J,K,KK)
 410              CONTINUE
                  DO 420 K=1,15
                     S(K)=SSWO(J,K,KK)
 420              CONTINUE
                     
                  DO 450 I=1,NPWAKS
                     
                     XLOC=ZERO
                     YLOC=ZERO
                     ZLOC=ZERO
                     DO 430 K=1,3
                        XLOC=XLOC+(XCTWO(I,K,1)-XCTWO(J,K,KK))
     %                       *DIRWO(J,1,K,KK)
                        YLOC=YLOC+(XCTWO(I,K,1)-XCTWO(J,K,KK))
     %                       *DIRWO(J,2,K,KK)
                        ZLOC=ZLOC+(XCTWO(I,K,1)-XCTWO(J,K,KK))
     %                       *DIRWO(J,3,K,KK)
 430                 CONTINUE
                     
                     IMR=IMR0
                     CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWSO(J,KK)
     %                    ,FS,FD,FSX,FSY,FDX,FDY,FDZ,0,IMR)

                     IF(IMR.EQ.2) THEN
                        DO 440 IXYZ=1,3
                           XMC(IXYZ)=XCTWO(I,IXYZ,KK)
 440                    CONTINUE
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
                     
 450              CONTINUE
                  CALL WRITE1(ISCWF,TEMP4,NPWAKS)
 460           CONTINUE
 470        CONTINUE
            DO 480 KK=1,NBLADE
               CALL WRITE1(90+KK,STRGTH1(1,KK),NPWAKS)
 480        CONTINUE
 490     CONTINUE
 500  CONTINUE
      DO 700 MM=MR,1,-1
         DO 605 I=1,NPWAKS
            WUSINF1(I,MM)=ZERO
 605    CONTINUE
 700  CONTINUE
      
      CALL CLEAR(CHDK1,NPWAKS)
      CALL CLEAR(CHDKDT1,NPWAKS)
      CALL CLEAR(WDK1,NPWAKS)

      IF(IHUB .NE. 0) CALL INFDISKSC2

      DO I = 1 , MR
         IF(RZP(I) .GE. 0.7) THEN
            MPDK = I
            GO TO 1100
         ENDIF
      ENDDO
 1100 CONTINUE
C      DO 730 I=1,NPWAKS
C         WUSINF1(I,MPDK)=WUSINF1(I,MPDK)
C     %        -(0.8/RULT/DELK/TANBUW)*WDK1(I)
C  730 CONTINUE

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
      DO K = 1 , NBLADE
         DO I = 1 , NPWAKS
            DO J = 1 , 3
               XCPW(I,J,K) = XCTWO(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END





