      SUBROUTINE CAVRHS_SP
************************************************************************
*     CAVity solution Right-Hand-Side                                  *
*     Adds the influence of                                            *
*        (i) The key blade wake                                        *
*        (ii) Other blades and wakes
*     to the right-hand-side of the cavity system of equations         *
*                                                                      *
*     Notes:  The different between this subroutine and cavrhs.f is    *
*             (1) we only consider the fully submerged panels.         *
*             (2) in the first time step, the strength of the dipoles  *
*                 in the same strip of the wake are identical.         *
*             (3) the effect of the image panels are included if IMG=1 *
*             (4) the effect of the other blades are considered only   *
*                 if NREV>1.                                           *
*                                                                      *
*     Author: Julie Young  11-24-99                                    *
************************************************************************
      USE CVRHS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'      

C      COMMON/CVRHS/RHS(NTZ)

C-----------------------------------------------------------------------
C     *****INFLUENCE OF THE KEY BLADE*****
C     Compute the influence of shedding vortices of the key blade
C-----------------------------------------------------------------------
      IREC=NTPOS(1)-1
      IF(IREC.LE.0)IREC=IREC+TWOPI/DELTAT

C....store last time step's shed vortices in wake to TEMP5..............
      NREAD=NWMINFW*MR
      CALL READ2(46,IREC,TEMP5,NREAD)

C....store key blade's wake influence coefficients to TEMP1.............
      REWIND 81
      REWIND 91
      IF(IMG.EQ.1) THEN
         REWIND 181
         REWIND 191
      END IF

      DO 10 L=1,NWMINFW
         DO 20 M=MR,1,-1

C..........IDX: t(k-1);  IDX1: t(k-2)...................................
            IDX=MR*(L-2)+M
            IDX1=MR*(L-1)+M
            CALL READ1(81,TEMP1,NPANEL)
            IF(IMG.EQ.1) THEN
               CALL READ1(181,TEMP2,NPANEL)
               DO I=1,NPANEL
                  TEMP1(I)=TEMP1(I)-TEMP2(I)
               END DO
            END IF

            IEQN=0
            
            IF(L.LE.MSW(M,IDXREV)) THEN

C.............RHS of Green's formula on the blade
               DO 30 MM=MR,1,-1
                  DO 40 NN=1,NC
                     I=INDEXB(NN,MM)
                     IF(ISUBM(I,IDXREV).EQ.1) THEN
                        IEQN=IEQN+1
                        IF(L.EQ.1) THEN
                           TERM1=(HALF*TEMP1(I)-WSUBIF(I,IDX1))*
     *                          TEMP5(IDX1)
                        ELSE
                           TERM1=TEMP1(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
                        ENDIF
                        RHS(IEQN)=RHS(IEQN)-TERM1
                     END IF
 40               CONTINUE
                  IEQN=IEQN+NNWC(MM)
 30            CONTINUE

C.............RHS of Green's formula on the hub
               IEQN=NBW
               IF(IHUB.NE.0) THEN
                  DO 530 NN=1,NHBX
                     DO 540 MM=1,MHBT
                        I=INDEXH(NN,MM)
                        IF(ISUBM(I,IDXREV).EQ.1) THEN
                           IEQN=IEQN+1                     
                           IF(L.EQ.1)THEN
                              TERM1=(TEMP1(I)*HALF-
     *                             WSUBIF(I,IDX1))*TEMP5(IDX1)
                           ELSE
                              TERM1=TEMP1(I)*HALF*(TEMP5(IDX)+
     *                             TEMP5(IDX1))
                           ENDIF
                           RHS(IEQN)=RHS(IEQN)-TERM1
                        END IF
 540                 CONTINUE
 530              CONTINUE
               ENDIF

            END IF

 20      CONTINUE
 10   CONTINUE

      DO 50 M=MR,1,-1
         DO 60 L=1,NWMINFW

C..........IDX: t(k-1);  IDX1: t(k-2)...................................
            IDX=MR*(L-2)+M
            IDX1=MR*(L-1)+M
            IF(L.LE.NSUB)THEN
               CALL READ1(91,TEMP4,NPWAKS)
               IF(IMG.EQ.1) THEN
                  CALL READ1(191,TEMP3,NPWAKS)
                  DO I=1,NPWAKS
                     TEMP4(I)=TEMP4(I)-TEMP3(I)
                  END DO
               END IF
            ELSE
               WRITE(*,*) 'Please set NSUB=NWMINFW in gwakescs.f and'
               WRITE(*,*) 're-run the code!'
               STOP
            END IF
            
C..........RHS of Green's formula on the wake..........................
            IEQN=0
            
            IF(L.LE.MSW(M,IDXREV)) THEN
               DO 70 MM=MR,1,-1
                  IEQN=IEQN+NPERM(MM)-NNWC(MM)
                  DO 80 NN=1,NNWC(MM)
                     IEQN=IEQN+1
                     I=(MR-MM)*NTRA+NN
                     IF(L.EQ.1)THEN
                        TERM1=(HALF*TEMP4(I))*TEMP5(IDX1)
                     ELSE
                        TERM1=TEMP4(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
                     ENDIF
                     RHS(IEQN)=RHS(IEQN)-TERM1
 80               CONTINUE
 70            CONTINUE
            END IF
 60      CONTINUE
 50   CONTINUE


      IF(NBLADE.EQ.1) GOTO 90

C-----------------------------------------------------------------------
C     *****INFLUENCE OF OTHER BLADES*****
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     1. Compute source influence from submerged panels of other blades
C-----------------------------------------------------------------------
      DO 100 KK=2,NBLADE
         IREC=NTPOS(KK)
         CALL READ2(47,IREC,STRGTH(1,KK),NPANEL)
 100  CONTINUE

C.....RHS of Green's formula on the blade...............................
      REWIND 42
      IF(IMG.EQ.1) REWIND 142
      DO 110 J=1,NPANEL
         CALL READ1(42,TEMP1,NPANEL)
         IF(IMG.EQ.1) CALL READ1(142,TEMP2,NPANEL)
         DO 120 KK=2,NBLADE
            CALL READ1(42,TEMP1,NPANEL)
            IF(IMG.EQ.1) THEN
               CALL READ1(142,TEMP2,NPANEL)
               DO I=1,NPANEL
                  TEMP1(I)=TEMP1(I)-TEMP2(I)
               END DO
            END IF

            IEQN=0
            IREC1=NTPOS(KK)
            IF(ISUBM(J,IREC1).EQ.1) THEN

C.............RHS of Green's formula on the blade
               DO 130 MM=MR,1,-1
                  DO 140 NN=1,NC
                     I=INDEXB(NN,MM)
                     IF(ISUBM(I,IDXREV).EQ.1) THEN
                        IEQN=IEQN+1
                        RHS(IEQN)=RHS(IEQN)+TEMP1(I)*STRGTH(J,KK)
                     END IF
 140              CONTINUE
                  IEQN=IEQN+NNWC(MM)
 130           CONTINUE

C.............RHS of Green's formula on the hub
               IEQN=NBW
               IF(IHUB.NE.0) THEN
                  DO NN=1,NHBX               
                     DO MM=1,MHBT
                        I=INDEXH(NN,MM)
                        IF(ISUBM(I,IDXREV).EQ.1) THEN
                           IEQN=IEQN+1
                           RHS(IEQN)=RHS(IEQN)+TEMP1(I)*STRGTH(J,KK)
                        END IF
                     END DO
                  END DO
               END IF

            END IF
 120     CONTINUE
 110  CONTINUE

C.... RHS of Green's formula in the wake................................
      REWIND 112
      IF(IMG.EQ.1) REWIND 212
      DO 150 J=1,NPANEL
         CALL READ1(112,TEMP4,NPWAKS)
         IF(IMG.EQ.1) CALL READ1(212,TEMP3,NPWAKS)
         DO 160 KK=2,NBLADE
            CALL READ1(112,TEMP4,NPWAKS)
            IF(IMG.EQ.1) THEN
               CALL READ1(212,TEMP3,NPWAKS)
               DO I=1,NPWAKS
                  TEMP4(I)=TEMP4(I)-TEMP3(I)
               END DO
            END IF

            IEQN=0
            IREC1=NTPOS(KK)
            IF(ISUBM(J,IREC1).EQ.1) THEN
               DO 170 MM=MR,1,-1
                  IEQN=IEQN+NPERM(MM)-NNWC(MM)
                  DO 180 NN=1,NNWC(MM)
                     IEQN=IEQN+1
                     I=(MR-MM)*NTRA+NN
                     RHS(IEQN)=RHS(IEQN)+TEMP4(I)*STRGTH(J,KK)
 180              CONTINUE
 170           CONTINUE
            END IF
 160     CONTINUE
 150  CONTINUE

C-----------------------------------------------------------------------
C     2. Compute dipole influence from submerged panels of other blades
C-----------------------------------------------------------------------
      DO 190 KK=2,NBLADE
         IREC=NTPOS(KK)
         CALL READ2(45,IREC,STRGTH(1,KK),NPANEL)
 190  CONTINUE

C.....RHS of Green's formula on the blade...............................
      REWIND 41
      IF(IMG.EQ.1) REWIND 141
      DO 200 J=1,NPANEL
         CALL READ1(41,TEMP1,NPANEL)
         IF(IMG.EQ.1) CALL READ1(141,TEMP2,NPANEL)
         DO 210 KK=2,NBLADE
            CALL READ1(41,TEMP1,NPANEL)
            IF(IMG.EQ.1) THEN
               CALL READ1(141,TEMP2,NPANEL)
               DO I=1,NPANEL
                  TEMP1(I)=TEMP1(I)-TEMP2(I)
               END DO
            END IF

            IEQN=0
            IREC1=NTPOS(KK)
            IF(ISUBM(J,IREC1).EQ.1) THEN

C.............RHS of Green's formula on the blade
               DO 220 MM=MR,1,-1
                  DO 230 NN=1,NC
                     I=INDEXB(NN,MM)
                     IF(ISUBM(I,IDXREV).EQ.1) THEN
                        IEQN=IEQN+1
                        RHS(IEQN)=RHS(IEQN)-TEMP1(I)*STRGTH(J,KK)
                     END IF
 230              CONTINUE
                  IEQN=IEQN+NNWC(MM)
 220           CONTINUE

C.............RHS of Green's formula on the hub
               IEQN=NBW
               IF(IHUB.NE.0) THEN
                  DO NN=1,NHBX               
                     DO MM=1,MHBT
                        I=INDEXH(NN,MM)
                        IF(ISUBM(I,IDXREV).EQ.1) THEN
                           IEQN=IEQN+1
                           RHS(IEQN)=RHS(IEQN)-TEMP1(I)*STRGTH(J,KK)
                        END IF
                     END DO
                  END DO
               END IF

            END IF
 210     CONTINUE
 200  CONTINUE

C.....RHS of Green's formula in the wake................................
      REWIND 111
      IF(IMG.EQ.1) REWIND 211
      DO 240 J=1,NPANEL
         CALL READ1(111,TEMP4,NPWAKS)
         IF(IMG.EQ.1) CALL READ1(211,TEMP3,NPWAKS)
         DO 250 KK=2,NBLADE
            CALL READ1(111,TEMP4,NPWAKS)
            IF(IMG.EQ.1) THEN
               CALL READ1(211,TEMP3,NPWAKS)
               DO I=1,NPWAKS
                  TEMP4(I)=TEMP4(I)-TEMP3(I)
               END DO
            END IF

            IEQN=0
            IREC1=NTPOS(KK)
            IF(ISUBM(J,IREC1).EQ.1) THEN
               DO 260 MM=MR,1,-1
                  IEQN=IEQN+NPERM(MM)-NNWC(MM)
                  DO 270 NN=1,NNWC(MM)
                     IEQN=IEQN+1
                     I=(MR-MM)*NTRA+NN
                     RHS(IEQN)=RHS(IEQN)-TEMP4(I)*STRGTH(J,KK)
 270              CONTINUE
 260           CONTINUE
            END IF
 250     CONTINUE
 240  CONTINUE

C-----------------------------------------------------------------------
C     3. Compute wake source influence from submerged panels of other 
C        blades
C-----------------------------------------------------------------------
      DO 280 KK=2,NBLADE
         IREC=NTPOS(KK)
         CALL READ2(48,IREC,STRGTH1(1,KK),NPWAKS)
 280  CONTINUE

C.....RHS of Green's formula on the blade...............................
      REWIND 110
      IF(IMG.EQ.1) REWIND 210 
      DO 290 M=MR,1,-1
         DO 300 L=1,NTRA
            J=(MR-M)*NTRA+L
            CALL READ1(110,TEMP1,NPANEL)
            IF(IMG.EQ.1) CALL READ1(210,TEMP2,NPANEL)
            DO 310 KK=2,NBLADE
               CALL READ1(110,TEMP1,NPANEL)
               IF(IMG.EQ.1) THEN
                  CALL READ1(210,TEMP2,NPANEL)
                  DO I=1,NPANEL
                     TEMP1(I)=TEMP1(I)-TEMP2(I)
                  END DO
               END IF

               IEQN=0
               IREC1=NTPOS(KK)
               IF(ICW(M,IREC1).EQ.1.AND.L.LE.
     *              (N1SUB+(MSW(M,IREC1)-1)*NWSUB1)) THEN

C................RHS of Green's formula on the blade
                  DO 320 MM=MR,1,-1
                     DO 330 NN=1,NC
                        I=INDEXB(NN,MM)
                        IF(ISUBM(I,IDXREV).EQ.1) THEN
                           IEQN=IEQN+1
                           RHS(IEQN)=RHS(IEQN)+TEMP1(I)*STRGTH1(J,KK)
                        END IF
 330                 CONTINUE
                     IEQN=IEQN+NNWC(MM)                     
 320              CONTINUE

C.............RHS of Green's formula on the hub
                  IEQN=NBW
                  IF(IHUB.NE.0) THEN
                     DO NN=1,NHBX               
                        DO MM=1,MHBT
                           I=INDEXH(NN,MM)
                           IF(ISUBM(I,IDXREV).EQ.1) THEN
                              IEQN=IEQN+1
                              RHS(IEQN)=RHS(IEQN)+TEMP1(I)*STRGTH(J,KK)
                           END IF
                        END DO
                     END DO
                  END IF

               END IF
 310        CONTINUE
 300     CONTINUE
 290  CONTINUE

C.....RHS of Green's formula in the wake................................
      REWIND 113
      IF(IMG.EQ.1) REWIND 213
      DO 340 M=MR,1,-1
         DO 350 L=1,NTRA
            J=(MR-M)*NTRA+L
            CALL READ1(113,TEMP4,NPWAKS)
            IF(IMG.EQ.1) CALL READ1(213,TEMP3,NPWAKS)
            DO 360 KK=2,NBLADE
               CALL READ1(113,TEMP4,NPWAKS)
               IF(IMG.EQ.1) THEN
                  CALL READ1(213,TEMP3,NPWAKS)
                  DO I=1,NPWAKS
                     TEMP4(I)=TEMP4(I)-TEMP3(I)
                  END DO
               END IF

               IEQN=0
               IREC1=NTPOS(KK)
               IF(ICW(M,IREC1).EQ.1.AND.L.LE.
     *              (N1SUB+(MSW(M,IREC1)-1)*NWSUB1)) THEN 
                  DO 370 MM=MR,1,-1
                     IEQN=IEQN+NPERM(MM)-NNWC(MM)
                     DO 380 NN=1,NNWC(MM)
                        IEQN=IEQN+1
                        I=(MR-MM)*NTRA+NN
                        RHS(IEQN)=RHS(IEQN)+TEMP4(I)*STRGTH1(J,KK)
 380                 CONTINUE
 370              CONTINUE
               END IF
 360        CONTINUE
 350     CONTINUE
 340  CONTINUE

C-----------------------------------------------------------------------
C     4. Compute wake dipole influence from submerged panels of other 
C        blades
C-----------------------------------------------------------------------
C.....RHS of Green's formula on the blade...............................
      DO 390 KK=2,NBLADE
         IO=80+KK
         REWIND IO
         IF(IMG.EQ.1) THEN
            IO1=180+KK
            REWIND IO1
         END IF

         IREC=NTPOS(KK)
         NREAD=NWMINFW*MR
         CALL READ2(46,IREC,TEMP5,NREAD)
         DO 400 L=1,NWMINFW
            DO 410 M=MR,1,-1
               IDX=MR*(L-1)+M
               IDX1=MR*L+M
               IF(L.EQ.NWMINFW) TEMP5(IDX1)=ZERO
               CALL READ1(IO,TEMP1,NPANEL)
               IF(IMG.EQ.1) THEN
                  CALL READ1(IO1,TEMP2,NPANEL)
                  DO I=1,NPANEL
                     TEMP1(I)=TEMP1(I)-TEMP2(I)
                  END DO
               END IF

               IEQN=0
               IF(L.LE.MSW(M,IREC)) THEN

C................RHS of Green's formula on the blade
                  DO 420 MM=MR,1,-1
                     DO 430 NN=1,NC
                        I=INDEXB(NN,MM)
                        IF(ISUBM(I,IDXREV).EQ.1) THEN
                           IEQN=IEQN+1
                           RHS(IEQN)=RHS(IEQN)-
     *                          TEMP1(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
                        END IF
 430                 CONTINUE
                     IEQN=IEQN+NNWC(MM)
 420              CONTINUE

C................RHS of Green's formula on the hub
                  IEQN=NBW
                  IF(IHUB.NE.0) THEN
                     DO NN=1,NHBX
                        DO MM=1,MHBT
                           I=INDEXH(NN,MM)
                           IF(ISUBM(I,IDXREV).EQ.1) THEN
                              IEQN=IEQN+1                     
                              RHS(IEQN)=RHS(IEQN)-TEMP1(I)*HALF*
     *                             (TEMP5(IDX)+TEMP5(IDX1))
                           END IF
                        END DO
                     END DO
                  ENDIF

               END IF
 410        CONTINUE
 400     CONTINUE
 390  CONTINUE

C.....RHS of Green's formula in the wake................................
      DO 440 KK=2,NBLADE
         IO=90+KK
         REWIND IO
         IF(IMG.EQ.1) THEN
            IO1=190+KK
            REWIND IO1
         END IF

         IREC=NTPOS(KK)
         NREAD=NWMINFW*MR
         CALL READ2(46,IREC,TEMP5,NREAD)
         DO 450 M=MR,1,-1
            DO 460 L=1,NWMINFW
               IDX=MR*(L-1)+M
               IDX1=MR*L+M
               IF(L.EQ.NWMINFW) TEMP5(IDX1)=ZERO
               IF(L.LE.NSUB)THEN
                  CALL READ1(IO,TEMP4,NPWAKS)
                  IF(IMG.EQ.1) THEN
                     CALL READ1(IO1,TEMP3,NPWAKS)
                     DO I=1,NPWAKS
                        TEMP4(I)=TEMP4(I)-TEMP3(I)
                     END DO
                  END IF
               ELSE
                WRITE(*,*) 'Please set NSUB=NWMINFW in gwakescs.f and'
                  WRITE(*,*) 're-run the code!'
                  STOP
               ENDIF
               IEQN=0
               IF(L.LE.MSW(M,IREC)) THEN
                  DO 470 MM=MR,1,-1
                     IEQN=IEQN+NPERM(MM)-NNWC(MM)
                     DO 480 NN=1,NNWC(MM)
                        IEQN=IEQN+1
                        I=(MR-MM)*NTRA+NN
                        RHS(IEQN)=RHS(IEQN)-
     *                       TEMP4(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
 480                 CONTINUE
 470              CONTINUE
               END IF
 460        CONTINUE
 450     CONTINUE
 440  CONTINUE

C-----------------------------------------------------------------------
C     End of other blade contribution to RHS
C-----------------------------------------------------------------------
 90   CONTINUE

      RETURN
      END


