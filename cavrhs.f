      SUBROUTINE CAVRHS
************************************************************************
*     CAVity solution Right-Hand-Side                                  *
*     Adds the influence of                                            *
*          (i) The other blades,hub,(duct, tunnel, or tip vortex)                                        *
*         (ii) The wakes of the other blades                           *
*        (iii) The key blade wake 
*         (iV) The wakes of duct                                     *
*     to the right-hand-side of the cavity system of equations         *
*  Date of last Revision         Revision                              *
*  ---------------------         --------                              *
*  JYHL102401                  Subroutine undergone major revision.    *
************************************************************************
      USE CVRHS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C      COMMON/CVRHS/RHS(NTZ)
      DIMENSION ISPE(2,MBZ),NSPT(MBZ)

      IF(ISTEADY.EQ.0) RETURN

C-----------------------------------------------------------------------
C     Define the panel number of the split panels on face and back
C-----------------------------------------------------------------------
      ISR=1
      IF(IFACE.EQ.2) ISR=2      
      DO 5 MM=1,MR
         ISPE(1,MM)=999
         ISPE(2,MM)=999
 5    CONTINUE

      DO 10 MM=MR,1,-1
         NSPT(MM)=0
         DO 20 II=1,ISR
            IF((IFACE.EQ.0).OR.(II.EQ.1.AND.IFACE.EQ.2)) THEN
               IDR=1
            ELSE
               IDR=2
            END IF
            IF(IDR.EQ.1) THEN
               IF(NSPP(MM,IDR).NE.0) ISPE(IDR,MM)=LCV(MM,IDR)-1
            ELSE IF(IDR.EQ.2) THEN
               IF(NSPP(MM,IDR).NE.0) ISPE(IDR,MM)=M0(MM,IDR)-1-
     *              JCV(MM,IDR)
            END IF
            NSPT(MM)=NSPT(MM)+NSPP(MM,IDR)
 20      CONTINUE
 10   CONTINUE
      
      IF(NBLADE.EQ.1)GOTO 480
C-----------------------------------------------------------------------
C     *****INFLUENCE OF OTHER BLADES*****
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     1. Compute source influence from the other blades
C-----------------------------------------------------------------------

C....store the source strengths of blade KK to STRGTH...................
      DO 30 KK=2,NBLADE
         IREC=NTPOS(KK)
         CALL READ2(47,IREC,STRGTH(1,KK),NPANEL)
 30   CONTINUE

      REWIND 42
      DO 40 J=1,NPANB+NPANH+NPAND+NPANTN
         CALL READ1(42,TEMP1,NPANEL)
         DO 50 KK=2,NBLADE
            CALL READ1(42,TEMP1,NPANEL)

C..........RHS of Green's formula on the blade..........................
            IEQN=0
            DO 60 MM=MR,1,-1
               DO 70 NN=1,NC-NSPT(MM)
                  IEQN=IEQN+1
                  I=INDEXB(NN,MM)
                  I1=I
                  IF(NSPT(MM).EQ.1) THEN
                     IF(NN.GE.ISPE(1,MM).OR.NN.GE.ISPE(2,MM)) I1=I1+1
                  ELSE IF(NSPT(MM).EQ.2) THEN
                     IF(NN.GE.ISPE(2,MM)) I1=I1+1
                     IF(NN+1.GE.ISPE(1,MM)) I1=I1+1
                  END IF
                  RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH(J,KK)
 70            CONTINUE
               IEQN=IEQN+NNWC(MM)
 60         CONTINUE

C..........RHS of Green's formula on the hub............................

            IEQN=NBW
            IF(IHUB.NE.0) THEN
               DO 80 NN=1,NHBX               
                  DO 90 MM=1,MHBT
                     IEQN=IEQN+1
                     I1=INDEXH(NN,MM)
                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH(J,KK)
 90               CONTINUE
 80            CONTINUE
            END IF

            IF(IDUCT.NE.0) THEN
               DO MM=1,MDUCT               
                  DO NN=1,NDUCT
                     IEQN=IEQN+1
                     I1=INDEXD(NN,MM)
                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH(J,KK)
                  ENDDO
               ENDDO
            END IF

            IF(ITUN .NE. 0) THEN
               DO NN=1,NAXT               
                  DO MM=1,MTUNEL
                     IEQN=IEQN+1
                     I1=INDEXTN(NN,MM)
                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH(J,KK)
                  ENDDO
               ENDDO
            END IF

cC..........RHS of Green's formula on the bulb & tip vortex
c
c          IF(IAN.EQ.2) THEN
c               DO NN=1,NTHX
c                  DO MM=1,MCVT
c                     IEQN=IEQN+1
c                     I1=INDEXT(NN,MM)
c                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH(J,KK)
c                  ENDDO
c               ENDDO
c
c               DO NN=1,NCVX
c                  DO MM=1,MCVT
c                     IEQN=IEQN+1
c                     I1=INDEXC(NN,MM)
c                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH(J,KK)
c                  ENDDO
c               ENDDO
c          ENDIF

 50      CONTINUE
 40   CONTINUE

C.... RHS of Green's formula in the wake................................

      REWIND 112
      DO 100 J=1,NPANEL
         CALL READ1(112,TEMP4,NPWAKS)
         DO 110 KK=2,NBLADE
            CALL READ1(112,TEMP4,NPWAKS)
            IEQN=0
            DO 120 MM=MR,1,-1
               IEQN=IEQN+NC-NSPT(MM)
               DO 130 NN=1,NNWC(MM)
                  IEQN=IEQN+1
                  I=(MR-MM)*NTRA+NN
                  RHS(IEQN)=RHS(IEQN)+TEMP4(I)*STRGTH(J,KK)
 130           CONTINUE
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE

C-----------------------------------------------------------------------
C     2. Compute dipole influence from the other blades
C-----------------------------------------------------------------------

      DO 140 KK=2,NBLADE
         IREC=NTPOS(KK)
         CALL READ2(45,IREC,STRGTH(1,KK),NPANEL)
 140  CONTINUE

C....store the dipole strengths of blade KK to STRGTH...................

      REWIND 41
      DO 150 J=1,NPANB + NPANH + NPAND + NPANTN
         CALL READ1(41,TEMP1,NPANEL)
         DO 160 KK=2,NBLADE
            CALL READ1(41,TEMP1,NPANEL)

C..........RHS of Green's formula on the blade..........................

            IEQN=0
            DO 170 MM=MR,1,-1
               DO 180 NN=1,NC-NSPT(MM)
                  IEQN=IEQN+1
                  I=INDEXB(NN,MM)
                  I1=I
                  IF(NSPT(MM).EQ.1) THEN
                     IF(NN.GE.ISPE(1,MM).OR.NN.GE.ISPE(2,MM)) I1=I1+1
                  ELSE IF(NSPT(MM).EQ.2) THEN
                     IF(NN.GE.ISPE(2,MM)) I1=I1+1
                     IF(NN+1.GE.ISPE(1,MM)) I1=I1+1
                  END IF
                  RHS(IEQN)=RHS(IEQN)-TEMP1(I1)*STRGTH(J,KK)
 180           CONTINUE
               IEQN=IEQN+NNWC(MM)
 170        CONTINUE
            
C..........RHS of Green's formula on the hub............................

            IEQN=NBW
            IF(IHUB.NE.0) THEN
               DO 190 NN=1,NHBX
                  DO 200 MM=1,MHBT
                     IEQN=IEQN+1
                     I1=INDEXH(NN,MM)
                     RHS(IEQN)=RHS(IEQN)-TEMP1(I1)*STRGTH(J,KK)
 200              CONTINUE
 190           CONTINUE
            END IF

C..........RHS of Green's formula on the Duct............................

            IF(IDUCT.NE.0) THEN
               DO MM=1,MDUCT
                  DO NN=1,NDUCT
                     IEQN=IEQN+1
                     I1=INDEXD(NN,MM)
                     RHS(IEQN)=RHS(IEQN)-TEMP1(I1)*STRGTH(J,KK)
                  ENDDO
               ENDDO
            END IF

C..........RHS of Green's formula on the Tunnel............................

            IF(ITUN .NE. 0) THEN
               DO NN=1,NAXT
                  DO MM=1,MTUNEL
                     IEQN=IEQN+1
                     I1=INDEXTN(NN,MM)
                     RHS(IEQN)=RHS(IEQN)-TEMP1(I1)*STRGTH(J,KK)
                  ENDDO
               ENDDO
            END IF

cC..........RHS of Green's formula on the bulb & tip vortex
c
c          IF(IAN.EQ.2) THEN
c               DO NN=1,NTHX
c                  DO MM=1,MCVT
c                     IEQN=IEQN+1
c                     I1=INDEXT(NN,MM)
c                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH(J,KK)
c                  ENDDO
c               ENDDO
c
c               DO NN=1,NCVX
c                  DO MM=1,MCVT
c                     IEQN=IEQN+1
c                     I1=INDEXC(NN,MM)
c                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH(J,KK)
c                  ENDDO
c               ENDDO
c          ENDIF

 160     CONTINUE
 150  CONTINUE

C....RHS of Green's formula in the wake.................................

      REWIND 111
      DO 210 J=1,NPANEL
         CALL READ1(111,TEMP4,NPWAKS)
         DO 220 KK=2,NBLADE
            CALL READ1(111,TEMP4,NPWAKS)
            IEQN=0
            DO 230 MM=MR,1,-1
               IEQN=IEQN+NC-NSPT(MM)
               DO 240 NN=1,NNWC(MM)
                  IEQN=IEQN+1
                  I=(MR-MM)*NTRA+NN
                  RHS(IEQN)=RHS(IEQN)-TEMP4(I)*STRGTH(J,KK)
 240           CONTINUE
 230        CONTINUE
 220     CONTINUE
 210  CONTINUE

C-----------------------------------------------------------------------
C     3. Compute wake source influence from the other blades
C-----------------------------------------------------------------------

      DO 250 KK=2,NBLADE
         IREC=NTPOS(KK)
         CALL READ2(48,IREC,STRGTH1(1,KK),NPWAKS)
 250  CONTINUE

C....store the wake source strengths of blade KK to STRGTH..............
      REWIND 110
      DO 260 J=1,NPWAKS
         CALL READ1(110,TEMP1,NPANEL)
         DO 270 KK=2,NBLADE
            CALL READ1(110,TEMP1,NPANEL)

C..........RHS of Green's formula on the blade..........................
            IEQN=0
            DO 280 MM=MR,1,-1
               DO 290 NN=1,NC-NSPT(MM)
                  IEQN=IEQN+1
                  I=INDEXB(NN,MM)
                  I1=I
                  IF(NSPT(MM).EQ.1) THEN
                     IF(NN.GE.ISPE(1,MM).OR.NN.GE.ISPE(2,MM)) I1=I1+1
                  ELSE IF(NSPT(MM).EQ.2) THEN
                     IF(NN.GE.ISPE(2,MM)) I1=I1+1
                     IF(NN+1.GE.ISPE(1,MM)) I1=I1+1
                  END IF
                  RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH1(J,KK)
 290           CONTINUE
               IEQN=IEQN+NNWC(MM)
 280        CONTINUE

C..........RHS of Green's formula on the hub............................
            IEQN=NBW
          IF(IHUB.NE.0) THEN
               DO 300 NN=1,NHBX
                  DO 310 MM=1,MHBT
                     IEQN=IEQN+1
                     I1=INDEXH(NN,MM)
                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH1(J,KK)
 310              CONTINUE
 300           CONTINUE
          ENDIF

          IF(IDUCT .NE. 0) THEN
               DO MM=1,MDUCT
                  DO NN=1,NDUCT
                     IEQN=IEQN+1
                     I1=INDEXD(NN,MM)
                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH1(J,KK)
                  ENDDO
               ENDDO
          ENDIF

          IF(ITUN .NE. 0) THEN
               DO NN=1,NAXT
                  DO MM=1,MTUNEL
                     IEQN=IEQN+1
                     I1=INDEXTN(NN,MM)
                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH1(J,KK)
                  ENDDO
               ENDDO
          ENDIF

c          IF(IAN.EQ.2) THEN
c               DO NN=1,NTHX
c                  DO MM=1,MCVT
c                     IEQN=IEQN+1
c                     I1=INDEXT(NN,MM)
c                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH1(J,KK)
c                  ENDDO
c               ENDDO
c               DO NN=1,NCVX
c                  DO MM=1,MCVT
c                     IEQN=IEQN+1
c                     I1=INDEXC(NN,MM)
c                     RHS(IEQN)=RHS(IEQN)+TEMP1(I1)*STRGTH1(J,KK)
c                  ENDDO
c               ENDDO
c          ENDIF

 270     CONTINUE
 260  CONTINUE

C....RHS of Green's formula in the wake.................................

      REWIND 113
      DO 320 J=1,NPWAKS
         CALL READ1(113,TEMP4,NPWAKS)
         DO 330 KK=2,NBLADE
            CALL READ1(113,TEMP4,NPWAKS)
            IEQN=0
            DO 340 MM=MR,1,-1
               IEQN=IEQN+NC-NSPT(MM)
               DO 350 NN=1,NNWC(MM)
                  IEQN=IEQN+1
                  I=(MR-MM)*NTRA+NN
                  RHS(IEQN)=RHS(IEQN)+TEMP4(I)*STRGTH1(J,KK)
 350           CONTINUE
 340        CONTINUE
 330     CONTINUE
 320  CONTINUE

C-----------------------------------------------------------------------
C     4. Compute wake dipole influence from the other blades
C-----------------------------------------------------------------------

C....store the Dphi(k) to TEMP5.........................................
C....store wake influence coefficients to TEMP1........................

      DO 360 KK=2,NBLADE
         IO=80+KK
         REWIND IO
         IREC=NTPOS(KK)
         NREAD=NWMINFW*MR
         CALL READ2(46,IREC,TEMP5,NREAD)
         DO 370 L=1,NWMINFW
            DO 380 M=MR,1,-1
               IDX=MR*(L-1)+M
               IDX1=MR*L+M
               IF(L.EQ.NWMINFW) THEN
                  TEMP5(IDX1)=ZERO
               END IF
               CALL READ1(IO,TEMP1,NPANEL)

C.............RHS of Green's formula on the blade.......................
               IEQN=0
               DO 390 MM=MR,1,-1
                  DO 400 NN=1,NC-NSPT(MM)
                     IEQN=IEQN+1
                     I=INDEXB(NN,MM)
                     I1=I
                     IF(NSPT(MM).EQ.1) THEN
                        IF(NN.GE.ISPE(1,MM).OR.NN.GE.ISPE(2,MM)) I1=I1+1
                     ELSE IF(NSPT(MM).EQ.2) THEN
                        IF(NN.GE.ISPE(2,MM)) I1=I1+1
                        IF(NN+1.GE.ISPE(1,MM)) I1=I1+1
                     END IF
                     RHS(IEQN)=RHS(IEQN)-
     *                    TEMP1(I1)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
 400              CONTINUE
                  IEQN=IEQN+NNWC(MM)
 390           CONTINUE

C.............RHS of Green's formula on the hub.........................

               IEQN=NBW
             IF(IHUB .NE. 0) THEN
               DO 410 NN=1,NHBX
                  DO 420 MM=1,MHBT
                     IEQN=IEQN+1
                     I1=INDEXH(NN,MM)
                     RHS(IEQN)=RHS(IEQN)-
     *                    TEMP1(I1)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
 420              CONTINUE
 410           CONTINUE
             ENDIF

             IF(IDUCT .NE. 0) THEN
                  DO MM=1,MDUCT
                     DO NN=1,NDUCT
                        IEQN=IEQN+1
                        I1=INDEXD(NN,MM)
                        RHS(IEQN)=RHS(IEQN)-
     *                       TEMP1(I1)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
                     ENDDO
                  ENDDO
             ENDIF

             IF(ITUN .NE. 0) THEN
                  DO NN=1,NAXT
                     DO MM=1,MTUNEL
                        IEQN=IEQN+1
                        I1=INDEXTN(NN,MM)
                        RHS(IEQN)=RHS(IEQN)-
     *                       TEMP1(I1)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
                     ENDDO
                  ENDDO
             ENDIF
               
c               IF(IAN.EQ.2) THEN
c                  DO NN=1,NTHX
c                     DO MM=1,MCVT
c                        IEQN=IEQN+1
c                        I1=INDEXT(NN,MM)
c                        RHS(IEQN)=RHS(IEQN)-
c     *                       TEMP1(I1)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
c                     ENDDO
c                  ENDDO
c                  DO NN=1,NCVX
c                     DO MM=1,MCVT
c                        IEQN=IEQN+1
c                        I1=INDEXC(NN,MM)
c                        RHS(IEQN)=RHS(IEQN)-
c     *                       TEMP1(I1)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
c                     ENDDO
c                  ENDDO
c             ENDIF
               
 380        CONTINUE
 370     CONTINUE
 360  CONTINUE

      DO 430 KK=2,NBLADE
         IO1=90+KK
         IO2=30+KK
         REWIND IO1
         REWIND IO2
         IREC=NTPOS(KK)
         NREAD=NWMINFW*MR
         CALL READ2(46,IREC,TEMP5,NREAD)
         DO 440 M=MR,1,-1
            DO 450 L=1,NWMINFW
               IDX=MR*(L-1)+M
               IDX1=MR*L+M
               IF(L.EQ.NWMINFW) THEN
                  TEMP5(IDX1)=ZERO
               END IF
               IF(L.LE.NSUB)THEN
                  CALL READ1(IO1,TEMP4,NPWAKS)
               ELSE
                  CALL READ1(IO2,TEMP4,NPWAKS)
               ENDIF
               IEQN=0
               DO 460 MM=MR,1,-1
                  IEQN=IEQN+NC-NSPT(MM)
                  DO 470 NN=1,NNWC(MM)
                     IEQN=IEQN+1
                     I=(MR-MM)*NTRA+NN
                     RHS(IEQN)=RHS(IEQN)-
     *                    TEMP4(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
 470              CONTINUE
 460           CONTINUE
 450        CONTINUE
 440     CONTINUE
 430  CONTINUE


      IF(IDUCT .NE. 0) THEN
C-----------------------------------------------------------------------
C     5. Compute wake dipole influence from duct
C-----------------------------------------------------------------------

C....store the Dphi(k) to TEMP5.........................................
C....store wake influence coefficients to TEMP1........................

         DO 1360 KK=2,NBLADE
            IO=500+KK
            REWIND IO
            IREC=NTPOS(KK)
            CALL READ2(49,IREC,TEMPD4,NPWAKED)
            DO 1370 L=1, NDWK
               DO 1380 M=1, MDUCT
                  IDX=MDUCT*(L-1)+M
                  IDX1=MDUCT*L+M
                  IF(L.EQ.NDWK) THEN
                     TEMP5(IDX1)=ZERO
                  END IF
                  CALL READ1(IO,TEMP1,NPANEL)

C.............RHS of Green's formula on the blade.......................
                  IEQN=0
                  DO 1390 MM=MR,1,-1
                     DO 1400 NN=1,NC-NSPT(MM)
                        IEQN=IEQN+1
                        I=INDEXB(NN,MM)
                        I1=I
                        IF(NSPT(MM).EQ.1) THEN
                           IF(NN.GE.ISPE(1,MM).OR. 
     %                             NN.GE.ISPE(2,MM)) I1=I1+1
                        ELSE IF(NSPT(MM).EQ.2) THEN
                           IF(NN.GE.ISPE(2,MM)) I1=I1+1
                           IF(NN+1.GE.ISPE(1,MM)) I1=I1+1
                        END IF
                        RHS(IEQN)=RHS(IEQN)-
     *                       TEMP1(I1)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
 1400                CONTINUE
                     IEQN=IEQN+NNWC(MM)
 1390             CONTINUE
                  
C.............RHS of Green's formula on the hub.........................

                  IEQN=NBW
                  IF(IHUB .NE. 0) THEN
                     DO 1410 NN=1,NHBX
                        DO 1420 MM=1,MHBT
                           IEQN=IEQN+1
                           I1=INDEXH(NN,MM)
                           RHS(IEQN)=RHS(IEQN)-
     *                       TEMP1(I1)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
 1420                   CONTINUE
 1410                CONTINUE
                  ENDIF
                  
                  IF(IDUCT .NE. 0) THEN
                     DO MM=1,MDUCT
                        DO NN=1,NDUCT
                           IEQN=IEQN+1
                           I1=INDEXD(NN,MM)
                           RHS(IEQN)=RHS(IEQN)-
     *                       TEMP1(I1)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
                        ENDDO
                     ENDDO
                  ENDIF

                  IF(ITUN .NE. 0) THEN
                     DO NN=1,NAXT
                        DO MM=1,MTUNEL
                           IEQN=IEQN+1
                           I1=INDEXTN(NN,MM)
                           RHS(IEQN)=RHS(IEQN)-
     *                       TEMP1(I1)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
                        ENDDO
                     ENDDO
                  ENDIF
               
c                  IF(IAN.EQ.2) THEN
c                     DO NN=1,NTHX
c                        DO MM=1,MCVT
c                           IEQN=IEQN+1
c                           I1=INDEXT(NN,MM)
c                           RHS(IEQN)=RHS(IEQN)-
c     *                       TEMP1(I1)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
c                        ENDDO
c                     ENDDO
c                     DO NN=1,NCVX
c                        DO MM=1,MCVT
c                           IEQN=IEQN+1
c                           I1=INDEXC(NN,MM)
c                           RHS(IEQN)=RHS(IEQN)-
c     *                       TEMP1(I1)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
c                        ENDDO
c                     ENDDO
c                  ENDIF
                  
 1380          CONTINUE
 1370       CONTINUE
 1360    CONTINUE
         
      ENDIF
      
C-----------------------------------------------------------------------
C     End of other blade contribution to RHS
C-----------------------------------------------------------------------
 480  CONTINUE

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
      REWIND 31

      DO 490 L=1,NWMINFW
         DO 500 M=MR,1,-1

C.........IDX: t(k-1);  IDX1: t(k-2)....................................
            IDX=MR*(L-2)+M
            IDX1=MR*(L-1)+M
            CALL READ1(81,TEMP1,NPANEL)

C.........RHS of Green's formula on the blade...........................
            IEQN=0
            DO 510 MM=MR,1,-1
               DO 520 NN=1,NC-NSPT(MM)
                  IEQN=IEQN+1
                  I=INDEXB(NN,MM)
                  I1=I
                  IF(NSPT(MM).EQ.1) THEN
                     IF(NN.GE.ISPE(1,MM).OR.NN.GE.ISPE(2,MM)) I1=I1+1
                  ELSE IF(NSPT(MM).EQ.2) THEN
                     IF(NN.GE.ISPE(2,MM)) I1=I1+1
                     IF(NN+1.GE.ISPE(1,MM)) I1=I1+1
                  END IF
                  IF(L.EQ.1)THEN
                     TERM1=(HALF*TEMP1(I1)-WSUBIF(I1,IDX1))*TEMP5(IDX1)
                  ELSE
                     TERM1=TEMP1(I1)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
                  ENDIF
                  RHS(IEQN)=RHS(IEQN)-TERM1
 520           CONTINUE
               IEQN=IEQN+NNWC(MM)
 510        CONTINUE

            IEQN=NBW

C..........RHS of Green's formula on the hub...........................       
          IF(IHUB .NE. 0) THEN
               DO 530 NN=1,NHBX
                  DO 540 MM=1,MHBT
                     IEQN=IEQN+1
                     I=INDEXH(NN,MM)
                     IF(L.EQ.1)THEN
                        TERM1=(TEMP1(I)*HALF-
     *                       WSUBIF(I,IDX1))*TEMP5(IDX1)
                     ELSE
                        TERM1=TEMP1(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
                     ENDIF
                     RHS(IEQN)=RHS(IEQN)-TERM1
 540              CONTINUE
 530           CONTINUE
          ENDIF

C..........RHS of Green's formula on the tunnel...........................       
          IF(IDUCT .NE. 0) THEN
               DO MM=1,MDUCT
                  DO NN=1,NDUCT
                     IEQN=IEQN+1
                     I=INDEXD(NN,MM)
                     IF(L.EQ.1)THEN
                        TERM1=(TEMP1(I)*HALF-
     *                       WSUBIF(I,IDX1))*TEMP5(IDX1)
                     ELSE
                        TERM1=TEMP1(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
                     ENDIF
                     RHS(IEQN)=RHS(IEQN)-TERM1
                  ENDDO
               ENDDO
          ENDIF

C..........RHS of Green's formula on the tunnel...........................       
          IF(ITUN .NE. 0) THEN
               DO NN=1,NAXT
                  DO MM=1,MTUNEL
                     IEQN=IEQN+1
                     I=INDEXTN(NN,MM)
                     IF(L.EQ.1)THEN
                        TERM1=(TEMP1(I)*HALF-
     *                       WSUBIF(I,IDX1))*TEMP5(IDX1)
                     ELSE
                        TERM1=TEMP1(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
                     ENDIF
                     RHS(IEQN)=RHS(IEQN)-TERM1
                  ENDDO
               ENDDO
          ENDIF

cC..........RHS of Green's formula on the bulb and tip vortex
c          IF(IAN.EQ.2) THEN
c               DO NN=1,NTHX
c                  DO MM=1,MCVT
c                     IEQN=IEQN+1
c                     I=INDEXT(NN,MM)
c                     IF(L.EQ.1)THEN
c                        TERM1=(TEMP1(I)*HALF-
c     *                       WSUBIF(I,IDX1))*TEMP5(IDX1)
c                     ELSE
c                        TERM1=TEMP1(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
c                     ENDIF
c                     RHS(IEQN)=RHS(IEQN)-TERM1
c                  ENDDO
c               ENDDO
c
c               DO NN=1,NCVX
c                  DO MM=1,MCVT
c                     IEQN=IEQN+1
c                     I=INDEXC(NN,MM)
c                     IF(L.EQ.1)THEN
c                        TERM1=(TEMP1(I)*HALF-
c     *                       WSUBIF(I,IDX1))*TEMP5(IDX1)
c                     ELSE
c                        TERM1=TEMP1(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
c                     ENDIF
c                     RHS(IEQN)=RHS(IEQN)-TERM1
c                  ENDDO
c               ENDDO
c          ENDIF

 500     CONTINUE
 490  CONTINUE

C....store Dphi(k) to TEMP5.............................................
C....store wake influence coefficients to TEMP4.........................

C....RHS of Green's formula in the wake.................................
      DO 550 M=MR,1,-1
         DO 560 L=1,NWMINFW
C..........IDX: t(k-1);  IDX1: t(k-2)...................................
            IDX=MR*(L-2)+M
            IDX1=MR*(L-1)+M
            IF(L.LE.NSUB)THEN
               CALL READ1(91,TEMP4,NPWAKS)
            ELSE
               CALL READ1(31,TEMP4,NPWAKS)
            ENDIF
            IEQN=0
            DO 570 MM=MR,1,-1
               IEQN=IEQN+NC-NSPT(MM)
               DO 580 NN=1,NNWC(MM)
                  IEQN=IEQN+1
                  I=(MR-MM)*NTRA+NN

                  IF(L.LE.NSUB.AND.NWDIR(MM).EQ.2) THEN
                     TEMP4(I)=WKFACE(I,M,L)
                  END IF

                  IF(L.EQ.1)THEN
                     TERM1=(HALF*TEMP4(I))*TEMP5(IDX1)
                  ELSE
                     TERM1=TEMP4(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
                  ENDIF

                  RHS(IEQN)=RHS(IEQN)-TERM1
 580           CONTINUE
 570        CONTINUE
 560     CONTINUE
 550  CONTINUE

C-----------------------------------------------------------------------
C     Add the influence of the rest of the wake
C      (assume with steady fully wetted strength)
C-----------------------------------------------------------------------

C....RHS of Green's formula on the blade................................
      IEQN=0
      DO 590 MM=MR,1,-1
         DO 600 NN=1,NC-NSPT(MM)
            IEQN=IEQN+1
            I=INDEXB(NN,MM)
            I1=I
            IF(NSPT(MM).EQ.1) THEN
               IF(NN.GE.ISPE(1,MM).OR.NN.GE.ISPE(2,MM)) I1=I1+1
            ELSE IF(NSPT(MM).EQ.2) THEN
               IF(NN.GE.ISPE(2,MM)) I1=I1+1
               IF(NN+1.GE.ISPE(1,MM)) I1=I1+1
            END IF
            DO 610 M=MR,1,-1
               RHS(IEQN)=RHS(IEQN)-WSTINF(I1,M)*DPHI(M,0)
 610        CONTINUE
 600     CONTINUE
         IEQN=IEQN+NNWC(MM)
 590  CONTINUE

      IEQN=NBW
C....RHS of Green's formula on the hub..................................
      IF(IHUB.NE.0)THEN
         DO 620 NN=1,NHBX
            DO 630 MM=1,MHBT
               IEQN=IEQN+1
               I1=INDEXH(NN,MM)
               DO 640 M=MR,1,-1
                  RHS(IEQN)=RHS(IEQN)-WSTINF(I1,M)*DPHI(M,0)
 640           CONTINUE
 630        CONTINUE
 620     CONTINUE
      ENDIF

C....RHS of Green's formula on the DUCT..................................
      IF(IDUCT .NE. 0)THEN
         DO MM=1,MDUCT
            DO NN=1,NDUCT
               IEQN=IEQN+1
               I1=INDEXD(NN,MM)
               DO M=MR,1,-1
                  RHS(IEQN)=RHS(IEQN)-WSTINF(I1,M)*DPHI(M,0)
               ENDDO
            ENDDO
         ENDDO
      ENDIF


C....RHS of Green's formula on the TUNNEL..................................
      IF(ITUN .NE. 0)THEN
         DO NN=1,NAXT
            DO MM=1,MTUNEL
               IEQN=IEQN+1
               I1=INDEXTN(NN,MM)
               DO M=MR,1,-1
                  RHS(IEQN)=RHS(IEQN)-WSTINF(I1,M)*DPHI(M,0)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

cC....RHS of Green's formula on the bulb and tip vortex
c      IF(IAN.EQ.2) THEN
c         DO NN=1,NTHX
c            DO MM=1,MCVT
c               IEQN=IEQN+1
c               I1=INDEXT(NN,MM)
c               DO M=MR,1,-1
c                  RHS(IEQN)=RHS(IEQN)-WSTINF(I1,M)*DPHI(M,0)
c             ENDDO
c            ENDDO
c       ENDDO
c         DO NN=1,NCVX
c            DO MM=1,MCVT
c               IEQN=IEQN+1
c               I1=INDEXC(NN,MM)
c               DO M=MR,1,-1
c                  RHS(IEQN)=RHS(IEQN)-WSTINF(I1,M)*DPHI(M,0)
c             ENDDO
c            ENDDO
c       ENDDO
c      ENDIF
      
C....RHS of Green's formula on the wake.................................
      IEQN=0
      DO 650 MM=MR,1,-1
         IEQN=IEQN+NC-NSPT(MM)
         DO 660 NN=1,NNWC(MM)
            IEQN=IEQN+1
            I=(MR-MM)*NTRA+NN
            DO 670 M=MR,1,-1
               RHS(IEQN)=RHS(IEQN)-WUSINF1(I,M)*DPHI(M,0)
 670        CONTINUE
 660     CONTINUE
 650  CONTINUE



      IF(IDUCT .EQ. 0) GO TO 2100

C-----------------------------------------------------------------------
C     *****INFLUENCE OF THE KEY DUCT WAKE*****
C     Compute the influence of shedding vortices of the key blade
C-----------------------------------------------------------------------

      IREC=NTPOS(1)-1
      IF(IREC.LE.0)IREC=IREC+TWOPI/DELTAT

C....store last time step's shed vortices in wake to TEMP5..............

      CALL READ2(49,IREC,TEMPD4,NPWAKED)

C....store key blade's wake influence coefficients to TEMP1.............
      REWIND 501
      REWIND 521

      DO 1490 L=1,NDWK
         DO 1500 M=1 , MDUCT

C.........IDX: t(k-1);  IDX1: t(k-2)....................................
            IDX=MDUCT*(L-2)+M
            IDX1=MDUCT*(L-1)+M
            CALL READ1(501,TEMP1,NPANEL)

C.........RHS of Green's formula on the blade...........................
            IEQN=0
            DO 1510 MM=MR,1,-1
               DO 1520 NN=1,NC-NSPT(MM)
                  IEQN=IEQN+1
                  I=INDEXB(NN,MM)
                  I1=I
                  IF(NSPT(MM).EQ.1) THEN
                     IF(NN.GE.ISPE(1,MM).OR.NN.GE.ISPE(2,MM)) I1=I1+1
                  ELSE IF(NSPT(MM).EQ.2) THEN
                     IF(NN.GE.ISPE(2,MM)) I1=I1+1
                     IF(NN+1.GE.ISPE(1,MM)) I1=I1+1
                  END IF
                  IF(L.EQ.1)THEN
                     TERM1=(HALF*TEMP1(I1)-WSUBIFD(I1,IDX1))
     %                             *TEMPD4(IDX1)
                  ELSE
                     TERM1=TEMP1(I1)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
                  ENDIF
                  RHS(IEQN)=RHS(IEQN)-TERM1
 1520           CONTINUE
               IEQN=IEQN+NNWC(MM)
 1510        CONTINUE

            IEQN=NBW

C..........RHS of Green's formula on the hub...........................       
          IF(IHUB .NE. 0) THEN
               DO 1530 NN=1,NHBX
                  DO 1540 MM=1,MHBT
                     IEQN=IEQN+1
                     I=INDEXH(NN,MM)
                     IF(L.EQ.1)THEN
                        TERM1=(TEMP1(I)*HALF-
     *                       WSUBIFD(I,IDX1))*TEMPD4(IDX1)
                     ELSE
                        TERM1=TEMP1(I)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
                     ENDIF
                     RHS(IEQN)=RHS(IEQN)-TERM1
 1540              CONTINUE
 1530           CONTINUE
          ENDIF

C..........RHS of Green's formula on the tunnel...........................       
          IF(IDUCT .NE. 0) THEN
               DO MM=1,MDUCT
                  DO NN=1,NDUCT
                     IEQN=IEQN+1
                     I=INDEXD(NN,MM)
                     IF(L.EQ.1)THEN
                        TERM1=(TEMP1(I)*HALF-
     *                       WSUBIFD(I,IDX1))*TEMPD4(IDX1)
                     ELSE
                        TERM1=TEMP1(I)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
                     ENDIF
                     RHS(IEQN)=RHS(IEQN)-TERM1
                  ENDDO
               ENDDO
          ENDIF

C..........RHS of Green's formula on the tunnel...........................       
          IF(ITUN .NE. 0) THEN
               DO NN=1,NAXT
                  DO MM=1,MTUNEL
                     IEQN=IEQN+1
                     I=INDEXTN(NN,MM)
                     IF(L.EQ.1)THEN
                        TERM1=(TEMP1(I)*HALF-
     *                       WSUBIFD(I,IDX1))*TEMPD4(IDX1)
                     ELSE
                        TERM1=TEMP1(I)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
                     ENDIF
                     RHS(IEQN)=RHS(IEQN)-TERM1
                  ENDDO
               ENDDO
          ENDIF

cC..........RHS of Green's formula on the bulb and tip vortex
c          IF(IAN.EQ.2) THEN
c               DO NN=1,NTHX
c                  DO MM=1,MCVT
c                     IEQN=IEQN+1
c                     I=INDEXT(NN,MM)
c                     IF(L.EQ.1)THEN
c                        TERM1=(TEMP1(I)*HALF-
c     *                       WSUBIFD(I,IDX1))*TEMPD4(IDX1)
c                     ELSE
c                        TERM1=TEMP1(I)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
c                     ENDIF
c                     RHS(IEQN)=RHS(IEQN)-TERM1
c                  ENDDO
c               ENDDO
c
c               DO NN=1,NCVX
c                  DO MM=1,MCVT
c                     IEQN=IEQN+1
c                     I=INDEXC(NN,MM)
c                     IF(L.EQ.1)THEN
c                        TERM1=(TEMP1(I)*HALF-
c     *                       WSUBIFD(I,IDX1))*TEMPD4(IDX1)
c                     ELSE
c                        TERM1=TEMP1(I)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
c                     ENDIF
c                     RHS(IEQN)=RHS(IEQN)-TERM1
c                  ENDDO
c               ENDDO
c          ENDIF

 1500     CONTINUE
 1490  CONTINUE


       DO 1560 L = 1, NDWK
          DO 1550 M= 1 , MDUCT

             IDX=MDUCT*(L-2)+M
             IDX1=MDUCT*(L-1)+M
             CALL READ1(521,TEMP4,NPWAKS)
             IEQN=0
             DO 1570 MM=MR, 1 , -1
                IEQN=IEQN+NC-NSPT(MM)
                DO 1580 NN=1,NNWC(MM)
                   IEQN=IEQN+1
                   I=(MR-MM)*NTRA+NN

                  IF(L.EQ.1)THEN
                     TERM1=(HALF*TEMP4(I))*TEMPD4(IDX1)
                  ELSE
                     TERM1=TEMP4(I)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
                  ENDIF

                  RHS(IEQN)=RHS(IEQN)-TERM1
 1580           CONTINUE
 1570        CONTINUE
 1550  CONTINUE
 1560     CONTINUE

 2100  CONTINUE

      RETURN
      END 



