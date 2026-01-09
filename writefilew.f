C ==============================================================
      SUBROUTINE WRITEFILEW
C ==============================================================
C
C   The format of 3-D plot has been changed to be compatible
C   with TECPLOT Version 10.0 or later
C
C ==============================================================
      USE TOGBFLOW
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C      COMMON /TOGBFLOW/ XCCC(NBHZP,MBPZ),YCCC(NBHZP,MBPZ),
C     %                  ZCCC(NBHZP,MBPZ),CPBMEAN(NSTEP,NBHZP,MBPZ)


CVV
      ALLOCATABLE :: YTMP(:,:),ZTMP(:,:)
      ALLOCATABLE :: POTMP(:,:)
CVV
      ALLOCATABLE :: VTOTSB(:,:),VTOTSH(:,:),VTOTSD(:,:),VTOTSTN(:,:)    !S.H.CHANG 03/10/2010
C      DIMENSION YTMP(400,100),ZTMP(400,100)
C      DIMENSION POTMP(400,100)

CVV
      ALLOCATE(YTMP(400,100),ZTMP(400,100))
      ALLOCATE(POTMP(400,100))
      ALLOCATE(VTOTSB(200,100),VTOTSH(400,100))   !S.H.CHANG 03/10/2010, S.N.KIM 02/20/2017
      IF(IDUCT.NE.0) ALLOCATE(VTOTSD(250,20))  !S.H.CHANG 03/10/2010
      IF(ITUN.NE.0) ALLOCATE(VTOTSTN(250,20))  !S.H.CHANG 03/10/2010

CVV

C-----------------------------------------------------------------------
C      Print time-averaged wetted circulation
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     Non-dimensional Circulation [100*DELP/(2*PI*R*UR)]
C     where UR/VS=SQRT(1+(.7*PI/ADVCO)**2)                      JY032898
C
C     Modified this section so that the averaged wetted circulation
C     is printed to file "wet.cir".                             JY010600
C-----------------------------------------------------------------------

      IF(NTSTEP.EQ.NTIME.AND.NTSTEP.NE.0) THEN

         UR=SQRT(1.+(.7*PI/ADVCO)**2.)

         WRITE(1700,6300)
         WRITE(1701,6300)

         DO 40 M=1,MR

            WRITE(15,6000) HRZP(1,M),NTPREV
            AVCIR=ZERO
            DO 35 NN=1, NTPREV
             WRITE(15,*) TT(NN),HUNTPI*DPHI(M,NN)/UR
             AVCIR=AVCIR+DPHI(M,NN)
 35            CONTINUE

            AVCIR1=(AVCIR/NTPREV)*HUNTPI
            AVCIR=AVCIR1/UR

            WRITE(1700,*) HRZP(1,M),AVCIR
            WRITE(1701,*) HRZP(1,M),AVCIR1

 40         CONTINUE

      END IF

C-----------------------------------------------------------------------
C       Print wetted pressures
C-----------------------------------------------------------------------
      IF(IDXREV.EQ.0) THEN

         WRITE(620,6101) 0., NCP,MRP

         IF(IHUB .NE. 0 .AND. IPHUB .EQ. 1) THEN
            WRITE(17,6100) 0.,NHBX,MHBT+1
         ENDIF

         IF(ITUN .NE. 0) THEN
            WRITE(621,6100) 0., NAXT+1,MTUNEL+1
         ENDIF

         IF(IDUCT .NE. 0) THEN
            WRITE(312,6100) 0.0,NDUCTP,MDUCTP
         ENDIF

      ELSE

         WRITE(620,6101) TT(IDXREV),NCP,MRP

         IF(IHUB .NE. 0 .AND. IPHUB .EQ. 1) THEN
            WRITE(17,6100) TT(IDXREV),NHBX,MHBT+1
         ENDIF

         IF(ITUN .NE. 0) THEN
            WRITE(621,6100) TT(IDXREV),NAXT+1,MTUNEL+1
         ENDIF

         IF(IDUCT .NE. 0) THEN
            WRITE(312,6100) TT(IDXREV),NDUCTP,MDUCTP
         ENDIF

      END IF

      IF(IDXREV.EQ.0) THEN
         DT1=ZERO
      ELSE
         DT1=-DELTAT*(IDXREV-1)
      END IF

      DO M=1,MR+1
         DO NN=1,NC+1
            YTMP(NN,M)=YB(NN,M)*COS(DT1)-ZB(NN,M)*SIN(DT1)
            ZTMP(NN,M)=YB(NN,M)*SIN(DT1)+ZB(NN,M)*COS(DT1)
         END DO
      END DO

C === Write 3D plot to file
C -- Pressure and potential on Blade

      DO N = 1, NC
         DO M = 1, MR
            L1 = INDEXB(N,M)
            POTMP(N,M) = POT(L1)
         ENDDO
      ENDDO

      WRITE(620,*) ((XB(N,M), N=1,NCP),M=1,MRP)
      WRITE(620,*) ((YTMP(N,M), N=1,NCP),M=1,MRP)
      WRITE(620,*) ((ZTMP(N,M), N=1,NCP),M=1,MRP)
      WRITE(620,*) ((-CPB(N,M)*ADVCO**2, N=1,NC), M=1,MR)
      WRITE(620,*) ((POTMP(N,M), N=1,NC), M=1,MR)
      WRITE(620,*) ((DPDTPRE(IDXREV,INDEXB(N,M)),N=1,NC), M=1,MR)
      WRITE(620,*) ((RER(N,M),N=1,NC), M=1,MR)
      WRITE(620,*) ((XCDF1(N,M),N=1,NC), M=1,MR)

      IF(ISCAV .EQ. 0) THEN
!Allen Du 01/11/2018 output the cp_min
         opt_cpmin_temp=-CPB(1,1)*ADVCO**2.
         DO M = 1, MR

            WRITE(690,6200) M, ITSTEP
!Allen Du 12/22/2017 output cp at the control points
            WRITE(7042,6200) M, ITSTEP

C---------Writing wetted pressure for pressure side.--------------------

            DO NN=1,NC/2
!Allen Du 12/22/2017 output cp at the control points
             LTEMP = indexb(nn,m)
             write(7042,*) xct(ltemp,1),-CPB(NN,M)*ADVCO**2.

             WRITE(690,*) SBP(NC/2-NN+1),-CPB(NN,M)*ADVCO**2.
!Allen Du 01/11/2018 output the cp_min
             if(opt_cpmin_temp<=-CPB(NN,M)*ADVCO**2.) then
                 opt_cpmin_temp=-CPB(NN,M)*ADVCO**2.
             endif
            ENDDO

C---------Writing wettedpressure for suction side.---------------------

            DO NN=1,NC/2
!Allen Du 12/22/2017 output cp at the control points
             LTEMP = indexb(nn+nc/2,m)
             write(7042,*) xct(ltemp,1),-CPB(NC/2+NN,M)*ADVCO**2.

             WRITE(690,*) SBP(NN),-CPB(NC/2+NN,M)*ADVCO**2.
!Allen Du 01/11/2018 output the cp_min
             if(opt_cpmin_temp<=-CPB(NC/2+NN,M)*ADVCO**2.) then
                 opt_cpmin_temp=-CPB(NC/2+NN,M)*ADVCO**2.
             endif
            ENDDO

            IF(IDXREV .NE. 0) THEN
             DO NN = 1 , NC/2
                L1 = INDEXB(NC/2+NN,M)
                L2 = INDEXB(NC/2-NN+1,M)

                CPMM = CPB(NC/2+NN,M)*SS(L1,1)-CPB(NC/2-NN+1,M)*SS(L2,1)
                CPBMEAN(IDXREV,NN,M) = -CPMM*ADVCO**2
             ENDDO
            ENDIF
!Allen Du 01/11/2018 output the cp_min
         write(7043,*)itstep,opt_cpmin_temp
         ENDDO
      ENDIF

C/S.N.KIM | Cp Blade Plotting for unsteady runs.
      IF (ISCAV.EQ.0.AND.IAN.EQ.2.AND.IVISC.EQ.1.AND.ISTEADY.NE.0) THEN
        DO M = 1, MR
          WRITE(1999+M,*) 'ZONE T= "TIME STEP =',ITSTEP,'"'  
          DO NN=1, NC/2
            L = INDEXB(NN,M)
            WRITE(1999+M,*) XCT(L,1),-CPB(NN,M)*ADVCO**2.
          ENDDO
          DO NN=1, NC/2
            L = INDEXB(NC/2+NN,M)
            WRITE(1999+M,*) XCT(L,1),-CPB(NC/2+NN,M)*ADVCO**2.
          ENDDO
        ENDDO
      ENDIF
C/S.N.KIM | Aug. 2018.

cc hongyang added
c      OPEN(691,FILE='wprs.plt',STATUS='UNKNOWN')
c      DO M=1,MR
c         WRITE(691,5122) M,ITSTEP
c         DO NN=1,NC/2
c            L = INDEXB(NN,M)
c            RTEMP = SQRT(XCT(L,2)**2+XCT(L,3)**2)
c            WRITE(691,5135) XCT(L,1),XCT(L,2),XCT(L,3),
c     &               RTEMP,-CPB(NN,M)*ADVCO**2
c         END DO
c         DO NN=1,NC/2
c            L = INDEXB(NC/2+NN,M)
c            RTEMP = SQRT(XCT(L,2)**2+XCT(L,3)**2)
c            WRITE(691,5135) XCT(L,1),XCT(L,2),XCT(L,3),
c     &                    RTEMP,-CPB(NC/2+NN,M)*ADVCO**2
c         END DO
c      END DO
c      CLOSE(691)
c 5135 FORMAT(1X,3(F10.6,1X),F6.3,1X,F12.6)
c 5122 FORMAT('ZONE T="M=',I3,' T=',I3,'"')
cc hongyang added


C -- Pressure and potential on HUB

      IF(IHUB .NE. 0 .AND. IPHUB .EQ. 1) THEN

         DO M=1,MHBT+1
            DO NN=1,NHBX
             YTMP(NN,M)=YH(NN,M)*COS(DT1)-ZH(NN,M)*SIN(DT1)
             ZTMP(NN,M)=YH(NN,M)*SIN(DT1)+ZH(NN,M)*COS(DT1)
            END DO
         END DO

         DO N = 1, NHBX-1
            DO M = 1, MHBT
             L1 = INDEXH(N,M)
             POTMP(N,M) = POT(L1)
            ENDDO
         ENDDO

         WRITE(17,*) ((XH(N,M),N=1,NHBX),M=1,MHBT+1)
         WRITE(17,*) ((YTMP(N,M),N=1,NHBX),M=1,MHBT+1)
         WRITE(17,*) ((ZTMP(N,M),N=1,NHBX),M=1,MHBT+1)
         WRITE(17,*) ((-CPH(N,M)*ADVCO**2,N=1,NHBX-1),M=1,MHBT)
         WRITE(17,*) ((POTMP(N,M),N=1,NHBX-1),M=1,MHBT)

C --- S.H.CHANG 03/10/2010
         WRITE(317,*) ((XB(N,M), N=1,NCP),M=1,MRP)
         WRITE(317,*) ((YB(N,M), N=1,NCP),M=1,MRP)
         WRITE(317,*) ((ZB(N,M), N=1,NCP),M=1,MRP)

         DO N = 1,NC
            DO M = 1,MR
               L1 = INDEXB(N,M)
               VTOTSB(N,M) = VTOTS(L1)
            END DO
         END DO

         WRITE(317,*) ((VTOTSB(N,M),N=1,NC),M=1,MR)
         WRITE(317,*) ((-CPB(N,M)*ADVCO**2,N=1,NC),M=1,MR)
         WRITE(317,*) ((POTMP(N,M), N=1,NC),M=1,MR)

         DO N = 1,NC
            DO M = 1,MR
               L1 = INDEXB(N,M)
               WRITE(318,8000) XCTP(L1,1,1),XCTP(L1,2,1),XCTP(L1,3,1),
     &                         UXTOT(N,M),UYTOT(N,M),UZTOT(N,M)
            END DO
         END DO

         WRITE(319,*) ((XH(N,M),N=1,NHBX+1),M=1,MHBT+1)
         WRITE(319,*) ((YH(N,M),N=1,NHBX+1),M=1,MHBT+1)
         WRITE(319,*) ((ZH(N,M),N=1,NHBX+1),M=1,MHBT+1)

         DO N = 1,NHBX
            DO M = 1,MHBT
               L1 = INDEXH(N,M)
               VTOTSH(N,M) = VTOTS(L1)
            END DO
         END DO

         WRITE(319,*) ((VTOTSH(N,M),N=1,NHBX),M=1,MHBT)
         WRITE(319,*) ((-CPH(N,M)*ADVCO**2,N=1,NHBX),M=1,MHBT)
         WRITE(319,*) ((POTMP(N,M),N=1,NHBX),M=1,MHBT)
         WRITE(319,*) ((DPDUH(N,M),N=1,NHBX),M=1,MHBT)
         WRITE(319,*) ((DPDVH(N,M),N=1,NHBX),M=1,MHBT)

         DO N = 1,NHBX
            DO M = 1,MHBT
               L1 = INDEXH(N,M)
               WRITE(320,8000) XCTP(L1,1,1),XCTP(L1,2,1),XCTP(L1,3,1),
     &                         UXHTOT(N,M),UYHTOT(N,M),UZHTOT(N,M)
            END DO
         END DO
C --- S.H.CHANG 03/10/2010

         DO M = 1 , MHBT
            WRITE(311,6200) M, ITSTEP
            SUM = 0.0
            DO N = 1 , NHBX
             NM = INDEXH(N,M)
             SUM = SUM + DELU(NM)
            ENDDO

            DUM1 = 0.0
            DO NN = 1 , NHBX
             IF(NN .EQ. 1) THEN
                NM = INDEXH(NN,M)
                DUM1 = DUM1 + 0.5 * DELU(NM) / SUM
             ELSE
                NM = INDEXH(NN,M)
                NM1 = INDEXH(NN-1,M)
                DUM1 = DUM1 + 0.5 * (DELU(NM) + DELU(NM1)) / SUM
             ENDIF

             WRITE(311,*) DUM1,-CPH(NN,M)*ADVCO**2
            ENDDO
         ENDDO

      ENDIF

C --- Pressure & Potential on TUNNEL

      IF(ITUN .NE. 0) THEN

         DO N = 1, NAXT
            DO M = 1, MTUNEL
             L1 = INDEXTN(N,M)
             POTMP(N,M) = POT(L1)
            ENDDO
         ENDDO

         WRITE(621,*) ((XTUN(N,M), N=1,NAXT+1),M=1,MTUNEL+1)
         WRITE(621,*) ((YTUN(N,M), N=1,NAXT+1),M=1,MTUNEL+1)
         WRITE(621,*) ((ZTUN(N,M), N=1,NAXT+1),M=1,MTUNEL+1)
         WRITE(621,*) ((-CPTN(N,M)*ADVCO**2, N=1,NAXT), M=1,MTUNEL)
         WRITE(621,*) ((POTMP(N,M), N=1,NAXT), M=1,MTUNEL)

C --- S.H.CHANG 03/10/2010
         WRITE(321,*) ((XTUN(N,M), N=1,NAXT+1),M=1,MTUNEL+1)
         WRITE(321,*) ((YTUN(N,M), N=1,NAXT+1),M=1,MTUNEL+1)
         WRITE(321,*) ((ZTUN(N,M), N=1,NAXT+1),M=1,MTUNEL+1)

         DO N = 1,NAXT
            DO M = 1,MTUNEL
               L1 = INDEXTN(N,M)
               VTOTSTN(N,M) = VTOTS(L1)
            END DO
         END DO

         WRITE(321,*) ((VTOTSTN(N,M),N=1,NAXT),M=1,MTUNEL)
         WRITE(321,*) ((-CPTN(N,M)*ADVCO**2,N=1,NAXT),M=1,MTUNEL)
         WRITE(321,*) ((POTMP(N,M),N=1,NAXT),M=1,MTUNEL)
         WRITE(321,*) ((DPDUTN(N,M),N=1,NAXT),M=1,MTUNEL)
         WRITE(321,*) ((DPDVTN(N,M),N=1,NAXT),M=1,MTUNEL)

         DO N = 1,NAXT
            DO M = 1,MTUNEL
               L1 = INDEXTN(N,M)
               WRITE(322,8000) XCTP(L1,1,1),XCTP(L1,2,1),XCTP(L1,3,1),
     &                         UXTNTOT(N,M),UYTNTOT(N,M),UZTNTOT(N,M)
            END DO
         END DO
C --- S.H.CHANG 03/10/2010

      ENDIF

C -- Pressure & Potential on DUCT

      IF(IDUCT .NE. 0) THEN

         DO N = 1, NDUCT
            DO M = 1, MDUCT
             L1 = INDEXD(N,M)
             POTMP(N,M) = POT(L1)
            ENDDO
         ENDDO

         WRITE(312,*) ((XD(N,M),N=1,NDUCTP),M=1,MDUCTP)
         WRITE(312,*) ((YD(N,M),N=1,NDUCTP),M=1,MDUCTP)
         WRITE(312,*) ((ZD(N,M),N=1,NDUCTP),M=1,MDUCTP)
         WRITE(312,*) ((-CPD(N,M)*ADVCO**2, N=1,NDUCT), M=1,MDUCT)
         WRITE(312,*) ((POTMP(N,M), N=1,NDUCT), M=1,MDUCT)

C --- S.H.CHANG 03/10/2010
         WRITE(315,*) ((XD(N,M),N=1,NDUCTP),M=1,MDUCTP)
         WRITE(315,*) ((YD(N,M),N=1,NDUCTP),M=1,MDUCTP)
         WRITE(315,*) ((ZD(N,M),N=1,NDUCTP),M=1,MDUCTP)

         DO N = 1,NDUCT
            DO M = 1,MDUCT
               L1 = INDEXD(N,M)
               VTOTSD(N,M) = VTOTS(L1)
            END DO
         END DO

         WRITE(315,*) ((VTOTSD(N,M),N=1,NDUCT),M=1,MDUCT)
         WRITE(315,*) ((-CPD(N,M)*ADVCO**2,N=1,NDUCT),M=1,MDUCT)
         WRITE(315,*) ((POTMP(N,M),N=1,NDUCT),M=1,MDUCT)
         WRITE(315,*) ((DPDUD(N,M), N=1,NDUCT), M=1,MDUCT)
         WRITE(315,*) ((DPDVD(N,M), N=1,NDUCT), M=1,MDUCT)

         DO N = 1,NDUCT
            DO M = 1,MDUCT
               L1 = INDEXD(N,M)
               WRITE(316,8000) XCTP(L1,1,1),XCTP(L1,2,1),XCTP(L1,3,1),
     &                         UXDTOT(N,M),UYDTOT(N,M),UZDTOT(N,M)
            END DO
         END DO
C --- S.H.CHANG 03/10/2010

         MDUCTTMP = MDUCT
         IF(IDUCT .EQ. 1 .AND. IDOPT .EQ. 1) MDUCTTMP=1

         DO M = 1 , MDUCTTMP

            WRITE(313,6200) M, ITSTEP
            SUM1 = 0.0
            SUM2 = 0.0

            DO N = 1 , NDUCTH
             N1 = NDUCTH - N + 1
             N2 = NDUCTH + N
             L1 = INDEXD(N1,M)
             L2 = INDEXD(N2,M)

             SUM1 = SUM1 + DELU(L1)
             SUM2 = SUM2 + DELU(L2)
            ENDDO

            DUM = 0.0
            DO N = 1 , NDUCTH
             N1 = NDUCTH - N + 1
             IF(N .EQ. 1) THEN
                L1 = INDEXD(N1,M)
                DUM = DUM + 0.5*DELU(L1) / SUM1
             ELSE
                L1 = INDEXD(N1,M)
                N10 = N1+1
                L10 = INDEXD(N10,M)
                DUM = DUM + 0.5*(DELU(L10) + DELU(L1)) /SUM1
             ENDIF

             WRITE(313,*) DUM, -CPD(N1,M)

            ENDDO

            DUM = 1.0 !SUM2
            DO N =  NDUCTH, 1 , -1
             N1 = NDUCTH + N
             IF(N .EQ. NDUCTH) THEN
                L1 = INDEXD(N1,M)
                DUM = DUM - 0.5*DELU(L1) / SUM2
             ELSE
                L1 = INDEXD(N1,M)
                N10 = N1 + 1
                L10 = INDEXD(N10,M)
                DUM = DUM - 0.5*(DELU(L10) + DELU(L1)) /SUM2
             ENDIF

             WRITE(313,*) DUM,-CPD(N1,M)

            ENDDO
         ENDDO

C -- Mean pressure

         IF(IDOPT .EQ. 0) THEN
            DO N = 1, NDUCT
             PMD = 0.0
             L1 = INDEXD(N,1)
             DO M = 1, MDUCT
                PMD = PMD -CPD(N,M)
             ENDDO
             PMD = PMD / REAL(MDUCT)
             WRITE(314,*) XCT(L1,1),PMD,-DCPMEAN(N)
            ENDDO
         ENDIF

      ENDIF

C-----------------------------------------------------------------------
C     printing out the wetted forces to *.wktkq.
C-----------------------------------------------------------------------

      IF(IDXREV.NE.0) WRITE(73,7000) TT(IDXREV),FBXP(1),FBXP(4),
     *           FBXPV(1),FBXPV(4)

C       Hong corrected unsteady forces on the duct
      IF(IDXREV.NE.0 .AND. IDUCT .NE. 0)
     *    WRITE(622,7000) TT(IDXREV),Nblade*FDXP(1),Nblade*FDXP(4),
     *           Nblade*FDXPV(1),Nblade*FDXPV(4)

C-----------------------------------------------------------------------
C     print CL,CD,CD/CL to screen for hydrofoils.               JY110300
C-----------------------------------------------------------------------
      IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
         WRITE(*,*) '-----------------------------------------'
         WRITE(*,*) 'WETTED:       CL        CD        CL/CD  '
         WRITE(*,'(8X,2(F11.6),F13.1)') CLIFT,CDRAG,CLIFT/CDRAG
         WRITE(*,*) '-----------------------------------------'
      END IF

 6000      FORMAT(1x,'ZONE T="r/R=',F6.3,'" I=',I3)
 6100      FORMAT(1x,'ZONE T="Time=',F8.1,'", I=',I5,', J=,',I5,
     *          'K=1, DATAPACKING=BLOCK,',
     *'VARLOCATION=(4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED)')
 6101      FORMAT(1x,'ZONE T="Time=',F8.1,'", I=',I5,', J=,',I5,
     *          'K=1, DATAPACKING=BLOCK,',
     *'VARLOCATION=(4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED
     *,7=CELLCENTERED,8=CELLCENTERED)')
 6200      FORMAT('ZONE T="M=',I3,' T=',I3,'"')
 6300      FORMAT('ZONE T="Fully wetted Mean Circulation"')

 7000      FORMAT(1X,F9.4,2X,F10.6,2X,F10.6,6X,F10.6,2X,F10.6)
 8000      FORMAT(6(1X,F12.6))

CVV
      DEALLOCATE(YTMP,ZTMP,POTMP)
CVV

      RETURN
      END


