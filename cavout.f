      SUBROUTINE CAVOUT
************************************************************************
*                                                                      *
*  Output various quantities from the cavity solution                  *
*  Author: Neal Fine  May 22, 1991                                     *
*  Date:      Revision                                                 *
*  -----      ---------------                                          *
*  121696CM   Changed file *.ht1 so that it outputs all cavity         *
*             heights for all spanwise sections for all timesteps      *
*  011097CM   File *.ht1 is obsolete.  Changed subroutine so that      *
*             a file is generated with the cavity panels on it         *
*  030597CM   Moved section of code that plots out the cavities to     *
*             a new subroutine, CAVPLT                                 *
*  041897CM   Added the Cp to file.cavp                                *
*  CM011898   Looking at the right hand side of matrix solution on a   *
*             contour plot.  That's why common statement was added     *
*  011998JY   Added the Cp(Suct) and Cp(Pres) to 2D contour plot       *
*             file *.cvp.                                              *
*  012398JY   Corrected the pressure calculation for the pressure side *
*             and changed the output format of *.prs so it can be      *
*             compare to that of HPUF3A.                               *
*  030198JY   Removed output files: *.pot, *.potw, *.src, *.qws, and   *
*             *.cht because they are in PROPLOT format.                *
*  080698JY   Completely redo subroutine vcavplt.f.                    *
*  092898JY   Completely redo subroutine cavplt2d.f.                   *
*  100598JY   Created a new subroutine, CAVVEC.F, which is being called*
*             in this subroutine to plot the cavity total velocity     *
*             vectors on the blade for the last cavitating revolution. *
*  JY030999  Modified subroutine to allow cavity to grow on both the   *
*            back and face of the foil.                                *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
*  JY061099  Moved the section to plot the expanded cavity to cavexp.f *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
C      DIMENSION YB1(NBPZ,MBPZ),ZB1(NBPZ,MBPZ)
      DIMENSION YTMP(400,100),ZTMP(400,100)
      DIMENSION POTMP(400,100)

C      SAVE
      CHARACTER*30 FNPLT,FNPRS,FNVEC
C.....initialize filenames for debugging................................

C      CHARACTER STRIP(40)*7
C      DATA STRIP/'STRIP01','STRIP02','STRIP03','STRIP04','STRIP05'
C     *          ,'STRIP06','STRIP07','STRIP08','STRIP09','STRIP10'
C     *          ,'STRIP11','STRIP12','STRIP13','STRIP14','STRIP15'
C     *          ,'STRIP16','STRIP17','STRIP18','STRIP19','STRIP20'
C     *          ,'STRIP21','STRIP22','STRIP23','STRIP24','STRIP25'
C     *          ,'STRIP26','STRIP27','STRIP28','STRIP29','STRIP30'
C     *          ,'STRIP31','STRIP32','STRIP33','STRIP34','STRIP35'
C     *          ,'STRIP36','STRIP37','STRIP38','STRIP39','STRIP40'/

      CALL CHRLEN(FN,LENCH)

C.....pressure..........................................................

C-----------------------------------------------------------------------
C     This purpose of this file is to check the convergence of the
C     cavity shape.  So I want it call at every time step.     JY100398
C-----------------------------------------------------------------------
C.....Only call CAVPLT2D at every timestep if ISP=0 (JY030900)
      IF(ISP.NE.1) THEN
         IF(ISC.EQ.0) THEN
            CALL CAVPLT2D
         ELSE
            CALL CAVPLT2D_SC
         END IF
      END IF

C---------Ok, here is where we add Cp to file.cavp----CM 041897---------
C---------We only want to add the last revolution to the file, so I-----
C---------think we should add an if statement here:  cm 051397----------
C-----------------------------------------------------------------------
C         This endif asks the question of weather or not we are in------
C         the last revolution                                  CM052598-
C-----------------------------------------------------------------------

      IF(NREV.EQ.NTREV)THEN

         IF(ISP.EQ.1) THEN
            CALL CAVPLT2D_SP
            CALL PLTDELP
         END IF

         WRITE(630,6101) TT(IDXREV), NCP,MRP

         IF(IDXREV.EQ.0) THEN
            DT1=ZERO
         ELSE
            DT1=-DELTAT*(IDXREV-1)
         END IF

         DO M=1,MR+1
            DO N=1,NC+1
               YTMP(N,M)=YB(N,M)*COS(DT1)-ZB(N,M)*SIN(DT1)
               ZTMP(N,M)=YB(N,M)*SIN(DT1)+ZB(N,M)*COS(DT1)
            END DO
         END DO

         DO N = 1, NC
            DO M = 1, MR
               L1 = INDEXB(N,M)
               POTMP(N,M) = POT(L1)
            ENDDO
         ENDDO

         WRITE(630,*) ((XB(N,M), N=1,NCP),M=1,MRP)
         WRITE(630,*) ((YTMP(N,M), N=1,NCP),M=1,MRP)
         WRITE(630,*) ((ZTMP(N,M), N=1,NCP),M=1,MRP)
         WRITE(630,*) ((-CPB(N,M)*ADVCO**2, N=1,NC), M=1,MR)
         WRITE(630,*) ((POTMP(N,M), N=1,NC), M=1,MR)
         WRITE(630,*) ((RER(N,M), N=1,NC), M=1,MR)
         WRITE(630,*) ((XCDF1(N,M), N=1,NC), M=1,MR)

cCyiran Added
c      OPEN(691,FILE='wprs_cav.plt',STATUS='UNKNOWN')
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
cCYIRANFINISHED


C -- Pressure and potential on HUB

      IF(IHUB .NE. 0 .AND. IPHUB .EQ. 1) THEN

           WRITE(18,6100) TT(IDXREV), NHBX,MHBT+1

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

         WRITE(18,*) ((XH(N,M), N=1,NHBX),M=1,MHBT+1)
         WRITE(18,*) ((YTMP(N,M), N=1,NHBX),M=1,MHBT+1)
         WRITE(18,*) ((ZTMP(N,M), N=1,NHBX),M=1,MHBT+1)
         WRITE(18,*) ((-CPH(N,M)*ADVCO**2, N=1,NHBX-1), M=1,MHBT)
         WRITE(18,*) ((POTMP(N,M), N=1,NHBX-1), M=1,MHBT)

         DO M = 1 , MHBT
            WRITE(310,6200) M, ITSTEP
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

             WRITE(310,*) DUM1,-CPH(NN,M)*ADVCO**2
            ENDDO
         ENDDO

         ENDIF

C --- Pressure & Potential on TUNNEL

      IF(ITUN .NE. 0) THEN

           WRITE(632,6100) TT(IDXREV), NAXT+1,MTUNEL+1

         DO N = 1, NAXT
            DO M = 1, MTUNEL
             L1 = INDEXTN(N,M)
             POTMP(N,M) = POT(L1)
            ENDDO
         ENDDO

         WRITE(632,*) ((XTUN(N,M), N=1,NAXT+1),M=1,MTUNEL+1)
         WRITE(632,*) ((YTUN(N,M), N=1,NAXT+1),M=1,MTUNEL+1)
         WRITE(632,*) ((ZTUN(N,M), N=1,NAXT+1),M=1,MTUNEL+1)
         WRITE(632,*) ((-CPTN(N,M)*ADVCO**2, N=1,NAXT), M=1,MTUNEL)
         WRITE(632,*) ((POTMP(N,M), N=1,NAXT), M=1,MTUNEL)

      ENDIF

C -- Pressure & Potential on DUCT
        IF(NREV.EQ.NTREV)THEN
           IF(IDUCT .NE. 0) THEN

              WRITE(312,6100) TT(IDXREV), NDUCTP,MDUCTP

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

                 DUM = 1.0 ! =SUM2/SUM2 (corrected by S.N.KIM 10/15/2016)
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

C-- Mean pressure

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

      ENDIF

      ENDIF

C-----------------------------------------------------------------------
C     IFILE=0: Only plot pressure, cavity, and velocity vectors if
C              in the last revolution.
C     IFILE=1: Plot pressure, cavity, and velocity vectors at every
C              time step, but overwrite the files at every new
C              revolution.
C-----------------------------------------------------------------------

      IPR=0
      IF(IFILE.EQ.0) THEN
         IF(NREV.EQ.NTREV) IPR=1
      ELSE IF(IFILE.EQ.1) THEN
         IPR=1
         IF(NREV.GT.1.AND.IDXREV.EQ.1) THEN

C..........open pressure file (*.prs)

            FNPRS=FN(1:LENCH)//'.prs'
            OPEN(680, FILE=FNPRS, STATUS='UNKNOWN')
            WRITE(680,6300)
            WRITE(680,6400)

C..........open 3-D cavity file (*.plt)

            CALL CHRLEN(FN,LENCH)
            FNPLT=FN(1:LENCH)//'.plt'
            OPEN(700, FILE=FNPLT, STATUS='UNKNOWN')
            WRITE(700,*) 'TITLE="3-D Cavity Plot"'
            WRITE(700,6500)
            CALL PROPLT2

C..........open velocity vector file (*.vec)

            IF(IAN .NE. 2) THEN
               CALL CHRLEN(FN,LENCH)
               FNVEC=FN(1:LENCH)//'.vec'
               OPEN(65,FILE=FNVEC,STATUS='UNKNOWN')
               WRITE(65,1003)
               WRITE(65,1004)
C               CALL PROPLT3
            ENDIF
         END IF
      END IF


 1003 FORMAT(1X,'TITLE="cavitating vector plot"')
 1004 FORMAT(1X,'VARIABLES="X","Y","Z","U","V","W","QC"')

 6100 FORMAT(1x,'ZONE T="Time=',F8.1,'", I=',I5,', J=',I5,
     *     ', K=1, DATAPACKING=BLOCK,',
     *     'VARLOCATION=(4=CELLCENTERED,5=CELLCENTERED)')
 6101 FORMAT(1x,'ZONE T="Time=',F8.1,'", I=',I5,', J=',I5,
     *     ', K=1, DATAPACKING=BLOCK,',
     *     'VARLOCATION=(4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED,7=CELLCENTERED)')
 6200 FORMAT('ZONE T="M=',I3,' T=',I3,'"')
 6300 FORMAT(1x,'TITLE="Cavitating Pressures on Blade (2D)"')
 6400 FORMAT(1x,'VARIABLES = "Xd/C", "-Cp"')
 6500 FORMAT(1X,'VARIABLES = "x", "y", "z"')

C 1003 FORMAT(1X,'TITLE="cavitating vector plot"')
C 1004 FORMAT(1X,'VARIABLES="X","Y","Z","U","V","W","QC"')

C-----------------------------------------------------------------------
C     Print cavitating pressure:
C     The pressures are printed in one column with the pressure side
C     first (T.E. to L.E.) then the suction side (L.E. to T.E.)
C-----------------------------------------------------------------------
      IF(IPR.EQ.1) THEN

 6600    FORMAT('ZONE T="M=',I3,' T=',I3,'"')

!Allen Du 01/11/2018 output the cp_min
         opt_cpmin_temp=-CPB(1,1)*ADVCO**2.

         DO M=1,MR
            WRITE(680,6600) M, ITSTEP

            DO N=1,NC/2
               WRITE(680,*) SBP(NC/2-N+1),-CPB(N,M)*ADVCO**2.
!S Kim 03/21/2022 output the cp_min | Allen forgot this part...
               if(opt_cpmin_temp<=-CPB(N,M)*ADVCO**2.) then
                  opt_cpmin_temp=-CPB(N,M)*ADVCO**2.
               endif
            END DO
            DO N=1,NC/2
               WRITE(680,*) SBP(N),-CPB(NC/2+N,M)*ADVCO**2.
!Allen Du 01/11/2018 output the cp_min
!S Kim 03/21/2022 changed -CPB(N,M) to -CPB(NC/2+N,M)...
               if(opt_cpmin_temp<=-CPB(NC/2+N,M)*ADVCO**2.) then
                  opt_cpmin_temp=-CPB(NC/2+N,M)*ADVCO**2.
               endif
            END DO
         END DO

C.......If IFILE=1, close file 680 at the end of every revolution.
         IF(IFILE.EQ.1.AND.IDXREV.EQ.NTPREV) CLOSE(680)
         IF(IFILE.EQ.1.AND.ISTEADY.EQ.0) CLOSE(680)
!Allen Du 01/11/2018 output the cp_min
         if(itstep.EQ.0) then
             WRITE(7043,7046)
         endif
         write(7043,*)itstep,opt_cpmin_temp
      END IF
 7046 FORMAT(1x,'ZONE T="cavitating"')

C-----------------------------------------------------------------------
C     Plot cavity
C-----------------------------------------------------------------------

      IF(IPR.EQ.1) THEN
         CALL HTPLT
         IF(ISP.NE.1) THEN
            IF(IAN .NE. 2) CALL VECPLT
         END IF
      END IF

      IF(ISC.EQ.0) THEN
         CALL VCAVPLT(IPR)
      ELSE
         CALL VCAVPLT_SC(IPR)
      END IF

      RETURN
      END

