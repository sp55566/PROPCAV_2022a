      SUBROUTINE SIMEQN
************************************************************************
*     SIMEQN: SIMultaneous EQuatioN solution                           *
*      --- Solve the linear system by a block iterative matrix solver  *
*                                                                      *
*  Date of last revision   Revision                                    *
*  ---------------------   --------                                    *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      character*29 :: potcheck,filename9 ! S.N.KIM tied to plot all the potentials.

C-----------------------------------------------------------------------
C     Set up parameters for the iterative matrix solver
C-----------------------------------------------------------------------
      NREAD=NWMIN*MR
      NREC = 360 / ndltat
      IUMAS=50
      IUMAU=51
      IUMAK=52

C.....Note: If ITMAX>100, You MUST set ITMX=ITMAX+1 in BLIC2.F(JY051800)
      ITMAX=50

      MR2=MRTIP
      NORD=NPANEL

C.....Tolerance of the matrix solution is set to be 0.000005............
      TOL=5.0E-06
      DIF=9999.0

C-----------------------------------------------------------------------
C     Solve the equation, and write down the solutions
C-----------------------------------------------------------------------
      IF(NTSTEP.EQ.0) THEN

         IF(IDUCT .NE. 0) THEN
            OPEN(996,FILE='ductpot.plt',STATUS='unknown')
            WRITE(996,2000)
            WRITE(996,2100)
         ENDIF

C.......Steady solution.................................................
         CALL STSOLN

         IF(ITUN .NE. 0) CALL TUNPOT
C.......assume steady solutions to be "the solution" at all time steps..
         NREC = 360 / ndltat
         DO 45 L=1,NWMIN
            DO 40 M=1,MR
               IDX = indexw2(l,m)
               TEMP4(IDX)=DELP(M)
 40         CONTINUE
 45      CONTINUE


         IF(IDUCT .NE. 0) THEN

            NPWAKED = MDUCT * NDWK  ! NDWK is defined as 50 in ductwakegeo.f

            WRITE(996,*) 'ZONE T="Steady Wetted Circulation"'

            DO M = 1, MDUCT
               DO N = 1, NDWK
                  L1 = INDEXWD(N,M)
                  TEMPD4(L1) = DELPD(M)
               ENDDO

               LL = INDEXD(1,M)
               ANGL = DANGLE(XCT(LL,3),XCT(LL,2))
               WRITE(996,*) ANGL*180./PI, HUNTPI*DELPD(M)
            ENDDO
         ENDIF

         NREAD=NWMIN*MR
         DO 50 I=1,NREC
            CALL WRITE2(45,I,POT,NPANEL)
            CALL WRITE2(47,I,DPDN,NPANEL)
            CALL WRITE2(46,I,TEMP4,NREAD)
            IF(IDUCT .NE. 0) CALL WRITE2(49,I,TEMPD4,NPWAKED)
 50      CONTINUE
      ELSE
C.......set W to be the infl. func. of key blade's first wake panel.....
         REWIND 81
         DO 55 M=MR,1,-1
            CALL READ1(81,TEMP1,NPANEL)
            DO 56 I=1,NPANEL
               W(I,M)=TEMP1(I)
 56         CONTINUE
 55      CONTINUE

C.......Unsteady solution...............................................

         CALL USSOLN

C.......Update the unsteady solutions PHI...............................
         IREC = itstep / ndltat + 1
         itmp = itstep / ndltat

         CALL WRITE2(45,IREC,POT,NPANEL)
         CALL WRITE2(47,IREC,DPDN,NPANEL)

C.......Read last time step DPhi........................................
         IREC1=IREC-1
         IF(IREC1.LE.0) THEN
            IREC1=IREC1+NREC
         END IF

         CALL READ2(46,IREC1,TEMP5,NREAD)
         IF(IDUCT .NE. 0) CALL READ2(49,IREC1,TEMPD4,NPWAKED)

C.......Generate current time step DPhi.................................

         IF(IAN .NE. 2) THEN
            DO 65 L=NWMIN-1,1,-1
               DO 60 M=MR,1,-1
                  idx0=indexw2(l,m)
                  idx1=indexw2(l+1,m)
                  TEMP5(IDX1)=TEMP5(IDX0)
 60            CONTINUE
 65         CONTINUE

            DO 70 M=MR,1,-1
C..............IAN 2 & 6 use different wake index than other IAN options. 
C..............We need to be careful on this! S.N.KIM Aug. 2019.
               IF(IAN.EQ.6) THEN
                  I1 = INDEXW2(1,M)
                  TEMP5(I1) = DELP(M)
               ELSE
                  TEMP5(M)=DELP(M)
               ENDIF
 70         CONTINUE

            IF(IDUCT .NE. 0) THEN

               WRITE(996,*) 'ZONE T="Unsteady Wetted Circulation"'

               DO N = NDWK-1, 1, -1
                  DO M = 1 , MDUCT
                     L1 = INDEXWD(N,M)
                     L2 = INDEXWD(N+1,M)
                     TEMPD4(L2) = TEMPD4(L1)
                  ENDDO
               ENDDO

               DO M = 1 , MDUCT
                 TEMPD4(M) = DELPD(M)
                 LL = INDEXD(1,M)
                 ANGL = DANGLE(XCT(LL,3),XCT(LL,2))
                 WRITE(996,*) ANGL*180./PI, HUNTPI*DELPD(M) 
              ENDDO
            ENDIF

         ELSE

            DO M = 1, MR
               DO L = NWMIN, 2, -1
                  I1 = INDEXW2(L,M)
                  TEMP5(I1) = TEMP5(I1-1)
               ENDDO
            ENDDO

            DO M = 1 , MR
               I1 = INDEXW2(1,M)
               TEMP5(I1) = DELP(M)
            ENDDO

         ENDIF

C.......Write current time step DPhi....................................
C../ Warning from S.N.KIM. The order of saving wake strength in the unit
C../ 46 is DIFFERENT from that of cavity calculation. Whenever reading the
C../ unit 46, careful attention is required on the type of calculation,
C../ i.e., fully-wetted or cavitating. cf) cavb.f
         CALL WRITE2(46,IREC,TEMP5,NREAD)
         IF(IDUCT .NE. 0) CALL WRITE2(49,IREC,TEMPD4,NPWAKED)

      END IF

 2000 FORMAT('TITLE="Circulation Distribution:G=DELP/(2*pi*R)"')
 2100 FORMAT('Variables="Angle(Deg)","100G"')

      RETURN
      END









