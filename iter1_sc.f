      SUBROUTINE ITER1_SC
************************************************************************
*                                                                      *
*  Subroutine ITER1 computes the first iteration cavity solution for   *
*  a given sigma. Also, this subroutine calls the other subroutines    *
*  that will compute subsequent iterations.                            *
*                                                                      *
*  Author: Neal Fine  July 15, 1991                                    *
*                                                                      *
*  Date of last Revision                     Revision                  *
*  ---------------------                 ----------------              *
*  JY100101                Copied from delr.f.  Modified for ISC=1.    *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION CLNGTH(MBZ,2)

      ITER=0
      ISR=1
      IF(IFACE.EQ.2) ISR=2
      CALL ARCL

C-----------------------------------------------------------------------
C     Initialize cavity length and detachment point.  Define inital
C     guess then use either solution from previous time step or 
C     revolution.    
C-----------------------------------------------------------------------
      DO 10 IDR=1,2
         IF(NTSTEP.EQ.1)THEN

C..........initial guess of cavity length...............................
            DO 30 M=1,MR

               CAVLTE=ARCLNG(NHP-NSR2,M,IDR)/ARCLNG(NHP,M,IDR)
               IF(NOCAV(M,IDXREV,IDR).EQ.1) THEN 
                  CAVL(M,IDR)=ZERO
               ELSE
                  DDUM=ARCLNG(NLEP(M,IDXREV,IDR)+1,M,IDR)/
     *                 ARCLNG(NHP,M,IDR)

                  IF((ICON.EQ.4).OR.(ICON.EQ.5).OR.
     *                 (ICON.EQ.6).OR.(ICON.EQ.8)) THEN
                     IF(IDR.EQ.1) CAVL(M,1)=CAVINIB
                     IF(IDR.EQ.2) CAVL(M,2)=CAVINIF
                  ELSE
                     CAVL(M,IDR)=CAVLTE-0.2
                  END IF
                  
                  IF(CAVL(M,IDR).LE.DDUM) THEN
                     CAVL(M,IDR)=MIN(DDUM+0.15,CAVLTE)
                  END IF
               END IF
               CAVLT(M)=1.2

 30         CONTINUE

         ELSE IF(NREV.GT.3) THEN

C..........Use the cavity length from the same time step in the ........
C..........previous revolution as the initial guess if we already.......
C..........completed two revolutions....................................
            DO 40 M=1,MR
               IF(NOCAV(M,IDXREV,IDR).EQ.1) THEN
                  CAVL(M,IDR)=ZERO
               ELSE
                  DUMM=ARCLNG(NLEP(M,IDXREV,IDR)+3,M,IDR)/
     *                 ARCLNG(NHP,M,IDR)
                  IF(CAVLSA(M,IDXREV,IDR).LT.(DUMM+.05)) THEN
                     CAVL(M,IDR)=DUMM+.05
                  ELSE
                     CAVL(M,IDR)=CAVLSA(M,IDXREV,IDR)
                  END IF
               END IF

               CAVLTE=ARCLNG(NHP-NSR2,M,IDR)/ARCLNG(NHP,M,IDR)
               IF(CAVL(M,IDR).GT.CAVLTE) CAVL(M,IDR)=CAVLTE
               CAVLT(M)=CAVLTSA(M,IDXREV)

 40         CONTINUE

         ELSE

C..........Use the cavity length from the previous time step as the
C..........initial guess if 1<NTSTEP<(2*NREV+1).........................
            DO 50 M=1,MR
               IF(NOCAV(M,IDXREV,IDR).EQ.1) THEN
                  CAVL(M,IDR)=ZERO
               ELSE
                  DUMM=ARCLNG(NLEP(M,IDXREV,IDR)+3,M,IDR)/
     *                 ARCLNG(NHP,M,IDR)
                  IF(CAVLP(M,IDR).LT.(DUMM+.08)) THEN
                     CAVL(M,IDR)=DUMM+.08                     
                  ELSE
                     CAVL(M,IDR)=CAVLP(M,IDR)
                  END IF
               END IF

               CAVLTE=ARCLNG(NHP-NSR2,M,IDR)/ARCLNG(NHP,M,IDR)
               IF(CAVL(M,IDR).GT.CAVLTE) CAVL(M,IDR)=CAVLTE
               CAVLT(M)=CAVLTP(M)

 50         CONTINUE

         END IF

 10   CONTINUE

C-----------------------------------------------------------------------
C     find the correct cavity planform
C-----------------------------------------------------------------------
      IF((ICON.EQ.4).OR.(ICON.EQ.6)) ITERMAX=1

 2001 FORMAT('Dtol Convergence achieved at iteration: ',I4)

      ICONVR=0

      DO 80 ITER=1, ITERMAX
         IF(ICONVR.EQ.0)THEN
C..........re-set the source strengths to the wetted source strengths...
            DO 90 J=1,NPANEL
               DPDNC(J)=DPDN(J)
 90         CONTINUE

C..........re-set the wake source strengths to zero. (JY060100)
            DO J=1,NPWAKS
               SORW(J)=ZERO
            END DO
            
C..........define the linear arclength on the suction foil and in the...
C..........wake must be repeated every iteration due to split panel.....
C..........modifications................................................
            CALL ARCL

C..........compute delta(r) for a the current cavity planform...........
            IF((ICON.EQ.4).OR.(ICON.EQ.6))THEN
               DO 100 IDR=1,2
                  DO 110 M=1,MR
                     CAVL(M,IDR)=CAVINI
                     IF(IFACE.EQ.2) THEN
                        CAVL(M,1)=CAVINIB
                        CAVL(M,2)=CAVINIF
                     END IF
                     CAVLTE=ARCLNG(NHP-NSR2,M,IDR)/ARCLNG(NHP,M,IDR)
                     IF(CAVL(M,IDR).GT.CAVLTE) CAVL(M,IDR)=CAVLTE

 110              CONTINUE
 100           CONTINUE
            END IF

C-----------------------------------------------------------------------
C     Let the cavity length of the last strip (MR) equal to the
C     previous strip (MR-1).                                    JY092799
C-----------------------------------------------------------------------
C..........Same for ICON=8 as well (JY110100)
            IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8) THEN
               DO M=MRTIP+1,MR
                  DO IDR=1,2
                     IF(NOCAV(M,IDXREV,IDR).EQ.0) THEN
                        NLEP(M,IDXREV,IDR)=NLEP(M-1,IDXREV,IDR)
                        CAVL(M,IDR)=CAVL(M-1,IDR)
                     END IF
                  END DO
                  CAVLT(M)=CAVLT(M-1)
               END DO
            END IF
C-----------------------------------------------------------------------

C..........Variable IT2 added for iteration of DPDVB within each 
C..........iteration for cavity length. (JY080600)
            IT2=0

            CALL DELR_SC

C..........Store cavity length (JY072001)
            DO M=1,MR
               DO IDR=1,2
                  CLNGTH(M,IDR)=CAVL(M,IDR)
               END DO
            END DO

C-----------------------------------------------------------------------
C     Let the cavity length and cavity height of the last strip (MR)
C     equal to the previous strip (MR-1).                       JY092799
C-----------------------------------------------------------------------
C..........Same for ICON=8 as well (JY110100)
            IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8) THEN
               DO M=MRTIP+1,MR
                  DO IDR=1,2
                     IF(NOCAV(M,IDXREV,IDR).EQ.0) THEN
                        NLEP(M,IDXREV,IDR)=NLEP(M-1,IDXREV,IDR)
                        CAVL(M,IDR)=CAVL(M-1,IDR)
                        DELTA(M,IDR)=DELTA(M-1,IDR)
                        DO N=1,NHP+NTRA
                           HT(N,M,IDR)=HT(N,M-1,IDR)
                        END DO
                     END IF
                  END DO
                  CAVLT(M)=CAVLT(M-1)
                  DELTAT1(M)=DELTAT1(M-1)
               END DO
            END IF
C-----------------------------------------------------------------------

C..........Do not allow any cavity to grow if NOCAV=1. (JY061700)
            DO IDR=1,2
               DO M=1,MR
                  IF(NOCAV(M,IDXREV,IDR).EQ.1) THEN
                     CAVL(M,IDR)=ZERO
                     DELTA(M,IDR)=ZERO
                  END IF
               END DO
            END DO

C..........Same for ICON=8 as well (JY110100)
            IF((ICON.EQ.5).OR.(ICON.EQ.6).OR.(ITER.LT.5).
     *           OR.(ICON.EQ.8))THEN
               MRTIP2=MR
            ELSE
               MRTIP2=MRTIP
            ENDIF

C..........Initialize convergence parameters
            MGT=0
            IF(ITER.GT.1) THEN
               DELMXP=DELMX
               MMAXP=MMAX
            END IF
            DELMX=ZERO

            DO M=1,MRTIP2
               DO IDR=1,2

C................Count total number of strips on the back and face 
C................that didn't satisfy convergence criterion.
                  IF(ABS(DELTA(M,IDR)).GT.DTOL.OR.
     *                 ABS(DELTAT1(M)).GT.DTOL) THEN
                     MGT=MGT+1
                  ENDIF
                  
C................Find maximum delta and correponding index.
                  IF(ABS(DELTA(M,IDR)).GT.DELMX.OR.
     *                 ABS(DELTAT1(M)).GT.DELMX) THEN
                     DELMX=MAX(ABS(DELTA(M,IDR)),ABS(DELTAT1(M)))
                     MMAX=M
                  END IF

               END DO
            END DO

C..........Avoid possible infinite loop: If all the strips converged
C..........except one, and the DELTA of that strip didn't change much in
C..........two consecative iteration, then stop the iteration.
            IF(MGT.EQ.1.AND.MMAX.EQ.MMAXP.AND.ITER.GT.1)THEN
               IF(ABS(DELMX-DELMXP).LE.DTOL) DELMX=ZERO
            END IF
            
            IF(DELMX.GT.DTOL)THEN
C.............If DELMAX>DTOL, then modify cavity length 
               DO M=1,MR

                  DO IDR=1,2

                     CAVLOLD=CAVL(M,IDR)
                     DELTAOLD=DELTA(M,IDR)
                     
C...................find the cavity length for the next iteration
                     IF(NOCAV(M,IDXREV,IDR).EQ.0.
     *                    AND.ITERMAX.GT.1) THEN

C...................special consideration is needed if
C...................cavity is negative at blade T.E.
                        CAVLTE=ARCLNG(NHP-NSR2,M,IDR)/
     *                       ARCLNG(NHP,M,IDR)
                        IF(CAVL(M,IDR).EQ.CAVLTE.AND.
     *                       DELTA(M,IDR).LT.0.) THEN
                           NDUM=NHP-NSR2-5
                           DUM1=ARCLNG(NDUM,M,IDR)/
     *                          ARCLNG(NHP,M,IDR)
                           CAVL(M,IDR)=DUM1
                        ELSE

                           IF(ITER.EQ.1) THEN
C.........................change the cavity length at each strip by 1.0%
                              CAVL(M,IDR)=CAVL(M,IDR)+
     *                             SIGN(0.01,DELTA(M,IDR))
                           ELSE
C.........................use Newton-Raphson extrapolation
                              CALL EXTNR(M,IDR)                              
                           ENDIF
                           
                        END IF

                     END IF

C...................Prevent negative cavity length
                     DUM=ARCLNG(NLEP(M,IDXREV,IDR)+1,M,IDR)/
     *                    ARCLNG(NHP,M,IDR)
                     IF(CAVL(M,IDR).LT.DUM) THEN
                        CAVL(M,IDR)=DUM
                        CAVLP(M,IDR)=DUM
                        DELTA(M,IDR)=ZERO
                        DELTAP(M,IDR)=ZERO
                     ELSE
                        CAVLP(M,IDR)=CAVLOLD
                        DELTAP(M,IDR)=DELTAOLD
                     END IF

C-----------------------------------------------------------------------
C     Let the cavity length and cavity height of the last strip (MR)
C     equal to the previous strip (MR-1).                       JY092799
C-----------------------------------------------------------------------
C..........Same for ICON=8 as well (JY110100)
                     IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8)THEN
                        DO M1=MRTIP+1,MR
                           IF(NOCAV(M1,IDXREV,IDR).EQ.0) THEN
                              CAVL(M1,IDR)=CAVL(M1-1,IDR)
                              DELTA(M1,IDR)=DELTA(M1-1,IDR)
                              CAVLP(M1,IDR)=CAVLP(M1-1,IDR)
                              DELTAP(M1,IDR)=DELTAP(M1-1,IDR)
                           END IF
                        END DO
                     END IF
C-----------------------------------------------------------------------

                  END DO

                  CAVLOLD=CAVLT(M)
                  DELTAOLD=DELTAT1(M)

C................find the cavity length for the next iteration
                  IF(ITERMAX.GT.1) THEN
                     IF(ITER.EQ.1) THEN

C......................change the cavity length at each strip by 1.0%
                        CAVLT(M)=CAVLT(M)+
     *                       SIGN(0.01,DELTAT1(M))
                     ELSE

C......................use Newton-Raphson extrapolation
                        CALL EXTNR_SC(M)
                        
                     ENDIF

C...................Prevent negative cavity length
                     DUM=ONE
                     IF(CAVLT(M).LT.DUM) THEN
                        CAVLT(M)=DUM
                        CAVLTP(M)=DUM
                        DELTAT1(M)=ZERO
                        DELTATP(M)=ZERO
                     END IF
                     
                  END IF

                  CAVLTP(M)=CAVLOLD
                  DELTATP(M)=DELTAOLD

C-----------------------------------------------------------------------
C     Let the cavity length and cavity height of the last strip (MR)
C     equal to the previous strip (MR-1).                       JY092799
C-----------------------------------------------------------------------
C..........Same for ICON=8 as well (JY110100)
                  IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8)THEN
                     DO M1=MRTIP+1,MR
                        CAVLT(M1)=CAVLT(M1-1)
                        DELTAT1(M1)=DELTAT1(M1-1)
                        CAVLTP(M1)=CAVLTP(M1-1)
                        DELTATP(M1)=DELTATP(M1-1)
                     END DO
                  END IF
C-----------------------------------------------------------------------

               END DO

            ELSE

C.............If DELMAX<DTOL, then the convergence criterion has........
C.............been met and we can exit the 80 LOOP......................
               ICONVR=1
               WRITE(*,2001) ITER

            ENDIF

         ENDIF
 80   CONTINUE

C....Print final cavity lengtho f this step to file. (JY072001)
      IF(NTSTEP.EQ.1) THEN
         OPEN(98,FILE='cavl.out',STATUS='UNKNOWN')
         write(98,1009) 
      END IF

 1009    FORMAT(
     %   '============================================================'/
     %   '                                                            '/
     %   '     Cavity Length Computed from the Blade Leading edge     '/
     %   '                                                            '/
     %   '============================================================')  

 1010 FORMAT('+ Time Step=',I5,'   Blade angle=',F4.0)
      WRITE(98,1010) NTSTEP, TT(IDXREV)
      WRITE(98,'(a4,1x,10(F6.3,1X))') 'Back',(CLNGTH(M,1),M=1,MR)
      WRITE(98,'(a4,1x,10(F6.3,1X))') 'Face',(CLNGTH(M,2),M=1,MR)

C....Save the converged cavity length for the next step (CAVLP) 
C....or next revolution (CAVLSA).
      DO M=1,MR
         DO IDR=1,2
            CAVLP(M,IDR)=CAVL(M,IDR)
            DELTAP(M,IDR)=DELTA(M,IDR)
            
            CAVLSA(M,IDXREV,IDR)=CAVL(M,IDR)
         END DO

         CAVLTP(M)=CAVLT(M)
         DELTATP(M)=DELTAT1(M)
         
         CAVLTSA(M,IDXREV)=CAVLT(M)
      END DO

C-----------------------------------------------------------------------
C     compute the time-derivative of the potentials on the blade
C-----------------------------------------------------------------------
      IF(ISTEADY.NE.0.AND.ISTEADY.NE.1) CALL DPOTDT

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine ITER1>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END




