      SUBROUTINE CAVB
************************************************************************
*                                                                      *
*  Subroutine CAVB computes the  CAVity on a propeller blade           *
*  Solves the entire Blade at once (as opposed to stripwise).          *
*                                                                      *
*  Author: Neal Fine  March 3, 1991                                    *
*                                                                      *
*  Date of last Revision                     Revision                  *
*  ---------------------                 ----------------              *
*  10/22/2001                 Undergone major clean up for v2.0        *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      CHARACTER*50 FNVEC
      DIMENSION DPHI0T(MBZ)
!s---- YE TIAN ---
!     write(*,*) 'cavb.f:21,isteady=',ISTEADY
!e---- YE TIAN ---

      CALL CHRLEN(FN,LENCH)
      IF(NTSTEP.EQ.1 .OR. (IAN.EQ.2.AND.ISTEADY.NE.0)) THEN

C-----------------------------------------------------------------------
C     Determine the arclength ds
C     Zero DPHI(M,0) for ISC=1.
C-----------------------------------------------------------------------

         IF(ISC.EQ.1) THEN
           DO M=1,MR
             DO N=1,NC
               DU=HALF*(XB(N+1,M)+XB(N+1,M+1)-XB(N,M)-XB(N,M+1))
               DV=HALF*(YB(N+1,M)+YB(N+1,M+1)-YB(N,M)-YB(N,M+1))
               DW=HALF*(ZB(N+1,M)+ZB(N+1,M+1)-ZB(N,M)-ZB(N,M+1))
               SUM=SQRT(DU*DU+DV*DV+DW*DW)
               DS(N,M)=SUM
             END DO
             DPHI(M,0)=ZERO
           END DO
         ELSE
C/s S.N.KIM | For the option IDRT=1, the arclength should be defined
c           | along the 'read-in' blade surface (XB, YB, and ZB). 
c           | Originally, there was only 'CALL DSCOMP' without if statement, and
c           | DSCOMP defines blade surface using THK, CAM, and etc. similar to GBLADE.  
           IF(IDRT.EQ.1) THEN
             CALL DSCOMP_IDRT
           ELSE
             CALL DSCOMP
           ENDIF
C/e S.N.KIM | Oct. 2018
         END IF

C-----------------------------------------------------------------------
C    compute the wake-foil and the foil-wake ic's for supercavitation
C-----------------------------------------------------------------------
       if(ian.eq.2) then
           CALL SCWKINF2
       else
           CALL SCWKINF ! As of Aug. 2018, FWA (IAN=6) also uses 'SCWKINF'. | S.N.KIM
       endif

C-----------------------------------------------------------------------
C     plot the wake geometry (JY081499)
C-----------------------------------------------------------------------

         CALL WAKEPLT2

C-----------------------------------------------------------------------
C       Add file to print cavitating vectors.                   JY092999
C-----------------------------------------------------------------------

         IF(IAN .NE. 2) THEN
            CALL CHRLEN(FN,LENCH)
            FNVEC=FN(1:LENCH)//'.vec'
            OPEN(65,FILE=FNVEC,STATUS='UNKNOWN')
            WRITE(65,1003)
            WRITE(65,1004)
C            CALL PROPLT3
         ENDIF

 1003    FORMAT(1X,'TITLE="cavitating vector plot"')
 1004    FORMAT(1X,'VARIABLES="X","Y","Z","U","V","W","QC"')

C-----------------------------------------------------------------------
C      Calculate steady influence coefficient if ISTEADY=0.     JY061300
C-----------------------------------------------------------------------
         IF(ISTEADY.EQ.0 .AND. ISP.NE.1) THEN
            CALL INFGENW
C-----------------------------------------------------------------------
C       If ISP=1 and IMG=1, then the influence coefficients are read in
C       subroutine INDPOT_IM.F because the IC's due to the image panels
C       need to be added.                                       JY011300
C-----------------------------------------------------------------------

         ELSE
            IF(ISP.NE.1 .OR. IMG.EQ.0) THEN

C-----------------------------------------------------------------------
C       read in influence coefficients for the key blade
C-----------------------------------------------------------------------
C.......set W to be the infl. func. of key blade's first wake panel on..
C.......the key blade...................................................

               REWIND 81
               DO 20 M=MR,1,-1
                  CALL READ1(81,TEMP1,NPANEL)
                  DO 10 I=1,NPANEL
                     W(I,M)=TEMP1(I)
 10               CONTINUE
 20            CONTINUE
               
C -- Read Inf. of key duct's first wake panel on key blade

               IF(IDUCT .NE. 0) THEN
                  REWIND 501
                  DO M = 1, MDUCT
                     CALL READ1(501,TEMP1,NPANEL)
                     DO I = 1, NPANEL
                        WD(I,M) = TEMP1(I)
                     ENDDO
                  ENDDO
               ENDIF


C.......set WK to be the infl. func. of key blade's first wake panel on.
C.......the key blade wake..............................................

               REWIND 91
               DO 40 M=MR,1,-1
                  CALL READ1(91,TEMP4,NPWAKS)
                  DO 30 I=1,NPWAKS
                     WK(I,M)=TEMP4(I)
 30               CONTINUE
                  DO 35 JJ=2,NSUB
                     CALL READ1(91,TEMP4,NPWAKS)
 35               CONTINUE
 40            CONTINUE

C --- Set WKD to be the infl. of key duct's first wake panel on key blade wake

               IF(IDUCT .NE. 0) THEN
                  REWIND 521
                  DO M = 1 , MDUCT
                     CALL READ1(521,TEMP4,NPWAKS)
                     DO I = 1 , NPWAKS
                        WKD(I,M) = TEMP4(I)
                     ENDDO
                  ENDDO
               ENDIF

C.......key blade influence coefficients AA (dipole) and BB (source)....
C.......file 42 contains source influence of all of the blades..........
C.......file 52 contains only the key blade dipole influence............
C.......coefficients....................................................

               REWIND 42
               REWIND 52
               
               IAAA1 = NPANB
               IAAA2 = IAAA1
               
               DO J = 1 , NPANB
                  CALL READ1(52,AA(1,J),NPANEL)
                  CALL READ1(42,BB(1,J),NPANEL)
                  IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=2
CSH--------------------------------------------     
                  DO KK=2,NBLADE
                     CALL READ1(42,TEMP1,NPANEL)
                  ENDDO
CSH--REPLACE, NBLADE=1-------------------------
                  IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=1
CSH--------------------------------------------
               ENDDO

               IF(IHUB .NE. 0) THEN
                  IAAA1 = IAAA2
                  IAAA2 = IAAA2 + NPANH
                  DO J = IAAA1 + 1, IAAA2
                     CALL READ1(52,AA(1,J),NPANEL)
                     CALL READ1(42,BB(1,J),NPANEL)
                     DO KK=2,NBLADE
                        CALL READ1(42,TEMP1,NPANEL)
                     ENDDO
                  ENDDO
               ENDIF

               IF(IDUCT .NE. 0) THEN            
                  IAAA1 = IAAA2
                  IAAA2 = IAAA2 + NPAND
                  DO J = IAAA1 + 1, IAAA2
                     CALL READ1(52,AA(1,J),NPANEL)
                     CALL READ1(42,BB(1,J),NPANEL)
                     DO KK=2,NBLADE
                        CALL READ1(42,TEMP1,NPANEL)
                     ENDDO
                  ENDDO
               ENDIF

               IF(ITUN .NE. 0) THEN
                  IAAA1 = IAAA2
                  IAAA2 = IAAA2 + NPANTN
                  DO J = IAAA1+1 , IAAA2
                     CALL READ1(52,AA(1,J),NPANEL)
                     CALL READ1(42,BB(1,J),NPANEL)
                     DO KK=2,NBLADE
                        CALL READ1(42,TEMP1,NPANEL)
                     ENDDO
                  ENDDO
               ENDIF
            
c               IF(IAN .EQ. 2) THEN
c                  IAAA1 = IAAA2
c                  DO J=IAAA1+1, NPANEL
c                     CALL READ1(52,AA(1,J),NPANEL)
c                     CALL READ1(42,BB(1,J),NPANEL)
c                  ENDDO
c               ENDIF
            
C.......the next four files contain supercavitation ic's (all blades)...

               REWIND 110
               REWIND 111
               REWIND 112
               REWIND 113
               
               DO 80 J=1,NPWAKS
                  CALL READ1(110,C(1,J),NPANEL)
                  CALL READ1(113,F(1,J),NPWAKS)
                  DO 70 KK=2,NBLADE
                     CALL READ1(110,TEMP1,NPANEL)
                     CALL READ1(113,TEMP4,NPWAKS)
 70               CONTINUE
 80            CONTINUE
               
               DO 100 J=1,NPANEL
                  CALL READ1(111,D(1,J),NPWAKS)
                  CALL READ1(112,E(1,J),NPWAKS)
                  DO 90 KK=2,NBLADE
                     CALL READ1(111,TEMP4,NPWAKS)
                     CALL READ1(112,TEMP4,NPWAKS)
 90               CONTINUE
 100           CONTINUE
      
C-----------------------------------------------------------------------
C       Add the influence of the hub dipole disk to D(I,J).     JY010600
C
C       Modified IF statements to accomodate new hub options.   JY052401
C-----------------------------------------------------------------------
               IF(IHUB.NE.0) THEN
                  DO NN=1,NHBX
                     DO MM=1,MHBT
                        J=INDEXH(NN,MM)
                        IF(NN.EQ.1) THEN
                           DO I=1,NPWAKS
                              D(I,J)=D(I,J)+CHDK1(I)/FLOAT(MHBT)
                           END DO
                        ELSE IF(NN.EQ.NHBX) THEN
                           DO I=1,NPWAKS
                              D(I,J)=D(I,J)+CHDKDT1(I)/FLOAT(MHBT)
                           END DO
                        END IF
                     END DO
                  END DO
               END IF
               
C.......initialize the wake source strength with zero...................

               DO 110 I=1,NPWAKS
                  TEMP4(I)=ZERO
 110           CONTINUE
               NREC = 360 / ndltat
               DO 120 N=1,NREC
                  CALL WRITE2(48,N,TEMP4,NPWAKS)
 120           CONTINUE
               
               IF(ISP.EQ.1) THEN
C........Initialize the blade source and dipole strengths and wake 
C........dipole strengths with zero if ISP=1 (JY030600)
                  DO I=1,NPANEL
                     TEMP1(I)=ZERO
                  END DO
                  NREAD=NWMIN*MR
                  DO I=1,NREAD
                     TEMP5(I)=ZERO
                  END DO
                  DO N=1,NREC
                     CALL WRITE2(45,N,TEMP1,NPANEL)
                     CALL WRITE2(47,N,TEMP1,NPANEL)
                     CALL WRITE2(46,N,TEMP5,NREAD)
                  END DO
               END IF
               
            END IF
         ENDIF
      ENDIF
C-----------------------------------------------------------------------
C     Cal. induced potentials due to image panels and add them to
C     the original induced coefficients if ISP=1.               JY011200
C-----------------------------------------------------------------------
      IF(ISP.EQ.1.AND.IMG.EQ.1) CALL INDPOT_IM
C-----------------------------------------------------------------------
C     compute the inflow velocities at the current time step
C-----------------------------------------------------------------------
C T....Temporary change:  If ISP=1, use velocity distribution in Page 56
C T....of Olofsson's dissertation. (JY090900)
      IF(ISP.EQ.1.AND.ICON.NE.7.AND.FN(1:LENCH).EQ.'m841b') THEN
         CALL M841BFLOW
      ELSE
         IF(ICON .EQ. 5 .AND. ITERGB .EQ. 1) THEN
            CALL INFLOWH
         ELSE
            CALL INFLOW
         ENDIF
C.....at the control points in the wake.................................
         CALL INFLOWK
      END IF
CT....end of temporary changes. (JY090900)
C....We do not need to iterate to find the correct cavity height for
C....surface piercing propeller. (JY112199)
      IF(ISP.EQ.1) THEN
         
         IF(ISEARCH.EQ.0) ITERMAX=1
         
         CALL ITER1_SP
         
         CALL CAVHTP

         IF(ISTEADY.NE.0.AND.ISTEADY.NE.1) THEN
            CALL DPOTDT_SP
            CALL DPOTDTW_SP
         END IF
         
         GO TO 185
      END IF
C-----------------------------------------------------------------------
C    compute the first iteration cavity shape and planform
C-----------------------------------------------------------------------

      IF(ISC.EQ.0) THEN
         CALL ITER1
      ELSE
         CALL ITER1_SC
      END IF
C-----------------------------------------------------------------------
C    Save the cavity height to HTP      
C-----------------------------------------------------------------------
      CALL CAVHTP
C-----------------------------------------------------------------------
C     compute the cavity volume and write it to *.vol
C-----------------------------------------------------------------------
      IF(ISC.EQ.0) THEN
         CALL CAVVOL
      ELSE
         CALL CAVVOL_SC
      END IF

 185  CONTINUE

C-----------------------------------------------------------------------
C     Calculate circumferentially averaged dipole strength at far wake
C     for ISC=1
C-----------------------------------------------------------------------
      IF(ISC.EQ.1.AND.ISP.EQ.0) THEN
         IF(IDXREV.EQ.1) THEN
            DO M=1,MR
               DPHI0T(M)=ZERO
            END DO
         END IF

         DO M=1,MR
            DPHI0T(M)=DPHI0T(M)+DELP(M)
         END DO
         
         IF(IDXREV.EQ.NTPREV) THEN
            DO M=1,MR
               DPHI(M,0)=DPHI0T(M)/FLOAT(NTPREV)
            END DO
         END IF
      END IF

C-----------------------------------------------------------------------
C     write various quantities to plotting files
C-----------------------------------------------------------------------
      CALL CAVOUT
 
C-----------------------------------------------------------------------
C     Calculate the forces
C-----------------------------------------------------------------------
      IF(ISP.EQ.1) THEN
         CALL FORCEV_SP
         WRITE(*,5080) NTSTEP
 5080    FORMAT(1X,'Time step ',I3,' .............  completed.')
      ELSE IF(ISP.EQ.0) THEN
         IF(ISC.EQ.1) THEN
            CALL FORCEV_SC
         ELSE
            CALL FORCEV
         END IF
      END IF

C,,,,,,,,,,,,,,,, Hong Added Viscous Run ,,,,,,,,,,,May 2006,,,,,,,,,,,,
c      IF(IVISC.EQ.1) THEN  
      IF(IVISC.EQ.1.AND.NREV.EQ.NTREV) THEN  
         CALL GEO2DC 
         write(*,*) 'geo2dc done' 
         CALL INFLOWBL2
         write(*,*) 'inflowbl2 done'
         CALL CAVBL2D_1
         write(*,*) 'cavbl2d_1 done'
       
c       DO N1 = 1, NC
c         DO M1 = 1, MR
c            IF(LVFLAG(M1).NE.1) CPB(N1,M1) = CPBL(N1,M1)
c         ENDDO
c       ENDDO

      CALL FORCE_BL       
       
        OPEN(810,FILE='bl.ktkq',STATUS='UNKNOWN')

        WRITE(810,5200)
 5200   FORMAT(1X,'TITLE="Pot. & vis. TOTAL KT & KQ (cavitating)"')
        WRITE(810,5220)
 5220   FORMAT(1X,'VARIABLES="ANGLE","KT","KQ","KTV","KQV"')
        
        WRITE(810,5260) NBLADE*FBXP(1),NBLADE*FBXP(4),
     &                  NBLADE*FBXPV(1),NBLADE*FBXPV(4)
 5260   FORMAT(5(1X,E14.7))
        
      ENDIF 
    
C................................


      IF(ICON .NE. 5) THEN
C.....Print the convergence history of the cavitating force coefficients
         IF(NTSTEP.EQ.1) THEN
            OPEN(997,FILE='force-his.cav',STATUS='UNKNOWN')
            
            IF(ISTEADY.EQ.0) THEN
               WRITE(997,*) 
     %           'TITLE="Cavitating 6-comp. Forces per Blade',
     %                   ' at each time step"'
               WRITE(997,*) 'VARIABLES="STEP","FX","FY","FZ",',
     %                                       '"MX","MY","MZ"'
            ELSE
               WRITE(997,*) 
     %           'TITLE="Cavitating 6-comp. Forces per Blade',
     %                   ' at each revolution"'
               WRITE(997,*) 'VARIABLES="ANGLE","FX","FY","FZ",',
     %                                        '"MX","MY","MZ"'
            END IF
         END IF

         IF(ISTEADY.EQ.0) THEN
C            WRITE(997,'(1X,I2,6(1X,E14.7))') NTSTEP,(FXV(I1),I1=1,6)
            WRITE(997,'(1X,I2,6(1X,E14.7))') NTSTEP,(FTOTAL(I1),I1=1,6)
         ELSE
            IF(IDXREV.EQ.1) 
     %           WRITE(997,*) 'ZONE T="Revolution No. = ',NREV,'"'
C            WRITE(997,'(7(1X,E14.7))') TT(IDXREV),(FXV(I1),I1=1,6)
            WRITE(997,'(7(1X,E14.7))') TT(IDXREV),(FTOTAL(I1),I1=1,6)
         END IF
         IF(NTSTEP.EQ.NCTIME) CLOSE(997)
      ENDIF

      IF(IDXREV.NE.0) THEN
         DO KK1=1,6
C            XKT(IDXREV,KK1)=FXP(KK1)
C            XKTV(IDXREV,KK1)=FXV(KK1)
            XKT(IDXREV,KK1)=FBXP(KK1)
            XKTV(IDXREV,KK1)=FBXPV(KK1)
         END DO
      END IF

C-----------------------------------------------------------------------
C     modify detachment line if SEARCH=1
C     NOTE: This routine should only be call AFTER CAVOUT!!!(JY062399)
C     NOTE: This routine should only be call AFTER FORCEV!!!(JY091899)
C
C     For surface piercing propellers, let's turn-off detachment
C     search routine.
C-----------------------------------------------------------------------
      IF(ISP.NE.1) CALL DETACH2

C-----------------------------------------------------------------------
C     write various quantities to scratch files for later use
C-----------------------------------------------------------------------
C.....write the potentials and the delp's to a scratch file ............
      NREAD=NWMINFW*MR
      NREC = 360 / ndltat

      IF(ISTEADY.EQ.0) THEN

         DO L=1,NWMIN
            DO M=1,MR
               IDX=MR*(L-1)+M
               TEMP4(IDX)=DELP(M)
            END DO
         END DO

         IF(IDUCT .NE. 0) THEN
            NPWAKED = MDUCT * NDWK
            DO M = 1, MDUCT
               DO N = 1, NDWK
                  L1 = INDEXWD(N,M)
                  TEMPD4(L1) = DELPD(M)
               ENDDO
            ENDDO
         ENDIF
         
         DO I=1,NREC
            CALL WRITE2(45,I,POT,NPANEL)
            CALL WRITE2(47,I,DPDNC,NPANEL)
            CALL WRITE2(48,I,SORW,NPANEL)
            CALL WRITE2(46,I,TEMP4,NREAD)
            IF(IDUCT .NE. 0) CALL WRITE2(49,I,TEMPD4,NPWAKED)
         END DO

      ELSE

         irec = itstep / ndltat + 1
         CALL WRITE2(45,IREC,POT,NPANEL)
         CALL WRITE2(47,IREC,DPDNC,NPANEL)
         CALL WRITE2(48,IREC,SORW,NPWAKS)
      
C.....Read last time step DPhi..........................................
         IREC1=IREC-1
         IF(IREC1.LE.0) THEN
            IREC1=IREC1+NREC
         END IF
         CALL READ2(46,IREC1,TEMP5,NREAD)
         IF(IDUCT .NE. 0) CALL READ2(49,IREC1,TEMPD4,NPWAKED)

C.....Write current time step DPhi......................................
         DO 140 L=NWMINFW-1,1,-1
            DO 130 M=MR,1,-1
               IDX0=MR*(L-1)+M
               IDX1=MR*L+M
               TEMP5(IDX1)=TEMP5(IDX0)
 130        CONTINUE
 140     CONTINUE               
         DO 150 M=MR,1,-1
            TEMP5(M)=DELP(M)
 150     CONTINUE

         IF(IDUCT .NE. 0) THEN
            DO N = NDWK-1, 1, -1
               DO M = 1 , MDUCT
                  L1 = INDEXWD(N,M)
                  L2 = INDEXWD(N+1,M)
                  TEMPD4(L2) = TEMPD4(L1)
               ENDDO
            ENDDO            
            DO M = 1 , MDUCT
               TEMPD4(M) = DELPD(M)
            ENDDO
         ENDIF
         
         CALL WRITE2(46,IREC,TEMP5,NREAD)
         IF(IDUCT .NE. 0) CALL WRITE2(49,IREC,TEMPD4,NPWAKED)
         
      END IF

      RETURN
C<<<<<<<<<<<<<<<<<<<<<<end of subroutine CAVB>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END







