C =======================================================================
      SUBROUTINE WRITEFILEC
C =======================================================================

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      
      CHARACTER*29 FNHIS,FNCIRC      
      DIMENSION FMAX(6),FVMAX(6)
      
C-----------------------------------------------------------------------
C     output the circulation history
C-----------------------------------------------------------------------
      CALL CHRLEN(FN,LENCH)

      IF(ICON .NE. 5 .AND. ISTEADY .NE. 0) THEN
         FNCIRC=FN(1:LENCH)//'.circ'
         OPEN(15,FILE=FNCIRC,STATUS='UNKNOWN',FORM='FORMATTED')
         WRITE(15,5000)
         WRITE(15,5020)
      ENDIF

 5000 FORMAT(1x,'TITLE = "Cavitating Circulation distribution"')
 5020 FORMAT(1x,'VARIABLES = "Angle","100*DELP/(2*PI*R*UR)"')

C-----------------------------------------------------------------------
C     Non-dimensional Circulation [100*DELP/(2*PI*R*UR)]
C     where UR/VS=SQRT(1+(.7*PI/ADVCO)**2)                      JY032898
C
C     Modified this section so that the averaged cavitating circulation
C     is printed to file "cav.cir".                             JY010600
C-----------------------------------------------------------------------
      IF(ISTEADY.NE.0) THEN

         WRITE(1700,*) 
     %        'ZONE T="Cavitating Mean Circulation"'
         WRITE(1701,*) 
     %        'ZONE T="Cavitating Mean Circulation"'         
         UR=SQRT(1.+(.7*PI/ADVCO)**2.)
         
         DO 310 M=1,MR
            WRITE(15,5115) HRZP(1,M), NTPREV
            
            AVCIR=ZERO
            DO 308 NN=1, NTPREV
               WRITE(15,*) TT(NN),HUNTPI*DPHI(M,NN)/UR
               AVCIR=AVCIR+DPHI(M,NN)
 308        CONTINUE

            AVCIR1=(AVCIR/FLOAT(NTPREV))*HUNTPI

            IF(ICON .EQ. 5) THEN
               WRITE(1700,*) HRZP(1,M),AVCIR1
               WRITE(1701,*) HRZP(1,M),AVCIR1
            ELSE
               AVCIR=AVCIR1/UR
               WRITE(1700,*) HRZP(1,M),AVCIR
               WRITE(1701,*) HRZP(1,M),AVCIR1
            ENDIF

 310     CONTINUE

 5115    FORMAT(1x,'ZONE T="r/R=',F6.3,'" I=',I3)

      END IF

C-----------------------------------------------------------------------
C     Write the forces
C-----------------------------------------------------------------------
      FNHIS=FN(1:LENCH)//'.ktkq'
      OPEN(710,FILE=FNHIS,STATUS='UNKNOWN')

      WRITE(710,5200)
 5200 FORMAT(1X,'TITLE="Pot. & vis. KT & KQ per blade (cavitating)"')
      WRITE(710,5220)
 5220 FORMAT(1X,'VARIABLES="ANGLE","KT","KQ","KTV","KQV"')
 5260 FORMAT(5(1X,E14.7))
 
C-----------------------------------------------------------------------
C     Write steady cavitating circulation and forces.          JY061300
C-----------------------------------------------------------------------
      IF(ISTEADY.EQ.0) THEN

         WRITE(1700,*) 
     %        'ZONE T="Steady cavitating Circulation"'
         WRITE(1701,*) 
     %        'ZONE T="Steady cavitating Circulation"'

         UR=SQRT(1.+(.7*PI/ADVCO)**2.)
         DO M=1,MR
            AVCIR1 = HUNTPI*DELP(M)/UR
            AVCIR = HUNTPI*DELP(M)
            WRITE(1700,*) HRZP(1,M),AVCIR1
            WRITE(1701,*) HRZP(1,M),AVCIR
         END DO

         CLOSE(1700)

         WRITE(710,5260) TT(1),XKT(1,1),XKT(1,4),XKTV(1,1),XKTV(1,4)
         CLOSE(710)
!Allen Du 01/06/2018 add the option of coupling with cavopt3d 
         opt_kt=XKTV(1,1)*nblade
         opt_kq=XKTV(1,4)*nblade

C-----------------------------------------------------------------------
C      New addition to calculate shaft forces for steady analysis.
C                                                               JY062701
C-----------------------------------------------------------------------
         IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8) THEN
            DO NN=1,NTPREV
               DO KK=1,6
                  XKTV(NN,KK)=FBXPV(KK)
               END DO
            END DO
            CALL PUF3HRM(NBLADE,NTPREV,IWET,XKTV,NSTEP)
         END IF
C-----------------------------------------------------------------------

      ELSE

C-----------------------------------------------------------------------
C     Print the cavitating forces for the last revolution.      JY100798
C-----------------------------------------------------------------------
         AVKT=0
         AVKQ=0
         AVKTV=0
         AVKQV=0
         
         DO 315 I=1,NTPREV
            WRITE(710,5260) TT(I),XKT(I,1),XKT(I,4),XKTV(I,1),XKTV(I,4)

            AVKT=AVKT+XKT(I,1)
            AVKQ=AVKQ+XKT(I,4)
            AVKTV=AVKTV+XKTV(I,1)
            AVKQV=AVKQV+XKTV(I,4) 
 315     CONTINUE

         WRITE(710,*) 'ZONE T="AVERAGEAD KT and KQ"'
         
         ANTPREV=FLOAT(NTPREV)
         AVKQ=AVKQ/ANTPREV
         AVKT=AVKT/ANTPREV
         AVKQV=AVKQV/ANTPREV
         AVKTV=AVKTV/ANTPREV
         
         WRITE(710,5260) TT(1),AVKT,AVKQ,AVKTV,AVKQV
         WRITE(710,5260) TT(NTPREV),AVKT,AVKQ,AVKTV,AVKQV

!Allen Du 01/06/2018 add the option of coupling with cavopt3d 
         opt_kt=AVKTV*nblade
         opt_kq=AVKQV*nblade
         
         CLOSE(710)

C-----------------------------------------------------------------------
C     Creat new force output file, force.cav, which will contain all 
C     six components of the viscous forces at all timesteps in the 
C     last revolution.                                          JY011600
C-----------------------------------------------------------------------
         OPEN(711,FILE='force.cav',STATUS='UNKNOWN')
         WRITE(711,*) 'TITLE="Cavitating 6-compt. Forces per Blade"'
         WRITE(711,*) 'VARIABLES="ANGLE","FX","FY","FZ","MX","MY","MZ"'
 7111    FORMAT(7(1X,E14.7))

         DO I=1,NTPREV
            WRITE(711,7111) TT(I),(XKTV(I,I1),I1=1,6)
         END DO
         CLOSE(711)

C-----------------------------------------------------------------------
C     Calculate the fully wetted blade and shaft harmonics for the last 
C     revolution and write the results to harmny.cav.           JY071399
C-----------------------------------------------------------------------
C....Same for ICON=8 as well (JY110100)
         IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8) 
     *        CALL PUF3HRM(NBLADE,NTPREV,IWET,XKTV,NSTEP)

C-----------------------------------------------------------------------
C     Print max. force coefficients at end of cavitating run.   JY030900
C-----------------------------------------------------------------------
         WRITE(*,*) 
         WRITE(*,5222)
 5222    FORMAT(1X,'------------------',1X,
     *      'Max cavitating force coefficients',1X,'------------------')
         DO L=1,6
            FMAX(L)=0
            FVMAX(L)=0
            DO I=1,NTPREV
               FMAX(L)=MAX(ABS(XKT(I,L)),FMAX(L))
               FVMAX(L)=MAX(ABS(XKTV(I,L)),FVMAX(L))
            END DO
         END DO
         WRITE(*,5225) 'FX','FY','FZ','MX','MY','MZ'
 5225    FORMAT(9X,4X,6(3X,'|',A2,'|',3X))
         WRITE(*,5223) (FMAX(L),L=1,6)
 5223    FORMAT(1X,'Potential',1X,6(1X,F9.5))
         WRITE(*,5224) (FVMAX(L),L=1,6)
 5224    FORMAT(1X,'Viscous  ',1X,6(1X,F9.5))
         WRITE(*,*)

      END IF

C-----------------------------------------------------------------------
C     print CL,CD,CD/CL to screen for hydrofoils.               JY110300
C-----------------------------------------------------------------------
      IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
         WRITE(*,*) '-----------------------------------------'
         WRITE(*,*) 'CAVITATING:   CL        CD        CL/CD  '
         WRITE(*,'(8X,2(F11.6),F11.1)') CLIFT,CDRAG,CLIFT/CDRAG
         WRITE(*,*) '-----------------------------------------'
      END IF

C-----------------------------------------------------------------------
C     close files
C-----------------------------------------------------------------------
      CLOSE(57)
      CLOSE(58)
      CLOSE(630)
      CLOSE(680)
      IF(IHUB.NE.0.AND.IPHUB.EQ.1) CLOSE(18)

      RETURN
      END
