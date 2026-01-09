      SUBROUTINE INFGEN
************************************************************************
*     INFGEN: INFluence coefficients matrix GENerator                  *
*     Compute total influence coefficients to establish the matrix     *
*     system                                                           *
*                                                                      *
*  Date of last revision      Revision                                 *
*  ---------------------      --------                                 *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C-----------------------------------------------------------------------
C     Initialization
C-----------------------------------------------------------------------

C.....FILE 50 (IUMAS) -- [A] for the steady problem.....................
C.....FILE 51 (IUMAU) -- [A] for the unsteady problem
C.....FILE 52 (IUMAK) -- [A] for the p.k. iterations

      IUMAS = 50
      IUMAU = 51
      IUMAK = 52

      IIII = NPANB + NPANH + 1

      IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=2

C-----------------------------------------------------------------------
C     Generate coefficients matrices
C-----------------------------------------------------------------------
      REWIND 41
      REWIND 81
      DO 100 M=MR,1,-1

C.......TEMP2 here is wake influence coeff. at L=1, KK=1................
         CALL READ1(81,TEMP2,NPANEL)
         DO 90 N=1,NC
            DO 10 I=1,NPANEL
               A(I)=ZERO
               B(I)=ZERO
 10         CONTINUE
C-----------------------------------------------------------------------
C         KEY blade's influence coefficients
C-----------------------------------------------------------------------
            CALL READ1(41,TEMP1,NPANEL)

            SUMINF=0
            DO 20 I=1,NPANEL
               A(I)=A(I)+TEMP1(I)
 20         CONTINUE

C.........Write matrix [A] without Morino's Kutta condition.............
          CALL WRITE1(IUMAK,TEMP1,NPANEL)

C.........Create [A] for the unsteady problem (with linear correction)..

          IF(N.EQ.1)THEN
             DO 25 I=1,NPANEL  
                TEMP1(I)=TEMP1(I)-HALF*TEMP2(I)-WSUBIF(I,M)
 25          CONTINUE
          ELSEIF(N.EQ.NC)THEN
             DO 30 I=1,NPANEL
                TEMP1(I)=TEMP1(I)+HALF*TEMP2(I)+WSUBIF(I,M)
 30          CONTINUE
          END IF
c          
c.........Tip seperation will be considered later
c.........
          CALL WRITE1(IUMAU,TEMP1,NPANEL)
C-----------------------------------------------------------------------
C         OTHER blades' influence coefficients
C----------------------------------------------------------------------- 

          DO 50 KK=2,NBLADE
             CALL READ1(41,TEMP1,NPANEL)
             DO 40 I=1,NPANEL
                A(I)=A(I)+TEMP1(I)
 40          CONTINUE
 50       CONTINUE

C.........Matrix for Morino's problem...................................
C.........for steady case (matrix is written on file 50, IUMAS).........
          DO 60 I=1,NPANEL

C...........Add the influence of the wake...............................
             IF(N.EQ.1) THEN
                A(I)=A(I)-WSTINF(I,M)
             ELSE IF(N.EQ.NC) THEN
                A(I)=A(I)+WSTINF(I,M)
             END IF
 60       CONTINUE

          CALL WRITE1(IUMAS,A,NPANEL)
 90    CONTINUE
 100  CONTINUE

C-----------------------------------------------------------------------
C     Compute influence coefficients due to the hub
C-----------------------------------------------------------------------

      IF(IHUB.NE.0) THEN
        DO 200 N=1,NHBX
          DO 190 M=1,MHBT
            DO 110 I=1,NPANEL
              A(I)=ZERO
              B(I)=ZERO
  110       CONTINUE   

C-----------------------------------------------------------------------
C           KEY blade's influence coefficients
C-----------------------------------------------------------------------
            CALL READ1(41,TEMP1,NPANEL)
            DO 120 I=1,NPANEL
              A(I)=A(I)+TEMP1(I)
  120       CONTINUE

C-----------------------------------------------------------------------
C     Add the influence of the hub disk to IUMAK & IUMAU.  This is 
C     necessary because otherwise the effect of the hub disk are not
C     included in the wetted nor cavitating solution.           
C-----------------------------------------------------------------------
            IF(N.EQ.1) THEN
               DO I=1,NPANEL
                  TEMP1(I)=TEMP1(I)+CHDK(I)/MHBT
               END DO
            ELSE IF(N.EQ.NHBX) THEN
               DO I=1,NPANEL
                  TEMP1(I)=TEMP1(I)+CHDKDT(I)/MHBT
               END DO
            END IF

C...........Write matrix [A] without Morino's Kutta condition...........
            CALL WRITE1(IUMAK,TEMP1,NPANEL)
            CALL WRITE1(IUMAU,TEMP1,NPANEL)
C-----------------------------------------------------------------------
C           OTHER blades' influence coefficients
C-----------------------------------------------------------------------

            DO 150 KK=2,NBLADE
              CALL READ1(41,TEMP1,NPANEL)
              DO 140 I=1,NPANEL
                A(I)=A(I)+TEMP1(I)
  140         CONTINUE
  150       CONTINUE

C...........Add the influence of the hub disk...........................

            IF(N.EQ.1) THEN
               DO I=1,NPANEL
                  A(I)=A(I)+CHDK(I)/MHBT
               END DO
            ELSE IF(N.EQ.NHBX) THEN
               DO I=1,NPANEL
                  A(I)=A(I)+CHDKDT(I)/MHBT
               END DO
            END IF
            CALL WRITE1(IUMAS,A,NPANEL)
 190     CONTINUE
 200  CONTINUE
      END IF


C-----------------------------------------------------------------------
C     Compute influence coefficients due to the Duct
C-----------------------------------------------------------------------

      IF(IDUCT .NE. 0) THEN

         REWIND 501

         DO 1200 M=1,MDUCT

            CALL READ1(501,TEMP2,NPANEL)

            DO 1190 N=1,NDUCT

               DO 1110 I=1,NPANEL
                  A(I)=ZERO
                  B(I)=ZERO
 1110          CONTINUE   
               
C-----------------------------------------------------------------------
C           KEY blade's influence coefficients
C-----------------------------------------------------------------------
               CALL READ1(41,TEMP1,NPANEL)
               DO 1120 I=1,NPANEL
                  A(I)=A(I)+TEMP1(I)
 1120          CONTINUE
               
               CALL WRITE1(IUMAK,TEMP1,NPANEL)

               IF(N.EQ.1)THEN
                  DO I=1,NPANEL  
                     TEMP1(I)=TEMP1(I)-HALF*TEMP2(I)-WSUBIFD(I,M)
                  ENDDO
               ELSEIF(N.EQ.NDUCT)THEN
                  DO I=1,NPANEL
                     TEMP1(I)=TEMP1(I)+HALF*TEMP2(I)+WSUBIFD(I,M)
                  ENDDO
               ENDIF

               CALL WRITE1(IUMAU,TEMP1,NPANEL)
               
               DO 1150 KK=2,NBLADE
                  CALL READ1(41,TEMP1,NPANEL)
                  DO 1140 I=1,NPANEL
                     A(I)=A(I)+TEMP1(I)
 1140             CONTINUE
 1150          CONTINUE
               
               DO I=1,NPANEL
                  IF(N.EQ.1) THEN
                     A(I)=A(I)-WSTINFD(I,M)
                  ELSE IF(N.EQ.NDUCT) THEN
                     A(I)=A(I)+WSTINFD(I,M)
                  END IF
               ENDDO
               CALL WRITE1(IUMAS,A,NPANEL)
 1190       CONTINUE
 1200    CONTINUE
      END IF
      
C-----------------------------------------------------------------------
C     Compute influence coefficients due to the tunnel
C-----------------------------------------------------------------------

      IF(ITUN.NE.0) THEN
         DO 2001 N=1,NAXT
            DO 1901 M=1,MTUNEL
               DO 1101 I=1,NPANEL
                  A(I)=ZERO
                  B(I)=ZERO
 1101          CONTINUE   
               
C-----------------------------------------------------------------------
C           KEY blade's influence coefficients
C-----------------------------------------------------------------------
               CALL READ1(41,TEMP1,NPANEL)
               DO 1201 I=1,NPANEL
                  A(I)=A(I)+TEMP1(I)
 1201          CONTINUE
               
C...........Write matrix [A] without Morino's Kutta condition...........
               CALL WRITE1(IUMAK,TEMP1,NPANEL)
               CALL WRITE1(IUMAU,TEMP1,NPANEL)
C-----------------------------------------------------------------------
C           OTHER blades' influence coefficients
C-----------------------------------------------------------------------
               DO 1501 KK=2,NBLADE
                  CALL READ1(41,TEMP1,NPANEL)
                  DO 1401 I=1,NPANEL
                     A(I)=A(I)+TEMP1(I)
 1401             CONTINUE
 1501          CONTINUE
               CALL WRITE1(IUMAS,A,NPANEL)
 1901       CONTINUE
 2001    CONTINUE
      END IF

C*********************************************************************
C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.      
C*********************************************************************
cC -- Begin Tip HSLEE(10/12/99)
c
c      if(ian.eq.2) then
c
c           do n = 1 , nthx
c              do m = 1 , mcvt
cc             Hong change the initialzation of A & B        09/13/07              
cc                 a = 0.0           
cc                 b = 0.0
c              DO I=1,NPANEL
c                 A(I)=ZERO
c                 B(I)=ZERO
c              ENDDO  
c
cC --- Read tip hub influence coeff
c      
c               call read1(41,temp1,npanel)
c      
c               do i = 1, npanel
c                   a(i) = a(i) + temp1(i)
c               enddo
c
c               call write1(iumak,temp1,npanel)
c               call write1(iumau,temp1,npanel)
c
c                 call write1(iumas,a,npanel)
c
c             enddo
c          enddo
c
c           do n = 1 , ncvx
c              do m = 1 , mcvt
cc             Hong change the initialzation of A & B        09/13/07              
cc                 a = 0.0           
cc                 b = 0.0
c              DO I=1,NPANEL
c                 A(I)=ZERO
c                 B(I)=ZERO
c              ENDDO  
c
cC --- Read tip hub influence coeff
c
c                 call read1(41,temp1,npanel)
c               do i = 1 , npanel
c                 a(i) = a(i) + temp1(i)
c               enddo
c              
c               call write1(iumak,temp1,npanel)
c               call write1(iumau,temp1,npanel)
c                 call write1(iumas,a,npanel)
c             enddo
c          enddo
c
c      endif
C*********************************************************************
C/e S.N.KIM | Aug. 2108.
C*********************************************************************

C -- End Tip (10/12/99)
CSH--REPLACE, NBLADE=1-------------------------
            IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=1
CSH--------------------------------------------
      RETURN
      END



