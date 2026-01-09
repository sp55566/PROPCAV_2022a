       SUBROUTINE INFGENW
************************************************************************
*      This subroutine calculates influence coefficients for the steady*
*      cavitating case.                                                *
*                                                                      *
*      Date      Comment or Revision                                   *
*      --------  -------------------                                   *
*      JY011100  Subroutine created.                                   *
*                                                                      *
************************************************************************

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       INCLUDE 'PUFCAVC.INC'

C-----------------------------------------------------------------------
C     W(I,M): 
C     Induced potential at blade panel I due to wake dipoles at strip
C     M of all blades (including far wake and wake sink disk). 
C-----------------------------------------------------------------------
       DO M=MR,1,-1
          DO I=1,NPANEL
             W(I,M)=WINF(I,M)
          END DO
       END DO

       IF(IDUCT .NE. 0) THEN
          DO M = 1 , MDUCT
             DO I = 1, NPANEL
                WD(I,M) = WSTINFD(I,M)
             ENDDO
          ENDDO
       ENDIF

C-----------------------------------------------------------------------
C     WK(I,M) & WK2(I,M):
C     Indueced potential at supercavitating wake panel I due to wake 
C     dipoles at strip M of all blades (including far wake and wake
C     sink disk). WK is used when calculating PHI+ and WK2 is used
C     when calculating PHI-.
C-----------------------------------------------------------------------
       DO M=MR,1,-1
          DO I=1,NPWAKS
             WK(I,M)=WUSINF1(I,M)
             WK2(I,M)=WUSINF1(I,M)
          END DO
       END DO

       DO KK=1,NBLADE
          IO1=90+KK
          IO2=30+KK
          REWIND IO1
          REWIND IO2
          DO M=MR,1,-1
             DO L=1,NWMINFW
                IF(L.LE.NSUB)THEN
                   CALL READ1(IO1,TEMP4,NPWAKS)
                ELSE
                   CALL READ1(IO2,TEMP4,NPWAKS)
                ENDIF
                DO I=1,NPWAKS
                   WK(I,M)=WK(I,M)+TEMP4(I)
                   IF(KK.EQ.1.AND.L.LE.NSUB) THEN
                      WK2(I,M)=WK2(I,M)+WKFACE(I,M,L)
                   ELSE
                      WK2(I,M)=WK2(I,M)+TEMP4(I)
                   END IF
                END DO
             END DO
          END DO
       END DO

C
C     Duct wake

       IF(IDUCT .NE. 0) THEN
          DO KK = 1 , NBLADE
             IO = 520+KK
             REWIND IO
             DO N = 1, NDWK
                DO M = 1, MDUCT
                   CALL READ1(IO,TEMP4,NPWAKS)
                   DO I = 1, NPWAKS
                      WKD(I,M) = WKD(I,M) + TEMP4(I)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF

       
CSH--Calculate the influence coeff from the image, so NBLADE=2
       IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=2
CSH--------------------------------------------           
C-----------------------------------------------------------------------
C     AA(I,J) & BB(I,J):
C     Induced potential at blade panel I due to dipoles and sources
C     in blade panel J of all blades, respectively.
C
C     D(I,J) & E(I,J):
C     Induced potential at wake panel I due to dipoles and sources
C     in blade panel J of all blades, respectively.
C-----------------------------------------------------------------------
       DO J=1,NPANEL
          DO I=1,NPANEL
             AA(I,J)=ZERO
             BB(I,J)=ZERO
          END DO
          DO I=1,NPWAKS
             D(I,J)=ZERO
             E(I,J)=ZERO
          END DO
       END DO
       
       REWIND 42
       REWIND 41
       REWIND 111
       REWIND 112

C -- For Blade

       III1 = 1
       III2 = NPANB

CSH--FOR, NBLADE=2-----------------------------
       IF((ICON.EQ.5).AND.(IHULL.EQ.1)) THEN
CSH--------------------------------------------
          DO J=III1,III2
             DO KK=1,NBLADE
                CALL READ1(41,TEMP1,NPANEL)
                CALL READ1(42,TEMP2,NPANEL)
                DO I=1,NPANEL
                   AA(I,J)=AA(I,J)+TEMP1(I)
                   BB(I,J)=BB(I,J)+TEMP2(I)
                END DO
             END DO
          END DO
C     
          NBLADE=1
          DO J=1,NPANB
             DO KK=1,NBLADE
                CALL READ1(111,TEMP4,NPWAKS)
                CALL READ1(112,TEMP5,NPWAKS)
                DO I=1,NPWAKS
                   D(I,J)=D(I,J)+TEMP4(I)
                   E(I,J)=E(I,J)+TEMP5(I)
             END DO
          END DO
       END DO
       
      ELSE

         DO J=III1,III2
            DO KK=1,NBLADE
               CALL READ1(41,TEMP1,NPANEL)
               CALL READ1(42,TEMP2,NPANEL)
               DO I=1,NPANEL
                  AA(I,J)=AA(I,J)+TEMP1(I)
                  BB(I,J)=BB(I,J)+TEMP2(I)
               END DO
               CALL READ1(111,TEMP4,NPWAKS)
               CALL READ1(112,TEMP5,NPWAKS)
               DO I=1,NPWAKS
                  D(I,J)=D(I,J)+TEMP4(I)
                  E(I,J)=E(I,J)+TEMP5(I)
               END DO
            END DO
         END DO
      ENDIF   

CSH--REPLACE, NBLADE=1-------------------------
      IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=1
CSH--------------------------------------------
C -- For  Hub

      IF(IHUB .NE. 0) THEN
         III1 = III2 + 1
         III2 = III2 + NPANH
         DO J=III1, III2
            DO KK=1,NBLADE
               CALL READ1(41,TEMP1,NPANEL)
               CALL READ1(42,TEMP2,NPANEL)
               DO I=1,NPANEL
                  AA(I,J)=AA(I,J)+TEMP1(I)
                  BB(I,J)=BB(I,J)+TEMP2(I)
               END DO
               CALL READ1(111,TEMP4,NPWAKS)
               CALL READ1(112,TEMP5,NPWAKS)
               DO I=1,NPWAKS
                  D(I,J)=D(I,J)+TEMP4(I)
                  E(I,J)=E(I,J)+TEMP5(I)
               END DO
            END DO
         END DO
      ENDIF
    
C -- For Duct
  
      IF(IDUCT .NE. 0) THEN
         III1 = III2 + 1
         III2 = III2 + NPAND
         DO J=III1, III2
            DO KK=1,NBLADE
               CALL READ1(41,TEMP1,NPANEL)
               CALL READ1(42,TEMP2,NPANEL)
               DO I=1,NPANEL
                  AA(I,J)=AA(I,J)+TEMP1(I)
                  BB(I,J)=BB(I,J)+TEMP2(I)
               END DO
               CALL READ1(111,TEMP4,NPWAKS)
               CALL READ1(112,TEMP5,NPWAKS)
               DO I=1,NPWAKS
                  D(I,J)=D(I,J)+TEMP4(I)
                  E(I,J)=E(I,J)+TEMP5(I)
               END DO
            END DO
         END DO
      ENDIF

C - For Tunnel

      IF(ITUN .NE. 0) THEN
         III1 = III2 + 1
         III2 = III2 + NPANTN
         DO J=III1, III2
            DO KK=1,NBLADE
               CALL READ1(41,TEMP1,NPANEL)
               CALL READ1(42,TEMP2,NPANEL)
               DO I=1,NPANEL
                  AA(I,J)=AA(I,J)+TEMP1(I)
                  BB(I,J)=BB(I,J)+TEMP2(I)
               END DO
               CALL READ1(111,TEMP4,NPWAKS)
               CALL READ1(112,TEMP5,NPWAKS)
               DO I=1,NPWAKS
                  D(I,J)=D(I,J)+TEMP4(I)
                  E(I,J)=E(I,J)+TEMP5(I)
               END DO
            END DO
         END DO
      ENDIF
      
C -- For Tip Vortex

      IF(IAN .EQ. 2) THEN
         III1 = III2 + 1
         III2 = NPANEL
         DO J=III1, III2
            CALL READ1(41,TEMP1,NPANEL)
            CALL READ1(42,TEMP2,NPANEL)
            DO I=1,NPANEL
               AA(I,J)=AA(I,J)+TEMP1(I)
               BB(I,J)=BB(I,J)+TEMP2(I)
            END DO
            DO KK = 1 , NBLADE
               CALL READ1(111,TEMP4,NPWAKS)
               CALL READ1(112,TEMP5,NPWAKS)
               DO I=1,NPWAKS
                  D(I,J)=D(I,J)+TEMP4(I)
                  E(I,J)=E(I,J)+TEMP5(I)
               END DO
            ENDDO
         END DO
      ENDIF
       

C-----------------------------------------------------------------------
C     Add the influence of the hub dipole disk to D(I,J).
C
C     Modified IF statements to accomodate new hub options.     JY052401
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
      
C-----------------------------------------------------------------------
C     C(I,J):
C     Induced potential at blade panel I due to sources in wake panel 
C     J of all blades.
C
C     F(I,J):
C     Induced potential at wake panel I due to sources in wake panel J 
C     of all blades.
C-----------------------------------------------------------------------

      DO J=1,NPWAKS
         DO I=1,NPANEL
            C(I,J)=ZERO
         END DO
         DO I=1,NPWAKS
            F(I,J)=ZERO
         END DO
      END DO
      
      REWIND 110
      REWIND 113
      
      DO J=1,NPWAKS
         DO KK=1,NBLADE
            CALL READ1(110,TEMP1,NPANEL)
            CALL READ1(113,TEMP4,NPWAKS)
            DO I=1,NPANEL
               C(I,J)=C(I,J)+TEMP1(I)
            END DO
            DO I=1,NPWAKS
               F(I,J)=F(I,J)+TEMP4(I)
            END DO
         END DO
      END DO

C-----------------------------------------------------------------------
C     Initialize the wake source strength with zero.
C-----------------------------------------------------------------------
      DO I=1,NPWAKS
         TEMP4(I)=ZERO
      END DO
      
       NREC = 360 / NDLTAT
       DO N=1,NREC
          CALL WRITE2(48,N,TEMP4,NPWAKS)
       END DO
       
       RETURN
       END
