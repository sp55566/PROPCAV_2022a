      SUBROUTINE DETACH
************************************************************************
*                                                                      *
*     This subroutine locates the DETACHment line for the wetted       *
*     pressure distribution.  Created 071797 by CM at the University   *
*     of Texas.                                                        *
*                                ----------            ----------      *
*                               |          |          |          |     *
*                               | PROPCAV  |--------->|  DETACH  |     *
*                               |          |          |          |---> *
*                                ------^---            ----------    | *
*                                      |                             | *
*                                      <-----------------------------v *
*     Date      Revision                                               *
*     --------  ---------                                              *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      ADVCO2=ADVCO*ADVCO

      IF(SIGMA.LE.ZERO) THEN
         SIGMA1=ZERO
      ELSE
         SIGMA1=SIGMA-.005
      END IF

      DO 10 IDR=1,2
         DO 20 M=1,MR
            NLEP(M,IDXREV,IDR)=0
            IF(ISEARCH.EQ.0) THEN
               IF(IDR.EQ.1) THEN
                  INDX=NHP
               ELSE
                  INDX=NH
               END IF
               IF(-CPB(INDX,M)*ADVCO2.LT.SIGMA1) THEN
                  NOCAV(M,IDXREV,IDR)=1
               END IF
            END IF
 20      CONTINUE
 10   CONTINUE

      IF(ISEARCH.EQ.1) THEN

         IF(ISC.EQ.0) THEN
            NLAST=NH-3
         ELSE
            NLAST=NHOLD
         END IF

C-----------------------------------------------------------------------
C     Search for detachment location based on fully wetted pressures on
C     both back and face of the blade.
C
C          Suction Side (IDR=1)
C          Search Direction (K=1)
C          ---->   
C                        *
C                 *
C             *
C          *
C         * <------Starting Point (Leading Edge)    
C          *
C             * 
C                 *
C         ---->          *
C         Search Direction (K=-1)
C         Pressure Side (IDR=2)
C-----------------------------------------------------------------------
         DO 30 IDR=1,2
            IF(IDR.EQ.1)THEN
               K=1
               ISF=0
               IF(ICON.EQ.8) SIGMA1=SIGMAB-.005
            ELSE IF(IDR.EQ.2) THEN
               K=-1
               ISF=1
               IF(ICON.EQ.8) SIGMA1=SIGMAF-.005
            ENDIF

            DO 45 M=1, MR
               DO 50 N=1, NC/2
                  INDX=K*N+NC/2+ISF
                  IF(-CPB(INDX,M)*ADVCO2.GE.SIGMA1) THEN
                     NLEP(M,IDXREV,IDR)=N-1
                     NOCAV(M,IDXREV,IDR)=0
                     GOTO 40
                  ENDIF               
 50            CONTINUE
               NOCAV(M,IDXREV,IDR)=1

 40            CONTINUE

C.............Do not let cavity detach later than 3 panels ahead of 
C.............blade T.E. (JY061400)
               IF(NLEP(M,IDXREV,IDR).GE.NLAST) 
     *              NOCAV(M,IDXREV,IDR)=1
            
               IF(IFACE.EQ.0.AND.IDR.EQ.2) NOCAV(M,IDXREV,IDR)=1
               IF(IFACE.EQ.1.AND.IDR.EQ.1) NOCAV(M,IDXREV,IDR)=1

 45         CONTINUE

 30      CONTINUE
         
      ENDIF

C-----------------------------------------------------------------------
C     Make sure that there's at least SIX panel between cavity L.E. on
C     back and cavity L.E. on face.
C-----------------------------------------------------------------------
      IF(IFACE.EQ.2) THEN
         DO M=1,MR
            IF(NOCAV(M,IDXREV,1).EQ.0.AND.
     *           NOCAV(M,IDXREV,2).EQ.0) THEN
               IWLE=NLEP(M,IDXREV,1)+NLEP(M,IDXREV,2)
               IF(IWLE.LT.6) THEN
                  IF(NLEP(M,IDXREV,1).GE.NLEP(M,IDXREV,2)) THEN
                     NLEP(M,IDXREV,2)=6-NLEP(M,IDXREV,1) 
                  ELSE
                     NLEP(M,IDXREV,1)=6-NLEP(M,IDXREV,2)
                  END IF

                  DO 55 IDR=1,2
                     DO N=NLEP(M,IDXREV,IDR)+1,NLAST
                        IF(IDR.EQ.1) THEN
                           INDX=NH+N
                        ELSE
                           INDX=NHP-N
                        END IF
                        IF(-CPB(INDX,M)*ADVCO2.GE.SIGMA1) THEN
                           NLEP(M,IDXREV,IDR)=N-1
                           GO TO 55
                        END IF
                     END DO
                     NOCAV(M,IDXREV,IDR)=1                     
 55               CONTINUE

               END IF
            END IF
         END DO
      END IF

C-----------------------------------------------------------------------
C     For ISC=1, set the detachment point equal to the trailing edge of
C     the original blade.                                       JY070401
C-----------------------------------------------------------------------
      IF(ISC.EQ.1) THEN
         DO M=1,MR
            DO IDR=1,2
               IF(NLEP(M,IDXREV,IDR).GE.NHOLD) THEN
                  NOCAV(M,IDXREV,IDR)=1
               END IF
               IF(NOCAV(M,IDXREV,IDR).EQ.1) THEN
                  NLEP(M,IDXREV,IDR)=NHOLD
               END IF
            END DO
         END DO
      END IF

      RETURN
C>>>>>>>>>>>>>>>>>>>>>>>END OF SUBROUTINE DETACH<<<<<<<<<<<<<<<<<<<<<<<<
      END













