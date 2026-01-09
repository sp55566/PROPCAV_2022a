      SUBROUTINE DETACH2
************************************************************************
*                                                                      *
*     This subroutine iterates to find the 2nd, (&3rd, etc).           *
*     DETACHment line after the first cavitating revolution.           *
*     Created by CM on 071797 at the University of                     *
*     of Texas.                                                        *
*                                ----------            ----------      *
*                               |          |          |          |     *
*                               |  CAVOUT  |--------->|  DETACH2 |     *
*                               |          |          |          |---> *
*                                ------^---            ----------    | *
*                                      |                             | *
*                                      <-----------------------------v *
*     Date      Revision                                               *
*     --------  ---------                                              *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

C-----------------------------------------------------------------------
C     Write the detachment line at the present time step (before
C     it is modified) to the file, .det. 
C-----------------------------------------------------------------------
      IF(ISC .NE.1) THEN 
         WRITE(720,5030) NTSTEP, MR
         DO 15 M=1, MR
            WRITE(720,5060) M,NLEP(M,IDXREV,1),CAVL(M,1), 
     *NOCAV(M,IDXREV,1),NLEP(M,IDXREV,2),CAVL(M,2),NOCAV(M,IDXREV,2)
 15      CONTINUE
      ENDIF

      ADVCO2=ADVCO*ADVCO

      IF(SIGMA.LE.ZERO) THEN
         SIGMA1=ZERO
      ELSE
         SIGMA1=SIGMA-.005
      END IF

      IF(ISEARCH.EQ.0) THEN
         DO IDR=1,2
            DO M=1,MR
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
            END DO
         END DO
         GO TO 1000
      END IF

      IF(ISC.EQ.0) THEN
         NLAST=NH-3
      ELSE
         NLAST=NHOLD
      END IF

      DO 10 IDR=1,2

C-----------------------------------------------------------------------
C     Search for NLE, the new location where the -CPB>SIGMA.
C-----------------------------------------------------------------------
         NLE=0
         IF(IDR.EQ.1)THEN
            K=1
            ISF=0
            IF(ICON.EQ.8) SIGMA1=SIGMAB-.005
         ELSE IF(IDR.EQ.2) THEN
            K=-1
            ISF=1
            IF(ICON.EQ.8) SIGMA1=SIGMAF-.005
         ENDIF

         DO 20 M=1,MR
            DO 30 N=1,NH
               INDX=K*N+NH+ISF
               IF(-CPB(INDX,M)*ADVCO2.GE.SIGMA1) THEN
                  NLE=N-1
                  NOCAV(M,IDXREV,IDR)=0
                  GOTO 40
               END IF
 30         CONTINUE
            NOCAV(M,IDXREV,IDR)=1
 40         CONTINUE

C-----------------------------------------------------------------------
C     Search criterian 1:
C       If we found no cavity at this step, but the cavitating pressure
C       says we should, then place new detachment point at NLE and
C       try grow the cavity from there.
C-----------------------------------------------------------------------
            IF((JCV(M,IDR).EQ.0).AND.(NLE.GE.0)) THEN
               NLEP(M,IDXREV,IDR)=NLE

C-----------------------------------------------------------------------
C     Search criterian 2:
C       If the initial cavity height is negative, then move the 
C       detachment point one to the right.  If the new detachment point
C       is still less than NLE, then move it one more to the right.
C-----------------------------------------------------------------------
            ELSE IF(HT(2,M,IDR).LT.ZERO.OR.HT(3,M,IDR).LT.ZERO) THEN
               NLEP(M,IDXREV,IDR)=NLEP(M,IDXREV,IDR)+1

C-----------------------------------------------------------------------
C     Search criterian 3:
C       If the cavitating pressure in front of the current detachment 
C       point is
C       greater than SIGMA, then move the detachment point one to the
C       left.  If the new detachment point is still more than NLE, then
C       move it one more to the left.
C-----------------------------------------------------------------------
            ELSE
               INDX=NHP+K*NLEP(M,IDXREV,IDR)-ISF
               IF((-CPB(INDX-K,M)*ADVCO2.GE.SIGMA1).AND.
     *              (NLEP(M,IDXREV,IDR).GT.0)) THEN                  
                  NLEP(M,IDXREV,IDR)=NLEP(M,IDXREV,IDR)-1
               END IF
            END IF

C-----------------------------------------------------------------------
C     The pressure at where the 1st cavitating panel should
C     be greater than SIGMA.
C-----------------------------------------------------------------------
            INDX=NHP+K*NLEP(M,IDXREV,IDR)-ISF
            IF(-CPB(INDX,M)*ADVCO2.LT.SIGMA1) 
     *           NOCAV(M,IDXREV,IDR)=1

C-----------------------------------------------------------------------
C     The cavity cannot detach after the blade T.E.
C-----------------------------------------------------------------------
            IF(NLEP(M,IDXREV,IDR).GE.NLAST) NOCAV(M,IDXREV,IDR)=1

c---- Reason for commenting out the following block:
C---- If we move to the bigger bubble by pushing NLEP toward T.E.,
C---- algorithm ignores cavity detached before around the L.E., which does not make sense.
C---- S.N.KIM 
cC-----------------------------------------------------------------------
cC     Check for the bigger cavity bubble.
cC-----------------------------------------------------------------------
c            IF(JCV(M,IDR).GT.0.AND.NOCAV(M,IDXREV,IDR).EQ.0) THEN
cC............1) before cavity
c                NCC1=0
c                DO N=NLE+1,NLEP(M,IDXREV,IDR)
c                   INDX=K*N+NH+ISF
c                   IF(-CPB(INDX,M)*ADVCO2.GE.SIGMA1) 
c     *                 NCC1=NCC1+1
c                END DO
c                IF(NCC1.GT.JCV(M,IDR)) THEN
c                   NLEP(M,IDXREV,IDR)=NLE
c                   GO TO 1200
c                END IF
c
cC............2) after cavity
c                NCC2=0
c                NLE1=NLAST
c                NLCAV=NLEP(M,IDXREV,IDR)+JCV(M,IDR)+NSPP(M,IDR)+1
c                DO N=NLCAV,NLAST
c                   INDX=K*N+NH+ISF
c                   IF(-CPB(INDX,M)*ADVCO2.GE.SIGMA1) THEN
c                      NCC2=NCC2+1
c                      NLE1=MIN(N-1,NLE1)
c                   END IF
c                END DO
c                IF(NCC2.GT.JCV(M,IDR)) THEN
c                   NLEP(M,IDXREV,IDR)=NLE1
c                END IF
c                
c 1200          CONTINUE
c            END IF

            IF(IFACE.EQ.0.AND.IDR.EQ.2) NOCAV(M,IDXREV,IDR)=1
            IF(IFACE.EQ.1.AND.IDR.EQ.1) NOCAV(M,IDXREV,IDR)=1

 20      CONTINUE
 10   CONTINUE

C-----------------------------------------------------------------------
C     Make sure that there's at least FOUR panel between cavity L.E. on
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

 1000 CONTINUE

C-----------------------------------------------------------------------
C     Write the detachment line at the present time step (after
C     it is modified) to the file, .det. 
C-----------------------------------------------------------------------
      IF(ISC .NE.1) THEN 
         WRITE(720,5040) NTSTEP, MR
         DO 16 M=1, MR
            WRITE(720,5060) M,NLEP(M,IDXREV,1),CAVL(M,1), 
     *NOCAV(M,IDXREV,1),NLEP(M,IDXREV,2),CAVL(M,2),NOCAV(M,IDXREV,2)
 16      CONTINUE
      ENDIF

 5030    FORMAT(1x,'ZONE T="Det Before ',I3,'" I=',I2)      
 5040    FORMAT(1x,'ZONE T="Det After ',I3,'" I=',I2)      
 5060    FORMAT(1X,I2,3X,I2,3X,F10.6,3X,I2,3X,I2,3X,F10.6,3X,I2)

      RETURN
C>>>>>>>>>>>>>>>>>>>>>>>END OF SUBROUTINE DETACH2<<<<<<<<<<<<<<<<<<<<<<<
      END







