       SUBROUTINE ITER1_SP

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       INCLUDE 'PUFCAVC.INC'
       DIMENSION NLEPOLD(MBZ)

C.....Iterate to find the correct cavity detachment location.
       DO ITER=1,ITERMAX
          
C........Calculate arclength
          CALL ARCL

C........Re-set the source strengths to the wetted source strengths
          DO J=1,NPANEL
             IF(ISUBM(J,IDXREV).EQ.1) THEN
                DPDNC(J)=DPDN(J)
             ELSE
                DPDNC(J)=ZERO
             END IF
          END DO

C........Re-set the wake source strengths to zero. 
          DO J=1,NPWAKS
             SORW(J)=ZERO
          END DO

C........Det. cavity detachment location
          IF(ITER.EQ.1.OR.NREV.EQ.1) THEN

C...........Initial guess of detachment location
             DO M=1,MR

C..............Face side
                IF(ICB(M,2,IDXREV).EQ.1) THEN
                   NLEP(M,IDXREV,2)=NH-IW(2,M,IDXREV)
                   NOCAV(M,IDXREV,2)=1
                END IF

C..............Back side
                IF(ICB(M,1,IDXREV).EQ.1) THEN
                   NMIN=IC(1,M,IDXREV)-NHP
                   NLEP(M,IDXREV,1)=NMIN
                   NOCAV(M,IDXREV,1)=0
                END IF

C..............Don't allow cavity to detach aft of actual blade
                IF(ISC.EQ.1) THEN
                   DO IDR=1,2
                      IF(NLEP(M,IDXREV,IDR).GE.NHOLD) THEN
                         NLEP(M,IDXREV,IDR)=NHOLD
                         NOCAV(M,IDXREV,IDR)=1
                      END IF
                   END DO
                END IF

             END DO
          ELSE

C...........Store detachment location of previous iteration
             DO M=1,MR
                NLEPOLD(M)=NLEP(M,IDXREV,1)
             END DO
             
C...........Det. new detachment location.
             CALL DETACH_SP

C...........Count the number of strips that needed changed in 
C...........cavity detachment points.
             NMOVE=0
             DO M=1,MRTIP
                IF(NLEPOLD(M).NE.NLEP(M,IDXREV,1)) NMOVE=NMOVE+1
             END DO

          END IF
          
          IF(ISC.EQ.0) THEN
             CALL DELR
          ELSE
             CALL DELR_SC
          END IF

C-----------------------------------------------------------------------
C     Let the cavity length and cavity height of the last strip (MR)
C     equal to the previous strip (MR-1).                       JY092799
C-----------------------------------------------------------------------
          DO M=MRTIP+1,MR
             IF(JCV(M,1).GT.0.AND.ICB(M,1,IDXREV).EQ.1) THEN
                HT(1,M,1)=ZERO
                NLAST=IC(2,M,IDXREV)-(NH+NLEP(M,IDXREV,1))+NNWC(M)
                DO N=2,NLAST+1
                   HT(N,M,1)=HT(N,MRTIP,1)
                END DO
             END IF
             IF(ICB(M,2,IDXREV).EQ.1.AND.
     *            IW(1,M,IDXREV).LT.N0(2).AND.
     *            IW(2,M,IDXREV).GE.N0(2)) THEN
                HT(1,M,2)=ZERO
                NLAST=NSR2
                DO N=2,NLAST+1
                   HT(N,M,2)=HT(N,MRTIP,2)
                END DO
             END IF
          END DO
C-----------------------------------------------------------------------

C........Exit do loop if the detachment location remains constant
C........for most of the strips.
          IF(ITER.GT.1.AND.NMOVE.LE.1) GO TO 1000

       END DO

 1000  CONTINUE

       RETURN
       END
