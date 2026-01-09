       SUBROUTINE DETACH_SP
************************************************************************
*                                                                      *
*  This subroutine locates the DETACHment line for the ISP=1.          *
*                                                                      *
*  Author: Julie Young  07-17-00                                       *
************************************************************************

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

       ADVCO2=ADVCO*ADVCO
       SIGMA1=SIGMA-0.002

       DO M=1,MR
C........Face side (Do not allow any cavity to grow on this side)
          IF(ICB(M,2,IDXREV).EQ.1) THEN
             NLEP(M,IDXREV,2)=NH-IW(2,M,IDXREV)
             NOCAV(M,IDXREV,2)=1
          END IF

C........Back side (Search for cavity detachment point)
          IF(ICB(M,1,IDXREV).EQ.1) THEN
             NMIN=IC(1,M,IDXREV)-NHP
             IF(ISEARCH.EQ.0) THEN
                NLEP(M,IDXREV,1)=NMIN
                NOCAV(M,IDXREV,1)=0
                GO TO 100
             END IF
     
             IF(IC(2,M,IDXREV).LT.N0(1).OR.ISC.EQ.0) THEN
                NMAX=IC(2,M,IDXREV)-NH
             ELSE
                NMAX=N0(1)-NHP
             END IF
             NLE=NMIN
             
             DO N=NMIN+1,NMAX
                INDX=NH+N
                IF(-CPB(INDX,M)*ADVCO2.GE.SIGMA1) THEN
                   NLE=N
                   NOCAV(M,IDXREV,1)=0
                   GO TO 20
                END IF
             END DO
             NOCAV(M,IDXREV,1)=1
 20          CONTINUE

             IF(HT(2,M,1).LT.ZERO.OR.HT(3,M,1).LT.ZERO) THEN
                NLEP(M,IDXREV,1)=MIN(NLEP(M,IDXREV,1)+1,NMAX)
             ELSE IF(JCV(M,1).EQ.0) THEN
                NLEP(M,IDXREV,1)=NLE
             ELSE IF(NLEP(M,IDXREV,1).GT.NMIN) THEN
                NDUM=NLEP(M,IDXREV,1)-1
                INDX=NHP+NDUM
                IF(-CPB(INDX,M)*ADVCO2.GE.SIGMA1) 
     *               NLEP(M,IDXREV,1)=NDUM
             END IF

             INDX=NHP+NLEP(M,IDXREV,1)
             IF(-CPB(INDX,M)*ADVCO.LT.SIGMA1) NOCAV(M,IDXREV,1)=1

             IF(NLEP(M,IDXREV,1).GE.NMAX) NOCAV(M,IDXREV,1)=1
             IF(NOCAV(M,IDXREV,1).EQ.1) THEN
                NLEP(M,IDXREV,1)=NMAX
             END IF

 100         CONTINUE
          END IF
       END DO

       RETURN
       END
