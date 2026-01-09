       SUBROUTINE ITER2
************************************************************************
*                                                                      *
*      This subroutine iterates within each cavitating iteration to    *
*      determine the correct DPDVB.                                    *
*                                                                      *
*      Author: Julie Young                                             *
*                                                                      *
*      Date        Revision/comments                                   *
*      --------    ---------------                                     *
*      JY080600    Subroutine created.                                 *
*                                                                      *
************************************************************************


       INCLUDE 'PUFCAV.INC'
       DIMENSION DPDVBOLD(NBZ,MBZ)

       IT2MAX=10

       IF(IT2MAX.EQ.0) GO TO 1000

       DO IT2=1,IT2MAX

C-----------------------------------------------------------------------
C      determine PHI1 for the dynamic boundary conditio
C-----------------------------------------------------------------------
          ISR=1
          IF(IFACE.EQ.2) ISR=2
          
          DO M=1,MR         
             DO I=1,ISR
                IF((IFACE.EQ.0).OR.(I.EQ.1.AND.IFACE.EQ.2)) THEN
                   IDR=1
                ELSE
                   IDR=2
                END IF
                
                IF(JCV(M,IDR).GT.0) THEN
                   CALL COMPPHI1(M,IDR)
                   IF(NSPP(M,IDR).EQ.1) CALL DPDVEXTRAP(M,IDR)
                END IF

             END DO
          END DO
          
C-----------------------------------------------------------------------
C     set up the linear system
C-----------------------------------------------------------------------
          CALL CAVSET
      
C-----------------------------------------------------------------------
C     solve the linear system
C-----------------------------------------------------------------------
          CALL CAVSOL

C-----------------------------------------------------------------------
C     store old DPDVB before it's recalculated.
C-----------------------------------------------------------------------
          DO M=1,MR
             DO N=1,NC
                DPDVBOLD(N,M)=DPDVB(N,M)
             END DO
          END DO

C-----------------------------------------------------------------------
C     compute the cavity pressure
C-----------------------------------------------------------------------
          CALL PRSDIF(1)

C-----------------------------------------------------------------------
C     check convergence of DPDVB
C-----------------------------------------------------------------------
          DMAX=ZERO
          DO M=1,MRTIP
             DO N=1,NC                
                DIFF=ABS(DPDVBOLD(N,M)-DPDVB(N,M))
                IF(DIFF.GT.DMAX) THEN
                   DMAX=DIFF
                   NN1=N
                   MM1=M
                END IF
             END DO
          END DO

C          WRITE(998,*) ITER,IT2,DMAX,NN1,MM1

          IF(DMAX.LE.1.E-4) GO TO 1000
       END DO

 1000  CONTINUE

       RETURN
       END
