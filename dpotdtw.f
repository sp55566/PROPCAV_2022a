       Subroutine DPOTDTW
************************************************************************
*      Compute the time derivative of the potential on the wake        *
*      using 2nd order central difference approximation.               *
*                                                                      *
*      Author: Julie Young                                             *
*                                                                      *
*      Date           Comments/Revision                                *
*      -----------    ----------------------------------------------   *
*      120798JY       subroutine created                               *
************************************************************************


       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       INCLUDE 'PUFCAVC.INC'

       IF(ISTEADY.NE.0.AND.ISTEADY.NE.1) THEN
          DTN=DELTAT*ADVCO/PI
          IF(NTSTEP.EQ.1) THEN
             DO 10  I=1,NPWAKS
                DO 15 J=1,NTPREV
                   POTWM(J,I)=0.0
                   DPDTPREW(J,I)=0.0
 15             CONTINUE
 10          CONTINUE
          END IF

          NSTART=NTPREV+3

          DO 20 M=1,MR
             DO 25 N=1,NTRA
                I=(MR-M)*NTRA+N
                POTWM(IDXREV,I)=POTW(I)
                IF(NTSTEP.GE.NSTART) THEN
                   IF(IDXREV.EQ.1) THEN
                      N1=NTPREV-1
                      POTP1=(POTWM(IDXREV,I)+POTWM(NTPREV,I)+
     *                     POTWM(NTPREV-1,I))/3.
                      POTM1=(POTWM(NTPREV-1,I)+POTWM(NTPREV-2,I)+
     *                     POTWM(NTPREV-3,I))/3.
                   ELSE IF(IDXREV.EQ.2) THEN
                      N1=NTPREV
                      POTP1=(POTWM(IDXREV,I)+POTWM(IDXREV-1,I)+
     *                     POTWM(NTPREV,I))/3.
                      POTM1=(POTWM(NTPREV,I)+POTWM(NTPREV-1,I)+
     *                     POTWM(NTPREV-2,I))/3.
                   ELSE IF(IDXREV.EQ.3) THEN
                      N1=IDXREV-2
                      POTP1=(POTWM(IDXREV,I)+POTWM(IDXREV-1,I)+
     *                     POTWM(IDXREV-2,I))/3.
                      POTM1=(POTWM(IDXREV-2,I)+POTWM(NTPREV,I)+
     *                     POTWM(NTPREV-1,I))/3.
                   ELSE IF(IDXREV.EQ.4) THEN
                      N1=IDXREV-2
                      POTP1=(POTWM(IDXREV,I)+POTWM(IDXREV-1,I)+
     *                     POTWM(IDXREV-2,I))/3.
                      POTM1=(POTWM(IDXREV-2,I)+POTWM(IDXREV-3,I)+
     *                     POTWM(NTPREV,I))/3.
                   ELSE IF(IDXREV.GE.5) THEN
                      N1=IDXREV-2
                      POTP1=(POTWM(IDXREV,I)+POTWM(IDXREV-1,I)+
     *                     POTWM(IDXREV-2,I))/3.
                      POTM1=(POTWM(IDXREV-2,I)+POTWM(IDXREV-3,I)+
     *                     POTWM(IDXREV-4,I))/3.
                   END IF

                   DPDTPREW(N1,I)=(POTP1-POTM1)/(2.*DTN)

                END IF
 25          CONTINUE
 20       CONTINUE
       END IF

       RETURN
       END
      
