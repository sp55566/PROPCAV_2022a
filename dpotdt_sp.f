       Subroutine DPOTDT_SP
************************************************************************
*      Compute the time derivative of the potential on the blade for   *
*      ISP=1.                                                          *
*                                                                      *
*      Author: Julie Young                                             *
*                                                                      *
*      Date           Comments/Revision                                *
*      -----------    ----------------------------------------------   *
*      JY072500       Subroutine created.                              *
*                                                                      *
************************************************************************
    
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

       double precision XX(200),YY(200),B1(3),XGloc1

       IF(IDXREV.NE.0.AND.ISTEADY.NE.0.AND.ISTEADY.NE.1) THEN
          DTN=DELTAT*ADVCO/PI
          IF(NTSTEP.EQ.1) THEN
             DO 10  I=1,NPANEL
                DO 15 J=1,NTPREV
                   POTM(J,I)=0.0
                   DPDTPRE(J,I)=0.0
 15             CONTINUE
 10          CONTINUE
          END IF  

          NSTART=1
          IF(NBLADE.GT.1) NSTART=2

          DO 20 I=1,NPANEL             
             POTM(IDXREV,I)=POT(I)
             IF(NREV.GT.NSTART.AND.IDXREV.EQ.NTPREV) THEN

C..............Find the first time that panel I is submerged -> ISTART
C..............Find the last time that panel I is submerged -> IEND
                ISTART=0
                IEND=0
                DO K=1,NTPREV
                   IF(K.EQ.1.AND.ISUBM(I,K).EQ.1) THEN
                      ISTART=K
                   ELSE IF(K.GT.1) THEN
                      IF(ISUBM(I,K).EQ.1.AND.ISUBM(I,K-1).EQ.0) THEN
                         ISTART=K
                      ELSE IF(ISUBM(I,K).EQ.0.AND.
     *                        ISUBM(I,K-1).EQ.1) THEN
                         IEND=K-1
                      END IF
                   END IF
                END DO

                IF(ISTART.EQ.0.OR.IEND.EQ.0) GO TO 20
                
C..............Total number of times panel I is submerged.
                ITOTAL=IEND-ISTART+1

C..............Store all the submerged solution in XX and YY
                DO K=1,5
                   XX(K)=DBLE(ISTART-6+K)*DBLE(DTN)
                   YY(K)=0D0
                END DO
                DO K=ISTART,IEND
                   XX(K-ISTART+6)=DBLE(K)*DBLE(DTN)
                   YY(K-ISTART+6)=DBLE(POTM(K,I))
                END DO
                DO K=1,5
                   XX(ITOTAL+5+K)=DBLE(IEND+K)*DBLE(DTN)
                   YY(ITOTAL+5+K)=0D0
                END DO
                ITOTAL=ITOTAL+10

C..............If statement added to prevent out-of-bound errors.
C..............(JY101200)
                IF(ITOTAL.GT.5000) THEN
                   WRITE(*,*)'PLEASE MAKE NMAX > ',ITOTAL,' IN MVLS1D!'
                   WRITE(*,*)'STOP PROGRAM'
                   STOP
                END IF  

C..............Use 2nd order moving least square method to calculate
C..............the time derivative of the potential.
                DO K=ISTART,IEND
                   XGloc1=DBLE(K)*DBLE(DTN)
                   CALL MVLS1D_2(ITOTAL,XX,YY,XGloc1,B1)
                   DPDTPRE(K,I)=SNGL(B1(2))
                END DO

             END IF
 20       CONTINUE
       END IF

       RETURN
       END
      

