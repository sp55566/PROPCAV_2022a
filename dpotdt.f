       Subroutine DPOTDT
************************************************************************
*      Compute the time derivative of the potential at the end of each *
*      revolution starting form the 2nd cavitating revolution using    *
*      2nd order central difference approximation.                     *
*                                                                      *
*      Author: Julie Young                                             *
*                                                                      *
*      Date           Comments/Revision                                *
*      -----------    ----------------------------------------------   *
*      102198JY       subroutine created                               *
*      052699JY       subroutine modified to include potentials on     *
*                     the hub.                                         *
************************************************************************


       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

C -- Call the unsteady version when IUNS=1  ---- Yiran Su 09/20/2017
       IF(IUNS.EQ.1) THEN
         CALL DPOTDT_UNS
         RETURN
       END IF

       IF(ISTEADY.NE.0.AND.ISTEADY.NE.1) THEN
          DTN=DELTAT*ADVCO/PI
          IF(NTSTEP.EQ.1) THEN
             DO I=1,NPANEL
                DO J=1,NTPREV
                   POTM(J,I)=0.0
                   DPDTPRE(J,I)=0.0
                END DO
             END DO
          END IF

          NSTART=NTPREV+3
          IF(IFACE.EQ.2.AND.IWET.EQ.0) NSTART=3*NTPREV+3
          IF(ISC.EQ.1) NSTART=3*NTPREV+3

          DO I=1,NPANEL
             POTM(IDXREV,I)=POT(I)
             IF(NTSTEP.GE.NSTART) THEN

cc Old Method Start - Using 3 Points
c                IF(IDXREV.EQ.1) THEN
c                   N1=NTPREV-1
c                   POTP1=(POTM(IDXREV,I)+POTM(NTPREV,I)+
c     *                  POTM(NTPREV-1,I))/3.
c                   POTM1=(POTM(NTPREV-1,I)+POTM(NTPREV-2,I)+
c     *                  POTM(NTPREV-3,I))/3.
c                ELSE IF(IDXREV.EQ.2) THEN
c                   N1=NTPREV
c                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
c     *                  POTM(NTPREV,I))/3.
c                   POTM1=(POTM(NTPREV,I)+POTM(NTPREV-1,I)+
c     *                  POTM(NTPREV-2,I))/3.
c                ELSE IF(IDXREV.EQ.3) THEN
c                   N1=IDXREV-2
c                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
c     *                  POTM(IDXREV-2,I))/3.
c                   POTM1=(POTM(IDXREV-2,I)+POTM(NTPREV,I)+
c     *                  POTM(NTPREV-1,I))/3.
c                ELSE IF(IDXREV.EQ.4) THEN
c                   N1=IDXREV-2
c                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
c     *                  POTM(IDXREV-2,I))/3.
c                   POTM1=(POTM(IDXREV-2,I)+POTM(IDXREV-3,I)+
c     *                  POTM(NTPREV,I))/3.
c                ELSE IF(IDXREV.GE.5) THEN
c                   N1=IDXREV-2
c                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
c     *                  POTM(IDXREV-2,I))/3.
c                   POTM1=(POTM(IDXREV-2,I)+POTM(IDXREV-3,I)+
c     *                  POTM(IDXREV-4,I))/3.
c                END IF
c
c                DPDTPRE(N1,I)=(POTP1-POTM1)/(2.*DTN)
cc Old Method End.

c New Method Start - Using 5 Points
                IF(IDXREV.EQ.1) THEN
                   N1=NTPREV-3
                   POTP1=(POTM(IDXREV,I)+POTM(NTPREV,I)+
     *                  POTM(NTPREV-1,I)+POTM(NTPREV-2,I)+
     *                  POTM(NTPREV-3,I))/5.
                   POTM1=(POTM(NTPREV-3,I)+POTM(NTPREV-4,I)+
     *                  POTM(NTPREV-5,I)+
     *                  POTM(NTPREV-6,I)+POTM(NTPREV-7,I))/5.
                ELSE IF(IDXREV.EQ.2) THEN
                   N1=NTPREV-2
                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
     *                  POTM(NTPREV,I)+POTM(NTPREV-1,I)+
     *                  POTM(NTPREV-2,I))/5.
                   POTM1=(POTM(NTPREV-2,I)+POTM(NTPREV-3,I)+
     *                  POTM(NTPREV-4,I)+
     *                  POTM(NTPREV-5,I)+POTM(NTPREV-6,I))/5.
                ELSE IF(IDXREV.EQ.3) THEN
                   N1=NTPREV-1
                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
     *                  POTM(IDXREV-2,I)+POTM(NTPREV,I)+
     *                  POTM(NTPREV-1,I))/5.
                   POTM1=(POTM(NTPREV-1,I)+POTM(NTPREV-2,I)+
     *                  POTM(NTPREV-3,I)+
     *                  POTM(NTPREV-4,I)+POTM(NTPREV-5,I))/5.
                ELSE IF(IDXREV.EQ.4) THEN
                   N1=NTPREV
                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
     *                  POTM(IDXREV-2,I)+POTM(IDXREV-3,I)+
     *                  POTM(NTPREV,I))/5.
                   POTM1=(POTM(NTPREV,I)+POTM(NTPREV-1,I)+
     *                  POTM(NTPREV-2,I)+
     *                  POTM(NTPREV-3,I)+POTM(NTPREV-4,I))/5.
                ELSE IF(IDXREV.EQ.5) THEN
                   N1=IDXREV-4
                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
     *                  POTM(IDXREV-2,I)+POTM(IDXREV-3,I)+
     *                  POTM(IDXREV-4,I))/5.
                   POTM1=(POTM(IDXREV-4,I)+POTM(NTPREV,I)+
     *                  POTM(NTPREV-1,I)+
     *                  POTM(NTPREV-2,I)+POTM(NTPREV-3,I))/5.
                ELSE IF(IDXREV.EQ.6) THEN
                   N1=IDXREV-4
                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
     *                  POTM(IDXREV-2,I)+POTM(IDXREV-3,I)+
     *                  POTM(IDXREV-4,I))/5.
                   POTM1=(POTM(IDXREV-4,I)+POTM(IDXREV-5,I)+
     *                  POTM(NTPREV,I)+
     *                  POTM(NTPREV-1,I)+POTM(NTPREV-2,I))/5.
                ELSE IF(IDXREV.EQ.7) THEN
                   N1=IDXREV-4
                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
     *                  POTM(IDXREV-2,I)+POTM(IDXREV-3,I)+
     *                  POTM(IDXREV-4,I))/5.
                   POTM1=(POTM(IDXREV-4,I)+POTM(IDXREV-5,I)+
     *                  POTM(IDXREV-6,I)+
     *                  POTM(NTPREV,I)+POTM(NTPREV-1,I))/5.
                ELSE IF(IDXREV.EQ.8) THEN
                   N1=IDXREV-4
                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
     *                  POTM(IDXREV-2,I)+POTM(IDXREV-3,I)+
     *                  POTM(IDXREV-4,I))/5.
                   POTM1=(POTM(IDXREV-4,I)+POTM(IDXREV-5,I)+
     *                  POTM(IDXREV-6,I)+
     *                  POTM(IDXREV-7,I)+POTM(NTPREV,I))/5.
                ELSE IF(IDXREV.GE.9) THEN
                   N1=IDXREV-4
                   POTP1=(POTM(IDXREV,I)+POTM(IDXREV-1,I)+
     *                  POTM(IDXREV-2,I)+POTM(IDXREV-3,I)+
     *                  POTM(IDXREV-4,I))/5.
                   POTM1=(POTM(IDXREV-4,I)+POTM(IDXREV-5,I)+
     *                  POTM(IDXREV-6,I)+
     *                  POTM(IDXREV-7,I)+POTM(IDXREV-8,I))/5.
                END IF

                DPDTPRE(N1,I)=(POTP1-POTM1)/(4.*DTN)
c New Method End.

             END IF

          END DO
       END IF

       RETURN
       END

      SUBROUTINE DPOTDT_UNS
************************************************************************
*      In unsteady application, we can no longer use the potential     *
*      the previous revolution in calculating the current dPOTdT.      *
*      A 2nd order backward approximation can be used in this case.    *
*                                                                      *
*      Author: Yiran Su                                                *
*                                                                      *
*      Date           Comments/Revision                                *
*      -----------    ----------------------------------------------   *
*      09/20/2017     subroutine created                               *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      DTN=DELTAT*ADVCO/PI

C Initialize at the first unsteady step
      IF(NTSTEP.EQ.1) THEN
        DO I=1,NPANEL
          DO J=1,NTPREV
            POTM(J,I)=0.0
            DPDTPRE(J,I)=0.0
          END DO
        END DO
      END IF

C when to start calculating dPOTdT
      NSTART=NTPREV/NBLADE

      DO I=1,NPANEL

        POTM(IDXREV,I)=POT(I)  ! Store data

        IF(NTSTEP_KEY.LE.NSTART) THEN ! Start at certain time step
          DPDTPRE(IDXREV,I) = 0.0
          CYCLE
        END IF

        J0 = IDXREV
        J1 = J0 - 1
        IF (J1.LE.0) J1=J1+NTPREV
        J2 = J0 - 2
        IF (J2.LE.0) J2=J2+NTPREV
        J3 = J0 - 3
        IF (J3.LE.0) J3=J3+NTPREV

        POTP0 = POTM(J0,I)
        POTP1 = (POTM(J0,I)+POTM(J1,I)+POTM(J2,I))/3.0
        POTP2 = (POTM(J1,I)+POTM(J2,I)+POTM(J3,I))/3.0

*Eq(1) failed:    DPDTPRE(J0,I) = (POTP0*3.0-POTP1*4.0+POTP2) / (2.*DTN)
*Eq(2) failed:    DPDTPRE(J0,I) = (POTP0-POTP1) / DTN
*Eq(3) Good!!!    DPDTPRE(J0,I) = (POTP1-POTP2) / DTN
        DPDTPRE(J0,I) = (POTP0-POTP2) / DTN / 2.0

C        DPDTPRE(J0,I) = (POTM(J0,I)-POTM(J1,I)) / DTN

      END DO
      END SUBROUTINE


      SUBROUTINE DPOTDT_UNS_OLD
************************************************************************
*      Old_Method_Stable_But_Not_Good_For_Unsteady                     *
*      In unsteady application, we can no longer use the potential     *
*      the previous revolution in calculating the current dPOTdT.      *
*      A 2nd order backward approximation can be used in this case.    *
*                                                                      *
*      Author: Yiran Su                                                *
*                                                                      *
*      Date           Comments/Revision                                *
*      -----------    ----------------------------------------------   *
*      09/20/2017     subroutine created                               *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      DTN=DELTAT*ADVCO/PI

C Initialize at the first unsteady step
      IF(NTSTEP.EQ.1) THEN
        DO I=1,NPANEL
          DO J=1,NTPREV
            POTM(J,I)=0.0
            DPDTPRE(J,I)=0.0
          END DO
        END DO
      END IF

C when to start calculating dPOTdT
      NSTART=NTPREV/NBLADE

      DO I=1,NPANEL

        POTM(IDXREV,I)=POT(I)  ! Store data


        J0 = IDXREV - 1
        IF (J0.LE.0) J0=J0+NTPREV
        J1 = J0 - 1
        IF (J1.LE.0) J1=J1+NTPREV
        J2 = J0 - 2
        IF (J2.LE.0) J2=J2+NTPREV
        J3 = J0 - 3
        IF (J3.LE.0) J3=J3+NTPREV

        IF(NTSTEP_KEY.LE.NSTART) THEN ! Start at certain time step
          DPDTPRE(J1,I) = 0.0
          CYCLE
        END IF

        POTP0 = (POTM(IDXREV,I)+POTM(J0,I)+POTM(J1,I))/3.0
        POTP1 = (POTM(J0,I)+POTM(J1,I)+POTM(J2,I))/3.0
        POTP2 = (POTM(J1,I)+POTM(J2,I)+POTM(J3,I))/3.0

        DPDTPRE(J1,I) = (POTP0-POTP2) / (2.*DTN)

      END DO
      END SUBROUTINE

