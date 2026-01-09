C ----------------------------------------------------
      SUBROUTINE KUTTA(ITMAX)
C ----------------------------------------------------
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
!     PARAMETER(MBZSQ=MBZ*MBZ,MBZCUB=4*MBZ-4)
!     DIMENSION XJ(MBZSQ),IPVTJ(MBZSQ),WORKJ(MBZSQ)
      DIMENSION XJ(MBZ*MBZ),IPVTJ(MBZ*MBZ),WORKJ(MBZ*MBZ)
      DIMENSION DPMOR(MBZ),DFCOR(MBZ)
      DIMENSION RSQ1(MBZ),DFCOR1(MBZ)
      DIMENSION XJSN(MBZ)

      integer MBZSQ, MBZCUB
      MBZSQ=MBZ*MBZ
      MBZCUB=4*MBZ-4

C      COMMON/BKUTTA/POTEMP(NTZ),WWK(NTZ,MBZ)

      NMT = NPANEL

      DO N = 1 , NMT
        POT(N) = POTEMP(N)
      ENDDO

      IT = 0

      CALL PRSDIF(2)
!YE TIAN 04/08/2012 add following line
      DELPonW=DELP
C/s S.N.KIM | Ductwake alignment | 1/30/2017
      IF(IDUCT.NE.0) DELPonDW=DELPD
C/e S.N.KIM | Ductwake alignment | 1/30/2017

      DO M = 1 , MR
        DPMOR(M) = DELP(M)
        DCPPRE(M) = DELCP(M)
        XJSN(M) = 1.0
      END DO

      ICOUNT = 0
      DO M = 1 , MR
        IF(ABS(DPMOR(M)) .LE. 1.E-5) ICOUNT = ICOUNT + 1
      END DO

      IF(ICOUNT .GT. 0) THEN
        IF(IPKT.EQ.1) THEN
          WRITE(*,*) 'It satisfies Kutta Condition w/o IPK '
        END IF
        RETURN
      END IF

      BETAM = -0.01

      DO M = 1 , MRTIP
        DO K = 1 , NMT
          POT(K) = POTEMP(K)+BETAM*DPMOR(M)*WWK(K,M)
        END DO

        CALL PRSDIF(2)
!YE TIAN 04/08/2012 add following line
        DELPonW=DELP
C/s S.N.KIM | Ductwake alignment | 1/30/2017
        IF(IDUCT.NE.0) DELPonDW=DELPD
C/e S.N.KIM | Ductwake alignment | 1/30/2017

        DO J = 1 , MRTIP
          JINDX = (M-1)*MRTIP + J
          XJDUM = 4.0*(DELCP(J)-DCPPRE(J))/(BETAM*DPMOR(M))
          XJ(JINDX) = XJDUM
        END DO
      END DO

      DO M = 1 , MR
        DELCP(M) = DCPPRE(M)
      END DO

      CALL SDECOMP(MRTIP,MRTIP,XJ,COND,IPVTJ,WORKJ)

      DO M = 1, MRTIP
         J=(M-1)*MRTIP + M
         IF(XJ(J) .EQ. 0.0) THEN
            XJ(J) = 1.0
            XJSN(M) = 0.0
         END IF
      END DO

      DO M = 1 , MR
        DFCOR(M) = XJSN(M)*DPMOR(M)*0.01*SIGN(1.,DELCP(M))
      ENDDO

      DO IT = 1, ITMAX

        DO K = 1 , NMT
          POT(K) = POTEMP(K)
          DO M = 1 , MR
            POT(K) = POT(K) + DFCOR(M)*WWK(K,M)
          END DO
        END DO

        CALL PRSDIF(2)
!YE TIAN 04/08/2012 add following line
        DELPonW=DELP
C/s S.N.KIM | Ductwake alignment | 1/30/2017
        IF(IDUCT.NE.0) DELPonDW=DELPD
C/e S.N.KIM | Ductwake alignment | 1/30/2017

        DO M = 1 , MR
          DELP(M) = DPMOR(M) + DFCOR(M)
        END DO

        DIF  = -9999.
        DO  M = 1 , MRTIP
          DIF = AMAX1(DIF,ABS(DELCP(M)))
        END DO

        IF(IT.EQ.ITMAX.AND.DIF.GT.TOLKT.AND.IPKT.EQ.1) THEN
          WRITE(*,1001)
 1001     FORMAT(/5X,'WARNING!!! IPK did NOT converge!')
          WRITE(*,1002) DIF
 1002     FORMAT(/5X,'Max difference in DELCP=',F10.4)
        END IF

        IF(DIF.LT.TOLKT) THEN
          IF(IPKT.EQ.1) THEN
            WRITE(*,1000) IT , DIF
 1000       FORMAT(/5X,'ITERATION FOR PRESS. KUTTA=',I5,
     &             /5X,'MAX DIFFERENCE IN DELCP=',F10.4)
            WRITE(*,*) (DELCP(M),M=1,MR)
          END IF
          GO TO 180
        END IF

        ! Terminate IPK when DELCP goes up. Yiran Su 10/17/2017
*        IF ((DIF.GT.DIF1).AND.(IT.GT.10)) THEN
*          WRITE(*,*) 'IPK terminated before convergence.'
*          GO TO 180
*        END IF
        DIF1=DIF
        ! if (NTSTEP_KEY.gt.175) WRITE(*,*) it,delcp(1:mrtip)

        CALL SDSOLVE(MRTIP,MRTIP,XJ,DELCP,IPVTJ)

        DO M=1,MRTIP
          DFCOR(M)=DFCOR(M)-DELCP(M)
        END DO

C........if MR2 ne. MR, apply square root extrapolations
C........03-08-90 modified

        IF(MRTIP.NE.MR) THEN
          MR2P=MRTIP+1
          DO M=1,MRTIP
            RSQ1(M)=HRZPSQ(1,M)
            DFCOR1(M)=DFCOR(M)
          END DO
          RSQ1(MR2P)=ONE
          DFCOR1(MR2P)=ZERO
          DO M=MR2P,MR
            CALL BQUAD(RSQ1,DFCOR1,HRZPSQ(1,M),IM,COA,COB,COC,YY,MR2P)
            DFCOR(M)=YY
          END DO
        END IF

      ENDDO

180   CONTINUE

!YE TIAN 04/08/2012 add following line
      DELPonW=DELP
C/s S.N.KIM | Ductwake alignment | 01/30/2017
      IF(IDUCT.NE.0) DELPonDW=DELPD
C/e S.N.KIM | Ductwake alignment | 01/30/2017

      CALL PRSDIF(1)

      RETURN
      END SUBROUTINE


