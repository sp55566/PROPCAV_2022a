      SUBROUTINE COMPPHI1_SC2(M,IDR)
************************************************************************
*  Define the potential (minus phi0) on SR+wake surface beneath the    *
*  cavity. phi=phi0t+phi2. Used in the dynamic boundary condition.     *
*                                                                      *
*  Date of latest revision      Revision                               *
*  -----------------------      --------                               *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      THETAT=-FLOAT(NTSTEP-1)*DELTAT
      ADVCO2=ADVCO*ADVCO
      PI2=PI*PI

C-----------------------------------------------------------------------
C     Calculate dphids for cavity on the SR
C-----------------------------------------------------------------------
      N0M=N0(IDR)
      N0M1=N0(IDR)-1

      IF(ISP.EQ.0) THEN
         NL1=NSR2
      ELSE
         IF(IDR.EQ.1) THEN
            NL1=MAX(0,IC(2,M,IDXREV)-N0(1)+1)
         ELSE
            NL1=MAX(0,N0(2)-IW(1,M,IDXREV))
         END IF
      END IF
            
      DO N=1,NL1
         IF(IDR.EQ.2)THEN
            N1=N0M-N
         ELSE IF(IDR.EQ.1) THEN
            N1=N0M1+N
         ENDIF
         J1=INDEXB(N1,M)
         RCP=SQRT(XCT(J1,2)*XCT(J1,2)+XCT(J1,3)*XCT(J1,3))
         THETA=THETAT+ACOS(XCT(J1,2)/RCP)

C.......This is the unsteady term in Bernoulli's equation...............
C         IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1) THEN
            DPDT=0.0
C         ELSE
C            IF(NREV.GT.4.AND.ISP.NE.1) THEN
C               DPDT=DPDTPRE(IDXREV,J1)
C            ELSE
C               DPDT=0.
C            END IF
C         END IF

C.......This is the gravity term in Bernoulli's equation................
         IF((ISTEADY.EQ.0).OR.(ISTEADY.EQ.3)) THEN
            TWOGY=0.
CSH---------------shree----------07/08/2002--------
          IF((ITERGB.EQ.1).OR.(IFLAP.EQ.1))THEN
             ZTEMP = RCP*COS(THETA)
             TWOGY=-ZTEMP*2./(FNRUDDER)
                     ENDIF
CSH---------------shree----------07/08/2002--------
         ELSE
            TWOGY=RCP*COS(THETA)/(FROUDE*ADVCO2)
         END IF
         
C.......The gravity term is zero if we assume infinite Froude number
C.......for the surface piercing propeller. (JY112199)
         IF(ISP.EQ.1) TWOGY=0.

C.......This is the total tangential velocites in the s and v ..........
C.......direction.......................................................
         US=-VOX1(J1)*UL(J1,1)-VOY1(J1)*UL(J1,2)-VOZ1(J1)*UL(J1,3)
         UT=VOX1(J1)*VL(J1,1)+VOY1(J1)*VL(J1,2)+VOZ1(J1)*VL(J1,3)
         IF(IDR.EQ.2) US=-US
         
C.......UW2 is the square of the inflow + rotation......................
         UW2=VINFSB(N1,M)

C.......If ITAN=1, the crossflow will be zeroed. (JY053000)
         ITAN=1
c         IF(ICON.NE.5.AND.ICON.NE.6) THEN
c            IF(ISP.EQ.1) THEN
c               ITAN=1
c            ELSE
c               IF(ITER.EQ.1.AND.IT2.EQ.0) ITAN=1
c            END IF
c         END IF

C.......DUM1 is the total tangential velocity in v direction............
         IF(ITAN.NE.1) DUM1=UT+DPDVB(N1,M)

         QC2=SIGMA/ADVCO2+UW2-TWOGY-TWO*DPDT

C.......DUM2 is the square of the total cavity velocity subtract........
C.......the square of the total tangential velocity in v direction......
         IF(ITAN.NE.1) DUM2=QC2-DUM1*DUM1

C.......DPHIDS is the total velocity in s direction.....................
C.......DPHIDV is the total velocity in v direction.....................

         IF(ITAN.EQ.1) THEN
            IF(QC2.LT.ZERO) THEN
               WRITE(*,*) 'COMPPHI1-A: negative inside square root!'
               WRITE(*,*) 'ITER=',ITER,'M=',M,'N=',N1
               WRITE(*,*) 'qc^2=',QC2
               IF(ISTEADY.NE.0.AND.ISTEADY.NE.1) THEN
                  WRITE(*,*) 'probably due to error in cal. DPDT'
                  WRITE(*,*) 'setting DPDT=ZERO'
                  QC2=SIGMA/ADVCO2+UW2-TWOGY
               END IF
               IF(QC2.LT.ZERO) THEN
                  WRITE(*,*) 'STOP! Negative inside square root!'
                  STOP
               END IF
            END IF
            DPHIDS2(N,M,IDR)=SQRT(QC2)
            IF(IDR.EQ.2) THEN
               DPHIDV2(N,M,IDR)=DPHIDS2(N,M,IDR)*SINPHI(J1)
            ELSE IF(IDR.EQ.1) THEN
               DPHIDV2(N,M,IDR)=-DPHIDS2(N,M,IDR)*SINPHI(J1)
            END IF
         ELSE

C.......Since DUM2 can be less than zero if the DPDVB has not converge,
C.......let's temporarily divide DUM1/2. (JY083100)
            IF(DUM2.LT.ZERO) THEN
               IF(N.GT.1) THEN
                  DUM1=DPHIDV(N-1,M,IDR)
                  IF(DUM1*DUM1.GT.QC2) THEN
                     DUM1=SIGN(SQRT(QC2),DUM1)/2.
                  END IF
               ELSE
                  DUM1=SIGN(SQRT(QC2),DUM1)/2.
               END IF
               DUM2=QC2-DUM1*DUM1
           END IF

            IF(DUM2.LT.ZERO) THEN
               WRITE(*,*) 'COMPPHI1-B: negative inside square root!'
               WRITE(*,*) 'ITER=',ITER,'M=',M,'N=',N1
               WRITE(*,*) 'qc^2=',QC2
               WRITE(*,*) 'qv^2=',DUM1*DUM1
               IF(ISTEADY.NE.0.AND.ISTEADY.NE.1) THEN
                  WRITE(*,*) 'probably due to error in cal. DPDT'
                  WRITE(*,*) 'setting DPDT=ZERO'
                  QC2=SIGMA/ADVCO2+UW2-TWOGY
                  DUM2=QC2-DUM1*DUM1
               END IF
               IF(DUM2.LT.ZERO) THEN                  
C                  IF(ISP.NE.1) CALL VECPLT2
                  WRITE(*,*) 'STOP! Negative inside square root!'
C                  WRITE(*,*) 'Plotted velocity vectors in check.vec!'
                  STOP
               END IF
            END IF
            IF(IDR.EQ.2) THEN
               DPHIDS2(N,M,IDR)=SINPHI(J1)*DUM1+COSPHI(J1)*SQRT(DUM2)
            ELSE IF(IDR.EQ.1) THEN
               DPHIDS2(N,M,IDR)=-SINPHI(J1)*DUM1+COSPHI(J1)*SQRT(DUM2)
            ENDIF
            DPHIDV2(N,M,IDR)=DUM1

         END IF

C.......QC is the perturbation velocity in s direction..................
         QCSR(N+1,M,IDR)=-US+DPHIDS2(N,M,IDR)

      END DO

C-----------------------------------------------------------------------
C     Calculate dphids for supercavity on the wake
C-----------------------------------------------------------------------
      N1=NSR2+1      
      DO 30 N=1,NWC(M,IDR)+NSPS(M,IDR)
         J=(MR-M)*NTRA+N
         RCP=SQRT(XCTW(J,2)*XCTW(J,2)+XCTW(J,3)*XCTW(J,3))
         THETA=THETAT+ACOS(XCTW(J,2)/RCP)
            
C.......This is the gravity term in Bernoulli's equation.............
         IF((ISTEADY.EQ.0).OR.(ISTEADY.EQ.3)) THEN
            TWOGY=0.
CSH---------------shree----------07/08/2002--------
          IF((ITERGB.EQ.1).OR.(IFLAP.EQ.1))THEN
             ZTEMP = RCP*COS(THETA)
             TWOGY=-ZTEMP*2.0/(FNRUDDER)
                     ENDIF
CSH---------------shree----------07/08/2002--------
         ELSE
            TWOGY=RCP*COS(THETA)/(FROUDE*ADVCO2)
         END IF

C.......The gravity term is zero if we assume infinite Froude number
C.......for the surface piercing propeller. (JY112199)
         IF(ISP.EQ.1) TWOGY=0.
   
C.......UW2 is the square of the inflow velocity.....................
         UW2=VOXW(J)*VOXW(J)+VOYW(J)*VOYW(J)+VOZW(J)*VOZW(J)
            
C.......This is the total tangential velocites in the s direction....
         US=-VOXW(J)*ULW(J,1)-VOYW(J)*ULW(J,2)-VOZW(J)*ULW(J,3)

C.......This is the unsteady term in Bernoulli's equation............
C.......Setting DPDT on the wake to be zero.
C         IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.ISP.NE.1) THEN
            DPDT=0.0
C         ELSE
C            IF(NREV.GT.2) THEN
C               DPDT=DPDTPREW(IDXREV,J)
C            ELSE
C               DPDT=0.
C            END IF
C         END IF

         QC2=SIGMA/ADVCO2+UW2-TWOGY-TWO*DPDT

         IF(QC2.LT.ZERO) THEN
C..........Set DP/DT=0 if QC2 is negative.
            QC2=SIGMA/ADVCO2+UW2-TWOGY
            IF(QC2.LT.ZERO) THEN
               WRITE(*,*) 'COMPPHI1-C: negative inside square root!'
               WRITE(*,*) 'ITER=',ITER,'M=',M
               WRITE(*,*) 'qc^2=',QC2
               STOP
            END IF
         END IF

         DPHIDS2(NSR2+N,M,IDR)=SQRT(QC2)
            
C.......QC is the perturbation velocity in s direction...............
         QCSR(N1+N,M,IDR)=-US+DPHIDS2(NSR2+N,M,IDR) 
         
 30   CONTINUE

C-----------------------------------------------------------------------
C     determine QC at the leading edge of the SR
C-----------------------------------------------------------------------
      IF(IDR.EQ.1) THEN
         K=1
         ISF=0
      ELSE
         K=-1
         ISF=1
      END IF
      AK1=FLOAT(K)

      IF(NL1.GE.4) THEN
         S0=AK1*(SPZ(N0M1+ISF+K*4,M)-SZ(N0M,M))
         S1=AK1*(SPZ(N0M1+ISF+K*4,M)-SPZ(N0M1+K*1+ISF,M))
         S2=AK1*(SPZ(N0M1+ISF+K*4,M)-SPZ(N0M1+K*2+ISF,M))
         S3=AK1*(SPZ(N0M1+ISF+K*4,M)-SPZ(N0M1+K*3+ISF,M))
         CALL CUBEXT(S0,S1,S2,S3,A1,A2,A3,A4)
         QCSR(1,M,IDR)=A1*QCSR(2,M,IDR)+A2*QCSR(3,M,IDR)+
     *        A3*QCSR(4,M,IDR)+A4*QCSR(5,M,IDR)
      ELSE IF(NL1.EQ.3) THEN
         S0=AK1*(SPZ(N0M1+ISF+K*3,M)-SZ(N0M,M))
         S1=AK1*(SPZ(N0M1+ISF+K*3,M)-SPZ(N0M1+K*1+ISF,M))
         S2=AK1*(SPZ(N0M1+ISF+K*3,M)-SPZ(N0M1+K*2+ISF,M))
         CALL QUADEXT(S0,S1,S2,A1,A2,A3)
         QCSR(1,M,IDR)=A1*QCSR(2,M,IDR)+A2*QCSR(3,M,IDR)+
     *        A3*QCSR(4,M,IDR)
      ELSE IF(NL1.EQ.2) THEN
         S0=AK1*(SPZ(N0M1+ISF+K*2,M)-SZ(N0M,M))
         S1=AK1*(SPZ(N0M1+ISF+K*2,M)-SPZ(N0M1+K*1+ISF,M))
         A1=S0/S1
         A2=ONE-S0/S1
         QCSR(1,M,IDR)=A1*QCSR(2,M,IDR)+A2*QCSR(3,M,IDR)
      ELSE IF(NL1.EQ.1) THEN
         A1=ONE
         QCSR(1,M,IDR)=A1*QCSR(2,M,IDR)
      END IF

      QINT=ZERO
      DO N=1,NL1
         N1=N0M1+ISF+K*N
         IF(N.EQ.1) THEN
            DDS=AK1*(SPZ(N1,M)-SZ(N0M,M))
         ELSE
            DDS=AK1*(SPZ(N1,M)-SPZ(N1-K*1,M))
         END IF            
         QINT=QINT+HALF*DDS*(QCSR(N,M,IDR)+QCSR(N+1,M,IDR))
         PHI2(N,M,IDR)=QINT
      END DO

      DO 50 N=1,NWC(M,IDR)+NSPS(M,IDR)
         IF(N.EQ.1)THEN
            IF(IDR.EQ.1) THEN
               DDS=SPZ(NCP,M)-SZ(NCP,M)+HALF*DS(NC,M)
            ELSE IF(IDR.EQ.2) THEN
               DDS=SPZ(NCP,M)-SZ(NCP,M)+HALF*DS(1,M)
            END IF
         ELSE
            DDS=SPZ(NC+N,M)-SPZ(NC-1+N,M)
         ENDIF
         N1=NSR2P+N
         QINT=QINT+HALF*DDS*(QCSR(N1-1,M,IDR)+QCSR(N1,M,IDR))
         PHI2(NSR2+N,M,IDR)=QINT
 50   CONTINUE

      RETURN
      END
      
