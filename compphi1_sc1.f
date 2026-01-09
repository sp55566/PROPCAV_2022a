      SUBROUTINE COMPPHI1_SC1(M,IDR)
************************************************************************
*  Define the potential (minus phi0) on blade surface beneath the      *
*  cavity. phi=phi0+phi1. Used in the dynamic boundary condition.      *
*                                                                      *
*  Date of latest revision      Revision                               *
*  -----------------------      --------                               *
*  JY093601            Copied from compphi1.f.  Modified for ISC=1.    *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      THETAT=-FLOAT(NTSTEP-1)*DELTAT
      ADVCO2=ADVCO*ADVCO
      PI2=PI*PI

C-----------------------------------------------------------------------
C     Calculate dphids for cavity on the blade
C-----------------------------------------------------------------------
      JCAV=JCV(M,IDR)
      M0M=M0(M,IDR)
      M0M1=M0(M,IDR)-1

      DO 20 N=1,JCAV+NSPP(M,IDR)
         IF(IDR.EQ.2)THEN
            N1=M0M-N
         ELSE IF(IDR.EQ.1) THEN
            N1=M0M1+N
         ENDIF
         J1=INDEXB(N1,M)
         RCP=SQRT(XCT(J1,2)*XCT(J1,2)+XCT(J1,3)*XCT(J1,3))
         THETA=THETAT+ACOS(XCT(J1,2)/RCP)
            
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

C.......This is the unsteady term in Bernoulli's equation...............
         IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1) THEN
            DPDT=0.0
         ELSE
            IF(ISP.EQ.0) THEN
               NMIN=4
            ELSE
               NMIN=2
               IF(NBLADE.GT.1) NMIN=3
            END IF
            IF(NREV.GT.NMIN) THEN
               DPDT=DPDTPRE(IDXREV,J1)
            ELSE
               DPDT=0.
            END IF
         END IF

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
            DPHIDS(N,M,IDR)=SQRT(QC2)
            IF(IDR.EQ.2) THEN
               DPHIDV(N,M,IDR)=DPHIDS(N,M,IDR)*SINPHI(J1)
            ELSE IF(IDR.EQ.1) THEN
               DPHIDV(N,M,IDR)=-DPHIDS(N,M,IDR)*SINPHI(J1)
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
               DPHIDS(N,M,IDR)=SINPHI(J1)*DUM1+COSPHI(J1)*SQRT(DUM2)
            ELSE IF(IDR.EQ.1) THEN
               DPHIDS(N,M,IDR)=-SINPHI(J1)*DUM1+COSPHI(J1)*SQRT(DUM2)
            ENDIF
            DPHIDV(N,M,IDR)=DUM1

         END IF

C.......QC is the perturbation velocity in s direction..................
         QC(N+1,M,IDR)=-US+DPHIDS(N,M,IDR)

 20   CONTINUE

C-----------------------------------------------------------------------
C     determine QC at the leading edge of the cavity
C-----------------------------------------------------------------------
      IF(IDR.EQ.1) THEN
         K=1
         ISF=0
      ELSE
         K=-1
         ISF=1
      END IF
      AK1=FLOAT(K)

      NDUM=JCAV
      IF(SOP(M,IDR).EQ.ONE) THEN
         IF(ISP.EQ.0) THEN
            NDUM=NDUM+NSR2
         ELSE
            IF(IDR.EQ.1.AND.IC(2,M,IDXREV).GE.N0(1))
     *           NDUM=NDUM+IC(2,M,IDXREV)-N0(1)+1
         END IF
      END IF

      AQC1=QC(2,M,IDR)
      IF(NDUM.GE.4) THEN
         IF(JCAV.GE.4) THEN
            AQC2=QC(3,M,IDR)
            AQC3=QC(4,M,IDR)
            AQC4=QC(5,M,IDR)
         ELSE IF(JCAV.EQ.3) THEN
            AQC2=QC(3,M,IDR)
            AQC3=QC(4,M,IDR)
            AQC4=QCSR(2,M,IDR)
         ELSE IF(JCAV.EQ.2) THEN
            AQC2=QC(3,M,IDR)
            AQC3=QCSR(2,M,IDR)
            AQC4=QCSR(3,M,IDR)
         ELSE 
            AQC2=QCSR(2,M,IDR)
            AQC3=QCSR(3,M,IDR)
            AQC4=QCSR(4,M,IDR)
         END IF
      ELSE IF(NDUM.EQ.3) THEN
         IF(JCAV.EQ.3) THEN
            AQC2=QC(3,M,IDR)
            AQC3=QC(4,M,IDR)
         ELSE IF(JCAV.EQ.2) THEN
            AQC2=QC(3,M,IDR)
            AQC3=QCSR(2,M,IDR)
         ELSE 
            AQC2=QCSR(2,M,IDR)
            AQC3=QCSR(3,M,IDR)
         END IF
      ELSE IF(NDUM.EQ.2) THEN
         IF(JCAV.EQ.2) THEN  
            AQC2=QC(3,M,IDR)
         ELSE
            AQC2=QCSR(2,M,IDR)
         END IF
      END IF
         
      IF(NDUM.GE.4) THEN
         S0=AK1*(SPZ(M0M1+ISF+K*4,M)-SZ(M0M,M))
         S1=AK1*(SPZ(M0M1+ISF+K*4,M)-SPZ(M0M1+K*1+ISF,M))
         S2=AK1*(SPZ(M0M1+ISF+K*4,M)-SPZ(M0M1+K*2+ISF,M))
         S3=AK1*(SPZ(M0M1+ISF+K*4,M)-SPZ(M0M1+K*3+ISF,M))
         CALL CUBEXT(S0,S1,S2,S3,A1,A2,A3,A4)
         QC(1,M,IDR)=A1*AQC1+A2*AQC2+A3*AQC3+A4*AQC4
      ELSE IF(NDUM.EQ.3) THEN
         S0=AK1*(SPZ(M0M1+ISF+K*3,M)-SZ(M0M,M))
         S1=AK1*(SPZ(M0M1+ISF+K*3,M)-SPZ(M0M1+K*1+ISF,M))
         S2=AK1*(SPZ(M0M1+ISF+K*3,M)-SPZ(M0M1+K*2+ISF,M))
         CALL QUADEXT(S0,S1,S2,A1,A2,A3)
         QC(1,M,IDR)=A1*AQC1+A2*AQC2+A3*AQC3
      ELSE IF(NDUM.EQ.2) THEN
         S0=AK1*(SPZ(M0M1+ISF+K*2,M)-SZ(M0M,M))
         S1=AK1*(SPZ(M0M1+ISF+K*2,M)-SPZ(M0M1+K*1+ISF,M))
         A1=S0/S1
         A2=ONE-S0/S1
         QC(1,M,IDR)=A1*AQC1+A2*AQC2
      ELSE IF(NDUM.EQ.1) THEN
         A1=ONE
         QC(1,M,IDR)=A1*AQC1
      END IF

C-----------------------------------------------------------------------
C     Determine phi-phi0 on the cavity.
C-----------------------------------------------------------------------
      QINT=ZERO
      
      DO 40 N=1,JCAV+NSPP(M,IDR)
         N1=M0M1+ISF+K*N
         IF(N.EQ.1) THEN
            DDS=AK1*(SPZ(N1,M)-SZ(M0M,M))
         ELSE
            DDS=AK1*(SPZ(N1,M)-SPZ(N1-K*1,M))
         END IF            
         QINT=QINT+HALF*DDS*(QC(N,M,IDR)+QC(N+1,M,IDR))
         PHI1(N,M,IDR)=QINT
 40   CONTINUE

C....Determine PSI0T (PHI-PHI0 @ blade TE) for SR
      IF(SOP(M,IDR).EQ.ONE) THEN
         DDS=AK1*(SZ(N1+K*1+ISF,M)-SPZ(N1,M))
         QINT=QINT+DDS*(QC(N,M,IDR)+3*QC(N+1,M,IDR))/4.
         PSI0T(M,IDR)=QINT
      END IF
      
      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine PHI1>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END






