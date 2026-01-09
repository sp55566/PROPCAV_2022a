      SUBROUTINE COMPPHI1(M,IDR)
************************************************************************
*  Define the potential (minus phi0) on blade surface beneath the      *
*  cavity. phi=phi0+phi1. Used in the dynamic boundary condition.      *
*                                                                      *
*  Date of latest revision      Revision                               *
*  -----------------------      --------                               *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      THETAT=-FLOAT(NTSTEP-1)*DELTAT
      ADVCO2=ADVCO*ADVCO
      PI2=PI*PI

      JCAV=JCV(M,IDR)
      M0M=M0(M,IDR)
      M0M1=M0(M,IDR)-1

C-----------------------------------------------------------------------
C     Calculate dphids for cavity on the blade
C-----------------------------------------------------------------------
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
          IF((ITERGB.EQ.1).AND.(IFLAP.EQ.1))THEN
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

C.......This is the pressure law equation...............................
         PLAW=ONE-FSJ(SPZ(N1,M)-SZ(NHP,M))

C.......This is the unsteady term in Bernoulli's equation...............
         IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1) THEN
            DPDT=0.0
         ELSE
            NMIN1=2
            IF(ISP.EQ.1.AND.NBLADE.GT.1) NMIN1=3
            IF(IFACE.EQ.2) NMIN1=4
            IF(NREV.GT.NMIN1) THEN
               DPDT=DPDTPRE(IDXREV,J1)
            ELSE
               DPDT=0.
            END IF
         END IF

C.......If ITAN=1, the crossflow will be zeroed. (JY053000)
C.......If IDRT=1, panel input method ignores the crossflow (SKIM072419)
C.......If IFACE=2, the crossflow will be zeroed. (SNKIM082719)
         ITAN=0
         IF(ICON.NE.5.AND.ICON.NE.6) THEN
            IF(ISP.EQ.1.OR.IFACE.EQ.2.OR.IDRT.EQ.1) THEN
               ITAN=1
            ELSE
C.............For the first iteration, it's better to assume to flow is
C.............tangent than using the wetted DPDV. (JY080300)
               IF(ITER.EQ.1.AND.IT2.EQ.0) ITAN=1
            END IF
C........Set flow to be tangent for fan-like SP-propeller (JY091100)
            IF(ICON.EQ.7) ITAN=1
         END IF

C.......DUM1 is the total tangential velocity in v direction............
         IF(ITAN.NE.1) DUM1=UT+DPDVB(N1,M)

C.......QC2 is the square of the total cavity velocity 
C.......Call special routine to calculate QC2 if ICON=8 (JY110100)
         IF(ICON.EQ.8) THEN
            CALL CALQC2_ICON8(N,N1,M,IDR,PLAW,ADVCO2,UW2,QC2)
         ELSE
            QC2=SIGMA*PLAW/ADVCO2+UW2-TWOGY-TWO*DPDT
         END IF

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

C.......Call special routine to calculate QC2 if ICON=8 (JY110100)
                  IF(ICON.EQ.8) THEN
                     CALL CALQC2_ICON8(N,N1,M,IDR,PLAW,ADVCO2,UW2,QC2)
                  ELSE
                     QC2=SIGMA*PLAW/ADVCO2+UW2-TWOGY
                  END IF

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

C.......Call special routine to calculate QC2 if ICON=8 (JY110100)
                  IF(ICON.EQ.8) THEN
                     CALL CALQC2_ICON8(N,N1,M,IDR,PLAW,ADVCO2,UW2,QC2)
                  ELSE
                     QC2=SIGMA*PLAW/ADVCO2+UW2-TWOGY
                  END IF

                  DUM2=QC2-DUM1*DUM1
               END IF
               IF(DUM2.LT.ZERO) THEN                  
                  IF(ISP.NE.1) THEN
                     IF(IAN .NE. 2) CALL VECPLT2
                  ENDIF
                  WRITE(*,*) 'STOP! Negative inside square root!'
                  WRITE(*,*) 'Plotted velocity vectors in check.vec!'
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
C     Calculate dphids for supercavity on the wake
C-----------------------------------------------------------------------
      IF(NWC(M,IDR).GT.0) THEN
         N1=JCAV+1      
         DO 30 N=1,NWC(M,IDR)+NSPS(M,IDR)
            J=(MR-M)*NTRA+N
            RCP=SQRT(XCTW(J,2)*XCTW(J,2)+XCTW(J,3)*XCTW(J,3))
            THETA=THETAT+ACOS(XCTW(J,2)/RCP)
            
C..........This is the gravity term in Bernoulli's equation.............
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

C..........The gravity term is zero if we assume infinite Froude number
C..........for the surface piercing propeller. (JY112199)
            IF(ISP.EQ.1) TWOGY=0.
   
C..........UW2 is the square of the inflow velocity.....................
            UW2=VOXW(J)*VOXW(J)+VOYW(J)*VOYW(J)+VOZW(J)*VOZW(J)
            
C..........This is the total tangential velocites in the s direction....
            US=-VOXW(J)*ULW(J,1)-VOYW(J)*ULW(J,2)-VOZW(J)*ULW(J,3)

C..........This is the unsteady term in Bernoulli's equation............
C..........Setting DPDT on the wake to be zero.
            IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1) THEN
               DPDT=0.0
            ELSE
               IF(NREV.GT.2) THEN
                  DPDT=DPDTPREW(IDXREV,J)
               ELSE
                  DPDT=0.
               END IF
            END IF

C..........DPHIDS is the total velocity in s direction..................
C..........If ICON=8, SIGMA=SIGMAB on supercavity (JY110100)
            IF(ICON.EQ.8) THEN
               QC2=SIGMAB/ADVCO2+UW2
            ELSE
               QC2=SIGMA/ADVCO2+UW2-TWOGY-TWO*DPDT
            END IF

            IF(QC2.LT.ZERO) THEN
C.............Set DP/DT=0 if QC2 is negative.

               IF(ICON.EQ.8) THEN
                  QC2=SIGMAB/ADVCO2+UW2
               ELSE
                  QC2=SIGMA/ADVCO2+UW2-TWOGY
               END IF
               IF(QC2.LT.ZERO) THEN
                  WRITE(*,*) 'COMPPHI1-C: negative inside square root!'
                  WRITE(*,*) 'ITER=',ITER,'M=',M
                  WRITE(*,*) 'qc^2=',QC2
                  STOP
               END IF
            END IF

            DPHIDS(JCAV+N,M,IDR)=SQRT(QC2)
            
C..........QC is the perturbation velocity in s direction...............
            QC(N1+N,M,IDR)=-US+DPHIDS(JCAV+N,M,IDR)       
 30      CONTINUE
      END IF

C-----------------------------------------------------------------------
C     determine QC at the leading edge
C-----------------------------------------------------------------------
      IF(IDR.EQ.1) THEN
         K=1
         ISF=0
      ELSE
         K=-1
         ISF=1
      END IF

      IF(JCAV.GT.0) THEN
         AK1=FLOAT(K)
         IF(JCAV.GE.4) THEN
            S0=AK1*(SPZ(M0M1+ISF+K*4,M)-SZ(M0M,M))
            S1=AK1*(SPZ(M0M1+ISF+K*4,M)-SPZ(M0M1+K*1+ISF,M))
            S2=AK1*(SPZ(M0M1+ISF+K*4,M)-SPZ(M0M1+K*2+ISF,M))
            S3=AK1*(SPZ(M0M1+ISF+K*4,M)-SPZ(M0M1+K*3+ISF,M))
            CALL CUBEXT(S0,S1,S2,S3,A1,A2,A3,A4)
            QC(1,M,IDR)=A1*QC(2,M,IDR)+A2*QC(3,M,IDR)+A3*QC(4,M,IDR)+
     *           A4*QC(5,M,IDR)
         ELSE IF(JCAV.EQ.3) THEN
            S0=AK1*(SPZ(M0M1+ISF+K*3,M)-SZ(M0M,M))
            S1=AK1*(SPZ(M0M1+ISF+K*3,M)-SPZ(M0M1+K*1+ISF,M))
            S2=AK1*(SPZ(M0M1+ISF+K*3,M)-SPZ(M0M1+K*2+ISF,M))
            CALL QUADEXT(S0,S1,S2,A1,A2,A3)
            QC(1,M,IDR)=A1*QC(2,M,IDR)+A2*QC(3,M,IDR)+A3*QC(4,M,IDR)
         ELSE IF(JCAV.EQ.2) THEN
            S0=AK1*(SPZ(M0M1+ISF+K*2,M)-SZ(M0M,M))
            S1=AK1*(SPZ(M0M1+ISF+K*2,M)-SPZ(M0M1+K*1+ISF,M))
            A1=S0/S1
            A2=ONE-S0/S1
            QC(1,M,IDR)=A1*QC(2,M,IDR)+A2*QC(3,M,IDR)
         ELSE IF(JCAV.EQ.1) THEN
            A1=ONE
            QC(1,M,IDR)=A1*QC(2,M,IDR)
         END IF

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
      
      IF(NWC(M,IDR).GT.0) THEN
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
            N1=JCAV+1+N
            QINT=QINT+HALF*DDS*(QC(N1-1,M,IDR)+QC(N1,M,IDR))
            PHI1(JCAV+N,M,IDR)=QINT
 50      CONTINUE
      END IF

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine PHI1>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END






