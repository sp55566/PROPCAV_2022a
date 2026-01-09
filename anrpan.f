
*      PROGRAM TEST
*      REAL,DIMENSION(3):: P1,P2,P3,P4,CP
*      REAL FA,FB
*      P1=(/1.0,1.0,0.0/)
*      P2=(/1.0,-1.0,0.0/)
*      P3=(/-1.0,-1.0,0.0/)
*      P4=(/-1.0,1.0,0.0/)
*      CP=(/0.1,0.1,0.0/)
*      CALL ANRPAN(P1,P2,P3,P4,CP,FA,FB)
*      WRITE(*,*) FA,FB
*      END


      SUBROUTINE ANRPAN(P1I,P2I,P3I,P4I,CPI,FAO,FBO)
C  ASPARA IS A EXTERNA L 3*3 VALUE
C  ASTOL IS A EXTERNA L VALUE TOLERANCE

      IMPLICIT NONE
      REAL,DIMENSION(3)::P1I,P2I,P3I,P4I,CPI
      REAL FAO,FBO
      DOUBLE PRECISION,DIMENSION(3):: P1,P2,P3,P4,P5,P6,P7,P8,P9
      DOUBLE PRECISION,DIMENSION(3):: CP,CP1,CP2
      DOUBLE PRECISION A0,B0,FA,FB,A1,B1,A2,B2
      INTEGER I,ST,OVERF
      DOUBLE PRECISION CL,LTH1,H,SGN,ATOL
      DOUBLE PRECISION,DIMENSION(3):: U,V,W

      ATOL=1.0D-6

      DO I=1,3
         P1(I)=DBLE(P1I(I))
         P2(I)=DBLE(P2I(I))
         P3(I)=DBLE(P3I(I))
         P4(I)=DBLE(P4I(I))
         CP(I)=DBLE(CPI(I))
      END DO

CPLUS    DETERMIN HEIGHT
      DO I=1,3
         P8(I)=(P1(I)+P4(I))/2;
         P6(I)=(P2(I)+P3(I))/2;
         P5(I)=(P1(I)+P2(I))/2;
         P7(I)=(P3(I)+P4(I))/2;
      END DO
C  CALCULATE CHARACTERISTIC LENGTH (CL)
      CALL PMREAL(1.0D0,P6,-1.0D0,P8,3,U)
      CALL LEN1REAL(U,3,CL)
      CALL PMREAL(1.0D0,P5,-1.0D0,P7,3,U)
      CALL LEN1REAL(U,3,LTH1)
      IF (LTH1.LT.CL) CL=LTH1
C  CALCULATE CONTROL POINT HEIGHT (H)
      CALL PMREAL(1.0D0,P6,-1.0D0,P5,3,U)
      CALL PMREAL(1.0D0,P7,-1.0D0,P6,3,V)
      CALL CROSSREAL(U,V,W)
      CALL LEN1REAL(W,3,LTH1)
      CALL PMREAL(1.0D0,CP,-1.0D0,P5,3,U)
      H=U(1)*W(1)+U(2)*W(2)+U(3)*W(3)

      ST=0
      IF ((DABS(H)/CL).LT.(1.0D-2)) THEN
         CALL GETPS(P1,P2,P3,P4,CP,H,OVERF,CL)
         IF ((OVERF.EQ.0).AND.((DABS(H)/CL).LT.(1.0D-3))) THEN
            ST=1
            SGN=-1.0D0
            IF (H.GE.0.0D0) SGN=1.0D0
            CALL PMREAL(1.0D0,CP,1.0D-2*CL*SGN,W,3,CP1)
            CALL PMREAL(1.0D0,CP,2.0D-2*CL*SGN,W,3,CP2)
         END IF
      END IF

C     RUNNING WITH 2 DIF CASES
      IF (ST.EQ.0) THEN
         CALL GETAB(P1,P2,P3,P4,CP,A0,B0)
         CALL RECURSIVE_AB(P1,P2,P3,P4,CP,A0,B0,FA,FB,ATOL)
      ELSEIF (ST.EQ.1) THEN
         CALL GETAB(P1,P2,P3,P4,CP1,A0,B0)
         CALL RECURSIVE_AB(P1,P2,P3,P4,CP1,A0,B0,A1,B1,ATOL)
         CALL GETAB(P1,P2,P3,P4,CP2,A0,B0)
         CALL RECURSIVE_AB(P1,P2,P3,P4,CP2,A0,B0,A2,B2,ATOL)
         FA=2.0D0*A1-A2
         FB=2.0D0*B1-B2
      ELSE
         WRITE(*,*) 'HOW DID WE GET HERE?'
      END IF

      FAO=SNGL(FA)
      FBO=SNGL(FB)
      RETURN
      END


      SUBROUTINE GETPS(P1,P2,P3,P4,CP,H,OVERF,CL)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(3):: P1,P2,P3,P4,CP
      DOUBLE PRECISION,DIMENSION(3):: SP1,SP2,U,V
      INTEGER OVERF
      DOUBLE PRECISION TM,TN,TM1,TN1,TOLNST,LTH1,LTH2,H,CL
      TM=0
      TN=0
      OVERF=0
      TOLNST=1.0D-5*CL
CPLUS    TM IS FIXED
 981     CALL GETREALLOC(P1,P2,P3,P4,TM,-1.0D0,SP1)
         CALL GETREALLOC(P1,P2,P3,P4,TM,1.0D0,SP2)
         CALL PMREAL(1.0D0,SP2,-1.0D0,SP1,3,U)
         CALL LEN1REAL(U,3,LTH1)
         CALL PMREAL(1.0D0,CP,-1.0D0,SP1,3,V)
         LTH2=U(1)*V(1)+U(2)*V(2)+U(3)*V(3)
         IF ((LTH2.LT.0.0D0).OR.(LTH2.GT.LTH1)) THEN
            OVERF=1
            RETURN
         END IF
         TN1=LTH2/LTH1*2.0D0-1.0D0
         CALL GETREALLOC(P1,P2,P3,P4,-1.0D0,TN1,SP1)
         CALL GETREALLOC(P1,P2,P3,P4,1.0D0,TN1,SP2)
         CALL PMREAL(1.0D0,SP2,-1.0D0,SP1,3,U)
         CALL LEN1REAL(U,3,LTH1)
         CALL PMREAL(1.0D0,CP,-1.0D0,SP1,3,V)
         LTH2=U(1)*V(1)+U(2)*V(2)+U(3)*V(3)
         IF ((LTH2.LT.0.0D0).OR.(LTH2.GT.LTH1)) THEN
            OVERF=1
            RETURN
         END IF
         TM1=LTH2/LTH1*2.0D0-1.0D0

         IF (((DABS(TM1-TM)).LT.TOLNST).AND.((DABS(TN1-TN)).LT.TOLNST))
     *                                                             THEN
            TM=TM1
            TN=TN1
            CALL GETREALLOC(P1,P2,P3,P4,TM,TN,SP1)
       H=DSQRT((SP1(1)-CP(1))**2+(SP1(2)-CP(2))**2+(SP1(3)-CP(3))**2)
            RETURN
         ELSE
            TM=TM1
            TN=TN1
            GOTO 981
         END IF
         RETURN
         END

*      SUBROUTINE GETPS(P1,P2,P3,P4,CP,H)
*      IMPLICIT NONE
*      DOUBLE PRECISION,DIMENSION(3):: P1,P2,P3,P4,CP
*      DOUBLE PRECISION,DIMENSION(3):: SP1,SP2,U,V
*      INTEGER CT
*      DOUBLE PRECISION TM,TN,TM1,TN1,TOL,LTH1,LTH2,H
*      TM=0.0D0
*      TN=0.0D0
*      CT=0
*      TOL=1.0D-10
*CPLUS    TM IS FIXED
* 981  CALL GETREALLOC(P1,P2,P3,P4,TM,-1.0D0,SP1)
*      CALL GETREALLOC(P1,P2,P3,P4,TM,1.0D0,SP2)
*      CALL PMREAL(1.0D0,SP2,-1.0D0,SP1,3,U)
*      CALL LEN1REAL(U,3,LTH1)
*      CALL PMREAL(1.0D0,CP,-1.0D0,SP1,3,V)
*      LTH2=U(1)*V(1)+U(2)*V(2)+U(3)*V(3)
*      TN1=LTH2/LTH1*2.0D0-1.0D0
*      IF (TN1.LT.(-1.0D0)) TN1=-1.0D0
*      IF (TN1.GT.(1.0D0)) TN1=1.0D0
*
*      CALL GETREALLOC(P1,P2,P3,P4,-1.0D0,TN1,SP1)
*      CALL GETREALLOC(P1,P2,P3,P4,1.0D0,TN1,SP2)
*      CALL PMREAL(1.0D0,SP2,-1.0D0,SP1,3,U)
*      CALL LEN1REAL(U,3,LTH1)
*      CALL PMREAL(1.0D0,CP,-1.0D0,SP1,3,V)
*      LTH2=U(1)*V(1)+U(2)*V(2)+U(3)*V(3)
*      TM1=LTH2/LTH1*2.0D0-1.0D0
*      IF (TM1.LT.(-1.0D0)) TM1=-1.0D0
*      IF (TM1.GT.(1.0D0)) TM1=1.0D0
*
*
*      IF (((DABS(TM1-TM)).LT.TOL).AND.((DABS(TN1-TN)).LT.TOL)) THEN
*        TM=TM1
*        TN=TN1
*        CALL GETREALLOC(P1,P2,P3,P4,TM,TN,SP1)
*        H=DSQRT((SP1(1)-CP(1))**2+(SP1(2)-CP(2))**2+(SP1(3)-CP(3))**2)
*        RETURN
*      END IF
*
*      CT=CT+1
*      IF (CT.GT.100) THEN
*        H=0.0D0
*        WRITE(*,*) "WARNING: NRPAN-GETPTS"
*        RETURN
*      END IF
*
*      TM=TM1
*      TN=TN1
*      GOTO 981
*      RETURN
*      END SUBROUTINE

      RECURSIVE SUBROUTINE RECURSIVE_AB(P1,P2,P3,P4,CP,
     *                     WHOLEA,WHOLEB,FA,FB,ATOL)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(3)::P1,P2,P3,P4,P5,P6,P7,P8,P9,CP,V1,V2
      DOUBLE PRECISION,DIMENSION(4)::A,B
      DOUBLE PRECISION FA,FB,WHOLEA,WHOLEB,TEMA,TEMB,ATOL,L1,L2
      INTEGER I

      DO I=1,3
         P8(I)=(P1(I)+P4(I))/2
         P6(I)=(P2(I)+P3(I))/2
         P5(I)=(P1(I)+P2(I))/2
         P7(I)=(P3(I)+P4(I))/2
         P9(I)=(P8(I)+P6(I))/2
         V1(I)=P6(I)-P8(I)
         V2(I)=P7(I)-P5(I)
      END DO
      L1=DSQRT(V1(1)**2+V1(2)**2+V1(3)**2)
      L2=DSQRT(V2(1)**2+V2(2)**2+V2(3)**2)

      IF ((L1/L2).GT.(2.0D0)) THEN
         CALL GETAB(P1,P5,P7,P4,CP,TEMA,TEMB)
         A(1)=TEMA
         B(1)=TEMB
         CALL GETAB(P2,P3,P7,P5,CP,TEMA,TEMB)
         A(2)=TEMA
         B(2)=TEMB
         TEMA=A(1)+A(2)
         TEMB=B(1)+B(2)
         IF ((DABS(TEMA-WHOLEA).GT.ATOL).OR.
     *                    (DABS(TEMB-WHOLEB).GT.ATOL)) THEN
            CALL RECURSIVE_AB(P1,P5,P7,P4,CP,A(1),B(1),TEMA,TEMB,ATOL)
            A(1)=TEMA
            B(1)=TEMB
            CALL RECURSIVE_AB(P2,P3,P7,P5,CP,A(2),B(2),TEMA,TEMB,ATOL)
            A(2)=TEMA
            B(2)=TEMB
         END IF
         FA=A(1)+A(2)
         FB=B(1)+B(2)

      ELSEIF ((L2/L1).GT.(2.0D0)) THEN

         CALL GETAB(P1,P2,P6,P8,CP,TEMA,TEMB)
         A(1)=TEMA
         B(1)=TEMB
         CALL GETAB(P3,P4,P8,P6,CP,TEMA,TEMB)
         A(2)=TEMA
         B(2)=TEMB
         TEMA=A(1)+A(2)
         TEMB=B(1)+B(2)
         IF ((DABS(TEMA-WHOLEA).GT.ATOL).OR.
     *                    (DABS(TEMB-WHOLEB).GT.ATOL)) THEN
            CALL RECURSIVE_AB(P1,P2,P6,P8,CP,A(1),B(1),TEMA,TEMB,ATOL)
            A(1)=TEMA
            B(1)=TEMB
            CALL RECURSIVE_AB(P3,P4,P8,P6,CP,A(2),B(2),TEMA,TEMB,ATOL)
            A(2)=TEMA
            B(2)=TEMB
         END IF
         FA=A(1)+A(2)
         FB=B(1)+B(2)


      ELSE
         CALL GETAB(P1,P5,P9,P8,CP,TEMA,TEMB)
         A(1)=TEMA
         B(1)=TEMB
         CALL GETAB(P5,P2,P6,P9,CP,TEMA,TEMB)
         A(2)=TEMA
         B(2)=TEMB
         CALL GETAB(P9,P6,P3,P7,CP,TEMA,TEMB)
         A(3)=TEMA
         B(3)=TEMB
         CALL GETAB(P8,P9,P7,P4,CP,TEMA,TEMB)
         A(4)=TEMA
         B(4)=TEMB

         TEMA=A(1)+A(2)+A(3)+A(4)
         TEMB=B(1)+B(2)+B(3)+B(4)
         IF ((DABS(TEMA-WHOLEA).GT.ATOL).OR.(DABS(TEMB-WHOLEB).GT.ATOL))
     *                                                      THEN
            CALL RECURSIVE_AB(P1,P5,P9,P8,CP,A(1),B(1),TEMA,TEMB,ATOL)
            A(1)=TEMA
            B(1)=TEMB
            CALL RECURSIVE_AB(P5,P2,P6,P9,CP,A(2),B(2),TEMA,TEMB,ATOL)
            A(2)=TEMA
            B(2)=TEMB
            CALL RECURSIVE_AB(P9,P6,P3,P7,CP,A(3),B(3),TEMA,TEMB,ATOL)
            A(3)=TEMA
            B(3)=TEMB
            CALL RECURSIVE_AB(P8,P9,P7,P4,CP,A(4),B(4),TEMA,TEMB,ATOL)
            A(4)=TEMA
            B(4)=TEMB
         END IF
         FA=A(1)+A(2)+A(3)+A(4)
         FB=B(1)+B(2)+B(3)+B(4)
      END IF
      RETURN
      END

      SUBROUTINE GETAB(P1,P2,P3,P4,CP,A,B)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(3)::P1,P2,P3,P4,CP,PP1
      DOUBLE PRECISION,DIMENSION(3)::V13,V24,NV,U,V,RV,PROD
      DOUBLE PRECISION,DIMENSION(3,3,3)::PP
      DOUBLE PRECISION A,B,TM,TN,RLEN,S,APARA(3),TP
      INTEGER I,J,K,IJ
      APARA(1)=16
      APARA(2)=4
      APARA(3)=1
      A=0.0d0
      B=0.0d0

C.....RAEL LOCATION
      DO I=1,3
         PP(1,1,I)=P1(I)
         PP(3,1,I)=P2(I)
         PP(3,3,I)=P3(I)
         PP(1,3,I)=P4(I)
         PP(2,1,I)=(P1(I)+P2(I))/2
         PP(2,3,I)=(P3(I)+P4(I))/2
         PP(1,2,I)=(P1(I)+P4(I))/2
         PP(3,2,I)=(P3(I)+P2(I))/2
         PP(2,2,I)=(PP(2,1,I)+PP(2,3,I))/2
      END DO

C     NORMAL VECTOR
      DO I=1,3
         V13(I)=PP(3,2,I)-PP(1,2,I)
         V24(I)=PP(2,3,I)-PP(2,1,I)
      END DO
      CALL CROSSREAL(V13,V24,NV)
      CALL LEN1REAL(NV,3,S)

      DO I=1,3
         DO J=1,3
            IJ=ABS(I-2)+ABS(J-2)+1

            TM=DBLE(I)-2.0d0
            TN=DBLE(J)-2.0d0

C  FIND THE REAL LOCATION
            RV(1)=CP(1)-PP(I,J,1)
            RV(2)=CP(2)-PP(I,J,2)
            RV(3)=CP(3)-PP(I,J,3)
            RLEN=DSQRT(RV(1)**2+RV(2)**2+RV(3)**2)

C  JACOBIAN
            DO K=1,3
               U(K)=P1(K)*(TN-1.0D0)+P2(K)*(1.0D0-TN)+
     *         P3(K)*(1.0D0+TN)+P4(K)*(-1.0D0-TN)
               V(K)=P1(K)*(TM-1.0D0)+P2(K)*(-1.0D0-TM)+
     *         P3(K)*(1.0D0+TM)+P4(K)*(1.0D0-TM)
            END DO
            CALL CROSSREAL(U,V,PROD)
            CALL LEN1REAL(PROD,3,S)
            TP=S/RLEN*APARA(IJ)
            A=A+(NV(1)*RV(1)+NV(2)*RV(2)+NV(3)*RV(3))*TP/RLEN**2
            B=B+TP
         END DO
      END DO
      A=A/9.0D0/16.0D0
      B=B/9.0D0/16.0D0

      RETURN
      END

      SUBROUTINE GETREALLOC(P1,P2,P3,P4,TM,TN,PP)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(3)::P1,P2,P3,P4,PP
      DOUBLE PRECISION TM,TN,F1,F2,F3,F4
      INTEGER I
      F1=(1-TM)*(1-TN)/4;
      F2=(1+TM)*(1-TN)/4;
      F3=(1+TM)*(1+TN)/4;
      F4=(1-TM)*(1+TN)/4;
      PP(1)=P1(1)*F1+P2(1)*F2+p3(1)*F3+P4(1)*F4;
      PP(2)=P1(2)*F1+P2(2)*F2+p3(2)*F3+P4(2)*F4;
      PP(3)=P1(3)*F1+P2(3)*F2+p3(3)*F3+P4(3)*F4;

      RETURN
      END
      SUBROUTINE LEN1REAL(X,N,L)
      IMPLICIT NONE
      DOUBLE PRECISION X,L,X1
      INTEGER N,I
      DIMENSION X(N)
      L=0.0D0;
      DO I=1,N
         L=L+X(I)**2
      END DO
      L=DSQRT(L)
      DO I=1,N
         X(I)=X(I)/L
      END DO
      RETURN
      END
      SUBROUTINE PMREAL(A,X1,B,X2,N,X3)
      IMPLICIT NONE
      DOUBLE PRECISION A,B,X1,X2,X3
      INTEGER N,I
      DIMENSION X1(N),X2(N),X3(N)
      DO I=1,N
         X3(I)=A*X1(I)+B*X2(I)
      END DO
      RETURN
      END
      SUBROUTINE CROSSREAL(U,V,W)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(3)::U,V,W
      W(1)= U(2)*V(3)-U(3)*V(2)
      W(2)=-U(1)*V(3)+U(3)*V(1)
      W(3)= U(1)*V(2)-U(2)*V(1)
      RETURN
      END
