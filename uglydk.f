C))))))))))))))))) End of PUF10.02 geometry library ((((((((((((((((((((
C
C
C
C***********************************************************************
C     1975 DUCK SERIES  J.E.KERWIN
C***********************************************************************
      SUBROUTINE UGLYDK(NIN,NCL,NCR,XIN,YIN,ESL,ESR,AE)
      real XIN(NIN),YIN(NIN),AE(4*(NIN-1))
      real H(NIN-1),D(NIN-1),AU(NIN-1),AM(NIN-1),AL(NIN-1)
      real S(NIN),X(NIN-1)
      DATA HALF/0.5/,TWO/2.0/,SIX/6.0/,RAD/1.745329252E-02/
      DATA ONE/1.0/
      NM1=NIN-1
      NM2=NM1-1
      NM3=NM2-1
      NEQ=NM2
      DO 1 N=1,NM1
      H(N)=XIN(N+1)-XIN(N)
 1    D(N)=(YIN(N+1)-YIN(N))/H(N)
      IF(NCL.EQ.2) NEQ=NEQ+1
      IF(NCR.EQ.2) NEQ=NEQ+1
      NSQ=NEQ**2
      J=1
      IF(NCL.LT.2) GO TO 6
      AM(1)=TWO*H(1)
      AU(1)=H(1)
      SLP=ESL*RAD
      S(1)=(D(1)-TAN(SLP))*SIX
      J=J+1
      AL(2)=H(1)
 6    DO 5 N=1,NM2
      IF(N.GT.1) AU(J-1)=H(N)
      AM(J)=TWO*(H(N)+H(N+1))
      IF(N.LT.NM2) AL(J+1)=H(N+1)
      IF(N.EQ.2.AND.NCL.EQ.1) AU(J-1)=AU(J-1)-H(N-1)**2/H(N)
      IF(N.EQ.1.AND.NCL.EQ.1) AM(J)=AM(J)+(ONE+H(N)/H(N+1))*H(N)
      IF(N.EQ.NM2.AND.NCR.EQ.1) AM(J)=AM(J)+(ONE+H(N+1)/H(N))*H(N+1)
      IF(N.EQ.NM3.AND.NCR.EQ.1) AL(J+1)=AL(J+1)-H(N+2)**2/H(N+1)
      S(J)=(D(N+1)-D(N))*SIX
      J=J+1
 5    CONTINUE
      IF(NCR.LT.2) GO TO 7
      AL(NEQ)=-H(NM1)
      AM(NEQ)=-TWO*H(NM1)
      AU(NEQ-1)=H(NM1)
      SLP=ESR*RAD
      S(J)=(D(NM1)+TAN(SLP))*SIX
 7    CONTINUE
      DO 4 K=2,NEQ
      AL(K)=AL(K)/AM(K-1)
      AM(K)=AM(K)-AL(K)*AU(K-1)
      S(K)=S(K)-AL(K)*S(K-1)
 4    CONTINUE
      X(NEQ)=S(NEQ)/AM(NEQ)
      DO 2 L=2,NEQ
      K=NEQ-L+1
      X(K)=(S(K)-AU(K)*X(K+1))/AM(K)
 2    CONTINUE
      DO 22 N=1,NEQ
 22   S(N)=X(N)
      HOLD=S(NEQ)
      IF(NCL.EQ.2) GO TO 8
      DO 9 N=1,NM2
      M=NM2-N+2
 9    S(M)=S(M-1)
      IF(NCL.EQ.0) S(1)=0.0
      BUG=H(1)/H(2)
      IF(NCL.EQ.1) S(1)=(1.0+BUG)*S(2)-BUG*S(3)
 8    IF(NCR.EQ.0) S(NIN)=0.0
      BUG=H(NM1)/H(NM2)
      IF(NCR.EQ.1) S(NIN)=(1.0+BUG)*S(NM1)-BUG*S(NM2)
      IF(NCR.EQ.2) S(NIN)=HOLD
      DO 10 N=1,NM1
      AE(N)=(S(N+1)-S(N))/(SIX*H(N))
      M=N+NM1
      AE(M)=HALF*S(N)
      M=M+NM1
      AE(M)=D(N)-H(N)*(TWO*S(N)+S(N+1))/SIX
      M=M+NM1
 10   AE(M)=YIN(N)
      RETURN
C))))))))))))))))))))) End of subroutine UGLYDK ((((((((((((((((((((((((
      END






