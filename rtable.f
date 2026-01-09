      SUBROUTINE RTABLE(RULT,RHULT,DCD,RTBL,XTBL,XULT,RW,MM,
     1                  XHBTE,XHBT,XHBFD)
C***********************************************************************
C     RTABLE: Radii TABLE of the transition wake geometry
C       --- Generate a table of radii for contracted transition wake
C
      DIMENSION RTBL(11,MM),XTBL(11),RW(MM)
      ZERO=0.0
      MM1=MM-1
      DO 10 M=1,MM
         RTBL(1,M)=RW(M)
10    CONTINUE
      DC=-TAN(1.745329E-02*DCD)
      H=RW(MM)-RULT
      IF(DCD.GT.1.0) THEN
         XULTUN=-3.0*H/DC
      ELSE
         XULTUN=XULT
      END IF
      XULTUN=AMAX1(XULTUN,0.3)
      XTBL(1)=ZERO
      RHUB=RW(1)
      DO 40 N=2,11
         XTBL(N)=0.1*(N-1)*XULT
         P=XTBL(N)/XULT
         Q=XTBL(N)/XULTUN
         X=XHBTE+XULT*P
         IF(X.LT.XHBFD) THEN
            PQ=(XHBFD-X)/XHBT
            RHWAKE=RNOSE(PQ,RHUB)
            IF(RHWAKE.LT.RHULT) THEN
               RHWAKE=RHULT
               GO TO 20
            END IF
         ELSE
            RHWAKE=RHULT
         END IF
20       RTBL(N,1)=RHWAKE
         RTBL(N,MM)=RTBL(1,MM)+H*(-3.0*Q+3.0*Q**2-Q**3)
         IF(Q.GT.1.0) THEN
            RTBL(N,MM)=RULT
         END IF
         DO 30 M=2,MM1
            D1=SQRT(RTBL(N,M-1)**2+RTBL(1,M)**2-RTBL(1,M-1)**2)
     *        -RTBL(N,M-1)
            IF(M.GT.MM/2+2) THEN
               D1=AMIN1( D1,(RTBL(N,M-1)/RTBL(N,MM))**2
     *           *(RTBL(N,MM)-RTBL(N,M-1)) )
            END IF 
            RTBL(N,M)=RTBL(N,M-1)+D1
30       CONTINUE
40    CONTINUE
      RETURN
C))))))))))))))))))))) End of subroutine RTABLE ((((((((((((((((((((((((
      END
