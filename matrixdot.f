C     THE FOLLOWING SUBROUTINES ARE COPIED FROM CAV2DBL BY HONG SUN
C
C***********************************************************************
      SUBROUTINE MAT2DOT(A,B,C,N1,NP1,N2,NP2,N3,NP3)
C
C     COMPUTE A DOT B = C
C     A - (N1,N2) actually; (NP1,NP2) physically
C     B - (N2,N3) actually; (NP2,NP3) physically
C     C - (N1,N3) actually; (NP1,NP3) physically

      DIMENSION A(NP1,NP2),B(NP2,NP3),C(NP1,NP3)

      SUM = 0.
      DO 10 I = 1, N1
          DO 20 J = 1, N3
            DO 30 K = 1,N2 
            SUM = SUM  + A(I,K)*B(K,J)
 30         CONTINUE
          C(I,J) = SUM
          SUM = 0.
 20      CONTINUE
 10     CONTINUE

      RETURN
      END

C***********************************************************************
      SUBROUTINE MAT1DOT(A,B,C,N1,NP1,N2,NP2)
C
C     COMPUTE A DOT B = C
C     A - (N1,N2) actually; (NP1,NP2) physically
C     B - (N2) actually; (NP2,1) physically
C     C - (N1) actually; (NP1,1) physically

      DIMENSION A(NP1,NP2),B(NP2),C(NP1)

      SUM = 0.
      DO 10 I = 1, N1
          DO 20 K = 1, N2 
          SUM = SUM  + A(I,K)*B(K)
 20       CONTINUE
        C(I) = SUM              
 10    CONTINUE
      return
      end      
