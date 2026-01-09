C
C
C
      SUBROUTINE ENPROD(U,V,W)
C-----------------------------------------------------------------------
C                                       _   _ 
C     EVALUATES THE INNER PRODUCT   W = U . V
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(3),V(3)
      W = U(1)*V(1)+U(2)*V(2)+U(3)*V(3)
      RETURN
      END
