C
C
C
      SUBROUTINE EXPROD(U,V,W)
C-----------------------------------------------------------------------
C                                      _   _   _
C     EVALUATES THE EXTERIOR PRODUCT   W = U x V
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(3),V(3),W(3)
      W(1)= U(2)*V(3)-U(3)*V(2)
      W(2)=-U(1)*V(3)+U(3)*V(1)
      W(3)= U(1)*V(2)-U(2)*V(1)
      RETURN
      END
