C
C
C
      FUNCTION GROW(X)
C***********************************************************************
C     GROW: GROWth function
C      --- Assumed growth fuction for induced velocities in the 
C          transition wake
C
      GROW=3.*X-3.*X*X+X**3
      IF(X.GT.1.0) GROW=1.0
      IF(X.LT.0.0) GROW=0.0
      RETURN
C))))))))))))))))))))))) End of function GROW  (((((((((((((((((((((((((
      END
