C
C
C
      FUNCTION RNOSE(X,RHUB)
C***********************************************************************
C     Shape of the hub nose
C
      DATA ZERO/0.0/
      RNOSE=RHUB*16./3.*(.5*X-3./8*X*X+X**4/16)
      IF(X.GT.1.0) THEN
         RNOSE=RHUB
      ELSE IF(X.LT.0.0) THEN
         RNOSE=ZERO
      ELSE 
         RNOSE=RHUB*16./3.*(.5*X-3./8*X*X+X**4/16)
      END IF 
      RETURN
C))))))))))))))))))))))) End of function RNOSE (((((((((((((((((((((((((
      END
