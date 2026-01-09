C
C
C
      FUNCTION RADLE(XTI0,XCHD0,ITHK)
C***********************************************************************
C     Leading edge radius
C
      DIMENSION RLEC(6)
      DATA RLEC/0.448,0.5,0.636,0.687,1.1019,0.5/
      RADLE=RLEC(ITHK)*(XTI0/XCHD0)**2
      RETURN
C)))))))))))))))))))))) End of function RADLE ((((((((((((((((((((((((((
      END
