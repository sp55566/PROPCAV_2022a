
      SUBROUTINE CHRLEN(A,LENCH) 
      CHARACTER A*(*)
      N =LEN(A)
      DO 10 I=1,N 
         IF(A(I:I).EQ.' ') THEN
            LENCH=I-1
            GO TO 20 
         END IF
10    CONTINUE    
20    CONTINUE
      RETURN
      END
