      FUNCTION INDEXW(N,M)
C***********************************************************************
C     INDEXW: INDEX of the Wake panels
C
      INCLUDE 'PUFCAV.INC'
c      SAVE
      INDEXW=0
      DO 10 M1=MR,M+1,-1
        INDEXW=INDEXW+NWPAN(M1)
   10 CONTINUE
      INDEXW=INDEXW+N
      RETURN
C)))))))))))))))))))))) End of function INDEXW (((((((((((((((((((((((((
      END
