      FUNCTION INDEXT(N,M)
C***********************************************************************
C     INDEXT: INDEX of the TIP hub panels

      INCLUDE 'PUFCAV.INC'

      indext = npanb + npanh + npand+ npantn + (n-1)*mcvt + m

      RETURN
      END
