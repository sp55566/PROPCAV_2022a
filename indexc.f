      FUNCTION INDEXC(N,M)
C***********************************************************************
C     INDEXC: INDEX of the TIP vortex cavity panels

      INCLUDE 'PUFCAV.INC'

      indexc = npanb + npanh + npand + npantn + npant 
     %         + (n-1)*mcvt + m

      RETURN
      END
