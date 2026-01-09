C ----------------------------
      FUNCTION INDEXD(N,M)
C ----------------------------
      INCLUDE 'PUFCAV.INC'

      INDEXD = NPANB + NPANH + (M - 1) * NDUCT + N

      RETURN
      END
