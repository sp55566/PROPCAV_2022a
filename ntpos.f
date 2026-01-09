
      FUNCTION NTPOS(KB)
C
C     This function figures out the position of blade KB in terms of 
C       the time step at certain time step
C
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
C      NTOT=TWOPI/DELTAT
C      NADV=DELK/DELTAT
      NTOT = 360 / NDLTAT
      NADV = 360 / NBLADE / NDLTAT
      KEYPOS = ITSTEP / NDLTAT + 1
C      KEYPOS=NINT(TSTEP/DELTAT)+1
      NTPOS=KEYPOS+(KB-1)*NADV
      IF(NTPOS.GT.NTOT) THEN
         NTPOS=NTPOS-NTOT
      END IF
      RETURN
C)))))))))))))))))))))) End of function NTPOS (((((((((((((((((((((((((
      END
