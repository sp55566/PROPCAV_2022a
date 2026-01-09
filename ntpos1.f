       FUNCTION NTPOS1(KB,IT)
C
C     This function figures out the position of blade KB in terms of 
C       the time step at certain time step IT
C
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

       NTOT = 360 / NDLTAT
       NADV = 360 / NBLADE / NDLTAT
       KEYPOS = IT / NDLTAT + 1
       
       NTPOS1=KEYPOS+(KB-1)*NADV
       IF(NTPOS1.GT.NTOT) THEN
          NTPOS1=NTPOS1-NTOT
       END IF
       
       RETURN
C)))))))))))))))))))))) End of function NTPOS1 (((((((((((((((((((((((((
       END
