C ====================================================
      SUBROUTINE ROTATE2(IP,ANG,X,Y)
C ====================================================

      TANG = ANG
      IF(IP .EQ. -1) TANG = -ANG
      
      XTMP = X
      YTMP = Y
      X = XTMP * COS(TANG) - YTMP * SIN(TANG)
      Y = XTMP * SIN(TANG) + YTMP * COS(TANG)
      
      RETURN
      END
    
         
