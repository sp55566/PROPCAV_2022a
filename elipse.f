C ===================================================
      SUBROUTINE ELIPSE(NTH,THK,RLE,SC,YTC)
C
C     Elliptical thickness form
C ===================================================
C 
      DIMENSION SC(*),YTC(*)
      DATA HALF,ONE/0.5,1.0/

      RLE=HALF*THK**2
      DO 10 N=1,NTH
         YTC(N)=THK*SQRT(SC(N)*(ONE-SC(N)))
10    CONTINUE
      RETURN
      END
