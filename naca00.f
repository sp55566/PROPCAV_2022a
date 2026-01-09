C
C
C
      SUBROUTINE NACA00(NTH,THK,RLE,SC,YTC)
C***********************************************************************
C     NACA 4-digit section
C
      DIMENSION SC(*),YTC(*)
c      DATA A,B,C,D,E/1.4845,-0.63,-1.758,1.4215,-.5075/
c Change 3rd coefficient to get finite TE thickness so panel method
c   with Morino Kutta condition doesn't bomb
      DATA A,B,C,D,E/1.4845,-0.63,-1.7685,1.4215,-.5075/

      RLE=1.1019*THK**2
C
C.....Section coordinates and slope at nodal points
      DO 10 N=1,NTH
         A1=SQRT(SC(N))
         A2=SC(N)
         A3=A2*A2
         A4=A3*A2
         A5=A3*A3
         YTC(N)=THK*(A*A1+B*A2+C*A3+D*A4+E*A5)
10    CONTINUE

c      YTC(NTH)=0.
      RETURN
C))))))))))))))))))))) End of subroutine NACA00 ((((((((((((((((((((((((
      END
