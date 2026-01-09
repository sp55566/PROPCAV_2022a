      FUNCTION FILLIN(X,AB,OR,NO)
C     *** FILLIN *** PARABOLIC INTERPOLATION
C     FIND Y(X) FROM TABLE OF
C     AB(N) AND OR(N) CONTAINING NO POINTS.
      DIMENSION AB(15),OR(15)
      ANTRA(X1,X2,X3,X,Y1,Y2,Y3)=Y1*(X-X2)*(X-X3)/((X1-X2)*(X1-X3))+
     1 Y2*(X-X1)*(X-X3)/((X2-X1)*(X2-X3))+Y3*(X-X1)*(X-X2)/((X3-X1)*
     2 (X3-X2))
      IF(X-AB(1))  1,3,2
 3    Y=OR(1)
      GO TO 99
 1    Y=ANTRA(AB(1),AB(2),AB(3),X,OR(1),OR(2),OR(3))
      GO TO 99
 2    IF(X-AB(2))1,6,5
 6    Y=OR(2)
      GO TO 99
 5    DO 7 I=3,NO
      M=I
      IF(X-AB(I))8,9,7
 9    Y=OR(I)
      GO TO 99
 7    CONTINUE
 8    Y=ANTRA(AB(M-2),AB(M-1),AB(M),X,OR(M-2),OR(M-1),OR(M))
 99   FILLIN=Y
      RETURN
      END

