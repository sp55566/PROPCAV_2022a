
      SUBROUTINE ARCLEN(S1,S2,X1,Y1,X2,Y2,X3,Y3,D1,D2)
C-----------------------------------------------------------------------
C     ARC-LENGTH OF THE CIRCLE PASSING THROUGTH THREE POINTS
C     Copied from CAV2D-BL                          Hong Sun 
C-----------------------------------------------------------------------
      DATA TOL/1.0E-8/

      PI = ATAN(1.0)*4.0
      DX=X3-X1
      DY=Y3-Y1
      D=SQRT(DX*DX+DY*DY)
      A=ABS( ATAN2(Y2-Y1,X2-X1)-ATAN2(Y3-Y2,X3-X2) )
      IF(A.GT.PI) A=A-PI
      IF(A.LT.TOL) THEN
                      S1=D1
                      S2=D2
                   ELSE
                      TR=D/SIN(A)
                      S1=TR*ASIN(D1/TR)
                      S2=TR*ASIN(D2/TR)
                   END IF
      RETURN
      END
