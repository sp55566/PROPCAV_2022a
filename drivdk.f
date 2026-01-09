C
C
C
      SUBROUTINE DRIVDK(NIN,NOUT,XIN,XOUT,DYDX,D2YDX,A)
      DIMENSION XIN(*),XOUT(*),DYDX(*),D2YDX(*),A(*)
      DATA TWO,THREE,SIX/2.0,3.0,6.0/
      NM1=NIN-1
      J=1
      DO 3 N=1,NOUT
      IF(XOUT(N).GE.XIN(2)) GO TO 4
      J=1
      GO TO 5
 4    IF(XOUT(N).LT.XIN(NM1)) GO TO 6
      J=NM1
      GO TO 5
 6    IF(XOUT(N).GE.XIN(J+1)) GO TO 7
 5    H1=XOUT(N)-XIN(J)
      H2=H1**2
      J2=J+NM1
      J3=J2+NM1
      DYDX(N)=THREE*A(J)*H2+TWO*A(J2)*H1+A(J3)
      D2YDX(N)=SIX*A(J)*H1+TWO*A(J2)
      GO TO 3
 7    J=J+1
      GO TO 6
 3    CONTINUE
      RETURN
C))))))))))))))))))))) End of subroutine DRIVDK ((((((((((((((((((((((((
      END
