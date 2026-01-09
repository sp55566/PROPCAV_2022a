C ========================================================
      SUBROUTINE EVALDK(NIN,NOUT,XIN,XOUT,YOUT,A)
C ========================================================
      DIMENSION XIN(NIN),XOUT(NOUT),YOUT(NOUT),A(4*(NIN-1))
      NM1=NIN-1
      MOUT=IABS(NOUT)
      IF(NOUT.GT.0) GO TO 1
      DEL=(XIN(NIN)-XIN(1))/(MOUT-1)
      DO 2 N=1,MOUT
 2    XOUT(N)=XIN(1)+(N-1)*DEL
 1    J=1
      DO 3 N=1,MOUT
      IF(XOUT(N).GE.XIN(2)) GO TO 4
      J=1
      GO TO 5
 4    IF(XOUT(N).LT.XIN(NM1)) GO TO 6
      J=NM1
      GO TO 5
 6    IF(XOUT(N).GE.XIN(J+1)) GO TO 7
    9 IF (XOUT(N).LT.XIN(J)) GO TO 8
 5    H1=XOUT(N)-XIN(J)
      H2=H1**2
      H3=H1*H2
      J2=J+NM1
      J3=J2+NM1
      J4=J3+NM1
      YOUT(N)=A(J)*H3+A(J2)*H2+A(J3)*H1+A(J4)
      GO TO 3
 7    J=J+1
      GO TO 6
    8 J=J-1
      GO TO 9
 3    CONTINUE
      RETURN
C))))))))))))))))))))) End of subroutine EVALDK ((((((((((((((((((((((((
      END


C ========================================================
      SUBROUTINE EVALDKs(NIN,NOUT,XIN,XOUT,YOUT,A)
C ========================================================
      REAL XOUT, YOUT
      DIMENSION XIN(NIN),A(4*(NIN-1))
      NM1=NIN-1
      MOUT=IABS(NOUT)
      IF(NOUT.GT.0) GO TO 1
      DEL=(XIN(NIN)-XIN(1))/(MOUT-1)
      DO 2 N=1,MOUT
 2    XOUT=XIN(1)+(N-1)*DEL
 1    J=1
      DO 3 N=1,MOUT
      IF(XOUT.GE.XIN(2)) GO TO 4
      J=1
      GO TO 5
 4    IF(XOUT.LT.XIN(NM1)) GO TO 6
      J=NM1
      GO TO 5
 6    IF(XOUT.GE.XIN(J+1)) GO TO 7
    9 IF (XOUT.LT.XIN(J)) GO TO 8
 5    H1=XOUT-XIN(J)
      H2=H1**2
      H3=H1*H2
      J2=J+NM1
      J3=J2+NM1
      J4=J3+NM1
      YOUT=A(J)*H3+A(J2)*H2+A(J3)*H1+A(J4)
      GO TO 3
 7    J=J+1
      GO TO 6
    8 J=J-1
      GO TO 9
 3    CONTINUE
      RETURN
C))))))))))))))))))))) End of subroutine EVALDK ((((((((((((((((((((((((
      END
