      SUBROUTINE RPAN_DBL(Xi,Yi,Zi,CHRLENSI,FSS,FDD,FSXI,FSYI,FDXI,FDYI,
     %                 FDZI,ID,IFLAG,XV1,YV1,SQ1,SIDE1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XV1(4),YV1(4),SQ1(15),SIDE1(4)
      DIMENSION XV(4),YV(4),S(15),SIDE(4)
      DIMENSION R(4),RR(4),RI(4),XRI(4),YRI(4),FE(4),
     *          B(5),XMXV(4),YMYV(4),N1(4)

C----------------------------------------------------------------------
C
C  DATA ARRAY B FOR 6D RATIONAL APPROXIMATION OF ARCTANGENT
C  SEE HART, ET AL, COMPUTER APPROXIMATIONS, WILEY, 1968, TABLE 5090
C
C----------------------------------------------------------------------
      DATA B/ 2.4091197D-01, 3.7851122D+00, 5.6770721D+00,
     *        5.6772854D+00, 5.6770747D+00/
C      DATA PI/ 3.1415927D+00 /, PI2/ 1.5707963D+00 /
C      DATA TWOPI/ 6.2831853D+00 /
      DATA ONE/ 1.D+00/, A3/ 3.D+00/, A5/ 5.D+00/
      DATA A7/ 7.D+00/, A9/ 9.D+00/, A11/ 11.D+00/
      DATA A14/ 14.D+00/, A35/ 35.D+00/
      DATA A49/ 49.D+00/, A63/ 63.D+00/, A99/ 99.D+00/
C      DATA ONE10/.1D+00/, ONE6/.1666667D+00/, ONE3/ .3333333D+00/
C      DATA FIVE3/ 1.666667D+00/, SEVEN3/2.333333D+00/, ZERO/ 0.0D0/
      DATA ZERO/ 0.0D0/
      DATA FT3/ 4.666667D+00/, TOL/ 1D-08 /
      DATA N1/ 2, 3, 4, 1 /

      pi = dacos(-1.d0)
      pi2 = 0.5d0*pi
      twopi = 2.d0*pi
      one10 = 0.1d0
      one6 = 1.d0/6.d0
      one3 = 1.d0/3.d0
      five3 = 5.d0/3.d0
      seven3 = 7.d0/3.d0

      x = xi
      y = yi
      z = zi
      chrlens = chrlensi
      do i = 1 , 4
      xv(i) = xv1(i)
      yv(i) = yv1(i)
      side(i) = side1(i)
      enddo
      do i = 1 , 15
        s(i) = sq1(i)
      enddo

      XMXC=X-S(6)
      YMYC=Y-S(2)
      XX=XMXC*XMXC
      YY=YMYC*YMYC
      ZZ=Z*Z
      RRC=XX+YY+ZZ
      IF (RRC.LT.100.d0*CHRLENS) THEN
         IF(IFLAG.EQ.1) THEN
            GO TO 11
         ELSE IF (IFLAG.EQ.0) THEN
            IFLAG=2
            RETURN
         END IF
      END IF
C----------------------------------------------------------------------
C
C   TWO-TERM MULTIPOLE EXPANSION INCLUDING SECOND MOMENTS
C
C----------------------------------------------------------------------
      R2=ONE/RRC
      R1=DSQRT(R2)
      R3=R1*R2
      R5=R3*R2
      ZR2=Z*R2
      XY=XMXC*YMYC
      SS1=S(1)*R1
      SS3=-(S(3)+S(10))*R3
      SS5=(XX*S(10)+XY*S(7)+YY*S(3))*R5
      FS=SS1+ONE3*SS3+SS5
      FDSUM=SS1+SS3+A5*SS5
      FD=ZR2*FDSUM

C----------------------------------------------------------------------
C    DELETE NEXT  14  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 8
      RSS3=R2*SS1
      SSX3=-XMXC*RSS3
      SSY3=-YMYC*RSS3
      SSX5=(XMXC*(S(3)+A3*S(10))+YMYC*S(7))*R5
      SSY5=(YMYC*(S(10)+A3*S(3))+XMXC*S(7))*R5
      A5R2=A5*R2
      RSS7=-A5R2*SS5
      SSX7=XMXC*RSS7
      SSY7=YMYC*RSS7
      FSX=SSX3+SSX5+SSX7
      FSY=SSY3+SSY5+SSY7
      FDX=ZR2*(A3*SSX3+A5*SSX5+A7*SSX7)
      FDY=ZR2*(A3*SSY3+A5*SSY5+A7*SSY7)
      ZZR4=ZR2*ZR2
      FDZ=R2*FDSUM-ZZR4*(A3*SS1+A5*SS3+A35*SS5)
  8   IF (RRC.GT.15.d0*CHRLENS) GO TO 99
C----------------------------------------------------------------------
C
C    THIRD AND FOURTH MOMENTS ADDED FOR RRC/AREA BETWEEN 40 AND 150
C
C----------------------------------------------------------------------
      S914=S(9)+S(14)
      S813=S(8)+S(13)
      S411=S(4)+S(11)
      S512=S(5)+S(12)
      S1215=S(12)+S(15)
      R7=R5*R2
      R9=R7*R2
      SS5=(-XMXC*S813-YMYC*S411+ONE10*(S512+S1215))*R5
      SS7=(FIVE3*((XMXC*XX*S(13)+YMYC*YY*S(4))+A3*XY*(XMXC*S(11)
     *+YMYC*S(8)))-XX*S1215-YY*S512-XY*S914)*R7
      SS9=(A7*(ONE6*(XX*XX*S(15)+YY*YY*S(5))+XX*YY*S(12))
     *+SEVEN3*XY*(XX*S(14)+YY*S(9)))*R9
      FS=FS+SS5+SS7+SS9
      FDSUM=A5*SS5+A7*SS7+A9*SS9
      FD=FD+ZR2*FDSUM

C----------------------------------------------------------------------
C    DELETE NEXT  20  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 99
      TXY=XY+XY
      SSX5=-S813*R5
      SSY5=-S411*R5
      RSS7=A5R2*SS5
      SSX7=(A5*(XX*S(13)+TXY*S(11)+YY*S(8))-S1215*(XMXC+XMXC)
     *-YMYC*S914)*R7-XMXC*RSS7
      SSY7=(A5*(YY*S(4)+XX*S(11)+TXY*S(8))-S512*(YMYC+YMYC)
     *-XMXC*S914)*R7-YMYC*RSS7
      RSS9=A7*SS7*R2
      SSX9=(FT3*XMXC*XX*S(15)+A14*XMXC*YY*S(12)+A49*YMYC*(XX*S(14)
     *+ONE3*YY*S(9)))*R9-XMXC*RSS9
      SSY9=(FT3*YMYC*YY*S(5)+A14*YMYC*XX*S(12)+A49*XMXC*(YY*S(9)
     *+ONE3*XX*S(14)))*R9-YMYC*RSS9
      RSS11=A9*SS9*R2
      SSX11=-XMXC*RSS11
      SSY11=-YMYC*RSS11
      FSX=FSX+SSX5+SSX7+SSX9+SSX11
      FSY=FSY+SSY5+SSY7+SSY9+SSY11
      FDX=FDX+ZR2*(A5*SSX5+A7*SSX7+A9*SSX9+A11*SSX11)
      FDY=FDY+ZR2*(A5*SSY5+A7*SSY7+A9*SSY9+A11*SSY11)
      FDZ=FDZ+R2*FDSUM-ZZR4*(A35*SS5+A63*SS7+A99*SS9)
      GO TO 99
C----------------------------------------------------------------------
C
C    NEAR-FIELD SECTION USES EXACT FORMULATION
C      SET Z=TOL IF Z.LT.TOL TO AVOID INDETERMINACY ON PANEL
C      ZVTX IS USED TO DETERMINE PROXIMITY TO VERTEX NORMALS
C      MFLAG=1 IF NEAR VERTEX NORMALS
C
C----------------------------------------------------------------------
  11  FD=0.0d0
      FS=0.0d0
      ABZ=DABS(Z)
      IF(ABZ.GT.TOL) GO TO 12
        Z=TOL
        ZZ=Z*Z
        ABZ=TOL
  12  ZVTX=1.005D0*ABZ
      MFLAG=0
C----------------------------------------------------------------------
C    DELETE NEXT 5 LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      FSX=0.0d0
      FSY=0.0d0
      FDX=0.0d0
      FDY=0.0d0
      FDZ=0.0d0
C----------------------------------------------------------------------
C
C    LOOP FOR CORNER FUNCTIONS
C
C----------------------------------------------------------------------
      DO N=1,4
        XMXV(N)=X-XV(N)
        YMYV(N)=Y-YV(N)
        XX=XMXV(N)*XMXV(N)
        YY=YMYV(N)*YMYV(N)
        FE(N)=ZZ+XX
        RR(N)=FE(N)+YY
        R(N)=DSQRT(RR(N))
        IF (R(N).LT.ZVTX) MFLAG=1
        RI(N)=ONE/R(N)
C----------------------------------------------------------------------
C    DELETE NEXT   3   LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 13
        XRI(N)=RI(N)*XMXV(N)
        YRI(N)=RI(N)*YMYV(N)
   13   CONTINUE
      END DO
C----------------------------------------------------------------------
C
C    LOOP FOR SIDE FUNCTIONS AND SUMS OVER FOUR SIDES
C
C----------------------------------------------------------------------
      DO N=1,4
        IF (SIDE(N).LT.TOL) GO TO 33
        SIDI=ONE/SIDE(N)
        CT=(XV(N1(N))-XV(N))*SIDI
        ST=(YV(N1(N))-YV(N))*SIDI
        V=XMXV(N)*ST-YMYV(N)*CT
        VV=V*V
        RADS=VV+ZZ
        U1=XMXV(N)*CT+YMYV(N)*ST
        U2=XMXV(N1(N))*CT+YMYV(N1(N))*ST
        RSUM=R(N)+R(N1(N))
        FLAG=RI(N)*RI(N1(N))*U1*U2
C----------------------------------------------------------------------
C
C       FLAG=1 ON EXTENSIONS, -1 ON SIDES
C         IN FOLLOWING SUBSECTIONS FS,FSX,FSY,FDZ ARE EVALUATED FROM
C           LAST FORM OF (3.9) IN NORMAL CASE
C         ELSE
C           FIRST FORM OF (3.9) IF NEAR SIDE OF PANEL
C
C----------------------------------------------------------------------
        IF (FLAG.GT.-.99D0) THEN
          RSP=RSUM+SIDE(N)
          RSM=RSUM-SIDE(N)
          FLN=DLOG(RSP/RSM)
        ELSE
          RU1=R(N)+U1
          RU2=R(N1(N))-U2
          RADI=ONE/RADS
          FLN=DLOG(RU1*RU2*RADI)
        ENDIF

        FS=FS+V*FLN

C----------------------------------------------------------------------
C    DELETE FOLLOWING LINES TO NEXT ENDIF TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 14
        IF (FLAG.GT.-.99D0) THEN
          FAC=V*(SIDE(N)+SIDE(N))/(RSP*RSM)
          FSX=FSX+FLN*ST-FAC*(XRI(N)+XRI(N1(N)))
          FSY=FSY-FLN*CT-FAC*(YRI(N)+YRI(N1(N)))
          FDZ=FDZ-FAC*(RI(N)+RI(N1(N)))
        ELSE
          RU1I=ONE/RU1
          RU2I=ONE/RU2
          FA=RU1I-RU2I
          FB=-(V+V)*RADI
          FSX=FSX+FLN*ST+V*(FA*CT+FB*ST+RU1I*XRI(N)+RU2I*XRI(N1(N)))
          FSY=FSY-FLN*CT+V*(FA*ST-FB*CT+RU1I*YRI(N)+RU2I*YRI(N1(N)))
          FDZ=FDZ+FB+V*(RU1I*RI(N)+RU2I*RI(N1(N)))
        ENDIF
C----------------------------------------------------------------------
C
C        IN FOLLOWING SUBSECTIONS FACTORS IN (2.15) ARE EVALUATED FROM
C          (2.7) IN NORMAL CASE
C        ELSE
C          (2.14) IF NEAR NORMAL TO A VERTEX
C
C----------------------------------------------------------------------
   14   IF (MFLAG.EQ.0) THEN
          S1=V*R(N)
          C1=ABZ*U1
          S2=V*R(N1(N))
          C2=ABZ*U2
        ELSE
          FH1=XMXV(N)*YMYV(N)
          FH2=XMXV(N1(N))*YMYV(N1(N))
          S1=FE(N)*ST-FH1*CT
          C1=ABZ*R(N)*CT
          S2=FE(N1(N))*ST-FH2*CT
          C2=ABZ*R(N1(N))*CT
        ENDIF
        S12=S1*C2-S2*C1
        C12=C1*C2+S1*S2
C----------------------------------------------------------------------
C
C    EVALUATE THIRD ARCTANGENT IN (2.15)
C         ANGLE (MODULO PI) BETWEEN -PI/4 AND PI/4
C       ELSE
C         USE INVERSE COTANGENT AND ADD/SUBTRACT PI/2
C
C----------------------------------------------------------------------
        IF (DABS(S12).LE.DABS(C12)) THEN
          U=S12/C12
          IF (C12.LT.ZERO) FD=FD+DSIGN(PI,S12)
        ELSE
          U=-C12/S12
          FD=FD+DSIGN(PI2,S12)
        ENDIF
        UU=U*U
        FD=FD+U*((B(1)*UU+B(2))*UU+B(3))/((UU+B(4))*UU+B(5))
C----------------------------------------------------------------------
C
C    FOLLOWING THREE SECTIONS EVALUATE FDX,FDY FOR
C       FIELD POINT NEAR NORMAL TO VERTEX (MFLAG=1)
C       NORMAL CASE (FLAG.LT.0.99)
C       NEAR EXTENSIONS OF SIDES (FLAG.GE.0.99)
C
C    DELETE ALL FOLLOWING LINES ABOVE LABEL 33 TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 33
        IF (MFLAG.EQ.0) GO TO 20
          FAC=C1/((C1*C1+S1*S1)*RR(N))
        IF(Z.LT.ZERO) FAC=-FAC   ! Fix sign of derivatives, Newman
          FDX=FDX+(RR(N)*V+FH1*U1)*FAC
          FDY=FDY-FE(N)*U1*FAC
          FAC=C2/((C2*C2+S2*S2)*RR(N1(N)))
        IF(Z.LT.ZERO) FAC=-FAC   ! Fix sign of derivatives, Newman
          FDX=FDX-(RR(N1(N))*V+FH2*U2)*FAC
          FDY=FDY+FE(N1(N))*U2*FAC
          GO TO 33
  20    IF (FLAG.LT.0.99D0) THEN
          U1V=U1*V
          FAC=Z/(C1*C1+S1*S1)
          FDX=FDX+(U1V*XRI(N)+R(N)*YMYV(N))*FAC
          FDY=FDY+(U1V*YRI(N)-R(N)*XMXV(N))*FAC
          U2V=U2*V
          FAC=Z/(C2*C2+S2*S2)
          FDX=FDX-(U2V*XRI(N1(N))+R(N1(N))*YMYV(N1(N)))*FAC
          FDY=FDY-(U2V*YRI(N1(N))-R(N1(N))*XMXV(N1(N)))*FAC
        ELSE
          ZS=Z*SIDE(N)
          USUM=U1+U2
          VRADS=V*RADS
          SFAC=VRADS*USUM
          SFS=-SFAC*ZS
          SFA=SFAC*C12
          CFAC=U2*R(N)+U1*R(N1(N))
          SFB=SFAC*CFAC
          CCF=C12*CFAC
          PA=(CCF+CCF)*VRADS-SFA*RSUM-SFB*ZZ*USUM
          PB=CCF*USUM*(VV+VV+RADS)-SFB*(S1+S1)*R(N1(N))
          PC=-SFA*U2-SFB*VV*R(N1(N))
          PD=-SFA*U1-SFB*VV*R(N)
          FAC=ZS/(CCF*CCF+SFS*SFS)
          FDX=FDX-(PA*CT+PB*ST+PC*XRI(N)+PD*XRI(N1(N)))*FAC
          FDY=FDY-(PA*ST-PB*CT+PC*YRI(N)+PD*YRI(N1(N)))*FAC
        ENDIF
  33    CONTINUE
      END DO
      IF (FD.LT.ZERO) FD=FD+TWOPI
      IF (Z.LT.ZERO) FD=-FD
      FS=FS-Z*FD

C----------------------------------------------------------------------
C    DELETE NEXT  3  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 99
      FSX=FSX-Z*FDX
      FSY=FSY-Z*FDY

99    CONTINUE

      fdd = fd
      fss = fs

      if(id .ne. 0) then
        fsxi = fsx
        fsyi = fsy
        fdxi = fdx
        fdyi = fdy
        fdzi = fdz
      endif

      RETURN
      END
