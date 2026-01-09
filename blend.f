C=======================================================================
C------------------------------------------------------------------
C     
      SUBROUTINE BLEND(N,XL,XR,DXL,DXR,XGRID)
C     
C     CREATES A BLENDED SPACING BETWEEN POINTS XL AND XR
C     
C     N    = NUMBER OF INTERVALS BETWEEN XL AND XR
C     DXL  = SPACING AT XL
C     DXR  = SPACING AT XR
C     XGRID(I) = GRID COORDINATES ( I=1, N+1 ) 
C     XGRID(1) = XL , XGRID(N+1) = XR
C     
C     SPYROS A. KINNAS  -- JAN. 26 1989
C     
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XGRID(*)
C     
      A = XL
      B = XR
      DX1 = DXL
      DXN = DXR
C     
      ALPHA=DX1/(B-A)
      BETA=DXN/(B-A)
      AA=ALPHA**2+ALPHA*BETA+BETA**2
      BB=3*(2*ALPHA+2*BETA+3*ALPHA**2+2*BETA**2+5*ALPHA*BETA)
      CC=3*(6*ALPHA+4*BETA+3)
      DETN=BB**2-4*AA*CC
C     
C-----DETERMINE INITIAL GUESSES FOR NCRIT 
C     FOR N > NCRIT SPACING BECOMES NEGATIVE 
C     
      NCRIT0=(BB+SQRT(DETN))/2./AA
      DET2=F(FLOAT(NCRIT0),A,B,DX1,DXN)
      SOL1=NCRIT0
      SOL2=NCRIT0+10
C     
C-----DETERMINE NCRIT BY USING SECANT METHOD
C     
 20   DET1=DET2
      DET2=F(SOL2,A,B,DX1,DXN)
      SOL=SOL2-DET2*(SOL1-SOL2)/(DET1-DET2)
      IF(ABS(SOL-SOL2).GT.1) THEN
       SOL1=SOL2
       SOL2=SOL
       GO TO 20
      ENDIF
C     
C     NCRIT=DMAX1(SOL,SOL2)
      NCRIT=AMAX1(SOL,SOL2)
C     
C-----DETERMINE SPACING
C     
      DELX=(B-A)/N
      CA=DX1/DELX
      CB=DXN/DELX
      CK=6./(N-2)*(1.-(CA+CB)/2.)
C     
      XGRID(1)=A
      DO 10 I=1,N
       FI=( CA*(N-I)+ CB*(I-1)+ (N-I)*(I-1)*CK )/(N-1)
       IF(FI.LE.0.) THEN
        WRITE(*,*) '  **** ERROR IN SUBR. BLEND !!!  SPACING < 0 **** '
        WRITE(*,*) '  ****       THE PROGRAM STOPS HERE !!       **** ' 
        WRITE(*,*) '  XL,  XR    : ',XL,XR
        WRITE(*,*) '  DXL, DXR   : ',DXL,DXR
        WRITE(*,*) '       N     : ',N
        WRITE(*,*) '     TRY AGAIN WITH ANOTHER N < NCRIT =', NCRIT
        STOP
       ENDIF
       DELXI=DELX*FI
       XGRID(I+1)=XGRID(I)+DELXI
 10   CONTINUE
C     
      RETURN
      END
C=======================================================================
C--------------------------------------------------------------
      FUNCTION F(X,A,B,DX1,DXN)
      CA=DX1/(B-A)*X
      CB=DXN/(B-A)*X
      CK=6./(X-2)*(1-(CA+CB)/2)
      F=(CB-CA+CK*(X+1))**2+4.*CK*((CA-CK)*X-CB)
      RETURN
      END      
C------------------------------------------------------------------
C     
      SUBROUTINE BLEND2(N,XL,XR,DXL,DXR,XGRID)
C     
C     CREATES A BLENDED SPACING BETWEEN POINTS XL AND XR
C     
C     PARABOLIC VARIATION OF DELTA X WITH I
C     
C     N    = NUMBER OF INTERVALS BETWEEN XL AND XR
C     DXL  = SPACING AT XL
C     DXR  = SPACING AT XR
C     XGRID(I) = GRID COORDINATES ( I=1, N+1 ) 
C     XGRID(1) = XL , XGRID(N+1) = XR
C     
C     SPYROS A. KINNAS  -- FEB. 1  1989
C     
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XGRID(*)
C     
      A = XL
      B = XR
      DX1 = DXL
      DXN = DXR
C     
C-----DETERMINE SPACING
C     
      DELX=(B-A)/FLOAT(N)
      CA=DX1/DELX
      CB=DXN/DELX
      AN=FLOAT(N)
      CK=30./(AN**3-4*AN**2+6*AN-4)*( AN-1-(2*AN-1)/6.*(CA+CB) )
C     
      XGRID(1)=A
      DO 10 I=1,N
       FI=( CA*(N-I)**2+ CB*(I-1)**2+ ((N-I)*(I-1))**2*CK )/(N-1)**2
       IF(FI.LE.0.) THEN
        WRITE(*,*) '  **** ERROR IN SUBR. BLEND !!!  SPACING < 0 **** '
        WRITE(*,*) '  ****       THE PROGRAM STOPS HERE !!       **** ' 
        WRITE(*,*) '  XL,  XR    : ',XL,XR
        WRITE(*,*) '  DXL, DXR   : ',DXL,DXR
        WRITE(*,*) '       N     : ',N
        WRITE(*,*) '     TRY AGAIN WITH ANOTHER N < NCRIT =  ***** '
        STOP
       ENDIF
       DELXI=DELX*FI
       XGRID(I+1)=XGRID(I)+DELXI
 10   CONTINUE
C     
      RETURN
      END
C=======================================================================



