C***********************************************************************
C     UTILIB: UTIlity LIBrary of PUF10.02
C           
C     This library includes the following subroutines:
C
C     --- Geometry:  Spacing:
C                        SPACE,ACOSH,ASINH
C                    Thickness forms:
C                        NACA66,NACA65A,NACA64A,NACA00,ELIPSE,RAE,RADLE
C                    Mean line:
C                        A8ML
C
C     --- Spline:    UGLYDK,EVALDK,DRIDK
C
C     --- Time:      TIME,WHEN (especially for MicroVax)
C
C
C
      SUBROUTINE SPACE(NT,ISPACE,RLE,SB)
C***********************************************************************
C     Chordwise spacing for 2D section
C
      DIMENSION SB(*),SC(1201)
      DATA PI/3.14159265/
      DATA ONE,HALF,TWO/1.0,0.5,2.0/
      NH=NT/2
      NHP=NH+1
      DELP=PI/NH
C-----------------------------------------------------------------------
C     ISPACE = 0: constant spacing
C              1: cosine spacing
C              2: COS*COSH spacing
C              3: SINH(COS) spacing
C              4: half cosine spacing
C-----------------------------------------------------------------------
      IF(ISPACE.EQ.0) THEN
         DS=1.0/FLOAT(NH)
         DO 10 N=1,NHP
            SC(N)=FLOAT(N-1)*DS
10       CONTINUE
      ELSE IF(ISPACE.EQ.4) THEN
         DTH=HALF*PI/NH
         DO 15 N=1,NHP
            SC(N)=ONE-COS( DTH*(N-1) )
15       CONTINUE
C-----------------------------------------------------------------------
C     ISPACE=5: blended spacing                           JY&HSLEE082800
C
C     smallest panel size at LE = 10% of the smallest size of cosine
C                                 spacing.
C     smallest panel size at TE = smallest size of cosine spacing.
C-----------------------------------------------------------------------
      ELSE IF(ISPACE.EQ.5) THEN
         DS=HALF*(ONE-COS(DELP))
         XL=0.0
         XR=1.0
         DXL=DS*.1
         DXR=DS
         CALL BLEND2(NH,XL,XR,DXL,DXR,SC)
C-----------------------------------------------------------------------
      ELSE
         DO 20 N=1,NHP
            SC(N)=HALF*(ONE-COS((N-1)*DELP))
20       CONTINUE
         IF(ISPACE.EQ.2) THEN
C-----------------------------------------------------------------------
C           COS*COSH spacing
C-----------------------------------------------------------------------
            ARG=HALF*PI*PI/(RLE*NT*NT)
            IF(ARG .LE. 1.04975) THEN
               CK=0.1
            ELSE
               CK=ACOSH(ARG)/PI
            END IF
            DO 30 N=1,NHP
               SC(N)=SC(N)*COSH(CK*(N-1)*DELP)/COSH(CK*PI)
30          CONTINUE
         ELSE IF(ISPACE.EQ.3) THEN
C-----------------------------------------------------------------------
C           SINH( COS ) spacing
C-----------------------------------------------------------------------
40          CONTINUE
            WRITE(*,'(A,$)') 
     *      ' PUF10>Enter fraction of first panel to l.e. radius:'
            READ(*,*) FRACT
            IMAX=100
            EPS=0.001
            CK1=HALF
            DO 50 KK=1,IMAX
               CK=0.5*ASINH(CK1*DELP*DELP/TWO/FRACT/RLE)
               ERR=ABS(CK-CK1)
               IF(ERR.LE.EPS) GO TO 60
               CK1=CK
50          CONTINUE
60          DEN=SINH(2.0*CK)
            DO 70 N=1,NHP
               THETA=DELP*(N-1)
               SC(N)=SINH( CK*( ONE-COS(THETA) ) ) / DEN
70          CONTINUE
            WRITE(*,'('' PUF10>  SC:'',/,(10(2X,F6.4)) )')
     *           (SC(N),N=1,NHP)
            WRITE(*,'(A,$)') 
     *      ' PUF10>Do you satisfied? (0 -- yes, 1 -- no):'
            READ(*,*) NOK
            IF(NOK.EQ.1) THEN
               GO TO  40
            END IF
         END IF
      END IF
      DO 80 I=1,NH
         SB(I)=SC(I+1)
80    CONTINUE
      RETURN
C))))))))))))))))))))) End of subroutine SPACE (((((((((((((((((((((((((
      END
