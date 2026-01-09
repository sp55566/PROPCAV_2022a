C ======================================================================
      SUBROUTINE NACA66(NH,THK,RLE,SC,YTC)
C     NSRDC modified NACA 66 with parabolic tail
C     NH = no. of SC points, YTC = half thickness at SC
C
C ======================================================================
      DIMENSION T66(19),PC(19),YTC(*),PSQ(19),CUBIC(72),SC(*)

      DATA PC/0.0,0.007596,0.030154,0.066987,0.116978,0.178606,
     A 0.25,0.328990,0.413176,0.50,0.586824,0.671010,0.75,0.821394,
     B 0.883022, 0.933013, 0.969846,0.992404,1.0/

      DATA T66/0.0000,0.1634,0.3216,0.4776,0.6270,0.7614,0.8726,
     A 0.9520, 0.9944, 0.9924,0.9404,0.8416,0.7057,0.5469,0.3825,
     B 0.2304,0.1072,0.0278,0.0000/

      DATA ZERO,HALF/0.0,0.5/

      PI = ACOS(-1.)
      RLE=0.448*THK**2
      RAD=PI/180.0

C.....Square root stretched coordinate for spline interpolation
      DO 10 N=1,19
         PSQ(N)=ACOS(1.0-2.0*PC(N))
10    CONTINUE

C.....Spline coefficients with leading edge radius incorporated
      A66LER=0.448
      SLE66=SQRT(2.0*A66LER)
      TLE66=ATAN(SLE66)/RAD
      CALL UGLYDK(19,2,1,PSQ,T66,TLE66,ZERO,CUBIC)

C.....Half thickness and slope at nodal points
      WHTC=HALF*THK
      DO 30 N=1,NH
         XSQ=ACOS(1.0-2.0*SC(N))
         CALL EVALDKs(19,1,PSQ,XSQ,YSPLN,CUBIC)
         YTC(N)=WHTC*YSPLN
30    CONTINUE
      RETURN
C))))))))))))))))))))) End of subroutine NACA66 ((((((((((((((((((((((((
      END
