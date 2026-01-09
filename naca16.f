C =================================================
      SUBROUTINE NACA16(NH,THK,RLE,SC,YTC)

C     NACA 16 
C     NH = no. of SC points, YTC = half thickness at SC
C
C =================================================

      DIMENSION T66(17),PC(17),YTC(*),PSQ(17),CUBIC(68),SC(*)

C----- N. FINE ---- 7/6/92 ----------------------------
C     NACA16 MODIFIED SECTION

      DATA PC/0.0,0.0125,0.025,0.05,0.075,0.1,0.15,0.2,0.3,0.4
     A ,0.5,0.6,0.7,0.8,0.9,0.95
     B ,1.0/
C
      DATA T66/0.0000,0.2153,0.301,.4183,.5053,.5763,.689,.7773
     A ,.903,.9757,1.0,.9723,.8783,.7,.4197,0.2357
     B ,0.0000/

      DATA ZERO,HALF/0.0,0.5/

      PI = ACOS(-1.)
      RLE=0.48889*THK**2
      RAD=PI/180.0
C
C.....Square root stretched coordinate for spline interpolation
      DO 10 N=1,17
         PSQ(N)=ACOS(1.-2.*PC(N))
10    CONTINUE
C
C.....Spline coefficients with leading edge radius incorporated

      A66LER=0.48889
      SLE66=0.5*SQRT(2.0*A66LER)
      TLE66=ATAN(SLE66)/RAD
      CALL UGLYDK(17,2,1,PSQ,T66,TLE66,ZERO,CUBIC)
C
C.....Half thickness and slope at nodal points
      HTC=HALF*THK
      DO 30 N=1,NH
         XSQ=ACOS(1.-2.*SC(N))
         CALL EVALDKs(17,1,PSQ,XSQ,YSPLN,CUBIC)
         YTC(N)=HTC*YSPLN
30    CONTINUE
      RETURN
C))))))))))))))))))))) End of subroutine NACA66 ((((((((((((((((((((((((
      END

C =================================================
      SUBROUTINE NAVSEA1(NH,THK,RLE,SC,YTC)

C     NACA 16
C     NH = no. of SC points, YTC = half thickness at SC
C
C =================================================

      DIMENSION T66(17),PC(17),YTC(*),PSQ(17),CUBIC(68),SC(*)

C----- N. FINE ---- 7/6/92 ----------------------------
C     NACA16 MODIFIED SECTION

      DATA PC/0.0,0.0125,0.025,0.05,0.075,0.1,0.15,0.2,0.3,0.4
     A ,0.5,0.6,0.7,0.8,0.9,0.95
     B ,1.0/
C
      DATA T66/0.0000,0.2153,0.301,.4183,.5053,.5763,.689,.7773
     A ,.903,.9757,1.0,.9600,.8400,.6400,.3600,0.1900
     B ,0.0000/

      DATA ZERO,HALF/0.0,0.5/

      PI = ACOS(-1.)
      RLE=0.48889*THK**2
      RAD=PI/180.0
C
C.....Square root stretched coordinate for spline interpolation
      DO 10 N=1,17
         PSQ(N)=ACOS(1.-2.*PC(N))
10    CONTINUE
C
C.....Spline coefficients with leading edge radius incorporated

      A66LER=0.48889
      SLE66=0.5*SQRT(2.0*A66LER)
c      SLE66=SQRT(2.0*A66LER)
      TLE66=ATAN(SLE66)/RAD
      CALL UGLYDK(17,2,1,PSQ,T66,TLE66,ZERO,CUBIC)
C
C.....Half thickness and slope at nodal points
      HTC=HALF*THK
      DO 30 N=1,NH
         XSQ=ACOS(1.-2.*SC(N))
         CALL EVALDKs(17,1,PSQ,XSQ,YSPLN,CUBIC)
         YTC(N)=HTC*YSPLN
30    CONTINUE
      RETURN
C))))))))))))))))))))) End of subroutine NACA66 ((((((((((((((((((((((((
      END


              


