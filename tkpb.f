C =====================================================================
      SUBROUTINE TKPB(NH,THK,RLE,SC,YTC)
C
C     PARABOLIC THICKNESS FORM.
C     NH = no. of SC points, YTC = half thickness at SC
C
C =====================================================================

      DIMENSION TRAE(26),PC(26),YTC(*),PSQ(26),CUBIC(104),SC(*)
      DATA PC/0.0,0.005,.0075,.0125,.025,   .05,.075,.1,.15,.20,
     *        0.25,.30,.35,.4,.45,    .5,.55,.6,.65,.7,
     *        0.75,.8,.85,.9,.95,     1.00/
      DATA ZERO,HALF,TWO,C687/0.0,0.5,2.0,0.687/

      PI = ACOS(-1.)
      RAD=PI/180.0
C
C.....Calculate thickness using parabolic function where the max. value  
C.....is located at PC=0.5
      DO I=1,26
         TRAE(I)=4.*(PC(I)-PC(I)**2.)
      END DO
C
C.....Assume RLE=0.687*THK**2
      RLE=0.687*THK**2
C
C.....Square root stretched coordinate for spline interpolation
      DO 10 N=1,26
         PSQ(N)=SQRT(PC(N))
10    CONTINUE
C
C.....Spline coefficients
      TRLE=ATAN(TWO*SQRT(TWO*C687))/RAD
      CALL UGLYDK(26,2,1,PSQ,TRAE,TRLE,ZERO,CUBIC)
C
C.....Half thickness at nodal point
      HTC=HALF*THK
      DO 20 N=1,NH
         XSQ=SQRT(SC(N))
         CALL EVALDKs(26,1,PSQ,XSQ,YSPLN,CUBIC)
         YTC(N)=HTC*YSPLN
20    CONTINUE
      RETURN
C))))))))))))))))))))) End of subroutine TKPB (((((((((((((((((((((((
      END

