C ======================================================================
      SUBROUTINE NACA65A(NH,THK,RLE,SC,YTC)
C
C     NACA 65A section from 006 to 100 ; theory of wing section pp 367
C     NH = no. of SC points, YTC = half thickness at SC
C
C ======================================================================
      DIMENSION TRAE(26),PC(26),YTC(*),PSQ(26),CUBIC(104),SC(*)
      DATA TRAE/0.0,.1547,.1877,.2393,.3270,
     *        .4377,.5303,.6080,.7313,.8247
     *       ,.8957,.9473,.9817,.9987,.9973
     *       ,.9750,.9310,.8673,.7880,.6957
     *       ,.5917,.4790,.3610,.2423,.1233, 0.0/
      DATA PC/0.0,0.005,.0075,.0125,.025,   .05,.075,.1,.15,.20,
     *        0.25,.30,.35,.4,.45,    .5,.55,.6,.65,.7,
     *        0.75,.8,.85,.9,.95,     1.00/
      DATA ZERO,HALF,TWO,C636/0.0,0.5,2.0,0.636/

      PI = ACOS(-1.)
      RAD=PI/180.0
C
C.....Assume RLE=0.636*THK**2
      RLE=0.636*THK**2
C
C.....Square root stretched coordinate for spline interpolation
      DO 10 N=1,26
         PSQ(N)=SQRT(PC(N))
10    CONTINUE
C
C.....Spline coefficients
      TRLE=ATAN(TWO*SQRT(TWO*C636))/RAD
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
C)))))))))))))))))))) End of subroutine NACA65A ((((((((((((((((((((((((
      END
