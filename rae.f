C =======================================================================
      SUBROUTINE RAE(NH,THK,SC,YTC)
C     RAE 010 foil sections from AGARD-AR-138 by Treadgold  et al.
C     NH = no. of SC points, YTC = half thickness at SC
C
C =======================================================================
      DIMENSION TRAE(18),PC(18),YTC(*),PSQ(18),CUBIC(68),SC(*)
      DATA TRAE/0.0,.2453,.3835,.5312,.7215,.9261,.9994,.9601,.9131
     *        ,.8534,.7062,.5361,.3577,.1789,.0894,.0447,.0179,0.0/
      DATA PC/0.0,0.01,0.025,0.05,0.1,0.2,0.3,
     *        0.4,0.45,0.5,0.6,0.7,0.8,0.9,0.95,
     *        0.975,0.99,1.00/
      DATA ZERO,HALF/0.0,0.5/

      PI = ACOS(-1.)
      RAD=PI/180.0
C
C.....Square root stretched coordinate for spline interpolation
      DO 10 N=1,18
         PSQ(N)=SQRT(PC(N))
10    CONTINUE
C
C.....Spline coefficients
      CALL UGLYDK(18,1,1,PSQ,TRAE,ZERO,ZERO,CUBIC)
C
C.....Half thickness at nodal point
      HTC=HALF*THK
      DO 20 N=1,NH
         XSQ=SQRT(SC(N))
         CALL EVALDKs(18,1,PSQ,XSQ,YSPLN,CUBIC)
         YTC(N)=HTC*YSPLN
20    CONTINUE
      RETURN
C))))))))))))))))))))))) End of subroutine RAE (((((((((((((((((((((((((
      END
