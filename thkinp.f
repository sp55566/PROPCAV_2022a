C =====================================================================
      SUBROUTINE THKINP(NH,TTI,SC,YTC)
C
C     NH = no. of SC points, YTC = half thickness at SC
C
C =====================================================================

      DIMENSION TTI(16),PC(17),YTC(*),PSQ(17),CUBIC(64),
     *          SC(*),TT(17)
      DATA PC/0.0,0.01,0.025,0.05,0.1,0.2,0.3,
     *        0.4,0.5,0.6,0.7,0.8,0.9,0.95,
     *        0.975,0.99,1.00/
      DATA ZERO,HALF/0.0,0.5/

      PI = ACOS(-1.)
      RAD=PI/180.0
C
C.....Square root stretched coordinate for spline interpolation
      DO 10 N=1,17
         PSQ(N)=SQRT(PC(N))
 10   CONTINUE
      DO 20 M=1,16
         TT(M+1)=TTI(M)
 20   CONTINUE
      TT(1)=0.0
      CALL UGLYDK(17,1,1,PSQ,TT,ZERO,ZERO,CUBIC)
C
C.....Half thickness and slope at nodal points
      DO 30 N=1,NH
         XSQ=SQRT(SC(N))
         CALL EVALDKs(17,1,PSQ,XSQ,YSPLN,CUBIC)
         YTC(N)=HALF*YSPLN
 30   CONTINUE
      RETURN
C))))))))))))))))))))) End of subroutine thkinp ((((((((((((((((((((((((
      END
