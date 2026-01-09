C
C
C
C
      SUBROUTINE CBPB(NH,CRT,SC,YCC,TCC)
C***********************************************************************
C     PARABOLIC CAMBER FORM
C     NH = no. of sc points, YCC = camber at SC, TCC = slope angle at SC
C
C     Date      Revisions or Comments
C     --------  ----------------------
C     JY042400  Modified subroutine to avoid out-of-bound error for 
C               variable D2YDX.
C
C***********************************************************************
      DIMENSION YCC(*),TCC(*),SC(*),PER(17),CAM(17),
     *          CUBIC(64),D2YDX(1200)

      DATA PER/0.0,0.01,0.025,0.05,0.1,0.2,0.3,
     *         0.4,0.5,0.6,0.7,0.8,0.9,0.95,
     *         0.975,0.99,1.0/
      DATA ZERO/0.0/
C
C.....Find camber using parabolic function with max. value located at
C.....PER=0.5
      DO I=1,17
         CAM(I)=4.*(PER(I)-PER(I)**2.)
      END DO
C
C.....Find spline cubic coefficients of camber line
      CALL UGLYDK(17,1,1,PER,CAM,ZERO,ZERO,CUBIC)
C
C.....Evaluate camber height and slope at nodal points
      CALL EVALDK(17,NH,PER,SC,YCC,CUBIC)
      CALL DRIVDK(17,NH,PER,SC,TCC,D2YDX,CUBIC)
      DO 10 N=1,NH
         YCC(N)=CRT*YCC(N)
         TCC(N)=ATAN(CRT*TCC(N))
 10   CONTINUE

      RETURN
C)))))))))))))))))))))) End of subroutine CBPB (((((((((((((((((((((((((
      END
