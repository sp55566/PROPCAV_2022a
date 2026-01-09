      SUBROUTINE CAMINP(NH,CAM1,SC,YCC,TCC)
C***********************************************************************
C     Copied from PSF10.3 on 11/27/98 by J.Young                       *
C
C     NH = no. of sc points, YCC = camber at SC, TCC = slope angle at SC
C
C     JY042400  Modified subroutine to avoid out-of-bound error for 
C               variable D2YDX.
C
C***********************************************************************
      DIMENSION YCC(*),TCC(*),SC(*),PER(17),CAM(17),CAM1(15),
     *          CUBIC(64),D2YDX(1000)
      DATA PER/0.0,0.01,0.025,0.05,0.1,0.2,0.3,
     *         0.4,0.5,0.6,0.7,0.8,0.9,0.95,
     *         0.975,0.99,1.0/
      DATA ZERO/0.0/
      DO 10 M=1,15
         CAM(M+1)=CAM1(M)
 10   CONTINUE
      CAM(1)=0.0
      CAM(17)=0.0
C
C.....Find spline cubic coefficients of camber line
      CALL UGLYDK(17,1,1,PER,CAM,ZERO,ZERO,CUBIC)

C.....Evaluate camber height and slope at nodal points
      CALL EVALDK(17,NH,PER,SC,YCC,CUBIC)
      CALL DRIVDK(17,NH,PER,SC,TCC,D2YDX,CUBIC)
      DO 20 N=1,NH
         TCC(N)=ATAN(TCC(N))
 20   CONTINUE

      RETURN
C)))))))))))))))))))))) End of subroutine caminp (((((((((((((((((((((((
      END
