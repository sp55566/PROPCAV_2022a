C     
C     
      SUBROUTINE A8ML(NH,CRT,SC,YCC,TCC)
************************************************************************
C     NACA A=.8 mean line
C     NH = no. of sc points, YCC = camber at SC, TCC = slope angle at SC
C     
C     Date      Revisions or Comments
C     --------  ----------------------
C     JY042400  Modified subroutine to avoid out-of-bound error for 
C     variable D2YDX.
C     
************************************************************************
      DIMENSION YCC(*),TCC(*),SC(*),PER(17),CAM(17),
     *  CUBIC(64),D2YDX(1200)

      DATA PER/0.0,0.01,0.025,0.05,0.1,0.2,0.3,
     *  0.4,0.5,0.6,0.7,0.8,0.9,0.95,
     *  0.975,0.99,1.0/
      DATA CAM/0.0,0.0755,0.1586,0.2711,0.4482,0.6993,
     *  0.8635,0.9614,1.0,0.9785,0.8891,0.7027,
     *  0.3586,0.1713,0.0823,0.0307,0.0/
      DATA ZERO/0.0/
C     
c.....Find spline cubic coefficients of camber line
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
C))))))))))))))))))))))End of subroutine A8ML (((((((((((((((((((((((((
      END

      SUBROUTINE NACA65(NH,CRT,SC,YCC,TCC)
************************************************************************
C     NACA 65 mean line
C     NH = no. of sc points, YCC = camber at SC, TCC = slope angle at SC
C
C     Date               Revisions or Comments
C     --------           ----------------------
C     SNKIM Aug/23/2017  added NACA 65 section meanline to test the 4710
C                        propeller in DTSNRDC report.
C
************************************************************************
      DIMENSION YCC(*),TCC(*),SC(*),PER(17),CAM(17),
     *  CUBIC(64),D2YDX(1200)

      DATA PER/0.0,0.0125,0.0250,0.0500,0.0750,0.1000,0.1500,
     *  0.2000,0.3000,0.4000,0.5000,0.6000,0.7000,0.8000,
     *  0.9000,0.9500,1.0000/
      DATA CAM/0.0,0.0494,0.0975,0.1900,0.2775,0.3600,
     *  0.5100,0.6400,0.8400,0.9600,1.0000,0.9600,
     *  0.8400,0.6400,0.3600,0.1900,0.0/
      DATA ZERO/0.0/
C
c.....Find spline cubic coefficients of camber line
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
C))))))))))))))))))))))End of subroutine A8ML (((((((((((((((((((((((((
      END

