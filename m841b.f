       SUBROUTINE M841B
************************************************************************
*     NX    = no. of radial input
*     NPTS  = no. of chord-wise input
*     DIA   = diameter of blade
*     XCHD  = chord length/diameter at each input radius
*     XSKEW = skew/diameter at each radius 
*     XRAKE = rake/diameter at each radius
*     XPI   = pitch in radiant at each radius
*     XC    = chord-wise locations (x/C) where the choordinates of the
*             pressure side and suctions sides are given
*     XR    = r/R at each radius
*     YPC   = Yp/C = coordinates of the pressure side at each radius
*     YSC   = Ys/C = coordinates of the suctions side at each radius
************************************************************************
       INCLUDE 'PARAM.INC'
       PARAMETER(NXM=10,NPTSM=16)
       COMMON/INT1/NPTS
       COMMON/REAL3/XC(NPTSM)
       COMMON/REAL2/YPC(NPTSM,NXM),YSC(NPTSM,NXM)
       COMMON/REAL1/DIA,XR1(NXM),C(NXM),SKEW1(NXM),RAKE1(NXM),RLE(NXM),
     *      PITCH1(NXM),YP(NPTSM*NXM),YS(NPTSM*NXM)
       COMMON/GINP/XR(NXMAX),XRSQ(NXMAX),XPI(NXMAX),XRAKE(NXMAX),
     *             XSKEW(NXMAX),XCHD(NXMAX),XCI(NXMAX),XTI(NXMAX),
     *             XVA(NXMAX),XVR(NXMAX),XVT(NXMAX),XUA(NXMAX),
     *             XUAU(NXMAX),XUT(NXMAX),XUTU(NXMAX)
       common/dat2/ nx, nblade

       CALL M841GEO

       PI=ATAN(1.)*4.
       RAD=DIA/2.

       DO I=1,NX
          DO N=1,NPTS
             IC=(I-1)*NPTS+N
             YPC(N,I)=YP(IC)/C(I)
             YSC(N,I)=YS(IC)/C(I)
          END DO
          XR(I)=XR1(I)
          XCHD(I)=C(I)/DIA
          XPI(I)=PITCH1(I)/180.*PI
          XSKEW(I)=SKEW1(I)/DIA
          XRAKE(I)=RAKE1(I)/DIA
       END DO  

       CALL M841BLADE

       RETURN
       END


       SUBROUTINE M841BLADE
************************************************************************
*      Generate blade geometry for SPP M841B.
************************************************************************
       use m_BLADEM
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       PARAMETER(NXM=10,NPTSM=16)
       COMMON/INT1/NPTS
       COMMON/REAL3/XC(NPTSM)
       COMMON/REAL2/YPC(NPTSM,NXM),YSC(NPTSM,NXM)
       DIMENSION XIM(NBHZP),ETAM(NBHZP),XIOLD(NBPZ,MBPZ),
     *      ETAOLD(NBPZ,MBPZ),XIMOLD(NBHZP,MBPZ),ETAMOLD(NBHZP,MBPZ)
       DIMENSION YP1(NPTSM,MBPZ),YS1(NPTSM,MBPZ),
     *      YTMP1(NXM),YTMP2(MBPZ),PSQ(NPTSM),XSQ(NBHZ),TMP1(NPTSM),
     *      YSPLN(NPHZ),CUBIC(4*NPTSM-4),XSQP(NBHZ),CUBTMP1(4*NXM-4)
!      COMMON/BLADEM/XBM(NBHZP,MBPZ),YBM(NBHZP,MBPZ),ZBM(NBHZP,MBPZ)
       COMMON /HUB5TMP/RHULT1
       CHARACTER*29 FNPLT  

       RHULT1=RHULT

C-----------------------------------------------------------------------
C      Check paramters for dimensioning errors
C-----------------------------------------------------------------------
       CALL CHECK_PARAM0
       CALL CHECK_PARAM

C-----------------------------------------------------------------------
C     Do not align the wake
C-----------------------------------------------------------------------
       DO N=1,NX
          XUA(N)=ZERO
          XUAU(N)=ZERO
          XUT(N)=ZERO
          XUTU(N)=ZERO
       END DO

C-----------------------------------------------------------------------
C     Read name of wake file (this part is moved from propcav.f), then
C     read and interpolate data from wake file by calling subroutine 
C     readwak.f.
C                                                               JY090799
C-----------------------------------------------------------------------
       WRITE(*,*) ' '
       WRITE(*,*) ' PROPCAV> ENTER INFLOW WAKE FILENAME: '
       READ(*,*) WKFILE
       WRITE(*,*) ' '
       
       CALL READWAK

C-----------------------------------------------------------------------
C     Calculate SPLINE cubic coefficients 
C-----------------------------------------------------------------------
       DO N=1,NX
          XRSQ(N)=1.0-SQRT(ABS(1.0-XR(N)))
       END DO
       CALL UGLYDK(NX,1,1,XR,XPI,0.0,0.0,PICUB)
       CALL UGLYDK(NX,1,1,XR,XRAKE,0.0,0.0,RKCUB)
       CALL UGLYDK(NX,1,1,XR,XSKEW,0.0,0.0,SKCUB)
       CALL UGLYDK(NX,1,1,XRSQ,XCHD,0.0,0.0,CHCUB)
     
       CALL UGLYDK(NX,1,1,XR,XVA,0.0,0.0,VACUB)
       CALL UGLYDK(NX,1,1,XR,XVR,0.0,0.0,VRCUB)
       CALL UGLYDK(NX,1,1,XR,XVT,0.0,0.0,VTCUB)
       CALL UGLYDK(NX,1,1,XR,XUA,0.0,0.0,UACUB)
       CALL UGLYDK(NX,1,1,XR,XUAU,0.0,0.0,UAUCUB)
       CALL UGLYDK(NX,1,1,XR,XUT,0.0,0.0,UTCUB)
       CALL UGLYDK(NX,1,1,XR,XUTU,0.0,0.0,UTUCUB)

C-----------------------------------------------------------------------
C     Parameters for generating the geometry
C-----------------------------------------------------------------------
       NH=NC/2
       NHP=NH+1
       NHM=NH-1
       NCP=NC+1
       MRP=MR+1
       NPANB=NC*MR
       DELK=TWOPI/NBLADE

C-----------------------------------------------------------------------
C     Create the blade geometry
C-----------------------------------------------------------------------
C.....Radial spacing: IRSPAC = 0:cosine; 1:constant; 2:half cosine
       IF(IRSPAC.EQ.0) THEN
          DTH=PI/MR
          DO 5 M=1,MRP
             RZ(M)=RHUB+0.5*(1.0-RHUB)*(1.0-COS(DTH*(M-1)))
 5        CONTINUE
       ELSE IF(IRSPAC.EQ.1) THEN 
          DELR=(1.0-RHUB)/MR
          DO 10 M=1,MRP
             RZ(M)=RHUB+DELR*(M-1)
 10       CONTINUE
       ELSE IF(IRSPAC.EQ.2) THEN
          DTH=HALFPI/MR
          DO 20 M=1,MRP
             RZ(M)=RHUB+(1.0-RHUB)*SIN(DTH*(M-1))
 20       CONTINUE
       END IF

       DO 30 M=1,MRP
          TERM1  =1.0 - RZ(M)
          IF(TERM1.LT.0.0)THEN
             WRITE(*,*) 'Notice: Radical adjusted in progeo'
             TERM1 = ABS(TERM1)
          ENDIF
          RZSQ(M)=1.0-SQRT(TERM1)
 30    CONTINUE

C.....Interpolate geometrical parameters in the radial direction
       CALL EVALDK(NX,MRP,XR,RZ,PITCH,PICUB)
       CALL EVALDK(NX,MRP,XR,RZ,RAKE,RKCUB)
       CALL EVALDK(NX,MRP,XR,RZ,SKEW,SKCUB)
       CALL EVALDK(NX,MRP,XRSQ,RZSQ,CHORD,CHCUB)

       DO N=1,NPTS

C........Pressure side coordinates
          DO I=1,NX
             YTMP1(I)=YPC(N,I)
          END DO
          CALL UGLYDK(NX,1,1,XR,YTMP1,0.0,0.0,CUBTMP1)
          CALL EVALDK(NX,MRP,XR,RZ,YTMP2,CUBTMP1)

          DO MM=1,MRP
             YP1(N,MM)=YTMP2(MM)
          END DO

C........Suction side coordinates
          DO I=1,NX
             YTMP1(I)=YSC(N,I)
          END DO
          CALL UGLYDK(NX,1,1,XR,YTMP1,0.0,0.0,CUBTMP1)
          CALL EVALDK(NX,MRP,XR,RZ,YTMP2,CUBTMP1)
          DO MM=1,MRP
             YS1(N,MM)=YTMP2(MM)
          END DO

       END DO
       
C-----------------------------------------------------------------------
C     Set up chordwise spacing
C-----------------------------------------------------------------------
       CALL SPACE(NC,ICSPAC,RLET,SB)

C.....Square root stretched coordinate for spline interpolation
       DO N=1,NPTS
          PSQ(N)=SQRT(XC(N))
       END DO

       DO N=1,NH
          XSQ(N)=SQRT(SB(N))
       END DO

       DO M=1,MRP
          TANP=TAN(PITCH(M))

C........Find chord length w.r.t. pressure side & camber
          NT1=10
          DO N=1,NT1
             N1=NPTS-NT1+N
             TMP1(N)=XC(N1)
             YSPLN(N)=YP1(N1,M)
             IF(YP1(N1,M).GT.0) NMARK=N1
          END DO
          CALL QUADFIT(NT1,TMP1,YSPLN,A0,A1,A2)
          DUM=A1*A1-4.*A2*A0
          IF(DUM.LT.0) THEN
             DUM1=YS1(NPTS,M)/TANP
 1           XCP=1.-DUM1
             XCM=1.-DUM1/2.
          ELSE
             X1=(-A1+SQRT(DUM))/(2.*A2)
             X2=(-A1-SQRT(DUM))/(2.*A2)
             IF(X1.GE.XC(NMARK).AND.X1.LE.XC(NMARK+1)) THEN
                XCP=X1
             ELSE
                XCP=X2
             END IF
             XCM=(XCP+1.)/2.
          END IF

C........Face side local coordinate
          DO N=1,NH
             N1=NHP-N
             XI(N1)=SB(N)*XCP
             XSQP(N)=SQRT(XI(N1))
          END DO

          DO N=1,NPTS
             TMP1(N)=YP1(N,M)
          END DO
          CALL UGLYDK(NPTS,1,1,PSQ,TMP1,ZERO,ZERO,CUBIC)
          CALL EVALDK(NPTS,NH,PSQ,XSQP,YSPLN,CUBIC)

          DO N=1,NH
             N1=NHP-N
             IF(N.LT.NH) THEN
                ETA(N1)=YSPLN(N)
             ELSE
                ETA(N1)=ZERO
             END IF
          END DO

C........Back side local coordinate
          DO N=1,NPTS
             TMP1(N)=YS1(N,M)
          END DO
          CALL UGLYDK(NPTS,1,1,PSQ,TMP1,ZERO,ZERO,CUBIC)
          CALL EVALDK(NPTS,NH,PSQ,XSQ,YSPLN,CUBIC)

          DO N=1,NH
             N1=NHP+N
             XI(N1)=SB(N)
             ETA(N1)=YSPLN(N)
          END DO

C........LE coordinate
          XI(NHP)=ZERO
          ETA(NHP)=ZERO

C........Camber coordinate
          XIM(1)=ZERO
          ETAM(1)=ZERO

          DO N=1,NH
             XIM(N+1)=SB(N)*XCM
             XSQP(N)=SQRT(XIM(N+1))
          END DO

          DO N=1,NPTS
             TMP1(N)=(YP1(N,M)+YS1(N,M))/2.
          END DO
          CALL UGLYDK(NPTS,1,1,PSQ,TMP1,ZERO,ZERO,CUBIC)
          CALL EVALDK(NPTS,NH,PSQ,XSQP,YSPLN,CUBIC)

          DO N=1,NH
             ETAM(N+1)=YSPLN(N)
          END DO

C........Set the origin of local coordinate at center of chord, and
C........Non-dimensionalize w.r.t. propeller radius
          DO N=1,NCP
             XI(N)=(XI(N)-0.5)*CHORD(M)*TWO
             ETA(N)=ETA(N)*CHORD(M)*TWO
             XIOLD(N,M)=XI(N)
             ETAOLD(N,M)=ETA(N)
             IF(N.LE.NHP) THEN
                XIM(N)=(XIM(N)-0.5)*CHORD(M)*TWO
                ETAM(N)=ETAM(N)*CHORD(M)*TWO
                XIMOLD(N,M)=XIM(N)
                ETAMOLD(N,M)=ETAM(N)
             END IF
          END DO

C........Non-dimensionalize RAKE (RAKE/D) & SKEW (DEG) & PITCH (P/D)
          COSP=COS(PITCH(M))
          SINP=SIN(PITCH(M))
          RAKE(M)=RAKE(M)+SKEW(M)*SINP
          SKEW(M)=SKEW(M)*180.*2.*COSP/RZ(M)/PI
          PITCH(M)=PI*RZ(M)*TANP

C........Generate key blade geometry
          DO N=1,NCP
             DX= XI(N)*SINP-ETA(N)*COSP
             XB(N,M)=RAKE(M)*TWO +DX
             THETA=SKEW(M)*RAD+(XI(N)*COSP+ETA(N)*SINP)/RZ(M)
             IF(M.EQ.1) THEN
                THR(N)=THETA
             END IF
             IF(N.EQ.1) THEN
                THT(M)=THETA
             END IF
             YB(N,M)=RZ(M)*COS(THETA)
             ZB(N,M)=RZ(M)*SIN(THETA)
          END DO

          DO N=1,NHP
             DX= XIM(N)*SINP-ETAM(N)*COSP
             XBM(N,M)=RAKE(M)*TWO +DX
             THETA=SKEW(M)*RAD+ (XIM(N)*COSP+ETAM(N)*SINP)/RZ(M)
             YBM(N,M)=RZ(M)*COS(THETA)
             ZBM(N,M)=RZ(M)*SIN(THETA)
          END DO

C.......Plot blade sectsions (JY071201)
         IF(M.EQ.1) THEN
            OPEN(55,FILE='bldsec.plt',STATUS='UNKNOWN')
            WRITE(55,*) 'VARIABLES="c/R","y/R"'
         END IF
         WRITE(55,*) 'ZONE T="M=',M,'"'
         DO N=1,NCP
            WRITE(55,*) XI(N)+CHORD(M),ETA(N)+RZ(M)
         END DO
         WRITE(55,*) XI(1)+CHORD(M),ETA(1)+RZ(M)
         IF(M.EQ.MRP) CLOSE(55)
       END DO

C-----------------------------------------------------------------------
C     Create geometry for the separated region.                 JY081601
C-----------------------------------------------------------------------
       IF(ISC.EQ.1) CALL SRGEO(XIOLD,ETAOLD,XIMOLD,ETAMOLD)

C-----------------------------------------------------------------------
C     Compute the normal vectors of the camber surface
C-----------------------------------------------------------------------
       DO 90 N=1,NH
          DO 80 M=1,MR
             XG(1,1,1)=XBM(N,M)
             XG(1,1,2)=YBM(N,M)
             XG(1,1,3)=ZBM(N,M)
             XG(1,2,1)=XBM(N,M+1)
             XG(1,2,2)=YBM(N,M+1)
             XG(1,2,3)=ZBM(N,M+1)
             XG(1,3,1)=XBM(N+1,M+1)
             XG(1,3,2)=YBM(N+1,M+1)
             XG(1,3,3)=ZBM(N+1,M+1)
             XG(1,4,1)=XBM(N+1,M)
             XG(1,4,2)=YBM(N+1,M)
             XG(1,4,3)=ZBM(N+1,M)
             
             CALL GEOM3D(1,XG,CHRLEPS,IER)
             IF(IER.EQ.0) THEN
                WRITE(*,'(A)') ' UNACCEPTABLE PANELS ON CAMBER SURFACE'
                STOP
             END IF
             XCON(N,M)=VEL(1,1)
             YCON(N,M)=VEL(1,2)
             ZCON(N,M)=VEL(1,3)

             IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XCON(N,M),YCON(N,M))

 80       CONTINUE
 90    CONTINUE

C-----------------------------------------------------------------------
C     Prepare parameters for hub and wake geometries
C-----------------------------------------------------------------------
       XHBLE=-XB(NC/2+1,1)
       XHBTE=XB(1,1)
       XHBNU=XHBU+XHBLE
       XHBND=XHBD+XHBTE
       XHBFD=XHBT+XHBND

C-----------------------------------------------------------------------
C     Generate wake geometry
C-----------------------------------------------------------------------
       CALL GWAKE
!      CALL GWAKES
       CALL GWAKES(5, 4)

       NPANW=0
       DO M=1,MR
          NPANW=NPANW+MIN(NSW(M),NSW(M+1))
          IF(NSW(M).EQ.NSW(M+1)) THEN
             NPANW=NPANW-1
          END IF 
       END DO

C-----------------------------------------------------------------------
C     Create the hub geometry
C
C     Added new option to IHUB (JY052401):
C        IHUB=0: hub off
C            =1: hub on (open ended far upstream, closed far downstream)
C            =2: hub on (closed far upstream, closed far downstream)
C            =3: hub on (open ended far upstream, open far downstream)
C            =4: hub on (closed far upstream, open far downstream)
C            =5: hub on (open ended far upstream, open far downstream)
C-----------------------------------------------------------------------
       IF(IHUB.EQ.0)THEN
          NHBX=0
       ELSE
          IF(IHUB.EQ.1) THEN
             CALL GHUB1
          ELSE IF(IHUB.EQ.2) THEN
             CALL GHUB2
          ELSE IF(IHUB.EQ.3) THEN
             CALL GHUB3
          ELSE IF(IHUB.EQ.4) THEN
             CALL GHUB4
          ELSE IF(IHUB.EQ.5) THEN
             RHULT=RHULT1
             NWPANEL=NSW(1)-1
             CALL GHUB5
          ELSE IF(IHUB.EQ.7) THEN
             CALL GHUB7
          END IF

C.......Next IF statement added to prevent out-of-bound errors(JY072400)
          IF(NHBX.GT.NHPZ) THEN
             WRITE(*,'(''!!! WARNING !!! NOT ENOUGH PANELS IN HUB'')')
             WRITE(*,'('' NHBX = '',I5,'' NHBZ = '',I5)') NHBX,NHBZ
             WRITE(*,'('' Increase NHBZ --> NHBX in PARAM.INC! '')')
             STOP
          ELSE IF(MHBT.GT.MHBZ) THEN
             WRITE(*,'(''!!! WARNING !!! NOT ENOUGH PANELS IN HUB'')')
             WRITE(*,'('' MHBT = '',I5,'' MHBZ = '',I5)') MHBT,MHBZ
             WRITE(*,'('' Increase MHBZ --> MHBT in PARAM.INC! '')')
             STOP
          END IF

       ENDIF

       NPANH=NHBX*MHBT
       NPANEL=NPANB+NPANH

C-----------------------------------------------------------------------
C      Plot the geometries of propeller blade (GBLADE), hub (GHUB)
C       and wake (GWAKE)
C-----------------------------------------------------------------------

       CALL CHRLEN(FN,LENCH)

       FNPLT=FN(1:LENCH)//'.plt'
       OPEN(700, FILE=FNPLT, STATUS='UNKNOWN')
       WRITE(700,5112)
 5112  FORMAT(1X,'VARIABLES = "x", "y", "z"')

       IF(ISP .NE. 0) THEN
          CALL TRANS_GEO1(0)
          CALL TRANS_GEO2(0)
       ENDIF

       CALL PROPLT
       CALL PROPLT_SR

       CALL WAKEPLT

C-----------------------------------------------------------------------
C     Calculate abscissae for plotting 
C     (Approximate positions of control points spanwise and chordwise)
C-----------------------------------------------------------------------
C
C.....Chordwise abscissa on the blade 
C.....  SBP:  from l.e. (0.0) to t.e. (1.0)
       SBP(1)=HALF*SB(1)
       DO I=2,NH
          SBP(I)=HALF*(SB(I)+SB(I-1))
       END DO
       
C.....Spanwise abscissa on the blade 
       DO M=1,MR
          RZP(M)=HALF*(RZ(M)+RZ(M+1))
       END DO

C-----------------------------------------------------------------------
C     Compute number of radii smaller than 0.95R (MRTIP)
C-----------------------------------------------------------------------
       DO M=1,MR
          RZPSQ(M)=1.0-SQRT(ABS(1.0-RZP(M)))
          IF(RZP(M).LT.rmrtip) THEN
             MRTIP=M
          END IF
       END DO
       
C-----------------------------------------------------------------------
C     No smoothing for hydrofoil case.                          JY112998
C-----------------------------------------------------------------------
       IF((ICON.EQ.5).OR.(ICON.EQ.6) .OR.(ICON.EQ.8)) MRTIP=MR

       RETURN
       END

       SUBROUTINE M841GEO
************************************************************************
*      This subroutine contains the data printed in Table 4.1 & 4.2    *
*      of (Olofsson, 96) for Propeller 841-B.                          *      
************************************************************************

       PARAMETER(NXM=10,NPTSM=16)
       COMMON/INT1/NPTS
       COMMON/REAL3/XC(NPTSM)
       COMMON/REAL1/DIA,XR1(NXM),C(NXM),SKEW1(NXM),RAKE1(NXM),RLE(NXM),
     *      PITCH1(NXM),YP(NPTSM*NXM),YS(NPTSM*NXM)
       common /dat2/ nx, nblade

       DATA NX,NPTS/10,16/
       DATA DIA/250./
       DATA XR1/0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.975,1.0/
       DATA C/75.1,87.2,98.0,105.5,105.5,97.3,75.8,53.7,40.5,26.7/
       DATA SKEW1/-9.94,-11.19,-11.13,-10.08,-5.71,1.85,14.38,25.24,
     *      31.95,39.05/
       DATA RAKE1/6.75,8.45,10.12,11.82,13.52,15.20,16.90,17.72,
     *      18.15,18.58/
       DATA RLE/0.117,0.113,0.095,0.060,0.034,0.021,0.012,0.008,
     *      0.007,0.005/
       DATA PITCH1/52.8,44.6,38.3,33.3,29.4,26.3,23.7,22.6,22.0,21.5/
       DATA XC/.0125,.025,.05,.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,
     *      .975,1.0/
       DATA YP/-0.24,-0.23,-0.18,-0.04, 0.30, 0.70, 1.13, 1.56,
     *          1.98, 2.41, 2.32, 1.73, 0.84,-0.29,-0.89,-1.48
     *      ,  -0.26,-0.25,-0.20,-0.04, 0.34, 0.79, 1.26, 1.73, 
     *          2.20, 2.67, 2.39, 1.64, 0.57,-0.72,-1.38,-2.03
     *      ,  -0.27,-0.27,-0.21,-0.05, 0.36, 0.85, 1.35, 1.86,
     *          2.36, 2.85, 2.40, 1.52, 0.31,-1.08,-1.78,-2.48
     *      ,  -0.28,-0.27,-0.21,-0.04, 0.38, 0.87, 1.38, 1.90, 
     *          2.41, 2.91, 2.36, 1.43, 0.16,-1.27,-1.98,-2.70
     *      ,  -0.26,-0.26,-0.20,-0.04, 0.35, 0.82, 1.31, 1.79, 
     *          2.28, 2.75, 2.27, 1.40, 0.22,-1.13,-1.80,-2.48
     *      ,  -0.23,-0.22,-0.17,-0.03, 0.30, 0.71, 1.13, 1.55, 
     *          1.97, 2.39, 2.05, 1.34, 0.34,-0.82,-1.41,-2.00
     *      ,  -0.17,-0.17,-0.13,-0.03, 0.22, 0.51, 0.82, 1.13, 
     *          1.44, 1.74, 1.58, 1.10, 0.41,-0.43,-0.86,-1.29
     *      ,  -0.12,-0.11,-0.09,-0.02, 0.15, 0.35, 0.56, 0.77, 
     *          0.98, 1.20, 1.10, 0.79, 0.22,-0.25,-0.54,-0.84
     *      ,  -0.09,-0.08,-0.07,-0.02, 0.11, 0.26, 0.42, 0.57, 
     *          0.73, 0.89, 0.83, 0.60, 0.26,-0.16,-0.38,-0.60
     *      ,  -0.06,-0.06,-0.04,-0.01, 0.07, 0.17, 0.27, 0.37, 
     *          0.47, 0.57, 0.56, 0.41, 0.19,-0.08,-0.22,-0.36/

       DATA YS/ 0.62, 0.95, 1.50, 2.48, 4.23, 5.90, 7.59, 9.28,
     *          9.90, 9.40, 7.60, 7.00, 7.00, 7.00, 6.60, 5.59
     *      ,   0.57, 0.88, 1.40, 2.31, 3.94, 5.50, 7.08, 8.66,
     *          9.50, 9.40, 7.80, 7.30, 7.50, 7.50, 7.20, 6.70
     *      ,   0.51, 0.78, 1.24, 2.05, 3.50, 4.90, 6.31, 7.71, 
     *          8.65, 8.80, 7.50, 7.20, 7.50, 7.65, 7.40, 7.00
     *      ,   0.43, 0.66, 1.05, 1.74, 2.97, 4.15, 5.35, 6.55, 
     *          7.40, 7.60, 6.40, 6.00, 6.50, 7.00, 6.80, 6.60
     *      ,   0.34, 0.53, 0.83, 1.38, 2.35, 3.29, 4.23, 5.18, 
     *          6.00, 6.30, 5.40, 4.80, 5.20, 5.80, 5.65, 5.40
     *      ,   0.26, 0.40, 0.62, 1.04, 1.77, 2.47, 3.18, 3.88, 
     *          4.70, 5.10, 4.70, 4.20, 4.40, 4.70, 4.50, 4.20
     *      ,   0.17, 0.27, 0.42, 0.70, 1.19, 1.66, 2.14, 2.61, 
     *          3.10, 3.30, 2.95, 2.60, 2.70, 2.90, 2.70, 2.50
     *      ,   0.12, 0.18, 0.29, 0.48, 0.81, 1.13, 1.46, 1.78, 
     *          2.05, 2.20, 2.05, 1.80, 1.95, 1.95, 1.80, 1.70
     *      ,   0.09, 0.14, 0.22, 0.36, 0.61, 0.85, 1.09, 1.33, 
     *          1.55, 1.60, 1.50, 1.40, 1.45, 1.35, 1.30, 1.20
     *      ,   0.06, 0.09, 0.14, 0.23, 0.40, 0.56, 0.72, 0.88, 
     *          1.05, 1.20, 1.15, 1.10, 1.00, 0.90, 0.80, 0.75/

       RETURN
       END





  
