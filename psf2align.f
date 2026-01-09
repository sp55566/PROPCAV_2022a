C************************ MIT PSF10 ************************************
C
C                PROPELLER STEADY FLOW ANALYSIS PROGRAM
C                        PANEL METHOD SOLUTION         
C
C                             Version 1.0
C
C       Copyright (c) Massachusetts Institute of Technology 1997
C
C                     Release Date: 1 January 1997
C     
C  ---- 113000 HSLEE -------------------------------------------
C    This version of PSF2ALIGN has ONE More Subroutine 
C                           at the END of this program 
C         SUBROUTINE LSQ(X,Y,N,M) : Least Square Curve Fitting
C             X(N) : X points
C             Y(N) : Y Points
C             N    : Number of Data points
C             M    : Degree of polynimial to  be fitted.
C  ---- 113000 HSLEE -------------------------------------------
C
C***********************************************************************
      SUBROUTINE PSF2ALIGN(KTHK,KCAM,DAMP,IHUB,IDUCT)
C***********************************************************************
C     PSF2ALIGN: PSF-2 type wake ALIGNment
C      --- Find induced velocities in the transition wake and 
C          utimate wake region
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      COMMON/HIMAGE/JHUB,JDUCT
      COMMON /DIMAGE/ DGAPSIZE,ENDFOIL

C -- 113000 HSLEE for the Least Square Curve Fitting
      Dimension xuan2(nxmax),xuau2(nxmax),xutn2(nxmax),xutu2(nxmax)
C -- mlsqp : The order of polynomial to be fitted. (It fixed to be 3)
      common /lsqorder/ mlsqp
      DIMENSION XRTMP(NXMAX),XPITMP(NXMAX),XRAKETMP(NXMAX),
     %          XSKEWTMP(NXMAX),XCHDTMP(NXMAX),XCITMP(NXMAX),
     %          XTITMP(NXMAX), XCSTMP(NXMAX,15),XTSTMP(NXMAX,16)
C ---- 113000 --------------------------------------

C-----------------------------------------------------------------------
C     Parameters of wake alignment procedure
C-----------------------------------------------------------------------

C --- 113000 HSLEE ----------
      MLSQP = 3
C ---------------------------

C      JHUB = 0
C      ICHECK=0
      JHUB = 1
      JDUCT = 1
      IF(IHUB .EQ. 0) JHUB = 0
      IF(IDUCT .EQ. 0) JDUCT = 0

      IRUN=0
      IPASS=0
      DO 1 N=1,NX
         XCAM(N)=XCI(N)
         XTHICK(N)=XTI(N)
         XUAN(N)=XUA(N)
         XUTN(N)=XUT(N)
 1    CONTINUE

      DO N = 1, NX
         XRTMP(N) = XR(N)
         XPITMP(N) = XPI(N)
         XRAKETMP(N) = XRAKE(N)
         XSKEWTMP(N) = XSKEW(N)
         XCHDTMP(N) = XCHD(N)
         XCITMP(N) = XCI(N)
         XTITMP(N) = XTI(N)
         DO K = 1, 15
            XCSTMP(N,K) = XCS(N,K)
         ENDDO
         DO K = 1, 16
            XTSTMP(N,K) = XTS(N,K)
         ENDDO
      ENDDO

C HSLEE (041300)
      rultsave = rult
      rhultsave = rhult
C HSLEE (041300)

C-----------------------------------------------------------------------
C     Wake alignment procedure
C-----------------------------------------------------------------------
      CALL BLINPT(KTHK,KCAM)
      CALL ADINPT
      CALL SETUP 
      CALL BLGRID

C --- 113000 HSLEE Change Alignment Routine ----------

      NALIGN=100

      IALIGN = 1

777   CONTINUE

      do n = 1 , nx
        xuan2(n) = xuan(n)
        xutn2(n) = xutn(n)
        xuau2(n) = xuau(n)
        xutu2(n) = xutu(n)
      enddo

      CALL SEPTV
      CALL PSFTWK
      CALL PSFUWK
      CALL SOURCE
      CALL KBLADE
      CALL OBLADE
      CALL ULTWAKE
      CALL SOLVE
      CALL LLFOR

      CALL ALIGN(DAMP)

        err = 0.0
        do n = 2 , nx-1
          err1 = abs(xuan2(n) - xuan(n))
          err2 = abs(xuau2(n) - xuau(n))
          err3 = abs(xutn2(n) - xutn(n))
          err4 = abs(xutu2(n) - xutu(n))
          err = amax1(err,err1,err2,err3,err4)
        enddo

        if(err .lt. 0.001 .or. ialign .gt. Nalign) go to 1000

        ialign = ialign + 1
        go to 777

1000    continue

        write(*,*)
      WRITE(*,*) ' ========================================'
      WRITE(*,*) '   Wake Alignment Routine was Completed. '
      WRITE(*,*) '   Number of Iteration = ', ialign       
      WRITE(*,*) '   Maximum Error       = ', err          
      WRITE(*,*) ' ========================================'
        write(*,*)

C 2    CONTINUE

C ------- 113000  --------------End of ALignment Routine ---------------
C
C-----------------------------------------------------------------------
C     Write ouput for PSF-10
C-----------------------------------------------------------------------
      DO 3 N=1,NX
         XUA(N)=XUAN(N)
         XUT(N)=XUTN(N)
 3    CONTINUE


      DO N = 1, NX
         XR(N) = XRTMP(N)
         XPI(N) = XPITMP(N)
         XRAKE(N) = XRAKETMP(N)
         XSKEW(N) = XSKEWTMP(N)
         XCHD(N) = XCHDTMP(N)
         XCI(N) = XCITMP(N)
         XTI(N) = XTITMP(N)
         DO K = 1, 15
            XCS(N,K) = XCSTMP(N,K)
         ENDDO
         DO K = 1, 16
            XTS(N,K) = XTSTMP(N,K)
         ENDDO
      ENDDO

      rult = rultsave
      rhult = rhultsave

      RETURN
      END



      SUBROUTINE BLINPT(KTHK,KCAM)
C***********************************************************************
C     BLINPT: BLade INPuT
C       --- Read blade input geometry, non-dimensionalize on prop radius
C           prepare factored matrices for MAP2 subroutine
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      COMMON/HIMAGE/JHUB,JDUCT
      COMMON/HIMAGE2/RRHUB
      DIMENSION XPCORR(NXMAX),CAM(15),THK(16)
      DATA PI/3.141593/
C---- NACA A=0.8 MEANLINE CAMBERLINE SHAPE:
      DATA CAM/.0755,.1586,.2712,.4482,.6993,.8635,.9615,1.0,.9786,
     *  .8892,.7027,.3588,.1713,.0823,.0307/
C---- NACA 66 (MOD) THICKNESS FORM:
C -- HSLEE (122199) 
C      DATA THK/.1870,.2932,.4132,.5814,.8,.9274,.9904,
C     *  0.9917,0.9256,0.7934,0.5950,0.3306,0.1736,0.0888,0.0360,0.0/
C -- HSLEE

      DATA THK/.1870,.2932,.4132,.5814,.8,.9274,.9904,
     *  0.9924,0.9281,0.7966,0.5981,0.3328,0.1747,0.0893,0.0364,0.0/

C-----------------------------------------------------------------------
C     Set NACA=0.8 load distribution (MLTYPE=1)
C         Modified NACA66 thickness form (MTHICK=1)
C         Viscous pitch correction (KPCORR=1)
C-----------------------------------------------------------------------
      MLTYPE=1
      MTHICK=1

C -- Kill viscous pitch correction to compare with PUF3a
C   HSLEE (121099)

      KPCORR = ivpitch

CT      KPCORR=1
CT      Kpcorr = 0

C --121099

      IF(KCAM .NE. 99) THEN
        DO 10 J=1,15
        DO 10 N=1,NX
          XCS(N,J)=CAM(J)*XCAM(N)
 10     CONTINUE
      ENDIF 

      IF(KTHK .NE. 99) THEN      
        DO 20 J=1,16
        DO 20 N=1,NX
          IF(XCHD(N).EQ.0.0) GO TO 30
          XTS(N,J)=THK(J)*XTHICK(N)/XCHD(N)
          GO TO 20
 30       XTS(N,J)=THK(J)*XTHICK(N)/0.1
 20       CONTINUE
      ENDIF

C-----------------------------------------------------------------------
C     Viscous pitch correction, after KERWIN and LEE 
C-----------------------------------------------------------------------
      DO 40 N=1,NX
      DELA=0.0
      IF(XCHD(N).LT.0.01) GO TO 50
C      DELA=1.9454*ABS(XCAM(N))*XTHICK(N)/XCHD(N)
      DELA=1.9454*XCAM(N)*XTHICK(N)/XCHD(N)
      IF(KPCORR.NE.1) DELA=0.0
 50   PHI=ATAN(XPI(N)/(PI*XR(N)))
      PHI=PHI-DELA
 40   XPCORR(N)=PI*XR(N)*TAN(PHI)

C-----------------------------------------------------------------------
C     Non-dimensionalize everything on propeller radius
C-----------------------------------------------------------------------
      RHUB=XR(1)
C      RRHUB = RHUB
      EAR=(2.*NBLADE/PI)*SIMPUN(XR,XCHD,NX)
      DO 100 N=1,NX
      XPI(N)=XPCORR(N)
      XRSQ(N)=1.0-SQRT(ABS(1.0-XR(N)))      
      XPI(N)=XPI(N)*2.0
      XRAKE(N)=XRAKE(N)*2.0
      XSKEW(N)=XSKEW(N)*1.745329E-02
      XCHD(N)=XCHD(N)*2.0
      XCAM(N)=XCAM(N)*XCHD(N)
 100  XTHICK(N)=XTHICK(N)*2.0
      DO 110 J=1,15
      DO 110 N=1,NX
 110  XCS(N,J)=XCS(N,J)*XCHD(N)
      DO 120 J=1,16
      DO 120 N=1,NX
      IF(XCHD(N).EQ.0.0) GO TO 130
      XTS(N,J)=XTS(N,J)*XCHD(N)
      GO TO 120
 130  XTS(N,J)=XTS(N,J)*0.1
 120  CONTINUE

      CALL MAP1(XR,XRFAC,IPVXR,NX)
      CALL MAP1(XRSQ,XRSFAC,IPVXRS,NX)
      DO 140 J=1,17
 140  PSQ(J)=SQRT(PER(J))

c      write(999,*) 'xpi,xskew,xrake,xchd'
c      do i = 1 , nx
c      write(999,1100) xpi(i),xskew(i),xrake(i),xchd(i)
c      enddo
c
c      write(999,*) ' THICKNESS'
c      do i = 1 , nx
c       write(999,1200) (xts(i,j), j = 1, 15)
c      enddo 
c      write(999,*) ' CAMBER'
c      do i = 1 , nx
c       write(999,1200) (xcs(i,j), j = 1, 15)
c      enddo 
c1100    format(4(1x,f14.7))
c1200    format(5(1x,f14.7))

      RETURN
C))))))))))))))))))) End of subroutiine BLINPT ((((((((((((((((((((((((
      END



      SUBROUTINE ADINPT
C***********************************************************************
C     ADINPT: ADministrative INPuT file 
C      --- Reads administrative data file and computes parameters
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION DUM(NXMAX)
C-----------------------------------------------------------------------
C     Set parameters
C-----------------------------------------------------------------------
      NF=2
      KSPACE=2
      DO 10 M=1,NBLADE
         MPAN(M)=8
         NPAN(M)=8
 10   CONTINUE
      DIAM=1.0
      VSKTS=10.0
      VSFPS=VSKTS*1.6889
      RPM=VSFPS/ADVCO/DIAM*60.0
      XTW=XULT
      NPUW=36
      KSEPTV=1
      DISPN=0.0
      DISPR=0.0
      KTHICK=1
      RHO=1.9905
      CDRAG=0.0070
      SFC=0.333
      IWTEST=0
      ICDRAG=0
      ISFC=0
      DO 60 N=1,NX
      XUAN(N)=0.0
      XUAU(N)=0.0
      XUTN(N)=0.0
 60   XUTU(N)=0.0
      IPASS=0
C      IRUN=IRUN+1
      DO 70 M=1,NX
 70   DUM(M)=XVA(M)*XR(M)
      VOLMNV=SIMPUN(XR,DUM,NX)*2.0/(1.0-RHUB**2)
      RADIUS=DIAM/2.
      RPS=RPM/60.
      OMEGA=RPS*2.*3.141593
      ADVCOM=ADVCO*VOLMNV
      DO 80 M=1,NX
 80   XVSTAR(M)=SQRT(XVA(M)**2+((3.141593*XR(M)/ADVCO)+XVT(M))**2)
      DO 100 N=1,NX
 100  XCDRAG(N)=CDRAG
      DO 120 N=1,NX
 120  XSFC(N)=SFC
      CONTINUE
      RETURN
C))))))))))))))))))) End of subroutiine ADINPT ((((((((((((((((((((((((
      END



      SUBROUTINE SETUP
C***********************************************************************
C     SETUP: SET UP index
C     Sets up indexing sums and checks for proper jumping
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
C      DATA IBGMAX/1200/
      DATA ICPMAX/144/
C      DATA ITWMAX/5400/
      NCPTOT=NPAN(1)*MPAN(1)
      IF(NCPTOT.LE.ICPMAX) GO TO 60
 60   DO 80 I=1,NBLADE
      MGRID(I)=MPAN(I)+1
 80   NGRID(I)=NPAN(I)+1
      NGPSUM(1)=0
      NCPSUM(1)=0
      NCVSUM(1)=0
      NTWSUM(1)=0
      DO 90 I=2,NBLADE
      NCVSUM(I)=NCVSUM(I-1)+MGRID(I-1)*NPAN(I-1)
      NGPSUM(I)=NGPSUM(I-1)+MGRID(I-1)*NGRID(I-1)
      NCPSUM(I)=NCPSUM(I-1)+MPAN(I-1)*NPAN(I-1)
 90   NTWSUM(I)=NTWSUM(I-1)+MGRID(I-1)*80
      NGPTOT=NGPSUM(NBLADE)+MGRID(NBLADE)*NGRID(NBLADE)
      NTWTOT=NTWSUM(NBLADE)+MGRID(NBLADE)*80
      RETURN
C))))))))))))))))))) End of subroutiine SETUP (((((((((((((((((((((((((
      END

      SUBROUTINE BLGRID
C***********************************************************************
C     BLGRID: BLade GRID
C      --- Subroutine to set up blade grid points on all blades
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION DUMX(16),DUM1(16),DUM2(16),DUM3(16)      
      DATA PI/3.141593/

      DELBL=2.0*PI/NBLADE
C-----------------------------------------------------------------------
C     Loop over blades
C-----------------------------------------------------------------------
      DO 100 JJ=1,NBLADE
      K=NBLADE-JJ+1
      TZERO=(K-1)*DELBL
      NCP=NPAN(K)
      MCP=MPAN(K)
      NGP=NGRID(K)
      MGP=MGRID(K)
      DR=(1.0-RHUB)/(4.0*MCP+2)
      RZ(1)=RHUB+DR
      DR=DR*4.0
      DO 10 M=1,MGP         
      IF(M.GT.1) RZ(M)=RZ(M-1)+DR
      RZSQ(M)=1.0-SQRT(1.0-RZ(M))
      IF(M.GT.1) R(M-1)=0.5*(RZ(M)+RZ(M-1))
 10   CONTINUE
C-----------------------------------------------------------------------
C     Splineinterpolation in radial direction
C-----------------------------------------------------------------------
      CALL MAP2(XR,XPI,NX,RZ,PITCH,MGP,XRFAC,IPVXR)
      CALL MAP2(XR,XRAKE,NX,RZ,RAKE,MGP,XRFAC,IPVXR)
      CALL MAP2(XR,XSKEW,NX,RZ,SKEW,MGP,XRFAC,IPVXR)
      CALL MAP2(XRSQ,XCHD,NX,RZSQ,CHORZ,MGP,XRSFAC,IPVXRS)
      IF(K.NE.1) GO TO 30
      CALL MAP1(R,RFAC,IPVR,MCP)
      DO 20 M=1,MCP
 20   RSQ(M)=1.0-SQRT(1.0-R(M))
      CALL MAP2(XRSQ,XCHD,NX,RSQ,CHORCP,MCP,XRSFAC,IPVXRS)
C-----------------------------------------------------------------------
C     Tip chord for equal area of parabolic segment
C-----------------------------------------------------------------------
 30   CHORZ(MGP)=AMAX1(CHORZ(MGP),7.0*CHORZ(MGP-1)/15.0)
      CALL SPACE2(XV,XC,FC,CPDC,NCP,K,KSPACE)
      DO 100 M=1,MGP
      CALL CAMFIT(RZ(M),XV,HC,NGP,SLPLE)
      PHI=ATAN(PITCH(M)/(2.0*PI*RZ(M)))
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      DO 100 N=1,NGP
      S=(XV(N)-0.5)*CHORZ(M)
      F=HC(N)
      T=SKEW(M)+(S*CPHI+F*SPHI)/RZ(M)+TZERO
      L=INDXBG(N,M,K)
C-----------------------------------------------------------------------
C     Cartesian coordinates of blade grid points
C-----------------------------------------------------------------------
      X(L)=RAKE(M)+S*SPHI-F*CPHI
      Y(L)=RZ(M)*COS(T)
      Z(L)=RZ(M)*SIN(T)
 100  CONTINUE

C-----------------------------------------------------------------------
C     Compute control point locations and normals on key blade
C-----------------------------------------------------------------------
      NCP=NPAN(1)
      MCP=MPAN(1)
      DO 110 M=1,MCP
      DO 110 N=1,NCP
C-----------------------------------------------------------------------
C     Control point location
C-----------------------------------------------------------------------
      LL=INDXCP(N,M,1)
      L1=INDXBG(N,M,1)
      L2=L1+1
      L3=INDXBG(N,M+1,1)
      L4=L3+1
      BUG=1.0-FC(N)
      XP(LL)=0.5*(FC(N)*(X(L2)+X(L4))+BUG*(X(L1)+X(L3)))
      YP(LL)=0.5*(FC(N)*(Y(L2)+Y(L4))+BUG*(Y(L1)+Y(L3)))
      ZP(LL)=0.5*(FC(N)*(Z(L2)+Z(L4))+BUG*(Z(L1)+Z(L3)))
C-----------------------------------------------------------------------
C     Normal vector
C-----------------------------------------------------------------------
      A1=X(L4)-X(L1)
      A2=Y(L4)-Y(L1)
      A3=Z(L4)-Z(L1)
      B1=X(L3)-X(L2)
      B2=Y(L3)-Y(L2)
      B3=Z(L3)-Z(L2)
      CALL VECXPR(A1,A2,A3,B1,B2,B3,C1,C2,C3)
      CALL VECNML(C1,C2,C3)
C-----------------------------------------------------------------------
C     Special treatment at the last pannel
C-----------------------------------------------------------------------
      DUMX(N)=XP(LL)
      DUM1(N)=C1
      DUM2(N)=C2
      DUM3(N)=C3
      IF(N.EQ.NCP) THEN         
        CALL SPLINEs(DUMX,DUM1,NCP-1,DUMX(NCP),C1,1)
        CALL SPLINEs(DUMX,DUM2,NCP-1,DUMX(NCP),C2,1)
        CALL SPLINEs(DUMX,DUM3,NCP-1,DUMX(NCP),C3,1)
        DDDD = SQRT(C1**2+C2**2+C3**2)
        C1 = C1 / DDDD
        C2 = C2 / DDDD
        C3 = C3 / DDDD
      END IF
      XON(LL)=C1
      YON(LL)=C2
      ZON(LL)=C3 
 110  CONTINUE
      RETURN
      END



      SUBROUTINE SEPTV
C***********************************************************************
C     This subroutine calculates geometry of separated tip
C     vortex on key blade only
C     DISPN = DISPlacement of collection point Normal to blade
C     DISPR = Radial DISPlacement
C
      INCLUDE 'WAKE.INC'
      DIMENSION DUM1(3,21),DUM2(20,21),DUM3(3,21)
      DATA PI/3.141593/
      RAD=PI/180.
      MCP=MPAN(1)
      NCP=NPAN(1)
      MGP=MGRID(1)
      NGP=NGRID(1)
C-----------------------------------------------------------------------
C     Calculate normal displacement
C-----------------------------------------------------------------------
      DO 20 N=1,NCP
      BUG=1.0-XV(N)
      DO 20 L=N,NGP
      LL=L-N+1
 20   DUM2(N,LL)=(XV(L)-XV(N))/BUG
C-----------------------------------------------------------------------
C     Calculate normal displacement vectors
C-----------------------------------------------------------------------
      DO 30 N=1,NGP
      LT=INDXBG(N,MGP,1)
      V1=X(LT)
      V2=Y(LT)
      V3=Z(LT)
      CALL VECNML(V1,V2,V3)
      DUM3(1,N)=V1
      DUM3(2,N)=V2
 30   DUM3(3,N)=V3
      DO 40 N=1,NCP
      L1=INDXBG(N,MGP,1)
      L2=INDXBG(N,MCP,1)
      L3=INDXBG(N+1,MGP,1)
      T1=X(L3)-X(L1)
      T2=Y(L3)-Y(L1)
      T3=Z(L3)-Z(L1)
      R1=X(L1)-X(L2)
      R2=Y(L1)-Y(L2)
      R3=Z(L1)-Z(L2)
      CALL VECXPR(T1,T2,T3,R1,R2,R3,EM1,EM2,EM3)
      CALL VECNML(EM1,EM2,EM3)
      DUM1(1,N)=EM1
      DUM1(2,N)=EM2
 40   DUM1(3,N)=EM3
      L1=INDXBG(NCP,MGP,1)
      L2=INDXBG(NGP,MCP,1)
      L3=INDXBG(NGP,MGP,1)
      T1=X(L3)-X(L1)
      T2=Y(L3)-Y(L1)
      T3=Z(L3)-Z(L1)
      R1=X(L3)-X(L2)
      R2=Y(L3)-Y(L2)
      R3=Z(L3)-Z(L2)
      CALL VECXPR(T1,T2,T3,R1,R2,R3,EM1,EM2,EM3)
      CALL VECNML(EM1,EM2,EM3)
      DUM1(1,NGP)=EM1
      DUM1(2,NGP)=EM2
      DUM1(3,NGP)=EM3
C-----------------------------------------------------------------------
C     Collection point coordinates
C-----------------------------------------------------------------------
      JJ=INDXBG(NGP,MGP,1)
      XCOLL=X(JJ)+DISPN*DUM1(1,NGP)+DISPR*DUM3(1,NGP)
      YCOLL=Y(JJ)+DISPN*DUM1(2,NGP)+DISPR*DUM3(2,NGP)
      ZCOLL=Z(JJ)+DISPN*DUM1(3,NGP)+DISPR*DUM3(3,NGP)
C-----------------------------------------------------------------------
C     Separated tip vortex locations
C-----------------------------------------------------------------------
      DO 50 N=1,NCP
      NTEM=NCP-N+1
      DO 60 L=1,NTEM
      KL=INDXBG(N+L-1,MGP,1)
      XTIP(N,L)=X(KL)+DUM2(N,L)*(DUM1(1,N)*DISPN+DUM3(1,N)*DISPR)
      YTIP(N,L)=Y(KL)+DUM2(N,L)*(DUM1(2,N)*DISPN+DUM3(2,N)*DISPR)
 60   ZTIP(N,L)=Z(KL)+DUM2(N,L)*(DUM1(3,N)*DISPN+DUM3(3,N)*DISPR)
      XTIP(N,NTEM+1)=XCOLL
      YTIP(N,NTEM+1)=YCOLL
      ZTIP(N,NTEM+1)=ZCOLL
 50   CONTINUE
      RETURN
C))))))))))))))))))) End of subroutiine SEPTV (((((((((((((((((((((((((
      END



      SUBROUTINE PSFTWK
C***********************************************************************
C     PSFTWK:PSF-2 Transition WaKe
C      --- Geometry of transition wake
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION RTBL(11,21),XTBL(11),RW(21),UASTAR(21),UAUSTR(21),
     *          UTSTAR(21),UTUSTR(21)
      DATA PI/3.141593/
      DTF=15.*ADVCO/180.0
      DTRADF=15.*1.745328E-02
      DTI=6.0*ADVCO/180.
      DTRADI=6.0*1.745328E-02
      DO 10 K=1,NBLADE
      MGP=MGRID(K)
      NGP=NGRID(K)
 20   DO 30 M=1,MGP
      L=INDXBG(NGP,M,K)
 30   RW(M)=SQRT(Y(L)**2+Z(L)**2)
      IF(K.EQ.1) RW(MGP)=SQRT(YCOLL**2+ZCOLL**2)
CHSLEE (041300)
      IF(DCD.EQ.0.0) RULT=RW(MGP)
      IF(DCD.EQ.0.0) RHULT=RW(1)
      IF(RHULT.LT.0.10) RHULT=0.10
      CALL MAP2(XR,XVA,NX,RW,VA,MGP,XRFAC,IPVXR)
      CALL MAP2(XR,XVT,NX,RW,VT,MGP,XRFAC,IPVXR)
      CALL MAP2(XR,XUAN,NX,RW,UASTAR,MGP,XRFAC,IPVXR)
      CALL MAP2(XR,XUAU,NX,RW,UAUSTR,MGP,XRFAC,IPVXR)
      CALL MAP2(XR,XUTN,NX,RW,UTSTAR,MGP,XRFAC,IPVXR)
      CALL MAP2(XR,XUTU,NX,RW,UTUSTR,MGP,XRFAC,IPVXR)
      CALL RWTABLE(RULT,DCD,RTBL,XTBL,XTW,RW,MGP)
      DO 10 M=1,MGP
C-----------------------------------------------------------------------
C     Initial coordinates of transition wake from blade grid
C-----------------------------------------------------------------------
      L1=INDXBG(NGP,M,K)
      L3=INDXTW(1,M,K)
      XW(L3)=X(L1)
      YW(L3)=Y(L1)
      ZW(L3)=Z(L1)
      IF(K.NE.1.OR.M.NE.MGP) GO TO 40
      XW(L3)=XCOLL
      YW(L3)=YCOLL
      ZW(L3)=ZCOLL
C-----------------------------------------------------------------------
C     Generate transition wake geometry
C-----------------------------------------------------------------------
 40   TWOLD=ATAN2(ZW(L3),YW(L3))
      UAINC=UAUSTR(M)-UASTAR(M)
      UTINC=UTUSTR(M)-UTSTAR(M)
      DO 50 N=2,80
      DT=DTF
      DTRAD=DTRADF
      IF(N.GE.4) GO TO 60
      DT=DTI
      DTRAD=DTRADI
 60   JJ=INDXTW(N-1,M,K)
      Q=(XW(JJ)-XW(L3))/XULT
      L=INDXTW(N,M,K)
      XW(L)=XW(L-1)+(VA(M)+UASTAR(M)+GROW(Q)*UAINC)*DT
      XWTE=XW(L)-XW(L3)
      RWK=RTBL(11,M)
      IF(XWTE.GE.XTW) GO TO 70
      DO 80 J=1,10
      IF(XWTE.LT.XTBL(J).OR.XWTE.GE.XTBL(J+1)) GO TO 80
      RWK=RTBL(J,M)+(RTBL(J+1,M)-RTBL(J,M))*(XWTE-XTBL(J))/
     *   (XTBL(J+1)-XTBL(J))
      GO TO 70
 80   CONTINUE
 70   TW=ATAN2(ZW(L-1),YW(L-1))
      IF(TW.LT.TWOLD) TW=TW+6.283185
      IF(RWK.LT.0.02) RWK=0.02
      TW=TW+DTRAD+(VT(M)+UTSTAR(M)+GROW(Q)*UTINC)*DT/RWK
      YW(L)=RWK*COS(TW)
      ZW(L)=RWK*SIN(TW)
      TWOLD=TW
      NSW(K,M)=N
      IF(XWTE.GE.XTW) GO TO 90
 50   CONTINUE
 90   CONTINUE
 10   CONTINUE
      IF(UWLAM.GT.0.01) GO TO 100
      UWLAM=(XVA(NX)+XUAU(NX))/((PI*RULT/ADVCO)+XUTU(NX)+XVT(NX))
 100  RETURN
C))))))))))))))))))) End of subroutiine PSFTWK ((((((((((((((((((((((((
      END



      SUBROUTINE RWTABLE(RULT,DCD,RTBL,XTBL,XTW,RW,MM)
C***********************************************************************
C     RWTABLE: TABLE of Radii for contracted transition Wake
C      --- Generates table of radii for contracted transition wake
C          RULT = ultimate slipstream radius   
C          RHULT= ultimate hub vortex radius
C          DCD  = contracion angle in degrees at blade tip
C          RTBL,XTBL = resulting  table at 11 axial positions
C          XTW = location of ultimate wake from tip
C          RW = initial radii
C
      DIMENSION RTBL(11,21),XTBL(11),RW(21)
      MM1=MM-1
      DO 10 M=1,MM
 10      RTBL(1,M)=RW(M)
      DC=-TAN(1.745329E-02*DCD)
      if(rult .ge. rw(mm)) then
       h = 0.0
      else
        H=RW(MM)-RULT
      endif

      IF(DCD.GT.1.0) XRW=3.0*(RULT-RW(MM))/DC
      IF(DCD.LE.1.0) XRW=XTW
 20   XTBL(1)=0.0
      IF(XRW.LT.0.1) XRW=0.1
      DO 30 N=2,11
      XTBL(N)=0.1*(N-1)*XTW
      P=XTBL(N)/XTW
      Q=XTBL(N)/XRW

C HSLEE (041300)
C      RTBL(N,1)=(RTBL(1,1)-RHULT)*(1.0-3.0*P**2+2.0*P**3)+RHULT

      rtbl(n,1) = rtbl(1,1)

      RTBL(N,MM)=RTBL(1,MM)+H*(-3.0*Q+3.0*Q**2-Q**3)
C      IF(Q.GT.1.0) RTBL(N,MM)=RULT
      IF(Q.GT.1.0) RTBL(N,MM)=rtbl(n-1,mm)
      DO 30 M=2,MM1
      D1=SQRT(RTBL(N,M-1)**2+RTBL(1,M)**2-RTBL(1,M-1)**2)-RTBL(N,M-1)
      IF(M.GT.MM/2+2) D1=AMIN1(D1,(RTBL(N,M-1)/RTBL(N,MM))**2*
     *   (RTBL(N,MM)-RTBL(N,M-1)))
      RTBL(N,M)=RTBL(N,M-1)+D1
 30   CONTINUE
      RETURN
C)))))))))))))))))) End of subroutiine RWTABLE ((((((((((((((((((((((((
      END



      SUBROUTINE PSFUWK
C***********************************************************************
C     PSFUWK: PSF-2 Ultimate WaKe
C      --- Geometry of ultimate wake
C 
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION DK(12)
      DATA PI/3.141593/
      MGP=MGRID(1)
      DTU=18.84956/NPUW
      DO 10 K=1,NBLADE
 10      DK(K)=6.283185*(K-1)/NBLADE
C-----------------------------------------------------------------------
C     Geometry of tip vortices from each blade
C-----------------------------------------------------------------------
      L1=INDXTW(1,MGP,1)
      XTE=XW(L1)
      L=INDXTW(NSW(1,MGP),MGP,1)
      THU=ATAN2(ZW(L),YW(L))
      DO 20 K=1,NBLADE
      YU(1,K)=RULT*COS(THU+DK(K))
      ZU(1,K)=RULT*SIN(THU+DK(K))
 20   XU(1,K)=XW(L)
      THU=THU+DTU
      VAVS=XVA(NX)
      VTVS=XVT(NX)
      UANEAR=XUAN(NX)
      UTNEAR=XUTN(NX)
      UAINC=XUAU(NX)-UANEAR
      UTINC=XUTU(NX)-UTNEAR
      DO 30 N=2,NPUW
      Q=(XU(N-1,1)-XTE)/XULT
      IF(Q.GE.1.0) GO TO 40
      P=GROW(Q)
      TANBW=(VAVS+UANEAR+P*UAINC)/((PI*RULT/ADVCO)+UTNEAR+VTVS+P*UTINC)
      DXU=TANBW*DTU*RULT
      GO TO 50
 40   DXU=UWLAM*DTU*RULT
 50   DO 60 K=1,NBLADE
      XU(N,K)=XU(N-1,K)+DXU
      YU(N,K)=RULT*COS(THU+DK(K))
 60   ZU(N,K)=RULT*SIN(THU+DK(K))
 30   THU=THU+DTU
C-----------------------------------------------------------------------
C     Coordinates of ultimate hub vortex
C-----------------------------------------------------------------------
      L=INDXTW(NSW(1,1),1,1)
      X1HV=XW(L)
      Y1HV=0.0
      Z1HV=0.0
      X2HV=XU(NPUW,1)
      Y2HV=0.0
      Z2HV=0.0
      RETURN
C))))))))))))))))))) End of subroutiine PSFUWK ((((((((((((((((((((((((
      END



      SUBROUTINE KBLADE
C***********************************************************************
C     KBLADE: Key BLADE 
C      --- Computes influence of key blade & key blade transition wake
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      MCP=MPAN(1)
      NCP=NPAN(1)
C-----------------------------------------------------------------------
C     Zero arrays before starting a matrix build-up
C-----------------------------------------------------------------------
      IF(IPASS.GE.2) GO TO 10
      DO 20 N=1,NCPTOT
 20   SOURCN(N)=0.0
 10   CONTINUE
      NCP2=NCPTOT**2
      DO 30 K=1,NCP2
 30   AA(K)=0.0
      IF(IPASS.GE.2) GO TO 40
      IF(IRUN.GE.2) GO TO 40
      DO 50 J=1,NCPTOT
      DO 50 I=1,NCPTOT
 50   AB(I,J)=0.0
 40   CONTINUE
      DO 60 J=1,MCP
      DO 60 L=1,NCP
      I=INDXCP(L,J,1)
 60   CALL HSKEY(I)
      RETURN
C))))))))))))))))))) End of subroutiine KBLADE ((((((((((((((((((((((((
      END



      SUBROUTINE OBLADE
C***********************************************************************
C     OBLADE: Other BLADEs 
C      --- Calculates influence of other blades and other blade transition
C          wakes
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      MCP=MPAN(1)
      NCP=NPAN(1)
      TIMEOV=0.0
      TIMEOW=0.0
      TIMEOS=0.0
C-----------------------------------------------------------------------
C     Set up interpolation matrices
C-----------------------------------------------------------------------
      IF(IPASS.GE.2.OR.IRUN.GE.2) GO TO 10
      CALL INTMAT
 10   CONTINUE
C-----------------------------------------------------------------------
C     Loop over other blades
C-----------------------------------------------------------------------
      IF(NBLADE.LE.1) RETURN

      DO 20 K=2,NBLADE
      MSO=MPAN(K)
      NSO=NPAN(K)
      ISO=MSO*NSO
      JUMPM=MPAN(1)/MSO
      JUMPN=NPAN(1)/NSO
      NSKB=JUMPM*JUMPN
      NSKBW=NCP*JUMPM
C-----------------------------------------------------------------------
C     Compute influence of other blade vortices
C-----------------------------------------------------------------------
      IF(IPASS.GE.2) GO TO 30
      IF(IRUN.GE.2) GO TO 30
      DO 40 J=1,ISO
      DO 40 I=1,NCPTOT
 40   AOBC(I,J)=0.0
C      DO 50 M=1,5
C     DO 50 N=1,5
      DO 50 M=1,MSO
      DO 50 N=1,NSO
      I=INDXCP(INDC(N),INDS(M),1)
 50   CALL HSOTHR(I,K)
      CALL COMPLT(ISO)
      DO 60 MO=1,MSO
      DO 60 NO=1,NSO
      L=(MO-1)*NSO+NO
      CALL JUMPB(NO,MO,JUMPN,JUMPM)
      DO 70 JJ=1,NSKB
      DO 70 I=1,NCPTOT
 70   AB(I,LKB(JJ))=AB(I,LKB(JJ))+AOBC(I,L)*DIST(JJ)
 60   CONTINUE
 30   CONTINUE
C-----------------------------------------------------------------------
C     Compute influence of other blade transition wake
C-----------------------------------------------------------------------
      DO 80 J=1,MSO
      DO 80 I=1,NCPTOT
 80   AOBC(I,J)=0.0
C      DO 90 M=1,5
C      DO 90 N=1,5
      DO 90 M=1,MSO
      DO 90 N=1,NSO
      I=INDXCP(INDC(N),INDS(M),1)
 90   CALL TWOTHR(I,K)
      CALL COMPLT(MSO)
      DO 100 MO=1,MSO
      CALL JUMPW(MO,JUMPM)
      DO 110 JJ=1,NSKBW
      IJG=(LKB(JJ)-1)*NCPTOT
      DO 110 I=1,NCPTOT
      JG=IJG+I
 110  AA(JG)=AA(JG)+AOBC(I,MO)*DIST(JJ)
 100  CONTINUE
C-----------------------------------------------------------------------
C     Compute influence of other blade sources
C-----------------------------------------------------------------------
      IF(IPASS.GE.2) GO TO 20
      IF(KTHICK.EQ.0) GO TO 20
      DO 120 I=1,NCPTOT
 120  AOBC(I,1)=0.0
C      DO 130 M=1,5
C      DO 130 N=1,5
      DO 130 M=1,MSO
      DO 130 N=1,NSO
      I=INDXCP(INDC(N),INDS(M),1)
 130  CALL SVOTHR(I,K)
      CALL COMPLT(1)
      DO 140 I=1,NCPTOT
 140  SOURCN(I)=SOURCN(I)+AOBC(I,1)
 20   CONTINUE
 180  CONTINUE

      RETURN
C))))))))))))))))))) End of subroutiine OBLADE ((((((((((((((((((((((((
      END



      SUBROUTINE HSKEY(I)
C***********************************************************************
C     HSKEY: HorseShoes on KEY blade 
C      --- Calculates normal velocities at control point due to
C          horseshoes on key blade and in key blade transition wake
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION DUM(20)
      NCP=NPAN(1)
      MCP=MPAN(1)
      NGP=NGRID(1)
      MGP=MGRID(1)
      SUM2=0.0
      MM1=MCP-1
C-----------------------------------------------------------------------
C     Compute influence of key blade singularities, excluding
C     tip panel
C-----------------------------------------------------------------------
      IF(IPASS.GE.2) GO TO 10
      DO 20 M=1,MM1
      IF(M.EQ.1) GO TO 30
      DO 40 N=1,NCP
 40   UTN(N,1)=UTN(N,2)
 30   SUM1=-SUM2
      SUM2=0.0
      DO 50 N=1,NCP
      JJ=INDXCP(N,M,1)
      L1=INDXBG(N,M,1)
      L2=INDXBG(N,M+1,1)
      CALL VORSGN(XP(I),YP(I),ZP(I),X(L1),Y(L1),Z(L1),X(L2),
     *            Y(L2),Z(L2),XON(I),YON(I),ZON(I),UBN(N),BUG,1)
      SOURCN(I)=SOURCN(I)+SB(JJ)*BUG*KTHICK
      IF(IRUN.GE.2) GO TO 50
      IF(M.GT.1) GO TO 60
      L3=INDXBG(N+1,M,1)
      CALL VORSGN(XP(I),YP(I),ZP(I),X(L1),Y(L1),Z(L1),X(L3),
     *            Y(L3),Z(L3),XON(I),YON(I),ZON(I),UTN(N,1),BUG,0)
      SUM1=SUM1-UTN(N,1)
 60   L4=INDXBG(N+1,M+1,1)
      CALL VORSGN(XP(I),YP(I),ZP(I),X(L2),Y(L2),Z(L2),X(L4),
     *            Y(L4),Z(L4),XON(I),YON(I),ZON(I),UTN(N,2),BUG,0)
      SUM2=SUM2+UTN(N,2)
 50   CONTINUE
      IF(IRUN.GE.2)GO TO 20
C-----------------------------------------------------------------------
C     Form horseshoes on blade
C-----------------------------------------------------------------------
      DO 70 N=1,NCP
      JJ=INDXCP(N,M,1)
      IF(N.GT.1) GO TO 80
      TEMP=SUM1+SUM2
      GO TO 70
 80   TEMP=TEMP+UTN(N-1,1)-UTN(N-1,2)
 70   AB(I,JJ)=UBN(N)+TEMP
 20   CONTINUE
 10   CONTINUE
C-----------------------------------------------------------------------
C     Compute influence of tip panel, including
C     separated tip vortex
C-----------------------------------------------------------------------
 90   SUM1=0.0
      SUM2=0.0
      DO 100 N=1,NCP
      L1=INDXBG(N,MCP,1)
      L2=INDXBG(N,MGP,1)
      JJ=INDXCP(N,MCP,1)
      CALL VORSGN(XP(I),YP(I),ZP(I),X(L1),Y(L1),Z(L1),X(L2),
     *            Y(L2),Z(L2),XON(I),YON(I),ZON(I),UBN(N),BUG,1)
      IF(IPASS.EQ.1)SOURCN(I)=SOURCN(I)+SB(JJ)*BUG*KTHICK
      L3=INDXBG(N+1,MCP,1)
      CALL VORSGN(XP(I),YP(I),ZP(I),X(L1),Y(L1),Z(L1),X(L3),
     *            Y(L3),Z(L3),XON(I),YON(I),ZON(I),UTN(N,1),BUG,0)
      SUM1=SUM1-UTN(N,1)
      BUG=0.0
      NTEM=NCP-N+1
      DO 110 LL=1,NTEM
      L1=LL+1
      CALL VORSGN(XP(I),YP(I),ZP(I),XTIP(N,LL),YTIP(N,LL),ZTIP(N,LL),
     *            XTIP(N,L1),YTIP(N,L1),ZTIP(N,L1),XON(I),YON(I),ZON(I),
     *            CAT,DOG,0)
 110  BUG=BUG+CAT
 100  DUM(N)=BUG
      DO 120 N=1,NCP
      JJ=INDXCP(N,MCP,1)
      JG=(JJ-1)*NCPTOT+I
      IF(N.GT.1) GO TO 130
      TEMP=SUM1
      GO TO 120
 130  TEMP=TEMP+UTN(N-1,1)
 120  AA(JG)=UBN(N)+TEMP+DUM(N)
C-----------------------------------------------------------------------
C     Compute influence of key blade transition wake
C-----------------------------------------------------------------------
      SUM2=0.0
      DO 150 M=1,MCP
      SUM1=-SUM2
      SUM2=0.0
      IF(M.GT.1) GO TO 160
      NSW1=NSW(1,M)-1
      DO 170 N=1,NSW1
      L1=INDXTW(N,M,1)
      L2=L1+1
      CALL VORSGN(XP(I),YP(I),ZP(I),XW(L1),YW(L1),ZW(L1),XW(L2),
     *            YW(L2),ZW(L2),XON(I),YON(I),ZON(I),BUG,CAT,0)
 170  SUM1=SUM1-BUG
 160  NSW2=NSW(1,M+1)-1
      DO 180 N=1,NSW2
      L1=INDXTW(N,M+1,1)
      L2=L1+1
      CALL VORSGN(XP(I),YP(I),ZP(I),XW(L1),YW(L1),ZW(L1),XW(L2),
     *            YW(L2),ZW(L2),XON(I),YON(I),ZON(I),BUG,CAT,0)
 180  SUM2=SUM2+BUG
      DO 150 NS=1,NCP
      JJ=INDXCP(NS,M,1)
      JG=(JJ-1)*NCPTOT+I
 150  AA(JG)=AA(JG)+SUM1+SUM2
      RETURN
C))))))))))))))))))) End of subroutiine HSKEY (((((((((((((((((((((((((
      END



      SUBROUTINE ULTWAKE
C***********************************************************************
C     ULTWAKE: ULTimate WAKE 
C      --- Computes influence of ultimate wake from all blades
C          on key blade control points
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      MCP=MPAN(1)
      NCP=NPAN(1)
      DO 10 I=1,NCPTOT
 10   AOBC(I,1)=0.0
C      DO 20 M=1,5
C      DO 20 N=1,5
      DO 20 M=1,MCP
      DO 20 N=1,NCP
      I=INDXCP(INDC(N),INDS(M),1)
      CALL UWKCP(I)
 20   CONTINUE
      CALL COMPLT(1)
      MDIV1=MCP/2
      IF(IPASS.LT.2) GO TO 30
C-----------------------------------------------------------------------
C     Search for max circulation
C-----------------------------------------------------------------------
      MDIV1=1
      GMAX=0.0
      DO 40 M=1,MCP
      SUMG=0.0
      DO 50 N=1,NCP
      I=INDXCP(N,M,1)
 50   SUMG=SUMG+GB(I)
      IF(SUMG.GE.GMAX) THEN
      MDIV1=M
      GMAX=SUMG
      END IF
 40   CONTINUE
 30   CONTINUE
      DO 70 N=1,NCP
      J1=INDXCP(N,MDIV1,1)
      IJG1=(J1-1)*NCPTOT
      DO 70 I=1,NCPTOT
      JG1=IJG1+I
 70   AA(JG1)=AA(JG1)+AOBC(I,1)
      RETURN
C)))))))))))))))))) End of subroutiine ULTWAKE ((((((((((((((((((((((((
      END



      SUBROUTINE UWKCP(I)
C***********************************************************************
C     UWKCP: Ultimate WaKe on selected Control Point 
C      --- Computes influence of ultimate wake from all blades on
C          selected control point
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
C-----------------------------------------------------------------------
C     Hub vortex from all blades
C-----------------------------------------------------------------------
      CALL VORSGN(XP(I),YP(I),ZP(I),X1HV,Y1HV,Z1HV,X2HV,
     *            Y2HV,Z2HV,XON(I),YON(I),ZON(I),BUG,CAT,0)
      AOBC(I,1)=-BUG*NBLADE
C-----------------------------------------------------------------------
C     Tip vortices from all blades
C-----------------------------------------------------------------------
      SUM2=0.0
      DO 10 K=1,NBLADE
      DO 10 N=2,NPUW
      CALL VORSGN(XP(I),YP(I),ZP(I),XU(N-1,K),YU(N-1,K),ZU(N-1,K),
     *            XU(N,K),YU(N,K),ZU(N,K),XON(I),YON(I),ZON(I),BUG,
     *            CAT,0)
 10   SUM2=SUM2+BUG
      AOBC(I,1)=AOBC(I,1)+SUM2
      RETURN
C))))))))))))))))))) End of subroutiine UWKCP (((((((((((((((((((((((((
      END



      SUBROUTINE HSOTHR(I,K)
C***********************************************************************
C     HSOTHR: HorseShoes ON OTHeR blades
C      --- Calculates normal velocities at control point due to
C          horseshoes on other blades
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      NSO=NPAN(K)
      MSO=MPAN(K)
      SUM2=0.0
C-----------------------------------------------------------------------
C     Compute influence of blade singularities
C-----------------------------------------------------------------------
      DO 10 M=1,MSO
      IF(M.EQ.1) GO TO 20
      DO 30 N=1,NSO
 30   UTN(N,1)=UTN(N,2)
 20   SUM1=-SUM2
      SUM2=0.0
      DO 40 N=1,NSO
      JJ=INDXCP(N,M,K)
      L1=INDXBG(N,M,K)
      L2=INDXBG(N,M+1,K)
      CALL VORSGN(XP(I),YP(I),ZP(I),X(L1),Y(L1),Z(L1),X(L2),
     *            Y(L2),Z(L2),XON(I),YON(I),ZON(I),UBN(N),BUG,0)
      IF(M.GT.1) GO TO 50
      L3=INDXBG(N+1,M,K)
      CALL VORSGN(XP(I),YP(I),ZP(I),X(L1),Y(L1),Z(L1),X(L3),
     *            Y(L3),Z(L3),XON(I),YON(I),ZON(I),UTN(N,1),BUG,0)
      SUM1=SUM1-UTN(N,1)
 50   L4=INDXBG(N+1,M+1,K)
      CALL VORSGN(XP(I),YP(I),ZP(I),X(L2),Y(L2),Z(L2),X(L4),
     *            Y(L4),Z(L4),XON(I),YON(I),ZON(I),UTN(N,2),BUG,0)
      SUM2=SUM2+UTN(N,2)
 40   CONTINUE
C-----------------------------------------------------------------------
C     Form horseshoes on blade
C-----------------------------------------------------------------------
      DO 60 N=1,NSO
      JK=(M-1)*NSO+N
      IF(N.GT.1) GO TO 70
      TEMP=SUM1+SUM2
      GO TO 60
 70   TEMP=TEMP+UTN(N-1,1)-UTN(N-1,2)
 60   AOBC(I,JK)=UBN(N)+TEMP
 10   CONTINUE
      RETURN
C))))))))))))))))))) End of subroutiine HSOTHR ((((((((((((((((((((((((
      END



      SUBROUTINE TWOTHR(I,K)
C***********************************************************************
C     TWOTHR: Transition Wake from OTHeR blades
C      --- Calculates normal velocity at control point due to
C          transition wake from other blades
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      MSO=MPAN(K)
      SUM2=0.0
      DO 10 M=1,MSO
      SUM1=-SUM2
      SUM2=0.0
      IF(M.GT.1) GO TO 20
      NSW1=NSW(K,M)-1
      DO 30 N=1,NSW1
      L1=INDXTW(N,M,K)
      L2=L1+1
      CALL VORSGN(XP(I),YP(I),ZP(I),XW(L1),YW(L1),ZW(L1),XW(L2),
     *            YW(L2),ZW(L2),XON(I),YON(I),ZON(I),BUG,CAT,0)
 30   SUM1=SUM1-BUG
 20   NSW2=NSW(K,M+1)-1
      DO 40 N=1,NSW2
      L1=INDXTW(N,M+1,K)
      L2=L1+1
      CALL VORSGN(XP(I),YP(I),ZP(I),XW(L1),YW(L1),ZW(L1),XW(L2),
     *            YW(L2),ZW(L2),XON(I),YON(I),ZON(I),BUG,CAT,0)
 40   SUM2=SUM2+BUG
 10   AOBC(I,M)=SUM1+SUM2
      RETURN
C))))))))))))))))))) End of subroutiine TWOTHR ((((((((((((((((((((((((
      END



      SUBROUTINE SOURCE
C***********************************************************************
C     SOURCE: SOURCE strengths
C      --- Calculates source strengths on all blades
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION DUM(21)
      DATA PI/3.141593/
      IF(IPASS.GE.2) RETURN
      MCP=MPAN(1)
      NCP=NPAN(1)
      CALL MAP2(XR,XVSTAR,NX,R,VSTAR,MCP,XRFAC,IPVXR)
      IF(IRUN.GE.2) GO TO 10
 10   DO 20 M=1,MCP
 20   SFACT(M)=VSTAR(M)/(2.0*PI)
      DO 30 M=1,MCP
      CALL THKFIT(R(M),XC,HT,NCP)
      DELR=RZ(M+1)-RZ(M)
      DO 40 N=1,NCP
      I=INDXCP(N,M,1)
      L1=INDXBG(N,M,1)
      L2=INDXBG(N,M+1,1)
      DUM(N)=VECDIS(X(L1),Y(L1),Z(L1),X(L2),Y(L2),Z(L2))
      CDELTH=DELR/DUM(N)
 40   SB(I)=SFACT(M)*HT(N)*CDELTH
C-----------------------------------------------------------------------
C     Correct source strengths for closure on each section
C-----------------------------------------------------------------------
      BUG=0.0
      DO 50 N=1,NCP
      I=INDXCP(N,M,1)
 50   BUG=BUG+SB(I)*DUM(N)
      CAT=BUG/NCP
      DO 60 N=1,NCP
      I=INDXCP(N,M,1)
      SB(I)=SB(I)-CAT/DUM(N)
 60   CONTINUE
 30   CONTINUE
      CALL BOUNDR(SB)
      RETURN
C))))))))))))))))))) End of subroutiine SOURCE ((((((((((((((((((((((((
      END



      SUBROUTINE BOUNDR(B)
C***********************************************************************
C     BOUNDR: BOUND singulaRity 
C      --- Calculates bound singularity strengths on other blades
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION B(*)
      IF(NBLADE.LE.1) RETURN
      MCP=MPAN(1)
      NCP=NPAN(1)
      DO 10 K=2,NBLADE
      MSO=MPAN(K)
      NSO=NPAN(K)
      JUMPM=MCP/MSO
      JUMPN=NCP/NSO
      NSKB=JUMPM*JUMPN
C-----------------------------------------------------------------------
C     Loop over coarse singularities
C-----------------------------------------------------------------------
      DO 10 MS=1,MSO
      DO 10 NS=1,NSO
      I=INDXCP(NS,MS,K)
      CALL JUMPB(NS,MS,JUMPN,JUMPM)
C-----------------------------------------------------------------------
C     Sum up contributing key blade singularities
C-----------------------------------------------------------------------
      B(I)=0.0
      DO 20 J=1,NSKB
 20   B(I)=B(I)+B(LKB(J))*DIST(J)
 10   CONTINUE
      RETURN
C))))))))))))))))))) End of subroutiine BOUNDR ((((((((((((((((((((((((
      END

 

      SUBROUTINE SVOTHR(I,K)
C***********************************************************************
C     SVOTHR: SolVe for OTHeR blades
C      --- Calculates influence of other blade sources on keyblade
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      MSO=MPAN(K)
      NSO=NPAN(K)
C-----------------------------------------------------------------------
C     Source strengths calculated in subroutine boundr
C-----------------------------------------------------------------------
      DO 10 M=1,MSO
      DO 10 N=1,NSO
      JJ=INDXCP(N,M,K)
      L1=INDXBG(N,M,K)
      L2=INDXBG(N,M+1,K)
      CALL VORSGN(XP(I),YP(I),ZP(I),X(L1),Y(L1),Z(L1),X(L2),
     *            Y(L2),Z(L2),XON(I),YON(I),ZON(I),CAT,BUG,1)
 10   AOBC(I,1)=AOBC(I,1)+SB(JJ)*BUG*KTHICK
      RETURN
C))))))))))))))))))) End of subroutiine SVOTHR ((((((((((((((((((((((((
      END



      SUBROUTINE INTMAT
C***********************************************************************
C     INTMAT: INTerpolation MATrices 
C      --- Computes interpolation matrices for use by COMPLT
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION DUM(8)
      MCP=MPAN(1)
      NCP=NPAN(1)

      do n = 1 , mcp
        inds(n) = n
        dum(n) = r(inds(n))
      enddo

       CALL MAPTBL(mcp,MCP,20,DUM,R,AS)

       do n = 1 , ncp
       indc(n) = n
       dum(n) = xc(indc(n))
       enddo

       CALL MAPTBL(ncp,NCP,20,DUM,XC,AC)

      RETURN
C))))))))))))))))))) End of subroutiine INTMAT ((((((((((((((((((((((((
      END



      SUBROUTINE COMPLT(NSING)
C***********************************************************************
C     COMPLT: COMPLeTe matrix
C      --- This subroutine does spanwise and chordwise interpolation
C          needed when calculating influence of other blades on
C          key blade.  AOBC array used for other blade, other wake,
C          and other source influence.
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION TEMP(20),DUM(8)
      NCP=NPAN(1)
      MCP=MPAN(1)
C-----------------------------------------------------------------------
C     Loop over number of singularities
C-----------------------------------------------------------------------
      DO 10 NS=1,NSING
C-----------------------------------------------------------------------
C     Chordwise interpolation
C-----------------------------------------------------------------------
C      DO 20 MC=1,5
      DO 20 MC=1,MCP
      MF=INDS(MC)
C      DO 30 NC=1,5
      DO 30 NC=1,NCP
      I=INDXCP(INDC(NC),MF,1)
 30   DUM(NC)=AOBC(I,NS)
      DO 40 NG=1,NCP
      I=INDXCP(NG,MF,1)
      TEMP(NG)=0.0
C      DO 50 JJ=1,5
      DO 50 JJ=1,8
 50   TEMP(NG)=TEMP(NG)+AC(NG,JJ)*DUM(JJ)
 40   AOBC(I,NS)=TEMP(NG)
 20   CONTINUE
C-----------------------------------------------------------------------
C     Spanwise interpolation
C-----------------------------------------------------------------------
      DO 60 NC=1,NCP
C      DO 70 MC=1,5
      DO 70 MC=1,MCP
      I=INDXCP(NC,INDS(MC),1)
 70   DUM(MC)=AOBC(I,NS)
      DO 80 MG=1,MCP
      I=INDXCP(NC,MG,1)
      TEMP(MG)=0.0
C      DO 90 JJ=1,5
      DO 90 JJ=1,8
 90   TEMP(MG)=TEMP(MG)+AS(MG,JJ)*DUM(JJ)
 80   AOBC(I,NS)=TEMP(MG)
 60   CONTINUE
 10   CONTINUE
      RETURN
C))))))))))))))))))) End of subroutiine COMPLT ((((((((((((((((((((((((
      END



      SUBROUTINE SOLVE
C***********************************************************************
C     SOLVE: SOLVE for bound circualtion strengths
C      --- This subroutine assembles rhs for equations and solves
C          for bound circualtion strengths
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION RP(20)
C-----------------------------------------------------------------------
C     Complete a matrix
C-----------------------------------------------------------------------
      DO 10 J=1,NCPTOT
      ILG=(J-1)*NCPTOT
      DO 10 I=1,NCPTOT
      LG=ILG+I
 10   AA(LG)=AA(LG)+AB(I,J)
C-----------------------------------------------------------------------
C     Apply weight factors to matrix
C-----------------------------------------------------------------------
      MCP=MPAN(1)
      NCP=NPAN(1)
      DO 20 I=1,NCPTOT
      DO 20 M=1,MCP
      DO 20 N=1,NCP
      JJ=INDXCP(N,M,1)
      LG=(JJ-1)*NCPTOT+I
 20   AA(LG)=AA(LG)*CPDC(N)
C-----------------------------------------------------------------------
C     Perform LU decomposition on a matrix
C-----------------------------------------------------------------------
      CALL FACTOR(AA,AA,IPIVA,RHS,NCPTOT,NCPTOT,IFLAG)
      IF(IFLAG.EQ.1) GO TO 30
 40   FORMAT(/2X,'A MATRIX SINGULAR',/)
      WRITE(6,40)
      STOP
 30   CONTINUE
C-----------------------------------------------------------------------
C     Compute right-hand side
C-----------------------------------------------------------------------
      DO 60 N=1,NCP
      DO 70 M=1,MCP
      I=INDXCP(N,M,1)
 70   RP(M)=SQRT(YP(I)**2+ZP(I)**2)
      CALL MAP2(XR,XVA,NX,RP,VA,MCP,XRFAC,IPVXR)
      CALL MAP2(XR,XVR,NX,RP,VR,MCP,XRFAC,IPVXR)
      CALL MAP2(XR,XVT,NX,RP,VT,MCP,XRFAC,IPVXR)
      DO 60 M=1,MCP
      L=INDXCP(N,M,1)
      V1=VA(M)
      VROT=(OMEGA*RP(M)*RADIUS/VSFPS)+VT(M)
      VRAD=VR(M)
      THET=ATAN2(ZP(L),YP(L))
      STHET=SIN(THET)
      CTHET=COS(THET)
      V2=-VROT*STHET+VRAD*CTHET
      V3=VROT*CTHET+VRAD*STHET
      VINFL=VECDPR(V1,V2,V3,XON(L),YON(L),ZON(L))
      RHS(L)=-VINFL-SOURCN(L)
 60   CONTINUE
C-----------------------------------------------------------------------
C     Perform back substitution
C-----------------------------------------------------------------------
      CALL SUBST(AA,RHS,GB,IPIVA,NCPTOT,NCPTOT)
C-----------------------------------------------------------------------
C     Unweight bounder strengths
C-----------------------------------------------------------------------
      DO 90 M=1,MCP
      DO 90 N=1,NCP
      I=INDXCP(N,M,1)
 90   GB(I)=GB(I)*CPDC(N)
      RETURN
C))))))))))))))))))) End of subroutiine SOLVE ((((((((((((((((((((((((
      END



      SUBROUTINE ALIGN(DAMP)
C***********************************************************************
C     ALIGN: ALIGNment
C      --- Aligns trailing vortex system with flow for a given
C          circulation distribution
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      common /lsqorder/ mlsqp
      DIMENSION CV(21,2),UANR(21),UAUR(21),UTNR(21),UTUR(21),
     *  XUANR(NXMAX),XUAUR(NXMAX),XUTNR(NXMAX),XUTUR(NXMAX)
      DIMENSION XNEAR(9),XFAR(9),XTE(9),UANEAR(9),UAFAR(9),
     *  UTNEAR(9),UTFAR(9)
      DIMENSION RCAL(9)
      DIMENSION INDSW(9)
      DATA PI/3.141593/
      DATA CORAD,PERULT/0.001,1.00/

      RAD=180./PI
C-----------------------------------------------------------------------
C     Select control point indices on trailers
C-----------------------------------------------------------------------
      MGP=MGRID(1)
      do i = 1 , 9
        indsw(i) = i
      enddo

      DO 10 I=1,9
      RCAL(I)=RZ(I)
10    CONTINUE

C-----------------------------------------------------------------------
C     Set bound singularity strengths on other blades
C-----------------------------------------------------------------------
      CALL BOUNDR(GB)
C-----------------------------------------------------------------------
C     Set chordwise vortex strengths
C-----------------------------------------------------------------------
      DO 20 K=1,NBLADE
         NGP=NGRID(K)
         MGP=MGRID(K)
         NCP=NGP-1
         MCP=MGP-1
         DO 30 J=1,NCP
         CV(J,1)=0.
 30      CV(J,2)=0.
         DO 40 M=1,MCP
         IF(M.EQ.1) GO TO 50
         DO 60 N=1,NCP
         CV(N,1)=CV(N,2)
 60      CV(N,2)=0.
 50      GSUM=0.
         DO 40 N=1,NCP
         I=INDXCP(N,M,K)
         GSUM=GSUM+GB(I)
         CV(N,1)=CV(N,1)-GSUM
         CV(N,2)=GSUM
         J=INDXCV(N,M,K)
         GC(J)=CV(N,1)
 40      CONTINUE
         DO 70 N=1,NCP
         I=INDXCV(N,MGP,K)
 70      GC(I)=CV(N,2)
 20   CONTINUE
      DO 80 K=1,NBLADE
      MGP=MGRID(K)
      NCP=NPAN(K)
      DO 80 M=1,MGP
      I=INDXCV(NCP,M,K)
 80   GAMTW(K,M)=GC(I)
      MGP=MGRID(1)
      NGP=NGRID(1)
      MCP=MGP-1
      NCP=NGP-1
      GAMULT=0.0
      DO 90 M=1,MCP
      SUM=0.0
      DO 100 N=1,NCP
      I=INDXCP(N,M,1)
 100  SUM=SUM+GB(I)
 90   GAMULT=AMAX1(SUM,GAMULT)
C-----------------------------------------------------------------------
C     Get directions from operator
C-----------------------------------------------------------------------
      VAMULT=1.0
      TBF=0.0
C      IF(ITUNL.EQ.0) GO TO 120
C      ICHECK=ICHECK+1
C 130  FORMAT(/2X,'ENTER APROP/ATUNNEL')
C      IF(ICHECK.EQ.1) THEN
C      WRITE(6,130)
C      READ(5,*)ALPHA
C      ELSE
C      WRITE(6,140)CT,ALPHA
C 140  FORMAT(/1X,'CT=',1X,F10.6,2X,'APROP/ATUNNEL=',1X,F10.6)
C      END IF
C      UBV=ALPHA*CT/(4.0*SQRT(1.0+CT))
C      VAMULT=(1.0-2.*UBV)/(1.0-UBV)
C      TBF=(VAMULT-1.0)*VOLMNV
C 120  CONTINUE
C-----------------------------------------------------------------------
C     Set ultimate wake
C-----------------------------------------------------------------------
      VAVS=VAMULT*XVA(NX)
      VTVS=XVT(NX)
      CALL UPITCH(NBLADE,RULT,GAMULT,PERULT,ADVCO,VAVS,VTVS,
     *            UWLAM,UATUW,UTTUW)
      UTHUW=-0.75*NBLADE*GAMULT*PERULT/RHULT
      UAHUW=0.75*GAMULT*PERULT*NBLADE/(RULT*UWLAM)
C-----------------------------------------------------------------------
C     Start loops in transition wake
C-----------------------------------------------------------------------
C      DO 160 NWI=1,NWIMAX
       DO 160 NWI = 1, 80
C-----------------------------------------------------------------------
C     Near wake induced velocities
C-----------------------------------------------------------------------
C      DO 180 I=1,5
      DO 180 I=1,9
      M=INDSW(I)
      LB=INDXBG(NGP,M,1)
      XTE(I)=X(LB)
      L1=INDXTW(1,M,1)
      L2=L1+1
      L3=L1+2
      XF=XW(L2)
      YF=YW(L2)
      ZF=ZW(L2)
      XNEAR(I)=XF
      CALL FPWAKE(XF,YF,ZF,UAI,VIY,VIZ,URI,UTI,1,1,1,1,M,2)
C      IF(I.NE.5) GO TO 190
      IF(I.NE.9) GO TO 190
      CALL SELFIN(XW(L1),YW(L1),ZW(L1),XW(L2),YW(L2),ZW(L2),XW(L3),
     *            YW(L3),ZW(L3),CORAD,UASI,VIY,VIZ,URSI,UTSI)
      UASI=UASI*GAMTW(1,M)
      UTSI=UTSI*GAMTW(1,M)
      URSI=URSI*GAMTW(1,M)
      GO TO 200
 190  UASI=0.0
      UTSI=0.0
      URSI=0.0
 200  UANEAR(I)=UAI+UASI
      UTNEAR(I)=UTI+UTSI
      URTOT=URI+URSI
 180  CONTINUE
C-----------------------------------------------------------------------
C     Far wake induced velocities
C-----------------------------------------------------------------------
C      DO 220 I=2,4

      if(xult .le. 2.5) then
         xxut = xult * 0.7
      else
         xxut = xult * 0.5
      endif

C -- Check the length of wake 

      x1max = -1000.0
      x1min = 1000.0

      DO I = 1 , 9
         J1 = indxtw(2,i,1)
         J2 = indxtw(nsw(1,I),i,1)

         xff1 = xw(j1) - xte(i)
         xff2 = xw(j2) - xte(i)

         x1max = amax1(xff1,x1max)
         x1min = amin1(xff2,x1min)
      ENDDO   

      if(xxut .le. x1max .or. xxut .gt. x1min) xxut = 0.8 * x1min
      
      DO 220 I=1,9

C --- 113000  HSLEE Calculate Panel index of the Induced velocity
C           calculation points having same distance from blade T.E.
C 
         do JJ = 2 , nsw(1,I)-1
            J2 = indxtw(jj,i,1)
            J3 = indxtw(jj+1,i,1)
            xff1 = xw(j2) - xte(I)
            xff2 = xw(j3) - xte(I)
            if(xff1 .le. xxut .and. xff2 .ge. xxut) then
               IW = JJ
               go to 1234
            endif
         enddo
         M=INDSW(I)
1234    continue

         M=INDSW(I)
         L1=INDXTW(IW-1,M,1)
         L2=L1+1
         L3=L1+2
         XF=XW(L2)
         YF=YW(L2)
         ZF=ZW(L2)
         XFAR(I)=XF
         CALL FPWAKE(XF,YF,ZF,UAI,VIY,VIZ,URI,UTI,1,0,1,1,M,IW)
         UASI=0.0
         UTSI=0.0
         URSI=0.0
         UAFAR(I)=UAI+UASI
         UTFAR(I)=UTI+UTSI
         URTOT=URI+URSI
 220  CONTINUE

      XFAR(1)=XTE(1)+XULT
      XFAR(9)=XTE(9)+XULT

C-----------------------------------------------------------------------
C     Get values just behind blade and in ultimate wake
C-----------------------------------------------------------------------
      DO 230 I=1,9
         QN=(XNEAR(I)-XTE(I))/XULT
         QF=(XFAR(I)-XTE(I))/XULT
         PN=GROW(QN)
         PF=GROW(QF)
         UAINC=(UAFAR(I)-UANEAR(I))/(PF-PN)
         UTINC=(UTFAR(I)-UTNEAR(I))/(PF-PN)
         UANR(I)=UANEAR(I)-UAINC*PN
         UTNR(I)=UTNEAR(I)-UTINC*PN
         UAUR(I)=UANR(I)+UAINC
         UTUR(I)=UTNR(I)+UTINC
230   CONTINUE

C --- 113000  HSLEE Smooth induced velocity in ultimate wake

        CALL LSQ(RCAL,UAUR,9,mlsqp)
        CALL LSQ(RCAL,UTUR,9,mlsqp)

C ---- 113000 ------------------------------------------------------

C-----------------------------------------------------------------------
C     Get new convection velocities at input radii
C-----------------------------------------------------------------------
      CALL SPLINE(RCAL,UANR,9,XR,XUANR,NX)
      CALL SPLINE(RCAL,UAUR,9,XR,XUAUR,NX)
      CALL SPLINE(RCAL,UTNR,9,XR,XUTNR,NX)
      CALL SPLINE(RCAL,UTUR,9,XR,XUTUR,NX)

C-----------------------------------------------------------------------
C     Correction for tunnel backflow
C-----------------------------------------------------------------------
      DO 240 M=1,NX
 240  XUAUR(M)=XUAUR(M)+TBF
C-----------------------------------------------------------------------
C     Set new values of wake convection velocities, including damping
C-----------------------------------------------------------------------
      BUG=1.0-DAMP
      ERR=0.0

      DO N = 2 , NX - 1
          ERR1=ABS(XUAN(N)-XUANR(N))
          ERR2=ABS(XUAU(N)-XUAUR(N))
          ERR3=ABS(XUTN(N)-XUTNR(N))
          ERR4=ABS(XUTU(N)-XUTUR(N))
          ERR=AMAX1(ERR,ERR1,ERR2,ERR3,ERR4)
      ENDDO

      DO 250 N=1,NX
         XUAN(N)=DAMP*XUAN(N)+BUG*XUANR(N)
         XUAU(N)=DAMP*XUAU(N)+BUG*XUAUR(N)
         XUTN(N)=DAMP*XUTN(N)+BUG*XUTNR(N)
         XUTU(N)=DAMP*XUTU(N)+BUG*XUTUR(N)
 250  CONTINUE

      XUAN(1) = 0.5*(XUAN(1) + XUAN(2))
      XUAU(1) = 0.5*(XUAU(1) + XUAU(2))
      XUTN(1) = 0.5*(XUTN(1) + XUTN(2))
      XUTU(1) = 0.5*(XUTU(1) + XUTU(2))

      CALL PSFTWK
      CALL PSFUWK
      IF(ERR.LT.0.005) GO TO 280
      IF(NWI.NE.80) GO TO 160

 160  CONTINUE
 280  CONTINUE

      RETURN
C))))))))))))))))))) End of subroutiine ALIGN ((((((((((((((((((((((((
      END



      SUBROUTINE FPWAKE(XF,YF,ZF,VIX,VIY,VIZ,VIR,VIT,IBLADE,ISEPTV,
     *                  ITWAKE,IUWAKE,MW,NW)
C***********************************************************************
C     FPWAKE: Field Point velocity in the WAKE
C      --- Calculates induced velocities at field points
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION VEL(3),VFP(3),VS(3),VSUM(3)
      ATHICK=KTHICK
      DO 10 J=1,3
 10   VFP(J)=0.0
      IF(IBLADE.NE.1) GO TO 20
C-----------------------------------------------------------------------
C     Bound singularities
C-----------------------------------------------------------------------
      DO 30 K=1,NBLADE
         NCP=NPAN(K)
         MCP=MPAN(K)
         DO 30 M=1,MCP
         DO 30 N=1,NCP
            L1=INDXBG(N,M,K)
            L2=INDXBG(N,M+1,K)
            I=INDXCP(N,M,K)
            II=INDXCP(N,M,1)
            CALL VORSEG(XF,YF,ZF,X(L1),Y(L1),Z(L1),X(L2),Y(L2),Z(L2),
     *                  VEL(1),VEL(2),VEL(3),VS(1),VS(2),VS(3),KTHICK)
            DO 30 J=1,3
            VFP(J)=VFP(J)+VEL(J)*GB(I)+VS(J)*SB(II)*ATHICK
30    CONTINUE
C-----------------------------------------------------------------------
C     Chordwise vortices
C-----------------------------------------------------------------------
      DO 40 K=1,NBLADE
      NCP=NPAN(K)
      MGP=MGRID(K)
      IF(K.EQ.1.AND.ISEPTV.EQ.1) MGP=MGP-1
      DO 40 M=1,MGP
      DO 40 N=1,NCP
         L1=INDXBG(N,M,K)
         L2=L1+1
         I=INDXCV(N,M,K)
         KDUMY=0
         CALL VORSEG(XF,YF,ZF,X(L1),Y(L1),Z(L1),X(L2),Y(L2),Z(L2),
     *               VEL(1),VEL(2),VEL(3),BUG,BUG,BUG,KDUMY)
         DO 40 J=1,3
            VFP(J)=VFP(J)+VEL(J)*GC(I)
 40   CONTINUE
 20   CONTINUE
      IF(ISEPTV.NE.1) GO TO 50
C-----------------------------------------------------------------------
C     Separated tip vortex
C-----------------------------------------------------------------------
      NCP=NPAN(1)
      MCP=MPAN(1)
      DO 60 N=1,NCP
      DO 70 J=1,3
 70   VSUM(J)=0.0
      I=INDXCP(N,MCP,1)
      NTEM=NCP-N+1
      DO 80 LL=1,NTEM
      L1=LL+1
      CALL VORSEG(XF,YF,ZF,XTIP(N,LL),YTIP(N,LL),ZTIP(N,LL),XTIP(N,L1),
     *     YTIP(N,L1),ZTIP(N,L1),VEL(1),VEL(2),VEL(3),BUG,BUG,BUG,KDUMY)
      DO 80 J=1,3
 80   VSUM(J)=VSUM(J)+VEL(J)
      DO 60 J=1,3
         VFP(J)=VFP(J)+VSUM(J)*GB(I)
 60   CONTINUE
 50   CONTINUE
      IF(ITWAKE.NE.1) GO TO 90
C-----------------------------------------------------------------------
C     Transition wake
C-----------------------------------------------------------------------
      NW1=NW-1
      DO 100 K=1,NBLADE
      MGP=MGRID(K)
      DO 100 M=1,MGP
      DO 110 J=1,3
 110  VSUM(J)=0.0
      NSW1=NSW(K,M)-1
      DO 120 N=1,NSW1
      IF(K.NE.1) GO TO 130
      IF(M.NE.MW) GO TO 130
      IF(N.EQ.NW1) GO TO 140
      IF(N.EQ.NW) GO TO 140
 130  I=INDXTW(N,M,K)
      CALL VORSEG(XF,YF,ZF,XW(I),YW(I),ZW(I),XW(I+1),YW(I+1),ZW(I+1),
     *            VEL(1),VEL(2),VEL(3),BUG,BUG,BUG,KDUMY)
      DO 120 J=1,3
         VSUM(J)=VSUM(J)+VEL(J)
 120  CONTINUE
 140  CONTINUE
      DO 100 J=1,3
      VFP(J)=VFP(J)+VSUM(J)*GAMTW(K,M)
 100  CONTINUE
 90   CONTINUE
      IF(IUWAKE.NE.1) GO TO 150
C-----------------------------------------------------------------------
C     Ultimate wake
C-----------------------------------------------------------------------
      CALL VORSEG(XF,YF,ZF,X1HV,Y1HV,Z1HV,X2HV,Y2HV,Z2HV,
     *            VEL(1),VEL(2),VEL(3),BUG,BUG,BUG,KDUMY)
      DO 160 J=1,3
      VFP(J)=VFP(J)-NBLADE*GAMULT*VEL(J)
 160  VSUM(J)=0.0
      DO 170 K=1,NBLADE
      DO 170 N=2,NPUW
      CALL VORSEG(XF,YF,ZF,XU(N-1,K),YU(N-1,K),ZU(N-1,K),XU(N,K),
     *     YU(N,K),ZU(N,K),VEL(1),VEL(2),VEL(3),BUG,BUG,BUG,KDUMY)
      DO 170 J=1,3
 170  VSUM(J)=VSUM(J)+VEL(J)
      DO 180 J=1,3
      VFP(J)=VFP(J)+GAMULT*VSUM(J)
 180  CONTINUE
 150  CONTINUE
      IF(YF .eq. 0.0) THEN
         THET = 0.0
      ELSE
         THET=ATAN2(ZF,YF)
      ENDIF
      CTHET=COS(THET)
      STHET=SIN(THET)
      VIX=VFP(1)
      VIY=VFP(2)
      VIZ=VFP(3)
      VIR=VIY*CTHET+VIZ*STHET
      VIT=-VIY*STHET+VIZ*CTHET
      RETURN
C))))))))))))))))))) End of subroutiine FPWAKE (((((((((((((((((((((((
      END



      SUBROUTINE UPITCH(NBLADE,RULT,GAMULT,PERULT,ADVCO,VAVS,VTVS,
     *                  TANBUW,UAVS,UTVS)
C***********************************************************************
C     UPITCH: Ulitimate PITCH
C      --- Revised version of LOUKAKIS tip vortex self-induced velocities
C          calculates ultimate wake pitch from LOUKAKIS theory.
C          tip vortex core radius=0.02 * wake radius
C          this version good for 1 to 12 blades,
C          tan(beta(w)) between 0.10 and 1.00
C
      DIMENSION UAS(98),UTS(98),UA(14),UT(14),TANB(14),
     *          UAC(14),UTC(14),UAST(70),UTST(70)
      DATA PI/3.141593/
      DATA TANB/.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9,1.0/
      DATA UAS/
     *  0.9762,0.7390,0.6212,0.5484,0.4962,0.4547,0.4195,
     *  0.3882,0.3597,0.3091,0.2655,0.2278,0.1955,0.1679,
     *  1.7155,1.2151,0.9663,0.8157,0.7129,0.6363,0.5756,
     *  0.5252,0.4820,0.4105,0.3526,0.3046,0.2645,0.2308,
     *  2.4785,1.7138,1.3333,1.1043,0.9496,0.8364,0.7488,
     *  0.6778,0.6185,0.5234,0.4493,0.3896,0.3405,0.2997,
     *  3.2445,2.2210,1.7093,1.4015,1.1946,1.0445,0.9294,
     *  0.8373,0.7614,0.6418,0.5507,0.4785,0.4198,0.3714,
     *  4.0244,2.7341,2.0902,1.7035,1.4441,1.2568,1.1140,
     *  1.0007,0.9079,0.7633,0.6548,0.5696,0.5011,0.4449,
     *  4.8164,3.2528,2.4750,2.0086,1.6966,1.4719,1.3103,
     *  1.1665,1.0566,0.8868,0.7605,0.6622,0.5836,0.5194,
     *  5.5837,3.7668,2.8597,2.3153,1.9509,1.6889,1.4904,
     *  1.3340,1.2070,1.0117,0.8675,0.7559,0.6670,0.5946/
      DATA UAST/
     *  6.3455,4.2815,3.2458,2.6235,2.2067,1.9072,1.6808,
     *  1.5027,1.3585,1.1376,0.9753,0.8503,0.7511,0.6704,
     *  7.1647,4.8103,3.6360,2.9329,2.4636,2.1266,1.8722,
     *  1.6724,1.5109,1.2643,1.0838,0.9453,0.8356,0.7467,
     *  7.9300,5.3260,4.0246,3.2434,2.7214,2.3470,2.0644,
     *  1.8429,1.6640,1.3916,1.1929,1.0408,0.9206,0.8233,
     *  8.7330,5.8572,4.4180,3.5551,2.9796,2.5677,2.2576,
     *  2.0138,1.8176,1.5194,1.3023,1.1366,1.0059,0.9002,
     *  9.5374,6.3779,4.8039,3.8631,3.2365,2.7787,2.4498,
     *  2.1851,1.9717,1.6477,1.4122,1.2328,1.0914,0.9773/
      DATA UTS/
     *  0.0455,0.0244,0.0035,-.0166,-.0355,-.0526,-.0678,
     *  -.0809,-.0919,-.1082,-.1180,-.1225,-.1233,-.1213,
     *  0.1302,0.1121,0.0936,0.0756,0.0586,0.0430,0.0289,
     *  0.0166,0.0060,-.0099,-.0198,-.0249,-.0263,-.0251,
     *  0.2130,0.1962,0.1792,0.1625,0.1467,0.1320,0.1187,
     *  0.1070,0.0969,0.0814,0.0716,0.0663,0.0645,0.0652,
     *  0.2935,0.2788,0.2630,0.2473,0.2322,0.2182,0.2056,
     *  0.1943,0.1846,0.1695,0.1597,0.1543,0.1522,0.1525,
     *  0.3750,0.3611,0.3460,0.3309,0.3164,0.3030,0.2908,
     *  0.2799,0.2704,0.2557,0.2460,0.2405,0.2381,0.2382,
     *  0.4431,0.4419,0.4287,0.4140,0.3999,0.3869,0.3750,
     *  0.3644,0.3552,0.3407,0.3310,0.3255,0.3230,0.3228,
     *  0.4803,0.5079,0.5055,0.4947,0.4821,0.4698,0.4584,
     *  0.4481,0.4391,0.4248,0.4153,0.4096,0.4070,0.4066/
      DATA UTST/
     *  0.6555,0.6167,0.5954,0.5793,0.5654,0.5528,0.5114,
     *  0.5313,0.5224,0.5084,0.4989,0.4932,0.4905,0.4899,
     *  0.6989,0.6822,0.6706,0.6586,0.6464,0.6347,0.6238,
     *  0.6140,0.6053,0.5915,0.5820,0.5763,0.5735,9.5727,
     *  0.7783,0.7654,0.7531,0.7408,0.7287,0.7171,0.7063,
     *  0.6965,0.6879,0.6742,0.6648,0.6590,0.6561,0.6552,
     *  0.8277,0.8418,0.8368,0.8250,0.8120,0.7996,0.7885,
     *  0.7786,0.7701,0.7565,0.7472,0.7414,0.7384,0.7375,
     *  1.0279,0.9718,0.9339,0.9100,0.8933,0.8802,0.8693,
     *  0.8598,0.8516,0.8385,0.8294,0.8236,0.8205,0.8195/
      ICOUNT=0
      NB=NBLADE
      IF(NB.LE.12) GO TO 10
      WRITE(6,20)
 20   FORMAT('  TOO MANY BLADES FOR UPITCH TABLES. ',/2X,
     *       'K=12 ASSUMED FOR ULTIMATE WAKE CALCULATION.')
      NB=12
 10   IF(NB.LE.7) THEN
      IK=(NB-1)*14
      DO 30 I=1,14
      UA(I)=UAS(I+IK)
 30   UT(I)=UTS(I+IK)
      ELSE
      IK=(NB-8)*14
      DO 40 I=1,14
      UA(I)=UAST(I+IK)
 40   UT(I)=UTST(I+IK)
      END IF
      CALL THEILC(TANB,UA,UAC,14)
      CALL THEILC(TANB,UT,UTC,14)
      GAMTIP=GAMULT*PERULT
      UTH=-NBLADE*GAMTIP/RULT
      FACT=2.0*PI*GAMTIP/RULT
      VROT=PI*RULT/ADVCO+VTVS
      TBOLD=0.15
 50   ICOUNT=ICOUNT+1
      CALL THEILP(TANB,UAC,14,TBOLD,UAD)
      CALL THEILP(TANB,UTC,14,TBOLD,UTD)
      UAVS=FACT*UAD
      UTVS=FACT*UTD+UTH
      TBNEW=(VAVS+UAVS)/(VROT+UTVS)
      ERR=ABS(TBOLD-TBNEW)/TBOLD
      IF(ERR.LT.0.002) GO TO 60
      IF(ICOUNT.GT.400) GO TO 70
      TBOLD=(TBNEW+TBOLD)/2.
      GO TO 50
 60   TANBUW=TBNEW
      RETURN
 70   WRITE(6,80)TBNEW
 80   FORMAT('  CONVERGENCE FAILURE IN UPITCH. TANBW=',F8.4)
      STOP
C))))))))))))))))))) End of subroutiine UPITCH (((((((((((((((((((((((
      END



      SUBROUTINE SELFIN(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,CORAD,VX,VY,VZ,VR,VT)
C***********************************************************************
C     SELFIN: SELF-INduced velocities
C      --- Calculates self-induced velocities of circular vortex segment
C          at point x2,y2,z2
C          equations from NASA CR-1911, S. GENE SADLER, 1971
C
      DATA PI/3.141593/
      IF(CORAD.LE.1.0E-10) GO TO 100
      EL1=VECDIS(X1,Y1,Z1,X2,Y2,Z2)
      EL2=VECDIS(X2,Y2,Z2,X3,Y3,Z3)
      EL1S=EL1*EL1
      EL2S=EL2*EL2
      AMX=(Y1-Y2)*(Z2-Z3)-(Z1-Z2)*(Y2-Y3)
      AMY=(Z1-Z2)*(X2-X3)-(X1-X2)*(Z2-Z3)
      AMZ=(X1-X2)*(Y2-Y3)-(Y1-Y2)*(X2-X3)
      AMXYZ=AMX**2+AMY**2+AMZ**2
      IF(AMXYZ.LE.1.0E-10) GO TO 100
      DEL=VECDIS(X1,Y1,Z1,X3,Y3,Z3)
      DELSQ=DEL*DEL
      DENOM=4.*EL1S*EL2S-(EL1S+EL2S-DELSQ)**2
      IF(DENOM.LE.0.0) GO TO 100
      DENOM=SQRT(DENOM)
      S=EL1*EL2*DEL/DENOM
      DF=SQRT(ABS(4.*S**2-EL1S))
      DG=SQRT(ABS(4.*S**2-EL2S))
      IF(EL1S.GT.(DELSQ+EL2S)) GO TO 10
      F=(2.*S-DF)/EL1
      GO TO 20
 10   F=(2.*S+DF)/EL1
 20   IF(EL2S.GT.(DELSQ+EL1S)) GO TO 30
      G=(2.*S-DG)/EL2
      GO TO 40
 30   G=(2.*S+DG)/EL2
 40   E1=8.*S*F/CORAD
      E2=8.*S*G/CORAD
      IF(E1.LE.0.0) GO TO 100
      IF(E2.LE.0.0) GO TO 100
      B=8.*PI*S*SQRT(AMXYZ)
      F=2.*PI*(ALOG(E1)+ALOG(E2)+0.5)/B
      VX=AMX*F
      VY=AMY*F
      VZ=AMZ*F
      THET=ATAN2(Z2,Y2)
      CTHET=COS(THET)
      STHET=SIN(THET)
      VR=VY*CTHET+VZ*STHET
      VT=-VY*STHET+VZ*CTHET
      RETURN
 100  VX=0.0
      VY=0.0
      VZ=0.0
      VR=0.0
      VZ=0.0
      RETURN
C))))))))))))))))))) End of subroutiine SELFIN (((((((((((((((((((((((
      END



      SUBROUTINE SPACE2(XV,XC,FC,CPDC,NT,KBLADE,KSPACE)
C***********************************************************************
C     SPACE2: SPACing
C      --- Subroutine to calculate vortex & control point locations
C
      DIMENSION XV(*),XC(*),FC(*),CPDC(*)
      DATA PI/3.141593/
      IF(KSPACE.NE.1) GO TO 10
C-----------------------------------------------------------------------
C     Super cosine spacing
C-----------------------------------------------------------------------
      DEL=PI/(2.0*NT)
      DO 20 N=1,NT
      THETA=(2.0*N-1.0)*DEL
      XV(N)=0.5*(1.0-COS(THETA))
      IF(KBLADE.GT.1) GO TO 20
      XC(N)=0.5*(1.0-COS(2.0*N*DEL))
      CPDC(N)=DEL*SIN(THETA)
 20   CONTINUE
      GO TO 30
 10   IF(KSPACE.NE.2) GO TO 40
C-----------------------------------------------------------------------
C     Uniform spacing
C-----------------------------------------------------------------------
      DEL=1.0/NT
      DO 50 N=1,NT
      XV(N)=(N-0.75)*DEL
      IF(KBLADE.GT.1) GO TO 50
      XC(N)=(N-0.25)*DEL
      CPDC(N)=DEL
 50   CONTINUE
 30   XV(NT+1)=1.0
      IF(KBLADE.GT.1) GO TO 60
      DO 70 N=1,NT
 70   FC(N)=(XC(N)-XV(N))/(XV(N+1)-XV(N))
 60   RETURN
 40   WRITE(6,80)KSPACE
 80   FORMAT(' KSPACE NOT VALID',I4)
      STOP
C))))))))))))))))))) End of subroutiine SPACE2 (((((((((((((((((((((((
      END



      SUBROUTINE LLFOR
C***********************************************************************
C     LLFOR: Lifting Line FORce
C      --- Estimate forces from lifting line theory 
C          this approach breaks down away from design j
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION XLL(NGPMAX),DTDX(NGPMAX),DQDX(NGPMAX),CLIFT(NGPMAX),
     *          DUM(NGPMAX)
      DATA PI/3.141593/
      MCP=MPAN(1)
      NCP=NPAN(1)
      DO 10 M=1,MCP
 10   XLL(M+1)=R(M)
 20   DO 30 M=1,MCP
      DUM(M)=0.0
      DO 30 N=1,NCP
      I=INDXCP(N,M,1)
 30   DUM(M)=DUM(M)+GB(I)
      DO 40 I=1,MCP
      CL=4.0*PI*DUM(I)/(VSTAR(I)*CHORCP(I))
      CLIFT(I)=CL
      PHI1=ATAN(PITCH(I)/(2.0*PI*RZ(I)))
      PHI2=ATAN(PITCH(I+1)/(2.0*PI*RZ(I+1)))
      PHI=0.5*(PHI1+PHI2)
      SPHI=SIN(PHI)
      CPHI=COS(PHI)
      FACT=0.5*RHO*CHORCP(I)*(VSTAR(I)*VSFPS*RADIUS)**2
      DTDX(I+1)=FACT*(CL*CPHI-CDRAG*SPHI)
 40   DQDX(I+1)=R(I)*FACT*RADIUS*(CL*SPHI+CDRAG*CPHI)
      DTDX(1)=0.
      DQDX(1)=0.
      MC2=MCP+2
      DTDX(MC2)=0.
      DQDX(MC2)=0.
      XLL(1)=RHUB
      XLL(MC2)=1.0
      THRUST=NBLADE*SIMPUN(XLL,DTDX,MC2)
      TORQUE=NBLADE*SIMPUN(XLL,DQDX,MC2)
      FACT=RHO*RPS**2*DIAM**4
      AKT=THRUST/FACT
      AKQ=TORQUE/(FACT*DIAM)
      ETA=ADVCOM*AKT/(AKQ*2.0*PI)
      CT=999.9999
      IF(ABS(ADVCOM).GT.0.0) CT=(8.*AKT)/(PI*ADVCOM**2)
      HORSEP=TORQUE*OMEGA/550.

      RETURN
C))))))))))))))))))) End of subroutiine LLFOR ((((((((((((((((((((((((
      END


      SUBROUTINE FACTOR(A,W,IPIVOT,D,N,NDIM,IFLAG)
C***********************************************************************
C     FACTOR:FACTORing
C      --- Factors square matrices. subroutine SUBST used for
C          back substituition to solve system of equations
C
      DIMENSION A(NDIM,NDIM),W(NDIM,NDIM),IPIVOT(*),D(*)
      IFLAG=1
      DO 10 I=1,N
      IPIVOT(I)=I
      ROWMAX=0.
      DO 9 J=1,N
      W(I,J)=A(I,J)
 9    ROWMAX=AMAX1(ROWMAX,ABS(W(I,J)))
      IF(ROWMAX.EQ.0.) GO TO 999
 10   D(I)=ROWMAX
      NM1=N-1
      IF(NM1.EQ.0) RETURN
      DO 20 K=1,NM1
      J=K
      KP1=K+1
      IP=IPIVOT(K)
      COLMAX=ABS(W(IP,K))/D(IP)
      DO 11 I=KP1,N
      IP=IPIVOT(I)
      AWIKOV=ABS(W(IP,K))/D(IP)
      IF(AWIKOV.LE.COLMAX) GO TO 11
      COLMAX=AWIKOV
      J=I
 11   CONTINUE
      IF(COLMAX.EQ.0.) GO TO 999
      IPK=IPIVOT(J)
      IPIVOT(J)=IPIVOT(K)
      IPIVOT(K)=IPK
      DO 20 I=KP1,N
      IP=IPIVOT(I)
      W(IP,K)=W(IP,K)/W(IPK,K)
      RATIO=-W(IP,K)
      DO 20 J=KP1,N
 20   W(IP,J)=RATIO*W(IPK,J)+W(IP,J)
      IF(W(IP,N).EQ.0.) GO TO 999
      RETURN
 999  IFLAG=2
      RETURN
C)))))))))))))))))))End of subroutiine FACTOR (((((((((((((((((((((((
      END



      SUBROUTINE SUBST(W,B,X,IPIVOT,N,NDIM)
C***********************************************************************
C     SUBST: SUBSTitution
C      --- Subroutine to perform back substitution after calling
C          subroutine FACTOR
C
      DIMENSION W(NDIM,NDIM),B(*),X(*),IPIVOT(*)
      IF(N.GT.1) GO TO 10
      X(1)=B(1)/W(1,1)
      RETURN
 10   IP=IPIVOT(1)
      X(1)=B(IP)
      DO 15 K=2,N
      IP=IPIVOT(K)
      KM1=K-1
      SUM=0.
      DO 14 J=1,KM1
 14   SUM=W(IP,J)*X(J)+SUM
 15   X(K)=B(IP)-SUM
      X(N)=X(N)/W(IP,N)
      K=N
      DO 20 NP1MK=2,N
      KP1=K
      K=K-1
      IP=IPIVOT(K)
      SUM=0.0
      DO 19 J=KP1,N
 19   SUM=W(IP,J)*X(J)+SUM
 20   X(K)=(X(K)-SUM)/W(IP,K)
      RETURN
C)))))))))))))))))))End of subroutiine SUBST ((((((((((((((((((((((((
      END



      FUNCTION INDXBG(N,M,K)
C***********************************************************************
C     INDXBG: INDeX for Blade Grid
C      --- Computes linear index for grid points on all blades
C          N = chordwise index
C          M = spanwise index
C          K = blade index
C
      INCLUDE 'WAKE.INC'
      INDXBG=NGPSUM(K)+(M-1)*NGRID(K)+N
      RETURN
C))))))))))))))))))) End of subroutiine INDXBG ((((((((((((((((((((((((
      END



      FUNCTION INDXCP(N,M,K)
C***********************************************************************
C     INDXCP: INDeX for Control Point
C      --- Computes linear index for control points
C          N = chordwise index
C          M = spanwise index
C          K = blade index
C
      INCLUDE 'WAKE.INC'
      INDXCP=NCPSUM(K)+(M-1)*NPAN(K)+N
      RETURN
C))))))))))))))))))) End of subroutiine INDXCP ((((((((((((((((((((((((
      END



      FUNCTION INDXCV(N,M,K)
C***********************************************************************
C     INDXCV: INDeX for Chordwise Vortices
C      --- Computes linear index for chordwise vortices
C          N = chordwise index
C          M = spanwise index
C          K = blade index
C
      INCLUDE 'WAKE.INC'
      INDXCV=NCVSUM(K)+(M-1)*NPAN(K)+N
      RETURN
C))))))))))))))))))) End of subroutiine INDXCV ((((((((((((((((((((((((
      END



      FUNCTION INDXTW(N,M,K)
C***********************************************************************
C     INDXTW: INDeX for Transition Wake element
C      --- Computes linear index for transition wake element
C
      INCLUDE 'WAKE.INC'
      INDXTW=NTWSUM(K)+(M-1)*80+N
      RETURN
C))))))))))))))))))) End of subroutiine INDXTW ((((((((((((((((((((((((
      END




      SUBROUTINE JUMPB(NO,MO,JUMPN,JUMPM)
C***********************************************************************
C     JUMPB: JUMP on other Blades
C      --- Calculates jump parameters for other blades
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      BUG=1./JUMPM
      DO 10 MJ=1,JUMPM
      MKB=(MO-1)*JUMPM+MJ
      DO 10 NJ=1,JUMPN
      NKB=(NO-1)*JUMPN+NJ
      JJ=(MJ-1)*JUMPN+NJ
      LKB(JJ)=INDXCP(NKB,MKB,1)
 10   DIST(JJ)=BUG
      RETURN
C))))))))))))))))))) End of subroutiine JUMPB ((((((((((((((((((((((((
      END



      SUBROUTINE JUMPW(MO,JUMPM)
C***********************************************************************
C     JUMPW: JUMP on other transition Wakes
C      --- Calculates jump parameters for other transition wakes
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      BUG=1./JUMPM
      NCP=NPAN(1)
      DO 10 MJ=1,JUMPM
      MKB=(MO-1)*JUMPM+MJ
      DO 10 N=1,NCP
      JJ=(MJ-1)*NCP+N
      LKB(JJ)=INDXCP(N,MKB,1)
 10   DIST(JJ)=BUG
      RETURN
C))))))))))))))))))) End of subroutiine JUMPW ((((((((((((((((((((((((
      END



      SUBROUTINE VECXPR(A1,A2,A3,B1,B2,B3,C1,C2,C3)
C***********************************************************************
C     VECXPR: VECtor cross PRoduct 
C      --- Calculates vector cross product A X B = C
C
      C1=A2*B3-A3*B2
      C2=A3*B1-A1*B3
      C3=A1*B2-A2*B1
      RETURN
C))))))))))))))))))) End of subroutiine VECXPR (((((((((((((((((((((((
      END



      FUNCTION VECDPR(A1,A2,A3,B1,B2,B3)
C***********************************************************************
C     VECDPR: VECtor Dot PRoduct 
C      --- Calculates vector dot product
C
      VECDPR=A1*B1+A2*B2+A3*B3
      RETURN
C))))))))))))))))))) End of subroutiine VECDPR (((((((((((((((((((((((
      END



      SUBROUTINE VECNML(A1,A2,A3)
C***********************************************************************
C     VECNML: VECtor to uNit Length
C      --- Normalizes vector to unit length
C
      DATA EPS/1.0E-13/
      C=SQRT(A1*A1+A2*A2+A3*A3)
      IF(C.LT.EPS) GO TO 10
      A1=A1/C
      A2=A2/C
      A3=A3/C
      RETURN
 10   A1=0.0
      A2=0.0
      A3=0.0
 20   FORMAT('VECTOR SET TO ZERO IN VECNML')
      WRITE(6,20)
      RETURN
C))))))))))))))))))) End of subroutiine VECNML (((((((((((((((((((((((
      END



      FUNCTION VECDIS(A1,A2,A3,B1,B2,B3)
C***********************************************************************
C     VECDIS: VECtor DIStance
C      --- Calculates distance between 2 points
C
      VECDIS=SQRT((A1-B1)**2+(A2-B2)**2+(A3-B3)**2)
      RETURN
C)))))))))))))))))))))) End of function VECDIS (((((((((((((((((((((((
      END


      SUBROUTINE THEILC(X,Y,C,N)
C***********************************************************************
C     THEILC: THEILheimer Coefficients
C      --- Computes theilheimer spline coefficients
C
      DIMENSION X(*),Y(*),C(*),B(44000),IPIV(201)
      DO 10 I=1,N
      DO 20 J=1,4
      K=(J-1)*N+I
      IF(J.EQ.1) GO TO 30
      B(K)=X(I)**(J-1)
      GO TO 20
 30   B(K)=1.0
 20   CONTINUE
      DO 10 J=5,N
      K=(J-1)*N+I
      TEST=X(I)-X(J-2)
      IF(TEST-0.00001) 40,40,50
 40   B(K)=0.0
      GO TO 10
 50   B(K)=TEST**3
 10   CONTINUE
      CALL FACTOR(B,B,IPIV,C,N,N,IFLAG)
      IF(IFLAG.EQ.1) GO TO 60
      WRITE(6,70) IFLAG
      STOP
 70   FORMAT(' THEILHEIMER MATRIX BAD',I5)
 60   CALL SUBST(B,Y,C,IPIV,N,N)
      RETURN
C))))))))))))))))))) End of subroutiine INFMAT ((((((((((((((((((((((( 
      END



      SUBROUTINE THEILP(X,C,N,XX,Y)
C***********************************************************************
C     THEILP: THEILheimer Polynomial spline
C      --- Evaluates theilheimer polynomial spline
C
      DIMENSION X(*),C(*)
      Y=C(1)+((C(4)*XX+C(3))*XX+C(2))*XX
      DO 10 I=5,N
      TEST=XX-X(I-2)
      IF(TEST-0.00001) 10,10,20
 20   Y=Y+C(I)*TEST**3
 10   CONTINUE
      RETURN
C))))))))))))))))))) End of subroutiine THEILP ((((((((((((((((((((((( 
      END



      SUBROUTINE THEILD(X,C,N,XX,YPRIME)
C***********************************************************************
C     THEILD: Derivative of THEILheimer 
C      --- Evaluates derivative of theilheimer polynomial
C
      DIMENSION X(*),C(*)
      YPRIME=C(2)+2.*C(3)*XX+3.*C(4)*XX**2
      DO 10 I=5,N
      TEST=XX-X(I-2)
      IF(TEST-0.00001) 10,10,20
 20   YPRIME=YPRIME+3.*C(I)*(XX-X(I-2))**2
 10   CONTINUE
      RETURN
C))))))))))))))))))) End of subroutiine THEILD ((((((((((((((((((((((( 
      END



      SUBROUTINE SPLINE(X,Y,NA,XX,YY,NB)
C***********************************************************************
C     SPLINE: SPLINE interpolation
C      --- Theilheimer spline interpolation
C 
      DIMENSION X(*),Y(*),XX(*),YY(*),C(201)

      CALL THEILC(X,Y,C,NA)
      DO 20 J=1,NB
      ABS=XX(J)
      CALL THEILP(X,C,NA,ABS,ORD)
20    YY(J)=ORD

      RETURN
      END

      SUBROUTINE SPLINEs(X,Y,NA,XX,YY,NB)
C***********************************************************************
C     SPLINE: SPLINE interpolation
C      --- Theilheimer spline interpolation
C 
      DIMENSION X(*),Y(*),C(201)
      real XX, YY

      CALL THEILC(X,Y,C,NA)
      DO 20 J=1,NB
      ABS=XX
      CALL THEILP(X,C,NA,ABS,ORD)
20    YY=ORD

      RETURN
      END


      SUBROUTINE MAP1(X,XFACT,IPIVX,NA)
C***********************************************************************
C     MAP1: MAPping
C      --- Replaces repetitive use of subroutine spline on
C          same base points
C
      DIMENSION X(*),XFACT(*),IPIVX(*),D(200)
      DO 10 I=1,NA
      DO 20 J=1,4
      K=(J-1)*NA+I
      IF(J.EQ.1) GO TO 30
      XFACT(K)=X(I)**(J-1)
      GO TO 20
 30   XFACT(K)=1.0
 20   CONTINUE
      DO 10 J=5,NA
      K=(J-1)*NA+I
      TEST=X(I)-X(J-2)
      IF(TEST-0.00001) 40,40,50
 40   XFACT(K)=0.0
      GO TO 10
 50   XFACT(K)=TEST**3
 10   CONTINUE
      CALL FACTOR(XFACT,XFACT,IPIVX,D,NA,NA,IFLAG)
      IF(IFLAG.LE.1) GO TO 60
 70   FORMAT(2X,'MATRIX BAD IN MAP1')
      WRITE(6,70)
      STOP
 60   RETURN
C))))))))))))))))))))) End of subroutiine MAP1 ((((((((((((((((((((((( 
      END



      SUBROUTINE MAP2(X,Y,NA,XX,YY,NB,XFACT,IPIVX)
C***********************************************************************
C     MAP2: MAPping
C      --- Replaces repetitive calls to spline using
C          same base points. MAP1 must be called before MAP2.
C
      DIMENSION X(*),Y(*),XX(*),YY(*),XFACT(*),IPIVX(*),C(200)
      CALL SUBST(XFACT,Y,C,IPIVX,NA,NA)
      DO 10 I=1,NB
 10   CALL THEILP(X,C,NA,XX(I),YY(I))
      RETURN
C))))))))))))))))))))) End of subroutiine MAP2 ((((((((((((((((((((((( 
      END



      SUBROUTINE MAPTBL(NIN,NOUT,NRD,XIN,XOUT,A)
C***********************************************************************
C     MAPTBL: MAtrix
C      --- Generates interpolation matrix
C
!     DIMENSION XIN(*),XOUT(*),A(NRD,1),XFACT(225),YD(40),IPIV(40)
      DIMENSION XIN(*),XOUT(*),A(NRD,NIN),XFACT(225),YD(40),IPIV(40)
      CALL MAP1(XIN,XFACT,IPIV,NIN)
      DO 10 N=1,NIN
      DO 20 L=1,NIN
 20   YD(L)=0.0
      YD(N)=1.0
 10   CALL MAP2(XIN,YD,NIN,XOUT,A(1,N),NOUT,XFACT,IPIV)
      RETURN
C))))))))))))))))))) End of subroutiine MAPTBL ((((((((((((((((((((((( 
      END


      FUNCTION SIMPUN(X,Y,N)
C***********************************************************************
C     SIMPUN: SIMPson 
C
      DIMENSION X(*),Y(*)
      IF(N.LT.2) THEN
         SIMPUN=0.0
         RETURN
      ELSE IF(N.EQ.2) THEN
         SIMPUN=(Y(1)+Y(2))*(X(2)-X(1))/2.
         RETURN          
      END IF
      M=N-1
 10   IF(M-2)20,30,40
 40   M=M-2
      GO TO 10
 20   S=(X(2)-X(1))/6.*(Y(1)*(3.-(X(2)-X(1))/(X(3)-X(1)))+Y(2)*
     *  (3.+(X(2)-X(1))/(X(3)-X(2)))-Y(3)*(((X(2)-X(1))**2)/
     *  ((X(3)-X(1))*(X(3)-X(2)))))
      L=3
      GO TO 50
 30   S=0.
      L=2
 50   M=N-1
      DO 60 K=L,M,2
         IF(ABS(X(K-1)-X(1)).GE.ABS(X(K)-X(1))) GO TO 70
         IF(ABS(X(K+1)-X(1)).GT.ABS(X(K)-X(1))) GO TO 80
 70      WRITE(6,90)K,X(K)
 90      FORMAT(1X,'NON MONOTONE X IN SIMPUN',I4,1PE12.4)
         S=0.0
         GO TO 100
 80      S=S+(X(K+1)-X(K-1))/6.*(Y(K-1)*(3.-(X(K+1)-X(K-1))/
     *     (X(K)-X(K-1)))+(Y(K)*(1.+(X(K+1)-X(K-1))/(X(K)-X(K-1))+
     *     (X(K)-X(K-1))/(X(K+1)-X(K)))+(Y(K+1)*(2.-(X(K)-X(K-1))/
     *     (X(K+1)-X(K))))))
 60   CONTINUE
 100  SIMPUN=S
      RETURN
C)))))))))))))))))))))) End of function SIMPUN ((((((((((((((((((((((( 
      END



      SUBROUTINE CAMFIT(RCAM,XVC,HCC,NP,SLPLE)
C***********************************************************************
C     CAMFIT: CAMber FITting
C      --- Computes camber offsets at desired radius
C
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION XVC(*),HCC(*),CHT(17),DUM(1),C(17)
      CHT(1)=0.0
      CHT(17)=0.0
      DUM(1)=RCAM
C---- GET COORDINATES AT STANDARD CHORDWISE LOCATIONS
      DO 10 J=1,15
 10   CALL MAP2(XR,XCS(1,J),NX,DUM,CHT(J+1),1,XRFAC,IPVXR)
      CALL THEILC(PER,CHT,C,17)
      DO 20 I=1,NP
 20   CALL THEILP(PER,C,17,XVC(I),HCC(I))
      CALL THEILD(PER,C,17,0.0,SLPLE)
      RETURN
C))))))))))))))))))) End of subroutiine CAMFIT ((((((((((((((((((((((( 
      END



C ======================================================================
      SUBROUTINE THKFIT(RTHK,XCT,HTT,NT)
C
C     THKFIT: THicKness FITting
C      --- Computes thickness values and thickness differences at
C          desired locations
C
C ======================================================================
      INCLUDE 'WAKE.INC'
      INCLUDE 'INPWAK.INC'
      DIMENSION XCT(*),HTT(*),DUM(1),XCSQ(21),TC(21),THT(17)
 
      THT(1)=0.0
      DUM(1)=RTHK
C---- GET VALUES AT STANDARD CHORDWISE LOCATIONS AT DESIRED RADIUS
      DO 10 J=1,16
 10   CALL MAP2(XR,XTS(1,J),NX,DUM,THT(J+1),1,XRFAC,IPVXR)
C---- INTERPOLATE OVER CHORD
      DO 20 J=1,NT
 20   XCSQ(J)=SQRT(XCT(J))
      CALL SPLINE(PSQ,THT,17,XCSQ,TC,NT)
      HTT(1)=TC(1)
      NTM1=NT-1
      DO 30 J=2,NTM1
 30   HTT(J)=TC(J)-TC(J-1)
      HTT(NT)=-TC(NTM1)
      RETURN
      END



      BLOCK DATA
      INCLUDE 'WAKE.INC'
C---- SELECTION OF INDICES FOR INTERPOLATION MATRICES
C      DATA INT2/4*0,4*2,5*3,3*4,4*5,6/
C      DATA INT3/4*0,2*3,2*4,2*5,6,2*7,3*8,2*9,2*10,11/
C      DATA INT4/4*0,4,5,6,7,7,8,9,2*10,11,2*12,13,14,2*15,16/
c---- CHORD FRACTIONS FOR SECTION INPUT DATA
      DATA PER/0.0,0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
     *  0.95,0.975,0.99,1.0/
      END


C ----------------------------------------
        SUBROUTINE LSQ(X,Y,N,M)
C ----------------------------------------
C -  X(N) : X points
C    Y(N) : Y Points
C    N    : Number of Data points
C    M    : Degree of polynimial to  be fitted.
C ----------------------------------------------
        dimension a(M+1,M+1),b(M+1),c(15,M+1)
        dimension x(N),y(N)
        dimension ipvt(M+1),work(M+1)

        DO I = 1 , N
          C(I,1) = 1.0
        ENDDO

        MP1 = M + 1
        DO J =  2, MP1
          DO I = 1 , N
             C(I,J) = C(I,J-1) * X(I)
          ENDDO
        ENDDO

        DO I = 1 , MP1
          DO J = 1 , MP1
            A(I,J) = 0.0
            DO K = 1 , N
               A(I,J) = A(I,J) + C(K,I)*C(K,J)
            ENDDO
          ENDDO
        ENDDO

        DO I = 1 , MP1
          B(I) = 0.0
          DO K = 1 , N
             B(I) = B(I) +C(K,I)*Y(K)
          ENDDO
        ENDDO

        CALL SDECOMP(MP1,MP1,A,COND,IPVT,WORK)
        CALL SDSOLVE(MP1,MP1,A,B,IPVT)


        DO I = 1 , N
          Y(I) = 0.0
          DO J = 1 , MP1
            Y(I) = Y(I) + B(J)*X(I)**(J-1)
          ENDDO
        ENDDO

        RETURN
        END


