
C
C     ***********************************************
C     *             Copyright  1989                 *
C     *   Massachusetts Institute of Technology     *
C     *           All Rights Reserved               *
C     *                                             *
C     *  Commercial use and/or reproduction without *
C     *  license prohibited. Licensing agent:       *
C     *                                             *
C     *       MIT Technology Licensing Office       *
C     *               (617) 253-6966                *
C     *                                             *
C     *  Academic and research use unrestricted by  *
C     *  verbal permission from:                    *
C     *                                             *
C     *  Mark Drela   (617) 253-0067                *
C     *                                             *
C     ***********************************************
C
      SUBROUTINE SETBL(XN,YN,DS,FL_ERR)
C................................................
C     Sets up the BL Newton system coefficients
C     for the current BL variables and the edge
C     velocities received from SETUP. The local
C     BL system coefficients are then
C     incorporated into the global Newton system.  
C...............................................
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
c---- wbrewer -- dimension xn,yn,ds -- this should be changes later
      DIMENSION XN(*), YN(*), DS(*)
c----
      REAL*8 USAV(IVX,2)
      REAL*8 U1_M(2*IVX), U2_M(2*IVX)
      REAL*8 D1_M(2*IVX), D2_M(2*IVX)
      REAL*8 ULE1_M(2*IVX), ULE2_M(2*IVX)
      REAL*8 UTE1_M(2*IVX), UTE2_M(2*IVX)


!s-- YE TIAN 06/17/2013----
      LOGICAL FL_ERR

      FL_ERR = .false.
!e-- YE TIAN 06/17/2013----
C
C---- set gas constant (= Cp/Cv)
ccc      GAMBL = GAMMA
ccc      GM1BL = GAMM1
      gambl = 1.4D0
      gm1bl = gambl - 1.D0          ! = 0.4
c
c--- set freestream Mach number
      minf = 0.D0      
C
C---- set parameters for compressibility correction
ccc      QINFBL = QINF
ccc      TKBL = TKLAM
      qinfbl = 1.D0
      tkbl = 0.D0
C
C---- stagnation density and 1/enthalpy
cccc      RSTBL = (1.D0 + 0.5D0*GM1BL*MINF**2) ** (1.D0/GM1BL)
cccc      HSTINV = GM1BL*(MINF/QINFBL)**2 / (1.D0 + 0.5D0*GM1BL*MINF**2)
      rstbl = 1.D0
      HSTINV = GM1BL*(MINF/QINFBL)**2 / (1.D0 + 0.5D0*GM1BL*MINF**2)  != 0
C
C---- Sutherland's const./To   (assumes stagnation conditions are at STP)
      HVRAT = 0.35D0
C
C---- set Reynolds number based on freestream density, velocity, viscosity
      HERAT = 1.D0 - 0.5D0*QINFBL**2*HSTINV
      REY = REINF * DSQRT((HERAT)**3) * (1.D0+HVRAT)/(HERAT+HVRAT) !=reinf
C
      lalfa = .true.
      retyp = 1
C---- set the CL used to define Reynolds number
      IF(LALFA) THEN
       CLREY = CL
      ELSE
       CLREY = CLSPEC
      ENDIF
C
      IF(LALFA .AND. RETYP.GT.1 .AND. CLREY.LT.0.05D0)
     & WRITE(*,*) 'SETBL warning:  CL < 0.05 ... Re(CL) may blow up.'
C
C---- set actual Reynolds number based on CL
      IF(RETYP.EQ.1) REYBL = REY                         !
      IF(RETYP.EQ.2) REYBL = REY/DSQRT(DABS(CLREY))
      IF(RETYP.EQ.3) REYBL = REY/DABS(CLREY)

      REYBL = REY
      AMCRIT = ACRIT
C
C---- save TE thickness
cccc      DWTE = WGAP(1)
      dwte=0.D0

C
      IF(.NOT.LBLINI) THEN
C----- initialize BL by marching with Ue (fudge at separation)
C       WRITE(*,*) ' '
C       WRITE(*,*) 'Initializing BL ...'
       CALL MRCHUE(FL_ERR)
       if (FL_ERR) then
!        write(*,*) 'floating number error.xbf.f:101'
         return
!        stop
       end if
       LBLINI = .TRUE.
      ENDIF
C
C      WRITE(*,*)
C
C---- march BL with current Ue and Ds to establish transition
c---- wbrewer -- 4/10/95 -- temporarily comment out MRCHDU
       if(ijump.ne.1)  CALL MRCHDU
C
      DO 5 IS=1, 2
        DO 6 IBL=2, NBL(IS)
          USAV(IBL,IS) = UEDG(IBL,IS)
    6   CONTINUE
    5 CONTINUE
C
      CALL UESET
C
      DO 7 IS=1, 2
c      write(71,*)' ZONE T=" IS =', IS , ' " ' 
        DO 8 IBL=2, NBL(IS)
          TEMP = USAV(IBL,IS)
          USAV(IBL,IS) = UEDG(IBL,IS)
          UEDG(IBL,IS) = TEMP
          DUE = UEDG(IBL,IS) - USAV(IBL,IS)
c      write(71,*)ibl,usav(ibl,is),uedg(ibl,is),due
    8   CONTINUE
    7 CONTINUE
C
      ILE1 = IPAN(2,1)
      ILE2 = IPAN(2,2)
      ITE1 = IPAN(IBLTE(1),1)
      ITE2 = IPAN(IBLTE(2),2)
C
      JVTE1 = ISYS(IBLTE(1),1)
      JVTE2 = ISYS(IBLTE(2),2)
C
      DULE1 = UEDG(2,1) - USAV(2,1)
      DULE2 = UEDG(2,2) - USAV(2,2)
C
C---- set LE and TE Ue sensitivities wrt all m values
      DO 10 JS=1, 2
        DO 110 JBL=2, NBL(JS)
          J  = IPAN(JBL,JS)
          JV = ISYS(JBL,JS)
          ULE1_M(JV) = -VTI(       2,1)*VTI(JBL,JS)*DIJ(ILE1,J)
          ULE2_M(JV) = -VTI(       2,2)*VTI(JBL,JS)*DIJ(ILE2,J)
          UTE1_M(JV) = -VTI(IBLTE(1),1)*VTI(JBL,JS)*DIJ(ITE1,J)
          UTE2_M(JV) = -VTI(IBLTE(2),2)*VTI(JBL,JS)*DIJ(ITE2,J)
  110   CONTINUE
   10 CONTINUE
C
      ULE1_A = UINV_A(2,1)
      ULE2_A = UINV_A(2,2)

C********************************************************************
C**** Go over each boundary layer/wake
      DO 2000 IS=1, 2
C
C---- there is no station "1" at similarity, so zero everything out
      DO 20 JS=1, 2
        DO 210 JBL=2, NBL(JS)
          JV = ISYS(JBL,JS)
          U1_M(JV) = 0.D0
          D1_M(JV) = 0.D0
  210   CONTINUE
   20 CONTINUE
      U1_A = 0.D0
      D1_A = 0.D0
C
      DUE1 = 0.D0
      DDS1 = 0.D0
c---- wbrewer -- initialize Cf
      CFTOTX = 0.D0
      CFTOTY = 0.D0
C
C---- similarity station pressure gradient parameter  x/u du/dx
      IBL = 2
      BULE = 1.D0
C
C---- set forced transition arc length position
      CALL XIFSET(IS)
C
      TRAN = .FALSE.
      TURB = .FALSE.
C
c---- wbrewer -- print Cf value side 2 (last call)
      OPEN(77,FILE='cf.dat',STATUS='unknown')
C**** Sweep downstream setting up BL equation linearizations
      DO 1000 IBL=2, NBL(IS)
C
      IV  = ISYS(IBL,IS)
C
c---- wbrewer -- 4/10/95 -- set the panel to add dtcav: jump in momemtum
c-    thickness; jpan=1 if it is the panel where the jump should occur
      if(ijump.eq.1.and.is.eq.1.and.ibl.eq.icte) then
         jpan=1
      else
         jpan=0
      endif
c
      SIMI = IBL.EQ.2
C ** Shige
      WAKE = IBL.GT.IBLTE(IS)
      ITYPCV = 0
C---- wbrewer -- if partial cavity (ICV = 2) set Cf to zero on cavity      
      IF(ICV.EQ.2) THEN
         IF ((IS.EQ.1).AND.(IBL.GE.ICLE).AND.(IBL.LE.ICTE)) ITYPCV = 1
      ENDIF
C---- wbrewer -- if supercavity (ICV = 3) set Cf to zero on pressure side      
      IF (ICV.EQ.3) THEN
       IF ((IS.EQ.2).AND.(IBL.GE.ICLE).AND.(IBL.LE.ICTE)) ITYPCV=1
       IF (IS.EQ.1) ITYPCV=1
      ENDIF
c---- wbrewer -- plot cf vs. x/c --------
c      write(*,*) ibl,cfa,is
      IF(IS.EQ.2) THEN 
       IF(NBL(IS)-IBL-NWvis+1.GE.1) THEN
        WRITE(77,*) X(NBL(IS)-IBL-NWvis+1)     !,nbl(is)-ibl-nw+1,nw
       ELSE 
        WRITE(77,*) X(NBL(1)+IBL-2),CFA        !,NBL(1)+IBL-2,ibl
       ENDIF
      ENDIF
c
c---- wbrewer -- sum Cf on foil for computing forces for sc foil ----

      IF(ICV.EQ.3.AND.IBL.LT.ICLE) THEN
        CFTOTX = CFTOTX 
     &         - CFA*YN(NBL(IS)-IBL-NWvis+1)*DS(NBL(IS)-IBL-NWvis+1)
        CFTOTY = CFTOTY  
     &         + CFA*XN(NBL(IS)-IBL-NWvis+1)*DS(NBL(IS)-IBL-NWvis+1)
      ENDIF
c      write(77,'(4I5,2f7.4)') NWvis,NBL(IS),IS,IBL,CFTOTX(IS),CFTOTY(IS)
C
      I = IPAN(IBL,IS)
C
C---- set primary variables for current station
      XSI = XSSI(IBL,IS)
      AMI = CTAU(IBL,IS)
      CTI = CTAU(IBL,IS)
      UEI = UEDG(IBL,IS)
      THI = THET(IBL,IS)
      MDI = MASS(IBL,IS)

      if(UEI.lt.0.0) then
!       write(*,*) 'UEI=', UEI, 'xbl.f:244'
!       stop
      end if
C
      DSI = MDI/UEI
C
      IF(WAKE) THEN
       IW = IBL - IBLTE(IS)
       DSWAKI = WGAP(IW)
      ELSE
       DSWAKI = 0.D0
      ENDIF
C
      D2_M2 =  1.D0/UEI
      D2_U2 = -DSI/UEI
C
      DO 30 JS=1, 2
        DO 310 JBL=2, NBL(JS)
          J  = IPAN(JBL,JS)
          JV = ISYS(JBL,JS)
          U2_M(JV) = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
          D2_M(JV) = D2_U2*U2_M(JV)
  310   CONTINUE
   30 CONTINUE
      D2_M(IV) = D2_M(IV) + D2_M2
C
      U2_A = UINV_A(IBL,IS)
      D2_A = D2_U2*U2_A
C
C---- "forced" changes due to mismatch between UEDG and USAV=UINV+dij*MASS
      DUE2 = UEDG(IBL,IS) - USAV(IBL,IS)
      DDS2 = D2_U2*DUE2
c      write(54,*)due2
C
C---- check for transition and set TRAN, XT, etc. if found
      IF((.NOT.SIMI) .AND. (.NOT.TURB)) CALL TRCHEK(IBL,IS)
C
C---- assemble 10x4 linearized system for dCtau, dTh, dDs, dUe, dXi
C     at the previous "1" station and the current "2" station
C
      IF(IBL.EQ.IBLTE(IS)+1) THEN
C
C----- define quantities at start of wake, adding TE base thickness to Dstar
       TTE = THET(IBLTE(1),1) + THET(IBLTE(2),2)
       DTE = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
       CTE = ( CTAU(IBLTE(1),1)*THET(IBLTE(1),1)
     &       + CTAU(IBLTE(2),2)*THET(IBLTE(2),2) ) / TTE
       CALL TESYS(CTE,TTE,DTE)
C
       TTE_TTE1 = 1.D0
       TTE_TTE2 = 1.D0
       DTE_MTE1 =              1.D0 / UEDG(IBLTE(1),1)
       DTE_UTE1 = -DSTR(IBLTE(1),1) / UEDG(IBLTE(1),1)
       DTE_MTE2 =              1.D0 / UEDG(IBLTE(2),2)
       DTE_UTE2 = -DSTR(IBLTE(2),2) / UEDG(IBLTE(2),2)
       CTE_CTE1 = THET(IBLTE(1),1)/TTE
       CTE_CTE2 = THET(IBLTE(2),2)/TTE
       CTE_TTE1 = (CTAU(IBLTE(1),1) - CTE)/TTE
       CTE_TTE2 = (CTAU(IBLTE(2),2) - CTE)/TTE
C
C----- re-define D1 sensitivities wrt m since D1 depends on both TE Ds values
       DO 35 JS=1, 2
         DO 350 JBL=2, NBL(JS)
           J  = IPAN(JBL,JS)
           JV = ISYS(JBL,JS)
           D1_M(JV) = DTE_UTE1*UTE1_M(JV) + DTE_UTE2*UTE2_M(JV)
  350    CONTINUE
   35  CONTINUE
       D1_M(JVTE1) = D1_M(JVTE1) + DTE_MTE1
       D1_M(JVTE2) = D1_M(JVTE2) + DTE_MTE2
C
C----- "forced" changes from  UEDG --- USAV=UINV+dij*MASS  mismatch
       DUE1 = 0.D0
       DDS1 = DTE_UTE1*(UEDG(IBLTE(1),1) - USAV(IBLTE(1),1))
     &      + DTE_UTE2*(UEDG(IBLTE(2),2) - USAV(IBLTE(2),2))
C
      ELSE
C
       CALL BLSYS
C
      ENDIF
C
C---- Save wall shear for plotting output
      TAU(IBL,IS) = 0.5D0*R2*U2*U2*CF2
C
C---- set XI sensitivities wrt LE Ue changes
      IF(IS.EQ.1) THEN
       XI_ULE1 =  SST_GO
       XI_ULE2 = -SST_GP
      ELSE
       XI_ULE1 = -SST_GO
       XI_ULE2 =  SST_GP
      ENDIF
C
C---- stuff BL system coefficients into main Jacobian matrix
C
      DO 40 JV=1, NSYS
        VM(1,JV,IV) = VS1(1,3)*D1_M(JV) + VS1(1,4)*U1_M(JV)
     &              + VS2(1,3)*D2_M(JV) + VS2(1,4)*U2_M(JV)
     &              + (VS1(1,5) + VS2(1,5))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   40 CONTINUE
C
      VB(1,1,IV) = VS1(1,1)
      VB(1,2,IV) = VS1(1,2)
C
      VA(1,1,IV) = VS2(1,1)
      VA(1,2,IV) = VS2(1,2)
C
      IF(LALFA) THEN
       VDEL(1,2,IV) = VSR(1)
      ELSE
       VDEL(1,2,IV) = 
     &       (VS1(1,4)*U1_A + VS1(1,3)*D1_A)
     &     + (VS2(1,4)*U2_A + VS2(1,3)*D2_A)
     &     + (VS1(1,5) + VS2(1,5))*(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF
C
      VDEL(1,1,IV) = VSREZ(1)
     &   + (VS1(1,4)*DUE1 + VS1(1,3)*DDS1)
     &   + (VS2(1,4)*DUE2 + VS2(1,3)*DDS2)
     &   + (VS1(1,5) + VS2(1,5))*(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      DO 50 JV=1, NSYS
        VM(2,JV,IV) = VS1(2,3)*D1_M(JV) + VS1(2,4)*U1_M(JV)
     &              + VS2(2,3)*D2_M(JV) + VS2(2,4)*U2_M(JV)
     &              + (VS1(2,5) + VS2(2,5))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   50 CONTINUE
C
      VB(2,1,IV)  = VS1(2,1)
      VB(2,2,IV)  = VS1(2,2)
C
      VA(2,1,IV) = VS2(2,1)
      VA(2,2,IV) = VS2(2,2)
C
      IF(LALFA) THEN
       VDEL(2,2,IV) = VSR(2)
      ELSE
       VDEL(2,2,IV) = 
     &       (VS1(2,4)*U1_A + VS1(2,3)*D1_A)
     &     + (VS2(2,4)*U2_A + VS2(2,3)*D2_A)
     &     + (VS1(2,5) + VS2(2,5))*(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF
C
      VDEL(2,1,IV) = VSREZ(2)
     &   + (VS1(2,4)*DUE1 + VS1(2,3)*DDS1)
     &   + (VS2(2,4)*DUE2 + VS2(2,3)*DDS2)
     &   + (VS1(2,5) + VS2(2,5))*(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      DO 60 JV=1, NSYS
        VM(3,JV,IV) = VS1(3,3)*D1_M(JV) + VS1(3,4)*U1_M(JV)
     &              + VS2(3,3)*D2_M(JV) + VS2(3,4)*U2_M(JV)
     &              + (VS1(3,5) + VS2(3,5))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   60 CONTINUE
C
      VB(3,1,IV) = VS1(3,1)
      VB(3,2,IV) = VS1(3,2)
C
      VA(3,1,IV) = VS2(3,1)
      VA(3,2,IV) = VS2(3,2)
C
      IF(LALFA) THEN
       VDEL(3,2,IV) = VSR(3)
      ELSE
       VDEL(3,2,IV) = 
     &       (VS1(3,4)*U1_A + VS1(3,3)*D1_A)
     &     + (VS2(3,4)*U2_A + VS2(3,3)*D2_A)
     &     + (VS1(3,5) + VS2(3,5))*(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF

C
      VDEL(3,1,IV) = VSREZ(3)
     &   + (VS1(3,4)*DUE1 + VS1(3,3)*DDS1)
     &   + (VS2(3,4)*DUE2 + VS2(3,3)*DDS2)
     &   + (VS1(3,5) + VS2(3,5))*(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      IF(IBL.EQ.IBLTE(IS)+1) THEN
C
C----- redefine coefficients for TTE, DTE, etc
       VZ(1,1)    = VS1(1,1)*CTE_CTE1
       VZ(1,2)    = VS1(1,1)*CTE_TTE1 + VS1(1,2)*TTE_TTE1
       VB(1,1,IV) = VS1(1,1)*CTE_CTE2
       VB(1,2,IV) = VS1(1,1)*CTE_TTE2 + VS1(1,2)*TTE_TTE2
C
       VZ(2,1)    = VS1(2,1)*CTE_CTE1
       VZ(2,2)    = VS1(2,1)*CTE_TTE1 + VS1(2,2)*TTE_TTE1
       VB(2,1,IV) = VS1(2,1)*CTE_CTE2
       VB(2,2,IV) = VS1(2,1)*CTE_TTE2 + VS1(2,2)*TTE_TTE2
C
       VZ(3,1)    = VS1(3,1)*CTE_CTE1
       VZ(3,2)    = VS1(3,1)*CTE_TTE1 + VS1(3,2)*TTE_TTE1
       VB(3,1,IV) = VS1(3,1)*CTE_CTE2
       VB(3,2,IV) = VS1(3,1)*CTE_TTE2 + VS1(3,2)*TTE_TTE2
C
      ENDIF
C
C---- turbulent intervals will follow if currently at transition interval
      IF(TRAN) TURB = .TRUE.
      TRAN = .FALSE.
C
      IF(IBL.EQ.IBLTE(IS)) THEN
C----- set "2" variables at TE to wake correlations for next station
C
       TURB = .TRUE.
       WAKE = .TRUE.
       CALL BLVAR(3)
      ENDIF
C
      DO 80 JS=1, 2
        DO 810 JBL=2, NBL(JS)
          JV = ISYS(JBL,JS)
          U1_M(JV) = U2_M(JV)
          D1_M(JV) = D2_M(JV)
  810   CONTINUE
   80 CONTINUE
C
      U1_A = U2_A
      D1_A = D2_A
C
      DUE1 = DUE2
      DDS1 = DDS2
C      
C---- set BL variables for next station
      DO 190 ICOM=1, NCOM
        COM1(ICOM) = COM2(ICOM)
  190 CONTINUE
C
C---- next streamwise station
 1000 CONTINUE
c---- wbrewer -- close cf.dat 
      CLOSE(77)
C
      IF(TFORCE(IS)) THEN
       WRITE(*,9100) IS,XTR(IS),ITRAN(IS)
 9100  FORMAT(1X,'Side',I2,' forced transition at x/c = ',F7.4,I5)
      ELSE
       WRITE(*,9200) IS,XTR(IS),ITRAN(IS)
 9200  FORMAT(1X,'Side',I2,'  free  transition at x/c = ',F7.4,I5)
      ENDIF
C
C---- next airfoil side
 2000 CONTINUE
C
      RETURN
      END


      SUBROUTINE IBLSYS
C............................................
C     Sets the BL Newton system line number
C     corresponding to each BL station.
C............................................
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
C
      IV = 0
      DO 10 IS=1, 2
        DO 110 IBL=2, NBL(IS)
          IV = IV+1
          ISYS(IBL,IS) = IV
  110   CONTINUE
   10 CONTINUE
C
      NSYS = IV
      IF(NSYS.GT.2*IVX) STOP '*** IBLSYS: BL system array overflow. ***'
C
      RETURN
      END


!     SUBROUTINE MRCHUE
!s-- YE TIAN 06/17/2013----
!--- added a flag to capture floating number error
!e-- YE TIAN 06/17/2013----
      SUBROUTINE MRCHUE(FL_ERR)
C...................................................
C     Marches the BLs and wake in direct mode using
C     the UEDG array. If separation is encountered,
C     a plausible value of Hk extrapolated from
C     upstream is prescribed instead.  Continuous
C     checking of transition onset is performed.
C...................................................
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      LOGICAL DIRECT
      REAL*8 MSQ         !!!  Bug  22 March 91
!s-- YE TIAN 06/17/2013----
      LOGICAL FL_ERR

      FL_ERR = .false.
!e-- YE TIAN 06/17/2013----
C
C---- shape parameters for separation criteria
      HLMAX = 3.8D0
      HTMAX = 2.5D0
C
      DO 2000 IS=1, 2
C
C      WRITE(*,*) '   side ', IS, ' ...'
C
C---- set forced transition arc length position
      CALL XIFSET(IS)
C
C---- initialize similarity station with Thwaites' formula
      IBL = 2
      XSI = XSSI(IBL,IS)
      UEI = UEDG(IBL,IS)
C      BULE = DLOG(UEDG(IBL+1,IS)/UEI) / DLOG(XSSI(IBL+1,IS)/XSI)
C      BULE = DMAX1( -.08d0 , BULE )
      BULE = 1.D0
      UCON = UEI/XSI**BULE
      TSQ = 0.45D0/(UCON*(5.D0*BULE+1.D0)*REYBL) * XSI**(1.D0-BULE)
      THI = DSQRT(TSQ)
      DSI = 2.2D0*THI
      AMI = 0.D0
C
C---- initialize Ctau for first turbulent station
      CTI = 0.03D0
C
      TRAN = .FALSE.
      TURB = .FALSE.
      ITRAN(IS) = IBLTE(IS)
C
C---- march downstream
      DO 1000 IBL=2, NBL(IS)
        IBM = IBL-1
C
        IW = IBL - IBLTE(IS)
C
        SIMI = IBL.EQ.2
        WAKE = IBL.GT.IBLTE(IS)
c        WAKE = (IBL.GT.IBLTE(IS)).OR.((IS.EQ.1).AND.(IBL.GE.ICLE)
c     *         .AND.(IBL.LE.ICTE))
      ITYPCV = 0
      IF ((IS.EQ.1).AND.(IBL.GE.ICLE).AND.(IBL.LE.ICTE)) THEN
        ITYPCV = 1
      ENDIF
C------ prescribed quantities
        XSI = XSSI(IBL,IS)
        UEI = UEDG(IBL,IS)

        if(UEI.lt.0.0) then
          FL_ERR = .true.
!         return
!         write(*,*) 'UEI=', UEI, 'xbf.f:602'
!         write(*,*) 'IBL=', IBL, 'IS=', IS
!         write(*,*) UEDG(IBL-1,IS), UEDG(IBL,IS), UEDG(IBL+1,IS)
!         do itmp = 1, NBL(IS)
!           write(*,*) itmp, UEDG(itmp,IS), UINV(itmp,IS)
!         end do
          return
!         stop
        end if
C
        IF(WAKE) THEN
         IW = IBL - IBLTE(IS)
         DSWAKI = WGAP(IW)
        ELSE
         DSWAKI = 0.D0
        ENDIF
C
C------ check for transition and set appropriate flags and things
        IF((.NOT.SIMI) .AND. (.NOT.TURB)) CALL TRCHEK(IBL,IS)
C
        DIRECT = .TRUE.
C
C------ Newton iteration loop for current station
        DO 100 ITBL=1, 25
C
C-------- assemble 10x3 linearized system for dCtau, dTh, dDs, dUe, dXi
C         at the previous "1" station and the current "2" station
C         (the "1" station coefficients will be ignored)
C
C
          IF(IBL.EQ.IBLTE(IS)+1) THEN
           TTE = THET(IBLTE(1),1) + THET(IBLTE(2),2)
           DTE = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
           CTE = ( CTAU(IBLTE(1),1)*THET(IBLTE(1),1)
     &           + CTAU(IBLTE(2),2)*THET(IBLTE(2),2) ) / TTE
           CALL TESYS(CTE,TTE,DTE)
          ELSE
           CALL BLSYS
          ENDIF
C
          IF(DIRECT) THEN
C
C--------- try direct mode (set dUe = 0 in currently empty 4th line)
           VS2(4,1) = 0.D0
           VS2(4,2) = 0.D0
           VS2(4,3) = 0.D0
           VS2(4,4) = 1.D0
           VSREZ(4) = 0.D0
C
      if(itbl.eq.1.or.itbl.eq.2)then
 123      format(4f12.5)      
      endif
C--------- solve Newton system for current "2" station
           CALL GAUSS(4,4,VS2,VSREZ,1)
C
C--------- determine max changes and underrelax if necessary
           DMAX = DMAX1( DABS(VSREZ(2)/THI),
     &                 DABS(VSREZ(3)/DSI)  )
           IF(IBL.LT.ITRAN(IS)) DMAX = DMAX1(DMAX,DABS(VSREZ(1)/AMCRIT))
           IF(IBL.GE.ITRAN(IS)) DMAX = DMAX1(DMAX,DABS(VSREZ(1)/CTI   ))
C
           RLX = 1.D0
           IF(DMAX.GT.0.3D0) then
             RLX = 0.3D0/DMAX
           END IF
           RLX = RLX*BETA_BL
C
C--------- see if direct mode is not applicable
           IF(IBL .NE. IBLTE(IS)+1) THEN
C
C---------- calculate resulting kinematic shape parameter Hk
            MSQ = UEI*UEI*HSTINV / (GM1BL*(1.D0-0.5D0*UEI*UEI*HSTINV))
            HTEST = (DSI + RLX*VSREZ(3)) / (THI + RLX*VSREZ(2))
           CALL HKIN( HTEST, MSQ, HKTEST, DUMMY, DUMMY)
C
C---------- decide whether to do direct or inverse problem based on Hk
            IF(IBL.LT.ITRAN(IS)) HMAX = HLMAX
            IF(IBL.GE.ITRAN(IS)) HMAX = HTMAX
            DIRECT = HKTEST.LT.HMAX
           ENDIF
C
           IF(DIRECT) THEN
C---------- update as usual
            IF(IBL.LT.ITRAN(IS)) AMI = AMI + RLX*VSREZ(1)
            IF(IBL.GE.ITRAN(IS)) CTI = CTI + RLX*VSREZ(1)
            THI = THI + RLX*VSREZ(2)
            DSI = DSI + RLX*VSREZ(3)
           ELSE
C---------- set prescribed Hk for inverse calculation at the current station
            IF(IBL.LT.ITRAN(IS)) THEN
C----------- laminar case: relatively slow increase in Hk downstream
             HTARG = HK1 + 0.03D0*(X2-X1)/T1
            ELSE IF(IBL.EQ.ITRAN(IS)) THEN
C----------- transition interval: weighted laminar and turbulent case
             HTARG = HK1 + (0.03D0*(XT-X1) - 0.15D0*(X2-XT))/T1
            ELSE IF(WAKE) THEN
C----------- turbulent wake case:
C-           asymptotic wake behavior with approximate Backward Euler
             CONST = 0.03D0*(X2-X1)/T1
             HK2 = HK1
             HK2 = HK2 - (HK2 +     CONST*(HK2-1.D0)**3 - HK1)
     &                /(1.D0 + 3.D0*CONST*(HK2-1.D0)**2)
             HK2 = HK2 - (HK2 +     CONST*(HK2-1.D0)**3 - HK1)
     &                /(1.D0 + 3.D0*CONST*(HK2-1.D0)**2)
             HK2 = HK2 - (HK2 +     CONST*(HK2-1.D0)**3 - HK1)
     &                /(1.D0 + 3.D0*CONST*(HK2-1.D0)**2)
             HTARG = HK2
            ELSE
C----------- turbulent case: relatively fast decrease in Hk downstream
             HTARG = HK1 - 0.15D0*(X2-X1)/T1
            ENDIF
C
C---------- limit specified Hk to something reasonable
            IF(WAKE) THEN
             HTARG = DMAX1( HTARG , 1.01d0 )
            ELSE
             HTARG = DMAX1( HTARG , HMAX )
            ENDIF
C         Commented out temporarily by Hong Sun
C            WRITE(*,1300) IBL, HTARG
C 1300       FORMAT(' MRCHUE: Inverse mode at', I4, '     Hk =', F8.3)
C
C---------- try again with prescribed Hk
            GO TO 100
C
           ENDIF
C
          ELSE
C
C-------- inverse mode (force Hk to prescribed value HTARG)
           VS2(4,1) = 0.D0
           VS2(4,2) = HK2_T2
           VS2(4,3) = HK2_D2
           VS2(4,4) = HK2_U2
           VSREZ(4) = HTARG - HK2
C
           CALL GAUSS(4,4,VS2,VSREZ,1)
C
           DMAX = DMAX1( DABS(VSREZ(2)/THI),
     &                 DABS(VSREZ(3)/DSI)  )
           IF(IBL.GE.ITRAN(IS)) DMAX = DMAX1( DMAX, DABS(VSREZ(1)/CTI))
C
           RLX = 1.D0
           IF(DMAX.GT.0.3D0) then
             RLX = 0.3D0/DMAX
           end if
           RLX = RLX*BETA_BL
C
C--------- update variables  (AMI is not updated since it will not change)
           IF(IBL.GE.ITRAN(IS)) CTI = CTI + RLX*VSREZ(1)
           THI = THI + RLX*VSREZ(2)
           DSI = DSI + RLX*VSREZ(3)
           UEI = UEI + RLX*VSREZ(4)
C
          ENDIF
C
C-------- eliminate absurd transients
          IF(IBL.LE.IBLTE(IS)) DSI = DMAX1( DSI , 1.05000D0*THI )
          IF(IBL.GT.IBLTE(IS)) DSI = DMAX1( DSI , 1.00005D0*THI )
          IF(IBL.GE.ITRAN(IS)) CTI = DMIN1( CTI , 0.20d0 )
C
          IF(DMAX.LE.1.0D-5) GO TO 110
C
  100   CONTINUE

C        Commented out temporarily  by Hong Sun
C        WRITE(*,1350) IBL, IS, DMAX 
C 1350   FORMAT(' MRCHUE: Convergence failed at',I4,'  side',I2,
C     &         '    Res =', E12.4)
C
C------ the current unconverged solution might still be reasonable...
        IF(DMAX .LE. 0.1D0) GO TO 110

        if(UEI.lt.0.0) then
!         write(*,*) 'UEI=', UEI, 'xbl.f:750'
!         stop
        end if
C
C------- the current solution is garbage --> extrapolate values instead
         IF(IBL.GT.3) THEN 
          IF(IBL.LE.IBLTE(IS)) THEN
           CTI = CTAU(IBM,IS)
           THI = THET(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5D0
           DSI = DSTR(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5D0
          ELSE IF(IBL.EQ.IBLTE(IS)+1) THEN
           CTI = CTE
           THI = TTE
           DSI = DTE
          ELSE
           CTI = CTAU(IBM,IS)
           THI = THET(IBM,IS)
           RATLEN = (XSSI(IBL,IS)-XSSI(IBM,IS)) / (10.D0*DSTR(IBM,IS))
           DSI = (DSTR(IBM,IS) + THI*RATLEN) / (1.D0 + RATLEN)
          ENDIF
          IF(IBL.EQ.ITRAN(IS)) CTI = 0.7D0*CQ1
          UEI = UEDG(IBL,IS)
          IF(IBL.GT.2 .AND. IBL.LT.NBL(IS))
     &     UEI = 0.5D0*(UEDG(IBL-1,IS) + UEDG(IBL+1,IS))
         ENDIF
         if(UEI.lt.0.0) then
           FL_ERR = .true.
!          write(*,*) 'UEI=', UEI, 'xbl.f:776'
!          write(*,*) 'IBL=', IBL, 'IS=', IS
!          write(*,*) UEDG(IBL-1,IS), UEDG(IBL,IS), UEDG(IBL+1,IS)
           return
!          do itmp = 1, NBL(IS)
!            write(*,*) itmp, UEDG(itmp,IS)
!          end do
!          stop
         end if
C
         AMPL2 = AMI
         S2  = CTI
         T2  = THI
         D2  = DSI - DSWAKI
         DW2 = DSWAKI
C
         U2 = UEI*(1.D0-TKBL) / (1.D0 - TKBL*(UEI/QINFBL)**2)
         if(U2.lt.0.0) write(*,*) UEI, TKBL, QINFBL
         U2_UEI = (1.D0 + TKBL*(2.D0*U2*UEI/QINFBL**2 - 1.D0))
     &          / (1.D0 - TKBL*(UEI/QINFBL)**2)
C
C------- set all other extrapolated values for current station
         IF(IBL.LT.ITRAN(IS)) CALL BLVAR(1)
         IF(IBL.GE.ITRAN(IS)) CALL BLVAR(2)
         IF(WAKE) CALL BLVAR(3)
C
C------ pick up here after the Newton iterations
  110   CONTINUE
C
C------ store primary variables
        IF(IBL.LT.ITRAN(IS)) CTAU(IBL,IS) = AMI
        IF(IBL.GE.ITRAN(IS)) CTAU(IBL,IS) = CTI
        THET(IBL,IS) = THI
        DSTR(IBL,IS) = DSI
        UEDG(IBL,IS) = UEI
        MASS(IBL,IS) = DSI*UEI
        TAU(IBL,IS)  = 0.5D0*R2*U2*U2*CF2
C
C------ set "1" variables to "2" variables for next streamwise station
        DO 310 ICOM=1, NCOM
          COM1(ICOM) = COM2(ICOM)
  310   CONTINUE
C
C------ turbulent intervals will follow transition interval or TE
        IF(TRAN .OR. IBL.EQ.IBLTE(IS)) TURB = .TRUE.
        TRAN = .FALSE.
C
        IF(IBL.EQ.IBLTE(IS)) THEN
         THI = THET(IBLTE(1),1) + THET(IBLTE(2),2)
         DSI = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
        ENDIF
C
 1000 CONTINUE
 2000 CONTINUE
C
      RETURN
      END !MRCHUE
  
 
      SUBROUTINE MRCHDU
C...................................................
C     Marches the BLs and wake in mixed mode using
C     the current Ue and Hk.  The calculated Ue
C     and Hk lie along a line quasi-normal to the
C     natural Ue-Hk characteristic line of the
C     current BL so that the Goldstein or Levy-Lees
C     singularity is never encountered.  Continuous
C     checking of transition onset is performed.
C...................................................
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      REAL*8 VTMP(4,5), VZTMP(4)
      real*8 SENNEW

      SENNEW = 0.0d0
C
C---- constant controlling how far Hk is allowed to deviate
C-    from the specified value.
      SENSWT = 1000.D0
C
      DO 2000 IS=1, 2
C
C---- set forced transition arc length position
      CALL XIFSET(IS)
C
C---- set leading edge pressure gradient parameter  x/u du/dx
      IBL = 2
      XSI = XSSI(IBL,IS)
      UEI = UEDG(IBL,IS)
CCC      BULE = DLOG(UEDG(IBL+1,IS)/UEI) / DLOG(XSSI(IBL+1,IS)/XSI)
CCC      BULE = DMAX1( -.08D0 , BULE )
      BULE = 1.D0
C
C---- old transition station
      ITROLD = ITRAN(IS)
C
      TRAN = .FALSE.
      TURB = .FALSE.
      ITRAN(IS) = IBLTE(IS)
C
C---- march downstream
      DO 1000 IBL=2, NBL(IS)
        IBM = IBL-1
C
        SIMI = IBL.EQ.2
        WAKE = IBL.GT.IBLTE(IS)
c        WAKE = (IBL.GT.IBLTE(IS)).OR.((IS.EQ.1).AND.(IBL.GE.ICLE)
c     *         .AND.(IBL.LE.ICTE))
        ITYPCV = 0
        IF ((IS.EQ.1).AND.(IBL.GE.ICLE).AND.(IBL.LE.ICTE)) THEN
          ITYPCV = 1
        ENDIF
C
C------ initialize current station to existing variables
        XSI = XSSI(IBL,IS)
        UEI = UEDG(IBL,IS)
        THI = THET(IBL,IS)
        DSI = DSTR(IBL,IS)
CCC        MDI = MASS(IBL,IS)
        AMI = CTAU(IBL,IS)
        CTI = CTAU(IBL,IS)
C
CCC        DSI = MDI/UEI
C
        IF(WAKE) THEN
         IW = IBL - IBLTE(IS)
         DSWAKI = WGAP(IW)
        ELSE
         DSWAKI = 0.D0
        ENDIF
C
C
C------ check for transition and set appropriate flags and things
        IF((.NOT.SIMI) .AND. (.NOT.TURB)) CALL TRCHEK(IBL,IS)
C
C------ reinitialize shear stress coefficient if transition point moved
        IF(TRAN .AND. IBL.LT.ITROLD) CTI = 0.03D0
        IF(TURB .AND. IBL.LT.ITROLD) CTI = CTAU(IBM,IS)
C
C------ Newton iteration loop for current station
        DO 100 ITBL=1, 25
C
C-------- assemble 10x3 linearized system for dCtau, dTh, dDs, dUe, dXi
C         at the previous "1" station and the current "2" station
C         (the "1" station coefficients will be ignored)
C
          IF(IBL.EQ.IBLTE(IS)+1) THEN
           TTE = THET(IBLTE(1),1) + THET(IBLTE(2),2)
           DTE = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
           CTE = ( CTAU(IBLTE(1),1)*THET(IBLTE(1),1)
     &           + CTAU(IBLTE(2),2)*THET(IBLTE(2),2) ) / TTE
           CALL TESYS(CTE,TTE,DTE)
          ELSE
           CALL BLSYS
          ENDIF
C
          IF(ITBL.EQ.1) THEN
C--------- set "current" Ue and Hk
           UEREF = U2
           HKREF = HK2
          ENDIF
C
          IF(SIMI .OR. IBL.EQ.IBLTE(IS)+1) THEN
C
C--------- for similarity station or first wake point, prescribe Ue
           VS2(4,1) = 0.D0
           VS2(4,2) = 0.D0
           VS2(4,3) = 0.D0
           VS2(4,4) = U2_UEI
           VSREZ(4) = UEREF - U2
C
          ELSE
C
C********* calculate Ue-Hk characteristic slope
C
           DO 20 K=1, 4
             VZTMP(K) = VSREZ(K)
             DO 201 L=1, 5
               VTMP(K,L) = VS2(K,L)
  201        CONTINUE
   20      CONTINUE
C
C--------- set unit dHk
           VTMP(4,1) = 0.D0
           VTMP(4,2) = HK2_T2
           VTMP(4,3) = HK2_D2
           VTMP(4,4) = HK2_U2*U2_UEI
           VZTMP(4)  = 1.D0
C
C--------- calculate dUe response
           CALL GAUSS(4,4,VTMP,VZTMP,1)
C
C--------- set  SENSWT * (normalized dUe/dHk)
           SENNEW = SENSWT * VZTMP(4) * HKREF/UEREF
           IF(ITBL.LE.1) THEN
            IF(SIMI) SENS = SENNEW
           ELSE IF(ITBL.LE.5) THEN
            SENS = SENNEW
           ELSE IF(ITBL.LE.15) THEN
            SENS = 0.5D0*(SENS + SENNEW)
           ENDIF
C
C--------- set prescribed Ue-Hk combination
           VS2(4,1) = 0.D0
           VS2(4,2) =  HK2_T2 * HKREF
           VS2(4,3) =  HK2_D2 * HKREF
           VS2(4,4) =( HK2_U2 * HKREF  +  SENS/UEREF )*U2_UEI
           VSREZ(4) = -(HKREF**2)*(HK2 / HKREF - 1.D0)
     &                     - SENS*(U2  / UEREF - 1.D0)
C
          ENDIF

C-------- solve Newton system for current "2" station
          CALL GAUSS(4,4,VS2,VSREZ,1)
C
C-------- determine max changes and underrelax if necessary
          DMAX = DMAX1( DABS(VSREZ(2)/THI),
     &                DABS(VSREZ(3)/DSI)  )
          IF(IBL.GE.ITRAN(IS)) DMAX = DMAX1(DMAX,DABS(VSREZ(1)/CTI))
C
          RLX = 1.D0
          IF(DMAX.GT.0.3D0) then
            RLX = 0.3D0/DMAX
          end if
          RLX = RLX*BETA_BL
C
C-------- update as usual  (AMI is not updated since it does not change)
          IF(IBL.GE.ITRAN(IS)) CTI = CTI + RLX*VSREZ(1)
          THI = THI + RLX*VSREZ(2)
          DSI = DSI + RLX*VSREZ(3)
          UEI = UEI + RLX*VSREZ(4)
C
C-------- eliminate absurd transients
          IF(IBL.LE.IBLTE(IS)) DSI = DMAX1( DSI , 1.05000D0*THI )
          IF(IBL.GT.IBLTE(IS)) DSI = DMAX1( DSI , 1.00005D0*THI )
          IF(IBL.GE.ITRAN(IS)) CTI = DMIN1( CTI , 0.20d0 )
C
          IF(DMAX.LE.5.0D-6) GO TO 110
C
  100   CONTINUE
C        Commented out temporarily   by Hong Sun
C        WRITE(*,1350) IBL, IS, DMAX 
C 1350   FORMAT(' MRCHDU: Convergence failed at',I4,'  side',I2,
C     &         '    Res =', E12.4)
C
C------ the current unconverged solution might still be reasonable...
        IF(DMAX .LE. 0.1D0) GO TO 110
C
C------- the current solution is garbage --> extrapolate values instead
         IF(IBL.GT.3) THEN
          IF(IBL.LE.IBLTE(IS)) THEN
           CTI = CTAU(IBM,IS)
           THI = THET(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5D0
           DSI = DSTR(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5D0
           UEI = UEDG(IBM,IS)
          ELSE IF(IBL.EQ.IBLTE(IS)+1) THEN
           CTI = CTE
           THI = TTE
           DSI = DTE
           UEI = UEDG(IBM,IS)
          ELSE
           CTI = CTAU(IBM,IS)
           THI = THET(IBM,IS)
           RATLEN = (XSSI(IBL,IS)-XSSI(IBM,IS)) / (10.D0*DSTR(IBM,IS))
           DSI = (DSTR(IBM,IS) + THI*RATLEN) / (1.D0 + RATLEN)
           UEI = UEDG(IBM,IS)
          ENDIF
          IF(IBL.EQ.ITRAN(IS)) CTI = 0.7D0*CQ1
         ENDIF
C
         AMPL2 = AMI
         S2  = CTI
         T2  = THI
         D2  = DSI - DSWAKI
         DW2 = DSWAKI
         U2 = UEI*(1.D0-TKBL) / (1.D0 - TKBL*(UEI/QINFBL)**2)
         U2_UEI = (1.D0 + TKBL*(2.D0*U2*UEI/QINFBL**2 - 1.D0))
     &          / (1.D0 - TKBL*(UEI/QINFBL)**2)
C
C------- set all other extrapolated values for current station
         IF(IBL.LT.ITRAN(IS)) CALL BLVAR(1)
         IF(IBL.GE.ITRAN(IS)) CALL BLVAR(2)
         IF(WAKE) CALL BLVAR(3)
C
C------ pick up here after the Newton iterations
  110   CONTINUE
C
        SENS = SENNEW
C
C------ store primary variables
        IF(IBL.LT.ITRAN(IS)) CTAU(IBL,IS) = AMI
        IF(IBL.GE.ITRAN(IS)) CTAU(IBL,IS) = CTI
        THET(IBL,IS) = THI
        DSTR(IBL,IS) = DSI
        UEDG(IBL,IS) = UEI
        MASS(IBL,IS) = DSI*UEI
        TAU(IBL,IS)  = 0.5D0*R2*U2*U2*CF2
C
C------ set "1" variables to "2" variables for next streamwise station
        DO 310 ICOM=1, NCOM
          COM1(ICOM) = COM2(ICOM)
  310   CONTINUE
C
C------ turbulent intervals will follow transition interval
        IF(TRAN .OR. IBL.EQ.IBLTE(IS)) TURB = .TRUE.
        TRAN = .FALSE.
C
 1000 CONTINUE
C
 2000 CONTINUE
C
      RETURN
      END   !MRCHDU
C
C  
C 
      SUBROUTINE XIFSET(IS)
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
C
      IF(XTR1(IS).GE.1.D0) THEN
       XIFORC = XSSI(IBLTE(IS),IS)
c       write(71,*) ' XIFORC1 =', XIFORC
       RETURN
      ENDIF

C
      IF(IS.EQ.1) THEN
C
C----- set approximate arc length value of forced transition point for SINVRT
       SFORCE = SLE + (    -SLE)*XTR1(IS)
C
C----- calculate actual arc length
!s--- YE TIAN 07/03/2013-------
! The sinvrt subroutine is declared as:
!     SUBROUTINE SINVRT(SI,XI,X,XS,S,N)
! where XI is a scalar.
!      CALL SINVRT(SFORCE,XTR1,X,XP,S,Nvis)
       CALL SINVRT(SFORCE,XTR1(IS),X,XP,S,Nvis)
!e--- YE TIAN 07/03/2013-------
C
C----- set BL coordinate value
       XIFORC = DMIN1( (SST - SFORCE) , XSSI(IBLTE(IS),IS) ) 
c       write(71,*) ' XIFORC2 =', XIFORC
      ELSE 
C
       SFORCE = SLE + (S(Nvis)-SLE)*XTR1(IS)

!s--- YE TIAN 07/03/2013-------
!      CALL SINVRT(SFORCE,XTR1,X,XP,S,Nvis)
       CALL SINVRT(SFORCE,XTR1(IS),X,XP,S,Nvis)
!e--- YE TIAN 07/03/2013-------
       XIFORC = DMIN1( (SFORCE - SST) , XSSI(IBLTE(IS),IS) )
c      write(71,*) ' XIFORC3 =', XIFORC
C
      ENDIF
C
      RETURN
      END
 
 
      SUBROUTINE TRCHEK(IBL,IS)
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
C
      X1 = XSSI(IBL-1,IS)
      X2 = XSSI(IBL  ,IS)
C
C---- calculate AMPL2 value
      CALL DAMPL( HK1, T1, RT1, AX, AX_HK1, AX_T1, AX_RT1 )
      AMPL2 = AMPL1 + AX*(X2-X1)
      AMI = AMPL2
C
C---- test for free or forced transition
      TRFREE = AMI.GE.AMCRIT
      TRFORC = XIFORC.GT.X1 .AND. XIFORC.LE.X2
C
C---- set transition interval flag
      TRAN = TRFORC .OR. TRFREE
C
      IF(.NOT.TRAN) THEN
       ITRAN(IS) = IBL+2
       RETURN
      ENDIF
C
C---- resolve if both forced and free transition
      IF(TRFREE .AND. TRFORC) THEN
       XT = (AMCRIT-AMPL1)/AX  +  X1
       TRFORC = XIFORC .LT. XT
       TRFREE = XIFORC .GE. XT
      ENDIF
C
      IF(TRFORC) THEN
C----- if forced transition, then XT is prescribed
       XT = XIFORC
       XT_A1 = 0.D0
       XT_X1 = 0.D0
       XT_T1 = 0.D0
       XT_D1 = 0.D0
       XT_U1 = 0.D0
       XT_X2 = 0.D0
       XT_T2 = 0.D0
       XT_D2 = 0.D0
       XT_U2 = 0.D0
      ELSE
C----- if free transition, XT is related to BL variables
       XT    =  (AMCRIT-AMPL1)/AX     + X1
       XT_AX = -(AMCRIT-AMPL1)/AX**2
C
       XT_A1 = -1.D0/AX
       XT_X1 = 1.D0
       XT_T1 = XT_AX*(AX_HK1*HK1_T1 + AX_T1 + AX_RT1*RT1_T1)
       XT_D1 = XT_AX*(AX_HK1*HK1_D1                        )
       XT_U1 = XT_AX*(AX_HK1*HK1_U1         + AX_RT1*RT1_U1)
       XT_X2 = 0.D0
       XT_T2 = 0.D0
       XT_D2 = 0.D0
       XT_U2 = 0.D0
      ENDIF
C
C---- save transition location
      ITRAN(IS) = IBL
      TFORCE(IS) = TRFORC
      XSSITR(IS) = XT
C
C---- save info for user output
      IF(TFORCE(IS)) THEN
       XTR(IS) = XTR1(IS)
      ELSE
C----- interpolate airfoil geometry to find transition x/c
       SB1 = SST - X1
       SB2 = SST - X2
       IF(IS.EQ.2) THEN
        SB1 = SST + X1
        SB2 = SST + X2
       ENDIF
       XB1 = SEVAL(SB1,X,XP,S,Nvis)
       XB2 = SEVAL(SB2,X,XP,S,Nvis)
       XTR(IS) = XB1 + (XB2-XB1)*(XT-X1)/(X2-X1)
      ENDIF
C
      RETURN
      END


      SUBROUTINE BLSYS
C..................................................................
C
C     Sets up the BL Newton system governing the current interval:
C
C     |       ||dA1|     |       ||dA2|       |     |
C     |  VS1  ||dT1|  +  |  VS2  ||dT2|   =   |VSREZ|
C     |       ||dD1|     |       ||dD2|       |     |
C              |dU1|              |dU2|
C              |dX1|              |dX2|
C
C        3x5    5x1         3x5    5x1          3x1
C
C     The system as shown corresponds to a laminar station
C     If TRAN, then  dS2  replaces  dA2
C     If TURB, then  dS1, dS2  replace  dA1, dA2
C
C..................................................................
C
      implicit real*8 (a-h,o-z)
      IMPLICIT REAL*8 (M)
      INCLUDE 'XBL.INC'
C
C---- set primary BL variables from current values
      X2 = XSI
      AMPL2 = AMI
      S2  = CTI
      T2  = THI
      D2  = DSI - DSWAKI
      DW2 = DSWAKI
C
!s--YE TIAN --- 06/03/2013---
!     write(*,*) 'UEI=',UEI, 'xbl.f:1229'
!     UEI = abs(UEI)
!e--YE TIAN --- 06/03/2013---

      U2 = UEI*(1.D0-TKBL) / (1.D0 - TKBL*(UEI/QINFBL)**2)
      U2_UEI = (1.D0 + TKBL*(2.D0*U2*UEI/QINFBL**2 - 1.D0))
     &       / (1.D0 - TKBL*(UEI/QINFBL)**2)
C
C---- calculate secondary BL variables and their sensitivities
      IF(WAKE) THEN
       CALL BLVAR(3)
      ELSE IF(TURB.OR.TRAN) THEN
       CALL BLVAR(2)
      ELSE
       CALL BLVAR(1)
      ENDIF
C
C---- for the similarity station, "1" and "2" variables are the same
      IF(SIMI) THEN
       DO 3 ICOM=1, NCOM
         COM1(ICOM) = COM2(ICOM)
    3  CONTINUE
      ENDIF
C
C---- set up appropriate finite difference system for current interval

      IF(TRAN) THEN
       CALL TRDIF
      ELSE IF(SIMI) THEN
       CALL BLDIF(0)
      ELSE IF(.NOT.TURB) THEN
       CALL BLDIF(1)
      ELSE IF(WAKE) THEN
       CALL BLDIF(3)
      ELSE IF(TURB) THEN
       CALL BLDIF(2)
      ENDIF

C
      IF(SIMI) THEN
C----- at similarity station, "1" variables are really "2" variables
       DO 10 K=1, 4
         DO 101 L=1, 5
           VS2(K,L) = VS1(K,L) + VS2(K,L)
           VS1(K,L) = 0.D0
  101    CONTINUE
   10  CONTINUE
      ENDIF
C
C---- change system over into incompressible dUe
      DO 20 K=1, 4
        VS1(K,4) = VS1(K,4)*U1_UEI
        VS2(K,4) = VS2(K,4)*U2_UEI
   20 CONTINUE
C
      RETURN
      END
 

      SUBROUTINE TESYS(CTE,TTE,DTE)
      implicit real*8 (a-h,o-z)

      INCLUDE 'XBL.INC'
C
      DO 55 K=1, 4
        VSREZ(K) = 0.D0
        DO 551 L=1, 5
          VS1(K,L) = 0.D0
          VS2(K,L) = 0.D0
  551   CONTINUE
   55 CONTINUE
C
C---- set primary BL variables from current values
      X2 = XSI
      AMPL2 = AMI
      S2  = CTI
      T2  = THI
      D2  = DSI - DSWAKI
      DW2 = DSWAKI
C
      U2 = UEI*(1.D0-TKBL) / (1.D0 - TKBL*(UEI/QINFBL)**2)
      U2_UEI = (1.D0 + TKBL*(2.D0*U2*UEI/QINFBL**2 - 1.D0))
     &       / (1.D0 - TKBL*(UEI/QINFBL)**2)
C
      CALL BLVAR(3)
C
      VS1(1,1) = -1.D0
      VS2(1,1) = 1.D0
      VSR(1)   = 0.D0
      VSREZ(1) = CTE - S2      
C
      VS1(2,2) = -1.D0
      VS2(2,2) = 1.D0
      VSR(2)   = 0.D0
      VSREZ(2) = TTE - T2
C
      VS1(3,3) = -1.D0
      VS2(3,3) = 1.D0
      VSR(3)   = 0.D0
      VSREZ(3) = DTE - D2 - DW2
C
      RETURN
      END

 
      SUBROUTINE BLVAR(ITYP)
C...................................................
C     Calculates all secondary "2" variables from
C     the primary "2" variables X2, U2, T2, D2, S2.
C     Also calculates the sensitivities of the
C     secondary variables wrt the primary variables.
C
C      ITYP = 1 :  laminar
C      ITYP = 2 :  turbulent
C      ITYP = 3 :  turbulent wake
C
C     ITYPCV = 0 : Cf not 0
C     ITYPCV = 1 : Cf = 0
C
C....................................................
      implicit real*8 (a-h,o-z)
 
      INCLUDE 'XBL.INC'
C
C---- shear coefficient constant (proportional to Cebeci & Smith 0.0168)
      CTCON = 0.015D0
C
C---- set edge Mach number ** 2
      M2    = U2*U2*HSTINV / (GM1BL*(1.D0 - 0.5D0*U2*U2*HSTINV))
      TR2   = 1.0D0 + 0.5D0*GM1BL*M2
      M2_U2 = 2.0D0*M2*TR2/U2
C
C---- set edge static density (isentropic relation)
      R2    = RSTBL*TR2**(-1.D0/GM1BL)
      R2_U2 = -R2/TR2 * 0.5D0*M2_U2
C
C---- set shape parameter
      H2    = D2/T2
      H2_D2 = 1.D0/T2
      H2_T2 = -H2/T2
C
C---- set edge static/stagnation enthalpy and molecular viscosity
      HERAT = 1.D0 - 0.5D0*U2*U2*HSTINV
      V2 = DSQRT((HERAT)**3) * (1.D0+HVRAT)/(HERAT+HVRAT)/REYBL
      V2_U2 = V2*U2*(1.D0/(HERAT+HVRAT)-1.5D0/HERAT)*HSTINV
      V2_RE = -V2/REYBL
C
C---- set kinematic shape parameter
      CALL HKIN( H2, M2, HK2, HK2_H2, HK2_M2 )
C
      IF(ITYP.EQ.3) HK2 = DMAX1(HK2,1.00005d0)
      IF(ITYP.NE.3) HK2 = DMAX1(HK2,1.05000d0)
C
      HK2_U2 =                HK2_M2*M2_U2
      HK2_T2 = HK2_H2*H2_T2
      HK2_D2 = HK2_H2*H2_D2
C
C---- set momentum thickness Reynolds number
      if (R2.lt.0.0) write(*,*) 'R2=',R2
      if (U2.lt.0.0) write(*,*) 'U2=',U2
      if (T2.lt.0.0) write(*,*) 'T2=',T2
      RT2    = R2*U2*T2/V2
      RT2_U2 = RT2*(1.D0/U2 + R2_U2/R2 - V2_U2/V2)
      RT2_T2 = RT2/T2
      RT2_RE = -RT2/V2 * V2_RE
C
C---- density thickness shape parameter     ( H** )
      CALL HCT( HK2, M2, HC2, HC2_HK2, HC2_M2 )
      HC2_U2 = HC2_HK2*HK2_U2 + HC2_M2*M2_U2
      HC2_T2 = HC2_HK2*HK2_T2
      HC2_D2 = HC2_HK2*HK2_D2
C
!     write(1237,*) 'xbl.f:1434'
!     write(1237,*) 'ITYP=', ITYP
C---- set KE thickness shape parameter from  H - H*  correlations
      IF(ITYP.EQ.1) THEN
       CALL HSL( HK2, RT2, M2, HS2, HS2_HK2, HS2_RT2, HS2_M2 )
      ELSE
       CALL HST( HK2, RT2, M2, HS2, HS2_HK2, HS2_RT2, HS2_M2 )
      ENDIF
C
      HS2_U2 = HS2_HK2*HK2_U2 + HS2_RT2*RT2_U2 + HS2_M2*M2_U2
      HS2_T2 = HS2_HK2*HK2_T2 + HS2_RT2*RT2_T2
      HS2_D2 = HS2_HK2*HK2_D2
      HS2_RE =                  HS2_RT2*RT2_RE
C
C---- normalized slip velocity  Us
      US2     = 0.5D0*HS2*( 3.D0 - 4.D0*(HK2-1.D0)/H2   )/3.D0
      US2_HS2 = 0.5D0  *  ( 3.D0 - 4.D0*(HK2-1.D0)/H2   )/3.D0
      US2_HK2 = 0.5D0*HS2*(      - 4.D0           /H2   )/3.D0
      US2_H2  = 0.5D0*HS2*(        4.D0*(HK2-1.D0)/H2**2)/3.D0
C
      US2_U2 = US2_HS2*HS2_U2 + US2_HK2*HK2_U2
      US2_T2 = US2_HS2*HS2_T2 + US2_HK2*HK2_T2 + US2_H2*H2_T2
      US2_D2 = US2_HS2*HS2_D2 + US2_HK2*HK2_D2 + US2_H2*H2_D2
      US2_RE = US2_HS2*HS2_RE
C
      IF(ITYP.LE.2 .AND. US2.GT.0.95D0) THEN
CCC       WRITE(*,*) 'BLVAR: Us clamped:', US2
       US2 = 0.98D0
       US2_U2 = 0.D0
       US2_T2 = 0.D0
       US2_D2 = 0.D0
       US2_RE = 0.D0
      ENDIF
C
      IF(ITYP.EQ.3 .AND. US2.GT.0.99995D0) THEN
CCC       WRITE(*,*) 'BLVAR: Wake Us clamped:', US2
       US2 = 0.99995D0
       US2_U2 = 0.D0
       US2_T2 = 0.D0
       US2_D2 = 0.D0
       US2_RE = 0.D0
      ENDIF
C
C---- equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2
      HKB = HK2 - 1.D0
      USB = 1.D0 - US2
!s--YE TIAN -- 06/10/2013----
      if (CTCON*HS2*HKB**3 / (USB*H2*HK2**2).lt.0.0) then
!       write(*,*) CTCON, HS2, HKB, USB, H2, HK2, 'xbl.f:1456'
!       HS2 = abs(HS2)
      end if
!e--YE TIAN -- 06/10/2013----
      CQ2     =
     &    DSQRT( CTCON*HS2*HKB**3 / (USB*H2*HK2**2) )
      CQ2_HS2 = CTCON  *  HKB**3 / (USB*H2*HK2**2)        * 0.5D0/CQ2
      CQ2_US2 = CTCON*HS2*HKB**3 / (USB*H2*HK2**2) / USB  * 0.5D0/CQ2
      CQ2_HK2 = CTCON*HS2*HKB**2 / (USB*H2*HK2**2) * 3.D0 * 0.5D0/CQ2
     &        - CTCON*HS2*HKB**3 / (USB*H2*HK2**3) * 2.D0 * 0.5D0/CQ2
      CQ2_H2  =-CTCON*HS2*HKB**3 / (USB*H2*HK2**2) / H2   * 0.5D0/CQ2
C
      CQ2_U2 = CQ2_HS2*HS2_U2 + CQ2_HK2*HK2_U2 + CQ2_US2*US2_U2
      CQ2_T2 = CQ2_HS2*HS2_T2 + CQ2_HK2*HK2_T2 + CQ2_US2*US2_T2
      CQ2_D2 = CQ2_HS2*HS2_D2 + CQ2_HK2*HK2_D2 + CQ2_US2*US2_D2
C
      CQ2_T2 = CQ2_T2 + CQ2_H2*H2_T2
      CQ2_D2 = CQ2_D2 + CQ2_H2*H2_D2
C
      CQ2_RE = CQ2_HS2*HS2_RE + CQ2_US2*US2_RE
C
      IF(ITYP.EQ.3) THEN
C----- 4x bigger  CtauEQ  for wake (eddy viscosity 4x bigger)
       CQ2    = CQ2   *2.D0
       CQ2_U2 = CQ2_U2*2.D0
       CQ2_T2 = CQ2_T2*2.D0
       CQ2_D2 = CQ2_D2*2.D0
       CQ2_RE = CQ2_RE*2.D0
      ENDIF
C
C---- skin friction coefficient  (zero in wake)
c      IF((ITYP.EQ.3).OR.(ITYPCV.EQ.1)) THEN
      IF(ITYP.EQ.3) THEN
       CF2     = 0.D0
       CF2_HK2 = 0.D0
       CF2_RT2 = 0.D0
       CF2_M2  = 0.D0
      ELSE IF(ITYP.EQ.1) THEN
       CALL CFL( HK2, RT2, M2, CF2, CF2_HK2, CF2_RT2, CF2_M2 )
      ELSE
       CALL CFT( HK2, RT2, M2, CF2, CF2_HK2, CF2_RT2, CF2_M2 )
      ENDIF
C
      CF2_U2 = CF2_HK2*HK2_U2 + CF2_RT2*RT2_U2 + CF2_M2*M2_U2
      CF2_T2 = CF2_HK2*HK2_T2 + CF2_RT2*RT2_T2
      CF2_D2 = CF2_HK2*HK2_D2
      CF2_RE =                  CF2_RT2*RT2_RE
C
C---- set similarity variables if not defined
      IF(SIMI) THEN
       HK1    = HK2
       HK1_T1 = HK2_T2
       HK1_D1 = HK2_D2
       HK1_U1 = HK2_U2
       RT1    = RT2
       RT1_T1 = RT2_T2
       RT1_U1 = RT2_U2
       RT1_RE = RT2_RE
       M1    = M2
       M1_U1 = M2_U2
      ENDIF
C
C---- define stuff for midpoint CF
      HKA = 0.5D0*(HK1 + HK2)
      RTA = 0.5D0*(RT1 + RT2)
      MA  = 0.5D0*(M1  + M2 )

C
C---- midpoint skin friction coefficient  (zero in wake)
      IF ((ITYP.EQ.3).OR.(ITYPCV.EQ.1)) THEN
c      IF(ITYP.EQ.3) THEN
       CFA     = 0.D0
       CFA_HKA = 0.D0
       CFA_RTA = 0.D0
       CFA_MA  = 0.D0
      ELSE IF(ITYP.EQ.1) THEN
       CFA_RTA = 1.0d0
       CALL CFL( HKA, RTA, MA, CFA, CFA_HKA, CFA_RTA, CFA_MA )
      ELSE
       CALL CFT( HKA, RTA, MA, CFA, CFA_HKA, CFA_RTA, CFA_MA )
      ENDIF

C
      CFA_U1 = 0.5D0*(CFA_HKA*HK1_U1 + CFA_MA*M1_U1 + CFA_RTA*RT1_U1)
      CFA_T1 = 0.5D0*(CFA_HKA*HK1_T1 +                CFA_RTA*RT1_T1)
      CFA_D1 = 0.5D0*(CFA_HKA*HK1_D1                                )
C
      CFA_U2 = 0.5D0*(CFA_HKA*HK2_U2 + CFA_MA*M2_U2 + CFA_RTA*RT2_U2)
      CFA_T2 = 0.5D0*(CFA_HKA*HK2_T2 +                CFA_RTA*RT2_T2)
      CFA_D2 = 0.5D0*(CFA_HKA*HK2_D2                                )
C
      CFA_RE = 0.5D0*CFA_RTA*(RT1_RE + RT2_RE)
C
C---- dissipation function    2 CD / H*
      IF(ITYP.EQ.1) THEN
       CALL DIL( HK2, RT2, DI2, DI2_HK2, DI2_RT2 )
C
       DI2_U2 = DI2_HK2*HK2_U2 + DI2_RT2*RT2_U2
       DI2_T2 = DI2_HK2*HK2_T2 + DI2_RT2*RT2_T2
       DI2_D2 = DI2_HK2*HK2_D2
       DI2_S2 = 0.D0
       DI2_RE =                  DI2_RT2*RT2_RE
      ELSE
       CALL DIT(     HS2,     US2,     CF2,     S2, DI2,
     &           DI2_HS2, DI2_US2, DI2_CF2, DI2_S2      )
C
       DI2_U2 = DI2_HS2*HS2_U2 + DI2_US2*US2_U2 + DI2_CF2*CF2_U2
       DI2_T2 = DI2_HS2*HS2_T2 + DI2_US2*US2_T2 + DI2_CF2*CF2_T2
       DI2_D2 = DI2_HS2*HS2_D2 + DI2_US2*US2_D2 + DI2_CF2*CF2_D2
       DI2_RE = DI2_HS2*HS2_RE + DI2_US2*US2_RE + DI2_CF2*CF2_RE
C
C----- add on CD contribution of inner shear layer
       IF(ITYP.EQ.3 .AND. DW2.GT.0.D0) THEN
        DKON = 0.03D0*0.75D0**3
        DDI = DKON*US2**3
        DDI_US2 = 3.0D0*DKON*US2**2
        DI2 = DI2 + DDI * DW2/DWTE
        DI2_U2 = DI2_U2 + DDI_US2*US2_U2 * DW2/DWTE
        DI2_T2 = DI2_T2 + DDI_US2*US2_T2 * DW2/DWTE
        DI2_D2 = DI2_D2 + DDI_US2*US2_D2 * DW2/DWTE
        DI2_RE = DI2_RE + DDI_US2*US2_RE * DW2/DWTE
       ENDIF
       IF(ITYP.EQ.3) THEN
C------ double dissipation for the wake (two wake halves)
        DI2    = DI2   *2.0D0
        DI2_U2 = DI2_U2*2.0D0
        DI2_T2 = DI2_T2*2.0D0
        DI2_D2 = DI2_D2*2.0D0
        DI2_RE = DI2_RE*2.0D0
       ENDIF
      ENDIF
C
C---- BL thickness (Delta) from simplified Green's correlation
      DE2     = (3.15D0 + 1.72D0/(HK2-1.0D0)   )*T2  +  D2
      DE2_HK2 = (       - 1.72D0/(HK2-1.0D0)**2)*T2
C
      DE2_U2 = DE2_HK2*HK2_U2
      DE2_T2 = DE2_HK2*HK2_T2 + (3.15D0 + 1.72D0/(HK2-1.D0))
      DE2_D2 = DE2_HK2*HK2_D2 + 1.0D0
C
      RETURN
      END
 
 
      SUBROUTINE TRDIF
C...............................................
C     Sets up the Newton system governing the
C     transition interval.  Equations governing
C     the  laminar  part  X1 < xi < XT  and
C     the turbulent part  XT < xi < X2
C     are simply summed.
C...............................................
      implicit real*8 (a-h,o-z)

      INCLUDE 'XBL.INC'
      REAL  BL1(4,5), BL2(4,5), BLREZ(4), BLR(4)
     &    , BT1(4,5), BT2(4,5), BTREZ(4), BTR(4)
C
C---- save variables and sensitivities for future restoration
      DO 5 ICOM=1, NCOM
        C1SAV(ICOM) = COM1(ICOM)
        C2SAV(ICOM) = COM2(ICOM)
    5 CONTINUE
C
C---- weighting factors for linear interpolation to transition point
      WF2    = (XT-X1)/(X2-X1)
      WF2_XT = 1.D0/(X2-X1)
C
      WF2_A1 = WF2_XT*XT_A1
      WF2_X1 = WF2_XT*XT_X1 + (WF2-1.D0)/(X2-X1)
      WF2_X2 = WF2_XT*XT_X2 -  WF2     /(X2-X1)
      WF2_T1 = WF2_XT*XT_T1
      WF2_T2 = WF2_XT*XT_T2
      WF2_D1 = WF2_XT*XT_D1
      WF2_D2 = WF2_XT*XT_D2
      WF2_U1 = WF2_XT*XT_U1
      WF2_U2 = WF2_XT*XT_U2
C
      WF1    = 1.D0 - WF2
      WF1_A1 = -WF2_A1
      WF1_X1 = -WF2_X1
      WF1_X2 = -WF2_X2
      WF1_T1 = -WF2_T1
      WF1_T2 = -WF2_T2
      WF1_D1 = -WF2_D1
      WF1_D2 = -WF2_D2
      WF1_U1 = -WF2_U1
      WF1_U2 = -WF2_U2
C
C
C**** FIRST,  do laminar part between X1 and XT
C
C-----interpolate primary variables to transition point
      TT    = T1*WF1    + T2*WF2
      TT_A1 = T1*WF1_A1 + T2*WF2_A1
      TT_X1 = T1*WF1_X1 + T2*WF2_X1
      TT_X2 = T1*WF1_X2 + T2*WF2_X2
      TT_T1 = T1*WF1_T1 + T2*WF2_T1 + WF1
      TT_T2 = T1*WF1_T2 + T2*WF2_T2 + WF2
      TT_D1 = T1*WF1_D1 + T2*WF2_D1
      TT_D2 = T1*WF1_D2 + T2*WF2_D2
      TT_U1 = T1*WF1_U1 + T2*WF2_U1
      TT_U2 = T1*WF1_U2 + T2*WF2_U2
C
      DT    = D1*WF1    + D2*WF2
      DT_A1 = D1*WF1_A1 + D2*WF2_A1
      DT_X1 = D1*WF1_X1 + D2*WF2_X1
      DT_X2 = D1*WF1_X2 + D2*WF2_X2
      DT_T1 = D1*WF1_T1 + D2*WF2_T1
      DT_T2 = D1*WF1_T2 + D2*WF2_T2
      DT_D1 = D1*WF1_D1 + D2*WF2_D1 + WF1
      DT_D2 = D1*WF1_D2 + D2*WF2_D2 + WF2
      DT_U1 = D1*WF1_U1 + D2*WF2_U1
      DT_U2 = D1*WF1_U2 + D2*WF2_U2
C
      UT    = U1*WF1    + U2*WF2
      UT_A1 = U1*WF1_A1 + U2*WF2_A1
      UT_X1 = U1*WF1_X1 + U2*WF2_X1
      UT_X2 = U1*WF1_X2 + U2*WF2_X2
      UT_T1 = U1*WF1_T1 + U2*WF2_T1
      UT_T2 = U1*WF1_T2 + U2*WF2_T2
      UT_D1 = U1*WF1_D1 + U2*WF2_D1
      UT_D2 = U1*WF1_D2 + U2*WF2_D2
      UT_U1 = U1*WF1_U1 + U2*WF2_U1 + WF1
      UT_U2 = U1*WF1_U2 + U2*WF2_U2 + WF2
C
C---- set "2" variables to primary "T" variables at XT
      X2 = XT
      U2 = UT
      T2 = TT
      D2 = DT
C
C---- calculate laminar secondary "T" variables
      CALL BLVAR(1)
C=
C=    at this point, all "2" variables are really "T" variables at XT
C=
C
C---- set up Newton system for dAm, dTh, dDs, dUe, dXi  at  X1 and XT
      CALL BLDIF(1)
C
C---- The current Newton system is in terms of "1" and "T" variables,
C-    so calculate its equivalent in terms of "1" and "2" variables.
C-    In other words, convert residual sensitivities wrt "T" variables
C-    into sensitivities wrt "1" and "2" variables.  The amplification
C-    equation is unnecessary here, so the K=1 row is left empty.
      DO 10 K=2, 3
        BLREZ(K) = VSREZ(K)
        BLR(K)   = VSR(K)
C
        BL1(K,1) = VS1(K,1)
     &           + VS2(K,2)*TT_A1
     &           + VS2(K,3)*DT_A1
     &           + VS2(K,4)*UT_A1
     &           + VS2(K,5)*XT_A1
        BL1(K,2) = VS1(K,2)
     &           + VS2(K,2)*TT_T1
     &           + VS2(K,3)*DT_T1
     &           + VS2(K,4)*UT_T1
     &           + VS2(K,5)*XT_T1
        BL1(K,3) = VS1(K,3)
     &           + VS2(K,2)*TT_D1
     &           + VS2(K,3)*DT_D1
     &           + VS2(K,4)*UT_D1
     &           + VS2(K,5)*XT_D1
        BL1(K,4) = VS1(K,4)
     &           + VS2(K,2)*TT_U1
     &           + VS2(K,3)*DT_U1
     &           + VS2(K,4)*UT_U1
     &           + VS2(K,5)*XT_U1
        BL1(K,5) = VS1(K,5)
     &           + VS2(K,2)*TT_X1
     &           + VS2(K,3)*DT_X1
     &           + VS2(K,4)*UT_X1
     &           + VS2(K,5)*XT_X1
C
        BL2(K,1) = 0.D0
        BL2(K,2) = VS2(K,2)*TT_T2
     &           + VS2(K,3)*DT_T2
     &           + VS2(K,4)*UT_T2
     &           + VS2(K,5)*XT_T2
        BL2(K,3) = VS2(K,2)*TT_D2
     &           + VS2(K,3)*DT_D2
     &           + VS2(K,4)*UT_D2
     &           + VS2(K,5)*XT_D2
        BL2(K,4) = VS2(K,2)*TT_U2
     &           + VS2(K,3)*DT_U2
     &           + VS2(K,4)*UT_U2
     &           + VS2(K,5)*XT_U2
        BL2(K,5) = VS2(K,2)*TT_X2
     &           + VS2(K,3)*DT_X2
     &           + VS2(K,4)*UT_X2
     &           + VS2(K,5)*XT_X2
C
   10 CONTINUE
C
C
C**** SECOND, set up turbulent part between XT and X2  ****
C
C---- calculate equilibrium shear coefficient CQT at transition point
      CALL BLVAR(2)
C
C---- set initial shear coefficient value ST at transition point
C-    ( note that CQ2, CQ2_T2, etc. are really "CQT", "CQT_TT", etc.)
C
      CTR     = 1.8D0*DEXP(-3.3D0/(HK2-1.0D0))
      CTR_HK2 = CTR * 3.3D0/(HK2-1.0D0)**2
C
CCC      CTR = 1.2D0
CCC      CTR = 0.7D0
CCC      CTR_HK2 = 0.0D0
C
      ST    = CTR*CQ2
      ST_TT = CTR*CQ2_T2 + CQ2*CTR_HK2*HK2_T2
      ST_DT = CTR*CQ2_D2 + CQ2*CTR_HK2*HK2_D2
      ST_UT = CTR*CQ2_U2 + CQ2*CTR_HK2*HK2_U2
      ST_RE = CTR*CQ2_RE
C
C---- calculate ST sensitivities wrt the actual "1" and "2" variables
      ST_A1 = ST_TT*TT_A1 + ST_DT*DT_A1 + ST_UT*UT_A1
      ST_X1 = ST_TT*TT_X1 + ST_DT*DT_X1 + ST_UT*UT_X1
      ST_X2 = ST_TT*TT_X2 + ST_DT*DT_X2 + ST_UT*UT_X2
      ST_T1 = ST_TT*TT_T1 + ST_DT*DT_T1 + ST_UT*UT_T1
      ST_T2 = ST_TT*TT_T2 + ST_DT*DT_T2 + ST_UT*UT_T2
      ST_D1 = ST_TT*TT_D1 + ST_DT*DT_D1 + ST_UT*UT_D1
      ST_D2 = ST_TT*TT_D2 + ST_DT*DT_D2 + ST_UT*UT_D2
      ST_U1 = ST_TT*TT_U1 + ST_DT*DT_U1 + ST_UT*UT_U1
      ST_U2 = ST_TT*TT_U2 + ST_DT*DT_U2 + ST_UT*UT_U2
C
      S2 = ST
C
C---- recalculate turbulent secondary "T" variables using proper ST
      CALL BLVAR(2)
C
C---- set "1" variables to "T" variables and reset "2" variables
C-    to their saved turbulent values
      DO 30 ICOM=1, NCOM
        COM1(ICOM) = COM2(ICOM)
        COM2(ICOM) = C2SAV(ICOM)
   30 CONTINUE
C
      CALL BLVAR(2)
C
C---- set up Newton system for dCt, dTh, dDs, dUe, dXi  at  XT and X2
      CALL BLDIF(2)
C
C---- convert sensitivities wrt "T" variables into sensitivities
C-    wrt "1" and "2" variables as done before for the laminar part
      DO 40 K=1, 3
        BTREZ(K) = VSREZ(K)
        BTR(K)   = VSR(K) + VS1(K,1)*ST_RE
C
        BT1(K,1) = VS1(K,1)*ST_A1
     &           + VS1(K,2)*TT_A1
     &           + VS1(K,3)*DT_A1
     &           + VS1(K,4)*UT_A1
     &           + VS1(K,5)*XT_A1
        BT1(K,2) = VS1(K,1)*ST_T1
     &           + VS1(K,2)*TT_T1
     &           + VS1(K,3)*DT_T1
     &           + VS1(K,4)*UT_T1
     &           + VS1(K,5)*XT_T1
        BT1(K,3) = VS1(K,1)*ST_D1
     &           + VS1(K,2)*TT_D1
     &           + VS1(K,3)*DT_D1
     &           + VS1(K,4)*UT_D1
     &           + VS1(K,5)*XT_D1
        BT1(K,4) = VS1(K,1)*ST_U1
     &           + VS1(K,2)*TT_U1
     &           + VS1(K,3)*DT_U1
     &           + VS1(K,4)*UT_U1
     &           + VS1(K,5)*XT_U1
        BT1(K,5) = VS1(K,1)*ST_X1
     &           + VS1(K,2)*TT_X1
     &           + VS1(K,3)*DT_X1
     &           + VS1(K,4)*UT_X1
     &           + VS1(K,5)*XT_X1
C
        BT2(K,1) = VS2(K,1)
        BT2(K,2) = VS2(K,2)
     &           + VS1(K,1)*ST_T2
     &           + VS1(K,2)*TT_T2
     &           + VS1(K,3)*DT_T2
     &           + VS1(K,4)*UT_T2
     &           + VS1(K,5)*XT_T2
        BT2(K,3) = VS2(K,3)
     &           + VS1(K,1)*ST_D2
     &           + VS1(K,2)*TT_D2
     &           + VS1(K,3)*DT_D2
     &           + VS1(K,4)*UT_D2
     &           + VS1(K,5)*XT_D2
        BT2(K,4) = VS2(K,4)
     &           + VS1(K,1)*ST_U2
     &           + VS1(K,2)*TT_U2
     &           + VS1(K,3)*DT_U2
     &           + VS1(K,4)*UT_U2
     &           + VS1(K,5)*XT_U2
        BT2(K,5) = VS2(K,5)
     &           + VS1(K,1)*ST_X2
     &           + VS1(K,2)*TT_X2
     &           + VS1(K,3)*DT_X2
     &           + VS1(K,4)*UT_X2
     &           + VS1(K,5)*XT_X2
C
   40 CONTINUE
C
C---- Add up laminar and turbulent parts to get final system
C-    in terms of honest-to-God "1" and "2" variables.
      VSREZ(1) =            BTREZ(1)
      VSREZ(2) = BLREZ(2) + BTREZ(2)
      VSREZ(3) = BLREZ(3) + BTREZ(3)
      VSR(1)   =            BTR(1)
      VSR(2)   = BLR(2)   + BTR(2)
      VSR(3)   = BLR(3)   + BTR(3)
      DO 60 L=1, 5
        VS1(1,L) =            BT1(1,L)
        VS2(1,L) =            BT2(1,L)
        VS1(2,L) = BL1(2,L) + BT1(2,L)
        VS2(2,L) = BL2(2,L) + BT2(2,L)
        VS1(3,L) = BL1(3,L) + BT1(3,L)
        VS2(3,L) = BL2(3,L) + BT2(3,L)
   60 CONTINUE
C
C---- To be sanitary, restore "1" quantities which got clobbered
C-    in all of the numerical gymnastics above.  The "2" variables
C-    were already restored for the XT-X2 differencing part.
      DO 70 ICOM=1, NCOM
        COM1(ICOM) = C1SAV(ICOM)
   70 CONTINUE
C
      RETURN
      END
 
 
      SUBROUTINE BLDIF(ITYP)
C............................................................
C     Sets up the Newton system coefficients and residuals
C
C        ITYP = 0 :  similarity station
C        ITYP = 1 :  laminar interval
C        ITYP = 2 :  turbulent interval
C        ITYP = 3 :  wake interval
C
C     This routine knows nothing about a transition interval,
C     which is taken care of by TRDIF.
C............................................................
      implicit real*8 (a-h,o-z)

      INCLUDE 'XBL.INC'
C
C---- shear coefficient lag constant (see M.D.'s thesis)
      SCC = 5.6D0
C
      IF(ITYP.EQ.0) THEN
C----- similarity logarithmic differences  (prescribed)
       XLOG = 1.D0
       ULOG = BULE
       TLOG = 0.5D0*(1.D0 - BULE)
       HLOG = 0.D0
       DDLOG = 0.D0
      ELSE
C----- usual logarithmic differences
       XLOG = DLOG(X2/X1)
!s--- YE TIAN ----
!      if((U2.lt.0.0).or.(U1.le.0.0)) write(*,*) U2,U1, 'xbl.f:1911'
       U2=ABS(U2)
       U1=ABS(U1)
!e--- YE TIAN ----
       ULOG = DLOG(U2/U1)
       TLOG = DLOG(T2/T1)
       HLOG = DLOG(HS2/HS1)
C       XLOG = 2.D0*(X2-X1)/(X2+X1)
C       ULOG = 2.D0*(U2-U1)/(U2+U1)
C       TLOG = 2.D0*(T2-T1)/(T2+T1)
C       HLOG = 2.D0*(HS2-HS1)/(HS2+HS1)
       DDLOG = 1.D0
      ENDIF
C
      DO 55 K=1, 4
        VSREZ(K) = 0.D0
        DO 551 L=1, 5
          VS1(K,L) = 0.D0
          VS2(K,L) = 0.D0
  551   CONTINUE
   55 CONTINUE
C
C---- set triggering constant for local upwinding
      HDCON = 5.D0
C
C---- less upwinding in the wake
      IF(WAKE) HDCON = 1.D0
C
C---- local upwinding is based on local change in   log (Hk-1)
      ARG = DABS((HK2-1.D0)/(HK1-1.D0))
      HL = DLOG(ARG)
      HL_HK1 = -1.D0/(HK1-1.D0)
      HL_HK2 =  1.D0/(HK2-1.D0)
C
C---- set local upwinding parameter UPW and linearize it
C
C       UPW = 0.5D0   Trapezoidal
C       UPW = 1.0D0   Backward Euler
C
      HLSQ = DMIN1( HL**2 , 15.0d0 )
      UPW = 1.D0 - 0.5D0*DEXP(-HLSQ*HDCON)
      UPW_HL = DEXP(-HLSQ*HDCON)*HL*HDCON
C
      UPW_HK1 = UPW_HL*HL_HK1
      UPW_HK2 = UPW_HL*HL_HK2
C
      UPW_U1 = UPW_HK1*HK1_U1
      UPW_T1 = UPW_HK1*HK1_T1
      UPW_D1 = UPW_HK1*HK1_D1
      UPW_U2 = UPW_HK2*HK2_U2
      UPW_T2 = UPW_HK2*HK2_T2
      UPW_D2 = UPW_HK2*HK2_D2
C
      IF(ITYP.EQ.0) THEN
C
C***** LE point -->  set zero amplification factor
       VS2(1,1) = 1.D0
       VSR(1)   = 0.D0
       VSREZ(1) = -AMPL2
C
      ELSE IF(ITYP.EQ.1) THEN
C
C***** laminar part -->  set amplification equation
       CALL DAMPL(HK1,T1,RT1,AX,AX_HK1,AX_T1,AX_RT1)
       REZC = AMPL2 - AMPL1 - AX*(X2-X1)
       Z_AX = -(X2-X1)
C
       VS1(1,1) =-1.D0
       VS1(1,2) = Z_AX*AX_HK1*HK1_T1 + Z_AX*AX_T1 + Z_AX*AX_RT1*RT1_T1
       VS1(1,3) = Z_AX*AX_HK1*HK1_D1
       VS1(1,4) = Z_AX*AX_HK1*HK1_U1              + Z_AX*AX_RT1*RT1_U1
       VS1(1,5) = AX
       VS2(1,1) = 1.D0
       VS2(1,2) = 0.D0
       VS2(1,3) = 0.D0
       VS2(1,4) = 0.D0
       VS2(1,5) =-AX
       VSR(1)   = 0.D0
       VSREZ(1) = -REZC
C
      ELSE
C
C***** turbulent part -->  set shear lag equation
C
       SA  = (1.D0-UPW)*S1  + UPW*S2
       CQA = (1.D0-UPW)*CQ1 + UPW*CQ2
       CFU = (1.D0-UPW)*CF1 + UPW*CF2
       HKA = (1.D0-UPW)*HK1 + UPW*HK2
C
       DEA = 0.5D0*(DE1 + DE2)
       DA  = 0.5D0*(D1  + D2 )
C
C----- set and linearize  equilibrium 1/Ue dUe/dx
       HR = (HKA - 1.D0) / (6.7D0*HKA)
       HR_HKA = 1.D0 / (6.7D0*HKA**2)
       UQ     = ( 0.5D0*CFU - HR**2 ) / (0.75D0*DA)
       UQ_CFA =   0.5D0               / (0.75D0*DA)
       UQ_HKA =      -2.0D0*HR*HR_HKA / (0.75D0*DA)
       UQ_DA  = -UQ/DA
C
       IF(ITYP.EQ.3) THEN
C------ quadruple it for wake  (just like equilibrium Ctau)
        UQ     = 4.D0*UQ
        UQ_CFA = 4.D0*UQ_CFA
        UQ_HKA = 4.D0*UQ_HKA
        UQ_DA  = 4.D0*UQ_DA
       ENDIF
C
       SLOG = DLOG(S2/S1)
       DXI = X2 - X1
C
       REZC = SCC*(CQA - SA)*DXI  -  DEA*2.D0*SLOG
     &      + DEA*2.D0*(UQ*DXI - ULOG)
C
       Z_UQ = DEA*2.D0*DXI
C
       Z_CFA = Z_UQ*UQ_CFA
       Z_HKA = Z_UQ*UQ_HKA
       Z_DA  = Z_UQ*UQ_DA
       Z_SL = -DEA*2.D0
       Z_UL = -DEA*2.D0
       Z_DXI = SCC*(CQA - SA) + DEA*2.D0*UQ
       Z_CQA = SCC*DXI
       Z_SA = -SCC*DXI
       Z_DEA = 2.D0*(UQ*DXI - ULOG - SLOG)
C
       Z_UPW = Z_CQA*(CQ2-CQ1) + Z_SA *(S2 -S1 )
     &       + Z_CFA*(CF2-CF1) + Z_HKA*(HK2-HK1)
       Z_DE1 = 0.5D0*Z_DEA
       Z_DE2 = 0.5D0*Z_DEA
       Z_D1  = 0.5D0*Z_DA
       Z_D2  = 0.5D0*Z_DA
       Z_U1  =                 - Z_UL/U1
       Z_U2  =                   Z_UL/U2
       Z_X1  = -Z_DXI
       Z_X2  =  Z_DXI
       Z_S1  = (1.D0-UPW)*Z_SA  - Z_SL/S1
       Z_S2  =       UPW *Z_SA  + Z_SL/S2
       Z_CQ1 = (1.D0-UPW)*Z_CQA
       Z_CQ2 =       UPW *Z_CQA
       Z_CF1 = (1.D0-UPW)*Z_CFA
       Z_CF2 =       UPW *Z_CFA
       Z_HK1 = (1.D0-UPW)*Z_HKA
       Z_HK2 =       UPW *Z_HKA
C
       VS1(1,1) = Z_S1
       VS1(1,2) =        Z_UPW*UPW_T1 + Z_DE1*DE1_T1
       VS1(1,3) = Z_D1 + Z_UPW*UPW_D1 + Z_DE1*DE1_D1
       VS1(1,4) = Z_U1 + Z_UPW*UPW_U1 + Z_DE1*DE1_U1
       VS1(1,5) = Z_X1
       VS2(1,1) = Z_S2
       VS2(1,2) =        Z_UPW*UPW_T2 + Z_DE2*DE2_T2
       VS2(1,3) = Z_D2 + Z_UPW*UPW_D2 + Z_DE2*DE2_D2
       VS2(1,4) = Z_U2 + Z_UPW*UPW_U2 + Z_DE2*DE2_U2
       VS2(1,5) = Z_X2
C
       VS1(1,2) = VS1(1,2) + Z_CQ1*CQ1_T1 + Z_CF1*CF1_T1 + Z_HK1*HK1_T1
       VS1(1,3) = VS1(1,3) + Z_CQ1*CQ1_D1 + Z_CF1*CF1_D1 + Z_HK1*HK1_D1
       VS1(1,4) = VS1(1,4) + Z_CQ1*CQ1_U1 + Z_CF1*CF1_U1 + Z_HK1*HK1_U1
C
       VS2(1,2) = VS2(1,2) + Z_CQ2*CQ2_T2 + Z_CF2*CF2_T2 + Z_HK2*HK2_T2
       VS2(1,3) = VS2(1,3) + Z_CQ2*CQ2_D2 + Z_CF2*CF2_D2 + Z_HK2*HK2_D2
       VS2(1,4) = VS2(1,4) + Z_CQ2*CQ2_U2 + Z_CF2*CF2_U2 + Z_HK2*HK2_U2
C
       VSR(1)   = Z_CQ1*CQ1_RE + Z_CQ2*CQ2_RE
     &          + Z_CF1*CF1_RE + Z_CF2*CF2_RE
       VSREZ(1) = -REZC
C
      ENDIF
C
C**** Set up momentum equation
      HA = 0.5D0*(H1 + H2)
      MA = 0.5D0*(M1 + M2)
      XA = 0.5D0*(X1 + X2)
      TA = 0.5D0*(T1 + T2)
      HWA = 0.5D0*(DW1/T1 + DW2/T2)
C
C---- set Cf term, using central value CFA for better accuracy in drag
      CFX     = 0.5D0*CFA*XA/TA  +  0.25D0*(CF1*X1/T1 + CF2*X2/T2)
      CFX_XA  = 0.5D0*CFA   /TA
      CFX_TA  = -.5D0*CFA*XA/TA**2
C
      CFX_X1  = 0.25D0*CF1   /T1     + CFX_XA*0.5D0
      CFX_X2  = 0.25D0*CF2   /T2     + CFX_XA*0.5D0
      CFX_T1  = -.25D0*CF1*X1/T1**2  + CFX_TA*0.5D0
      CFX_T2  = -.25D0*CF2*X2/T2**2  + CFX_TA*0.5D0
      CFX_CF1 = 0.25D0*    X1/T1
      CFX_CF2 = 0.25D0*    X2/T2
      CFX_CFA = 0.50D0*    XA/TA
C
      BTMP = HA + 2.D0 - MA + HWA
C
c---- wbrewer -- 4/10/95 -- add cavity drag in terms of momentum thickness
c-    dtcav = delta theta due to cavity drag
      if(jpan.eq.1) then
         REZT  = TLOG + BTMP*ULOG - XLOG*0.5D0*CFX - DTCAV/TA
      else
         REZT  = TLOG + BTMP*ULOG - XLOG*0.5D0*CFX 
      endif
c
      Z_CFX = -XLOG*0.5D0
      Z_HA  =  ULOG
      Z_HWA =  ULOG
c    ----------------------------- dianr
      Z_MA  = -ULOG*2.D0
      Z_XL  =-DDLOG * 0.5D0*CFX
      Z_UL  = DDLOG * BTMP
      Z_TL  = DDLOG
C
      Z_CFA = Z_CFX*CFX_CFA
      Z_CF1 = Z_CFX*CFX_CF1
      Z_CF2 = Z_CFX*CFX_CF2
C
c---- wbrewer -- 4/10/95 -- add contribution due to dtcav
      if(jpan.eq.1) then
         Z_T1 = -Z_TL/T1 + Z_CFX*CFX_T1 + Z_HWA*0.5D0*(-DW1/T1**2) 
     &           - 0.5D0 + DTCAV/TA**2
         Z_T2 =  Z_TL/T2 + Z_CFX*CFX_T2 + Z_HWA*0.5D0*(-DW2/T2**2)
     &           - 0.5D0 + DTCAV/TA**2
      else
         Z_T1 = -Z_TL/T1 + Z_CFX*CFX_T1 + Z_HWA*0.5D0*(-DW1/T1**2)
         Z_T2 =  Z_TL/T2 + Z_CFX*CFX_T2 + Z_HWA*0.5D0*(-DW2/T2**2)
      endif
c
      Z_X1 = -Z_XL/X1 + Z_CFX*CFX_X1
      Z_X2 =  Z_XL/X2 + Z_CFX*CFX_X2
      Z_U1 = -Z_UL/U1
      Z_U2 =  Z_UL/U2
C
      VS1(2,2) = 0.5D0*Z_HA*H1_T1 + Z_CFA*CFA_T1 + Z_CF1*CF1_T1 + Z_T1
      VS1(2,3) = 0.5D0*Z_HA*H1_D1 + Z_CFA*CFA_D1 + Z_CF1*CF1_D1
      VS1(2,4) = 0.5D0*Z_MA*M1_U1 + Z_CFA*CFA_U1 + Z_CF1*CF1_U1 + Z_U1
      VS1(2,5) =                                                Z_X1
      VS2(2,2) = 0.5D0*Z_HA*H2_T2 + Z_CFA*CFA_T2 + Z_CF2*CF2_T2 + Z_T2
      VS2(2,3) = 0.5D0*Z_HA*H2_D2 + Z_CFA*CFA_D2 + Z_CF2*CF2_D2
      VS2(2,4) = 0.5D0*Z_MA*M2_U2 + Z_CFA*CFA_U2 + Z_CF2*CF2_U2 + Z_U2
      VS2(2,5) =                                                Z_X2
 123      format(4f12.5)      
 124      format(3f15.5)
C
      VSR(2)   = Z_CFA*CFA_RE + Z_CF1*CF1_RE + Z_CF2*CF2_RE
      VSREZ(2) = -REZT
C
C**** Set up shape parameter equation
C
      XOT1 = X1/T1
      XOT2 = X2/T2
C
      HA  = 0.5D0*(H1  + H2 )
      HSA = 0.5D0*(HS1 + HS2)
      HCA = 0.5D0*(HC1 + HC2)
      HWA = 0.5D0*(DW1/T1 + DW2/T2)
C
      DIX = (1.D0-UPW)*DI1*XOT1 + UPW*DI2*XOT2
      CFX = (1.D0-UPW)*CF1*XOT1 + UPW*CF2*XOT2
      DIX_UPW = DI2*XOT2 - DI1*XOT1
      CFX_UPW = CF2*XOT2 - CF1*XOT1
C
      BTMP = 2.D0*HCA/HSA + 1.D0 - HA - HWA
C
      REZH  = HLOG + BTMP*ULOG + XLOG*(0.5D0*CFX-DIX)
      Z_CFX =  XLOG*0.5D0
      Z_DIX = -XLOG
      Z_HCA = 2.D0*ULOG/HSA
      Z_HA  = -ULOG
      Z_HWA = -ULOG
      Z_XL  = DDLOG * (0.5D0*CFX-DIX)
      Z_UL  = DDLOG * BTMP
      Z_HL  = DDLOG
C
      Z_UPW = Z_CFX*CFX_UPW + Z_DIX*DIX_UPW
C
      Z_HS1 = -HCA*ULOG/HSA**2 - Z_HL/HS1
      Z_HS2 = -HCA*ULOG/HSA**2 + Z_HL/HS2
C
      Z_CF1 = (1.D0-UPW)*Z_CFX*XOT1
      Z_CF2 =       UPW *Z_CFX*XOT2
      Z_DI1 = (1.D0-UPW)*Z_DIX*XOT1
      Z_DI2 =       UPW *Z_DIX*XOT2
C
      Z_T1 = (1.D0-UPW)*(Z_CFX*CF1 + Z_DIX*DI1)*(-XOT1/T1)
      Z_T2 =       UPW *(Z_CFX*CF2 + Z_DIX*DI2)*(-XOT2/T2)
      Z_X1 = (1.D0-UPW)*(Z_CFX*CF1 + Z_DIX*DI1)/ T1        - Z_XL/X1
      Z_X2 =       UPW *(Z_CFX*CF2 + Z_DIX*DI2)/ T2        + Z_XL/X2
      Z_U1 =                                               - Z_UL/U1
      Z_U2 =                                                 Z_UL/U2
C
      Z_T1 = Z_T1 + Z_HWA*0.5D0*(-DW1/T1**2)
      Z_T2 = Z_T2 + Z_HWA*0.5D0*(-DW2/T2**2)
C
      VS1(3,1) =                               Z_DI1*DI1_S1
      VS1(3,2) = Z_HS1*HS1_T1 + Z_CF1*CF1_T1 + Z_DI1*DI1_T1 + Z_T1
      VS1(3,3) = Z_HS1*HS1_D1 + Z_CF1*CF1_D1 + Z_DI1*DI1_D1
      VS1(3,4) = Z_HS1*HS1_U1 + Z_CF1*CF1_U1 + Z_DI1*DI1_U1 + Z_U1
      VS1(3,5) =                                              Z_X1
      VS2(3,1) =                               Z_DI2*DI2_S2
      VS2(3,2) = Z_HS2*HS2_T2 + Z_CF2*CF2_T2 + Z_DI2*DI2_T2 + Z_T2
      VS2(3,3) = Z_HS2*HS2_D2 + Z_CF2*CF2_D2 + Z_DI2*DI2_D2
      VS2(3,4) = Z_HS2*HS2_U2 + Z_CF2*CF2_U2 + Z_DI2*DI2_U2 + Z_U2
      VS2(3,5) =                                              Z_X2
C
      VS1(3,2) = VS1(3,2) +0.5D0*(Z_HCA*HC1_T1+Z_HA*H1_T1) +Z_UPW*UPW_T1
      VS1(3,3) = VS1(3,3) +0.5D0*(Z_HCA*HC1_D1+Z_HA*H1_D1) +Z_UPW*UPW_D1
      VS1(3,4) = VS1(3,4) +0.5D0*(Z_HCA*HC1_U1           ) +Z_UPW*UPW_U1
      VS2(3,2) = VS2(3,2) +0.5D0*(Z_HCA*HC2_T2+Z_HA*H2_T2) +Z_UPW*UPW_T2
      VS2(3,3) = VS2(3,3) +0.5D0*(Z_HCA*HC2_D2+Z_HA*H2_D2) +Z_UPW*UPW_D2
      VS2(3,4) = VS2(3,4) +0.5D0*(Z_HCA*HC2_U2           ) +Z_UPW*UPW_U2
C
      VSR(3)   = Z_HS1*HS1_RE + Z_CF1*CF1_RE + Z_DI1*DI1_RE
     &         + Z_HS2*HS2_RE + Z_CF2*CF2_RE + Z_DI2*DI2_RE
     &                        + Z_CFA*CFA_RE
      VSREZ(3) = -REZH
C
      RETURN
      END



c      SUBROUTINE DAMPL( HK, TH, RT, AX, AX_HK, AX_TH, AX_RT )
cC........................................................
cC
cC     ORIGINAL VERSION
cC
cC     output:   AX   spatial amplification rate
cC
cC     input :   HK     kinematic shape parameter
cC               TH     momentum thickness
cC               RT     momentum thickness Reynolds number
cC........................................................
c      implicit real*8 (a-h,o-z)
c      IMPLICIT REAL (A-H,M,O-Z)
cC
cC---- log10(Critical Rth) - H   correlation for Falkner-Skan profiles
c      HMI = 1.D0/(HK-1.D0)
c      GRCRIT = (1.415D0*HMI-0.489D0)*DTANH(20.D0*HMI-12.9D0) + 
C               3.295D0*HMI + 0.440D0
cC
c      IF(DLOG10(RT) .LT. GRCRIT) THEN
cC
cC----- no amplification for Rtheta < Rcrit
c       AX    = 0.D0
c       AX_HK = 0.D0
c       AX_TH = 0.D0
c       AX_RT = 0.D0
cC
c      ELSE
cC
cC----- Amplification envelope slope correlation for Falkner-Skan
c       TANHHK = DTANH(1.5D0*(HK-3.1D0))
c       ARG = 2.4D0*HK - 3.7D0 + 2.5D0*TANHHK
c       DADR    = 0.01D0*DSQRT(ARG*ARG + 0.25D0)
c       DADR_HK = 0.01D0/DSQRT(ARG*ARG + 0.25D0) * ARG
c     &         * (2.4D0 + 3.75D0*(1.D0-TANHHK**2))
cC
cC----- convert amplification rate in Rtheta to rate in xi
c       TFS    =  (6.54D0*HK - 14.07D0   )/HK**2
c       TFS_HK = -(6.54D0    - 28.14D0/HK)/HK**2
c       BUH    = ( 0.058D0*(HK-4.D0)**2/(HK-1.D0) - 0.068D0 )/TFS
c       BUH_HK = ( 0.058D0*(HK-4.D0) /  (HK-1.D0)
c     &            * (2.D0 - (HK-4.D0)/(HK-1.D0))        )/TFS
c     &        - BUH/TFS * TFS_HK
c       AX    = 0.5D0*(BUH+1.D0)*TFS*DADR/TH
c       AX_HK = 0.5D0     *     TFS*DADR/TH *  BUH_HK
c     &       + 0.5D0*(BUH+1.D0)  *  DADR/TH *  TFS_HK
c     &       + 0.5D0*(BUH+1.D0)*TFS   /  TH * DADR_HK
c       AX_TH =-AX/TH
c       AX_RT = 0.D0
cC
c      ENDIF
cC
c      IF(AX.LE.0.D0) THEN
c       AX    = 0.00001D0
c       AX_HK = 0.D0
c       AX_TH = 0.D0
c       AX_RT = 0.D0
c      ENDIF
cC
c      RETURN
c      END ! DAMPL

 
      SUBROUTINE DAMPL( HK, TH, RT, AX, AX_HK, AX_TH, AX_RT )
C........................................................
C
C     NEW VERSION.   March 1991
C
C     output:   AX   spatial amplification rate
C
C     input :   HK     kinematic shape parameter
C               TH     momentum thickness
C               RT     momentum thickness Reynolds number
C........................................................
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA DGR / 0.04D0 /
C
C---- log10(Critical Rth) - H   correlation for Falkner-Skan profiles
      HMI = 1.D0/(HK - 1.D0)
      GRCRIT = 2.492D0*HMI**(0.43D0)
     &       + 0.7D0*(1.D0+DTANH(14.D0*HMI-9.24D0))
C
!s-- YE TIAN 6/10/2013-----
      if (RT.lt.0.0) then
!       write(*,*) 'RT=',RT,'xbl.f:2319'
!       RT = abs(RT)
      end if
!e-- YE TIAN 6/10/2013-----
      GR = DLOG10(RT)
      GR_RT = 1.D0 / (2.3026D0*RT)
C
      IF(GR .LT. GRCRIT-DGR) THEN
C
C----- no amplification for Rtheta < Rcrit
       AX    = 0.D0
       AX_HK = 0.D0
       AX_TH = 0.D0
       AX_RT = 0.D0
C
      ELSE
C
C----- Set steep cubic ramp used to turn on AX smoothly as Rtheta 
C-     exceeds Rcrit (previously, this was done discontinuously).
C-     The ramp goes between  -DGR < log10(Rtheta/Rcrit) < DGR
C
       RNORM = (GR - (GRCRIT-DGR)) / (2.D0*DGR)
       IF(RNORM .GE. 1.D0) THEN
        RFAC = 1.D0
        RFAC_RT = 0.D0
       ELSE
        RFAC    =  3.D0*RNORM**2 - 2.D0*RNORM**3
        RFAC_RT = (6.D0*RNORM    - 6.D0*RNORM**2)/(2.D0*DGR) * GR_RT
       ENDIF
C
C----- Amplification envelope slope correlation for Falkner-Skan
       HMI_HK = -HMI**2
C
       ARG    = 3.87D0*HMI    - 2.52D0
       ARG_HK = 3.87D0*HMI_HK
C
       EX    = DEXP(-ARG**2)
       EX_HK = EX * (-2.D0*ARG*ARG_HK)
C
       DADR    = 0.028D0*(HK-1.D0) - 0.0345D0*EX
       DADR_HK = 0.028D0           - 0.0345D0*EX_HK
C
C----- new m(H) correlation    1 March 91
       AF = -0.05D0 + 2.7D0*HMI -  5.5D0*HMI**2 + 3.D0*HMI**3
       AF_HMI =       2.7D0     - 11.D0*HMI     + 9.D0*HMI**2
       AF_HK = AF_HMI*HMI_HK
C
       AX    = (AF   *DADR/TH                ) * RFAC
       AX_HK = (AF_HK*DADR/TH + AF*DADR_HK/TH) * RFAC
       AX_TH = -AX/TH
       AX_RT = (AF   *DADR/TH                ) * RFAC_RT
C
      ENDIF
C
      IF(AX.LE.0.D0) THEN
       AX    = 0.00001D0
       AX_HK = 0.D0
       AX_TH = 0.D0
       AX_RT = 0.D0
      ENDIF
C
      RETURN
      END ! DAMPL

 
 
      SUBROUTINE HKIN( H, MSQ, HK, HK_H, HK_MSQ )
           implicit real*8 (a-h,o-z) 
      REAL*8 MSQ
C
C---- calculate kinematic shape parameter (assuming air)
C     (from Whitfield )
      HK     = (H - 0.29D0*MSQ) / (1.D0 + 0.113D0*MSQ)
      HK_H   =  1.D0            / (1.D0 + 0.113D0*MSQ)
      HK_MSQ = (-.29D0 - 0.113D0*HK) / (1.D0 + 0.113D0*MSQ)
C
      RETURN
      END
 

C======================================================================
C---- new correlations
C
      SUBROUTINE DIL( HK, RT, DI, DI_HK, DI_RT )
           implicit real*8 (a-h,o-z) 
C
C---- Laminar dissipation function  ( 2 CD/H* )     (from Falkner-Skan)
      IF(HK.LT.4.D0) THEN
       DI    = ( 0.00205D0  *  (4.D0-HK)**(5.5D0) + 0.207D0 ) / RT
       DI_HK = ( -.00205D0*5.5D0*(4.D0-HK)**(4.5D0)         ) / RT
      ELSE
       HKB = HK - 4.D0
       DEN = 1.D0 + 0.02D0*HKB**2
       DI    = ( -.0016D0  *  HKB**2  /DEN   + 0.207D0           ) / RT
       DI_HK = ( -.0016D0  *  2.D0 * HKB * 
     &         ( 1.D0/DEN - 0.02D0*HKB**2/DEN**2 ) ) / RT
      ENDIF
      DI_RT = -DI/RT
C
      RETURN
      END


      SUBROUTINE HSL( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
           implicit real*8 (a-h,o-z) 
      REAL*8 MSQ
C
C---- Laminar HS correlation
      IF(HK.LT.4.35D0) THEN
       TMP = HK - 4.35D0
       HS    = 0.0111D0*TMP**2/(HK+1.D0)
     &       - 0.0278D0*TMP**3/(HK+1.D0)  + 1.528D0
     &       - 0.0002D0*(TMP*HK)**2
       HS_HK = 0.0111D0*(2.D0*TMP    - TMP**2/(HK+1.D0))/(HK+1.D0)
     &       - 0.0278D0*(3.D0*TMP**2 - TMP**3/(HK+1.D0))/(HK+1.D0)
     &       - 0.0002D0*2.D0*TMP*HK * (TMP + HK)
      ELSE
       HS    = 0.015D0*     (HK-4.35D0)**2/HK + 1.528D0
       HS_HK = 0.015D0*2.D0*(HK-4.35D0)   /HK
     &       - 0.015D0*     (HK-4.35D0)**2/HK**2
      ENDIF
C
      HS_RT  = 0.D0
      HS_MSQ = 0.D0
C
      RETURN
      END


      SUBROUTINE CFL( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )
      implicit real*8 (a-h,o-z) 
      REAL*8 HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ
C
C---- Laminar skin friction function  ( Cf )    ( from Falkner-Skan )
      IF(HK.LT.5.5D0) THEN
       TMP = (5.5D0-HK)**3 / (HK+1.D0)
       CF    = ( 0.0727D0*TMP                      - 0.07D0       )/RT
       CF_HK = (-.0727D0*TMP*3.D0/(5.5D0-HK)-0.0727D0*TMP/(HK+1.D0))/RT
      ELSE
       TMP = 1.D0 - 1.D0/(HK-4.5D0)
       CF    = ( 0.015D0*TMP**2      - 0.07D0  ) / RT
       CF_HK = ( 0.015D0*TMP*2.D0/(HK-4.5D0)**2 ) / RT
      ENDIF
      CF_RT = -CF/RT
      CF_MSQ = 0.D0
C
      RETURN
      END
C
C=====================================================================

 
C=====================================================================
C---- old correlations
C
c     SUBROUTINE DIL( HK, RT, DI, DI_HK, DI_RT )
c
c---- Laminar dissipation function  ( 2 CD/H* )     (from Falkner-Skan)
c     IF(HK.LT.4.D0) THEN
c      DI    = ( 0.00205D0  *  (4.D0-HK)**(5.5D0) + 0.207D0 ) / RT
c      DI_HK = ( -.00205D0*5.5D0*(4.D0-HK)**(4.5D0)         ) / RT
c     ELSE
c      HKB = HK - 4.D0
c      DEN = 1.D0 + 0.02D0*HKB**2
c      DI    = ( -.003D0  *  HKB**2  /DEN    + 0.207D0             ) / RT
c      DI_HK = ( -.003D0*2.D0*HKB*(1.0D0/DEN - 0.02D0*HKB**2/DEN**2) ) / RT
c     ENDIF
c     DI_RT = -DI/RT
c
c     RETURN
c     END
c
c
c     SUBROUTINE HSL( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
c     REAL MSQ
c
c---- Laminar HS correlation    ( from Falkner-Skan )
c     IF(HK.LT.4.D0) THEN
c      HS    = 0.076D0*(HK-4.D0)**2/HK + 1.515D0
c      HS_HK = 0.076D0*(1.D0-16.D0/HK**2)
c     ELSE
c      HS    = 0.040D0*(HK-4.D0)**2/HK + 1.515D0
c      HS_HK = 0.040D0*(1.D0-16.D0/HK**2)
c     ENDIF
c
c     HS_RT  = 0.D0
c     HS_MSQ = 0.D0
c
c     RETURN
c     END
c
c
c     SUBROUTINE CFL( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )
c     REAL MSQ
c
c---- Laminar skin friction function  ( Cf )    ( from Falkner-Skan )
c     IF(HK .LT. 7.4D0) THEN
c      TMP = (7.4D0-HK)**2   / (HK-1.D0)
c      CF    = ( 0.03954D0*TMP                         - 0.134D0 ) / RT
c      CF_HK = ( -.03954D0*TMP * (2.0D0/(7.4D0-HK) + 1.D0/(HK-1.D0)) ) / RT
c     ELSE
c      TMP = 1.D0 - 1.4D0/(HK-6.D0)
c      CF    = ( 0.044D0*TMP**2      - 0.134D0 ) / RT
c      CF_HK = ( 0.088D0*TMP*1.4D0/(HK-6.D0)**2 ) / RT
c     ENDIF
c     CF_RT = -CF/RT
c     CF_MSQ = 0.D0
c
c     RETURN
c     END
c
C===================================================================== 


      SUBROUTINE DIT( HS, US, CF, ST, DI, DI_HS, DI_US, DI_CF, DI_ST )
           implicit real*8 (a-h,o-z) 
C
C---- Turbulent dissipation function  ( 2 CD/H* )
      DI    =  ( 0.5D0*CF*US + ST*ST*(1.D0-US) ) * 2.D0/HS
      DI_HS = -( 0.5D0*CF*US + ST*ST*(1.D0-US) ) * 2.D0/HS**2
      DI_US =  ( 0.5D0*CF    - ST*ST           ) * 2.D0/HS
      DI_CF =  ( 0.5D0   *US                   ) * 2.D0/HS
      DI_ST =  (             2.D0*ST*(1.D0-US) ) * 2.D0/HS
C
      RETURN
      END
 
 
 
      SUBROUTINE HST( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
           implicit real*8 (a-h,o-z) 
      REAL*8 MSQ
C
C---- Turbulent HS correlation   ( from Swafford profile family )
!s--YE TIAN ----
!     RT = abs(RT)
!e--YE TIAN ----
      
      GRT = DLOG(RT)
      SRT = DSQRT(RT)
      HO    = 3.D0 + 400.D0/RT
      HO_RT =      - 400.D0/RT**2
      IF(RT.LT.400.D0) THEN
       HO    = 4.D0
       HO_RT = 0.D0
      ENDIF
C
      IF(HK.LT.HO) THEN
       HEX = (HO-HK)**(1.6D0)
       RTMP = 0.165D0 - 1.6D0/SRT
       HT    = 1.505D0 + 4.D0/RT + RTMP*HEX/HK
       HT_HK = RTMP*HEX/HK*(-1.6D0/(HO-HK) - 1.D0/HK)
       HT_RT = -4.D0/RT**2 + HEX/HK*0.8D0/SRT/RT
     &             + RTMP*HEX/HK*1.6D0/(HO-HK)*HO_RT
      ELSE
C===========================================================
C----- old correlation
c      HDIF = HK - HO
c      RTMP = HK - HO + 4.D0/GRT
c      HTMP = 0.04D0/HK + 0.007D0*GRT/RTMP**2
c      HT    = 1.505D0 + 4.D0/RT + HDIF**2 * HTMP
c      HT_HK = 2.D0*HDIF*HTMP + HDIF**2*(-.04D0/HK**2 - 0.014D0*GRT/RTMP**3)
c      HT_RT = -4.D0/RT**2
c    &  + HDIF**2*0.014D0/RTMP**2*(1.D0 + 2.D0*GRT/RTMP*4.D0/GRT**2)/RT
c    &  + 2.D0*HDIF*HTMP*(-HO_RT)
c
C===========================================================
C----- new correlation  25 Apr 90
       HDIF = HK - HO 
       RTMP = HK - HO + 4.D0/GRT
       HTMP = 0.015D0/HK + 0.007D0*GRT/RTMP**2
       HT    = 1.505D0 + 4.D0/RT + HDIF**2 * HTMP
       HT_HK = 2.D0*HDIF*HTMP + HDIF**2*( -.015D0/HK**2
     &         - 0.014D0*GRT/RTMP**3 )
       HT_RT = -4.D0/RT**2
     &  + HDIF**2*0.014D0/RTMP**2*(1.D0+2.D0*GRT/RTMP*4.D0/GRT**2)/RT
     &  + 2.D0*HDIF*HTMP*(-HO_RT)
C===========================================================
      ENDIF
C
C---- Whitfield's minor additional compressibility correction
      FM = 1.D0 + 0.014D0*MSQ
      HT     = ( HT + 0.028D0*MSQ ) / FM
      HT_HK  = ( HT_HK            ) / FM
      HT_RT  = ( HT_RT            ) / FM
      HT_MSQ = 0.028D0/FM  -  0.014D0*HT/FM
C
C---- fudge HS slightly to make sure   HS -> 2   as   HK -> 1
      HTF    = 0.485D0/9.D0 * (HK-4.D0)**2/HK  +  1.515D0
      HTF_HK = 0.485D0/9.D0 * (1.D0-16.D0/HK**2)
      ARG = DMAX1( 10.D0*(1.D0 - HK) , -15.D0 )
      HXX = DEXP(ARG)
      HXX_HK = -10.D0*HXX
C
      HS     = (1.D0-HXX)*HT     +  HXX*HTF
      HS_HK  = (1.D0-HXX)*HT_HK  +  HXX*HTF_HK
     &       + (        -HT     +      HTF    )*HXX_HK
      HS_RT  = (1.D0-HXX)*HT_RT
      HS_MSQ = (1.D0-HXX)*HT_MSQ
C
      RETURN
      END
 
 
      SUBROUTINE CFT( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )
      implicit real*8 (a-h,o-z) 
      REAL*8 MSQ
C
C---- Turbulent skin friction function  ( Cf )    (Swafford)
      FC = DSQRT(1.D0 + 0.2D0*MSQ)
      GRT = DLOG(RT/FC)
      GRT = DMAX1(GRT,3.0d0)
      GEX = -1.74D0 - 0.31D0*HK
      ARG = 4.D0 - HK/0.875D0
      ARG = DMIN1( 10.0d0, ARG )
      ARG = DMAX1(-10.0d0, ARG )
      CFO =  0.3D0*DEXP(-1.33D0*HK) * (GRT/2.3026D0)**GEX
      CF     = ( CFO  +  1.1D-4*(DTANH(ARG)-1.D0) ) / FC
      CF_HK  = (-1.33D0*CFO - 0.31D0*DLOG(GRT/2.3026D0)*CFO
     &         - 1.1D-4/DCOSH(ARG)**2 / 0.875D0    ) / FC
      CF_RT  = GEX*CFO/(FC*GRT) / RT
      CF_MSQ = GEX*CFO/(FC*GRT) * (-0.1D0/FC**2)  -  0.1D0*CF/FC**2
C
      RETURN
      END


 
      SUBROUTINE HCT( HK, MSQ, HC, HC_HK, HC_MSQ )
           implicit real*8 (a-h,o-z) 
      REAL*8 MSQ
C
C---- density shape parameter    (from Whitfield)
      HC     = MSQ * (0.064D0/(HK-0.8D0) + 0.251D0)
      HC_HK  = MSQ * (-.064D0/(HK-0.8D0)**2     )
      HC_MSQ =        0.064D0/(HK-0.8D0) + 0.251D0
C
      RETURN
      END
 
 

      SUBROUTINE BLSOLV
C     ............................................................
C      Custom solver for coupled viscous-inviscid Newton system:
C
C        A  |  |  .  |  |  .  |    d       R       S
C        B  A  |  .  |  |  .  |    d       R       S
C        |  B  A  .  |  |  .  |    d       R       S
C        .  .  .  .  |  |  .  |    d   =   R - dRe S
C        |  |  |  B  A  |  .  |    d       R       S
C        |  Z  |  |  B  A  .  |    d       R       S
C        .  .  .  .  .  .  .  |    d       R       S
C        |  |  |  |  |  |  B  A    d       R       S
C
C       A, B, Z  3x3  blocks containing linearized BL equation coefficients
C       |        3x1  vectors containing mass defect influence 
C                     coefficients on Ue
C       d        3x1  unknown vectors (Newton deltas for Ctau, Theta, m)
C       R        3x1  residual vectors
C       S        3x1  Re influence vectors
C     ........................................................
           implicit real*8 (a-h,o-z) 
      INCLUDE 'XFOIL.INC'
C
      IVTE1 = ISYS(IBLTE(1),1)
C
      DO 1000 IV=1, NSYS
C
        IVP = IV + 1
C
C====== Invert VA(IV) block
C
C------ normalize first row
        PIVOT = 1.D0 / VA(1,1,IV)
        VA(1,2,IV) = VA(1,2,IV) * PIVOT
        DO 10 L=IV, NSYS
          VM(1,L,IV) = VM(1,L,IV)*PIVOT
   10   CONTINUE
        VDEL(1,1,IV) = VDEL(1,1,IV)*PIVOT
        VDEL(1,2,IV) = VDEL(1,2,IV)*PIVOT
C
C------ eliminate lower first column in VA block
        DO 15 K=2, 3
          VTMP = VA(K,1,IV)
          VA(K,2,IV) = VA(K,2,IV) - VTMP*VA(1,2,IV)
          DO 150 L=IV, NSYS
            VM(K,L,IV) = VM(K,L,IV) - VTMP*VM(1,L,IV)
  150     CONTINUE
          VDEL(K,1,IV) = VDEL(K,1,IV) - VTMP*VDEL(1,1,IV)
          VDEL(K,2,IV) = VDEL(K,2,IV) - VTMP*VDEL(1,2,IV)
   15   CONTINUE
C
C
C------ normalize second row
        PIVOT = 1.D0 / VA(2,2,IV)
        DO 20 L=IV, NSYS
          VM(2,L,IV) = VM(2,L,IV)*PIVOT
   20   CONTINUE
        VDEL(2,1,IV) = VDEL(2,1,IV)*PIVOT
        VDEL(2,2,IV) = VDEL(2,2,IV)*PIVOT
C
C------ eliminate lower second column in VA block
        K = 3
        VTMP = VA(K,2,IV)
        DO 250 L=IV, NSYS
          VM(K,L,IV) = VM(K,L,IV) - VTMP*VM(2,L,IV)
  250   CONTINUE
        VDEL(K,1,IV) = VDEL(K,1,IV) - VTMP*VDEL(2,1,IV)
        VDEL(K,2,IV) = VDEL(K,2,IV) - VTMP*VDEL(2,2,IV)
C
C
C------ normalize third row
        PIVOT = 1.D0/VM(3,IV,IV)
        DO 350 L=IVP, NSYS
          VM(3,L,IV) = VM(3,L,IV)*PIVOT
  350   CONTINUE
        VDEL(3,1,IV) = VDEL(3,1,IV)*PIVOT
        VDEL(3,2,IV) = VDEL(3,2,IV)*PIVOT
C
C
C------ eliminate upper third column in VA block
        VTMP1 = VM(1,IV,IV)
        VTMP2 = VM(2,IV,IV)
        DO 450 L=IVP, NSYS
          VM(1,L,IV) = VM(1,L,IV) - VTMP1*VM(3,L,IV)
          VM(2,L,IV) = VM(2,L,IV) - VTMP2*VM(3,L,IV)
  450   CONTINUE
        VDEL(1,1,IV) = VDEL(1,1,IV) - VTMP1*VDEL(3,1,IV)
        VDEL(2,1,IV) = VDEL(2,1,IV) - VTMP2*VDEL(3,1,IV)
        VDEL(1,2,IV) = VDEL(1,2,IV) - VTMP1*VDEL(3,2,IV)
        VDEL(2,2,IV) = VDEL(2,2,IV) - VTMP2*VDEL(3,2,IV)
C
C------ eliminate upper second column in VA block
        VTMP = VA(1,2,IV)
        DO 460 L=IVP, NSYS
          VM(1,L,IV) = VM(1,L,IV) - VTMP*VM(2,L,IV)
  460   CONTINUE
        VDEL(1,1,IV) = VDEL(1,1,IV) - VTMP*VDEL(2,1,IV)
        VDEL(1,2,IV) = VDEL(1,2,IV) - VTMP*VDEL(2,2,IV)
C
C
        IF(IV.EQ.NSYS) GO TO 1000
C
C====== Eliminate VB(IV+1) block, rows  1 -> 3
        DO 50 K=1, 3
          VTMP1 = VB(K, 1,IVP)
          VTMP2 = VB(K, 2,IVP)
          VTMP3 = VM(K,IV,IVP)
          DO 510 L=IVP, NSYS
            VM(K,L,IVP) = VM(K,L,IVP)
     &        - (  VTMP1*VM(1,L,IV)
     &           + VTMP2*VM(2,L,IV)
     &           + VTMP3*VM(3,L,IV) )
  510     CONTINUE
          VDEL(K,1,IVP) = VDEL(K,1,IVP)
     &        - (  VTMP1*VDEL(1,1,IV)
     &           + VTMP2*VDEL(2,1,IV)
     &           + VTMP3*VDEL(3,1,IV) )
          VDEL(K,2,IVP) = VDEL(K,2,IVP)
     &        - (  VTMP1*VDEL(1,2,IV)
     &           + VTMP2*VDEL(2,2,IV)
     &           + VTMP3*VDEL(3,2,IV) )
   50   CONTINUE
C
        IF(IV.EQ.IVTE1) THEN
C------- eliminate VZ block
         IVZ = ISYS(IBLTE(2)+1,2)
C
         DO 55 K=1, 3
           VTMP1 = VZ(K,1)
           VTMP2 = VZ(K,2)
           DO 515 L=IVP, NSYS
             VM(K,L,IVZ) = VM(K,L,IVZ)
     &         - (  VTMP1*VM(1,L,IV)
     &            + VTMP2*VM(2,L,IV) )
  515      CONTINUE
           VDEL(K,1,IVZ) = VDEL(K,1,IVZ)
     &         - (  VTMP1*VDEL(1,1,IV)
     &            + VTMP2*VDEL(2,1,IV) )
           VDEL(K,2,IVZ) = VDEL(K,2,IVZ)
     &         - (  VTMP1*VDEL(1,2,IV)
     &            + VTMP2*VDEL(2,2,IV) )
   55    CONTINUE
        ENDIF
C
        IF(IVP.EQ.NSYS) GO TO 1000
C
C====== Eliminate lower VM column
        DO 60 KV=IV+2, NSYS
          VTMP1 = VM(1,IV,KV)
          VTMP2 = VM(2,IV,KV)
          VTMP3 = VM(3,IV,KV)
C
          IF(DABS(VTMP1).GT.VACCEL) THEN
          DO 610 L=IVP, NSYS
            VM(1,L,KV) = VM(1,L,KV) - VTMP1*VM(3,L,IV)
  610     CONTINUE
          VDEL(1,1,KV) = VDEL(1,1,KV) - VTMP1*VDEL(3,1,IV)
          VDEL(1,2,KV) = VDEL(1,2,KV) - VTMP1*VDEL(3,2,IV)
          ENDIF
C
          IF(DABS(VTMP2).GT.VACCEL) THEN
          DO 620 L=IVP, NSYS
            VM(2,L,KV) = VM(2,L,KV) - VTMP2*VM(3,L,IV)
  620     CONTINUE
          VDEL(2,1,KV) = VDEL(2,1,KV) - VTMP2*VDEL(3,1,IV)
          VDEL(2,2,KV) = VDEL(2,2,KV) - VTMP2*VDEL(3,2,IV)
          ENDIF
C
          IF(DABS(VTMP3).GT.VACCEL) THEN
          DO 630 L=IVP, NSYS
            VM(3,L,KV) = VM(3,L,KV) - VTMP3*VM(3,L,IV)
  630     CONTINUE
          VDEL(3,1,KV) = VDEL(3,1,KV) - VTMP3*VDEL(3,1,IV)
          VDEL(3,2,KV) = VDEL(3,2,KV) - VTMP3*VDEL(3,2,IV)
          ENDIF
C
   60   CONTINUE
C
 1000 CONTINUE
C
C
C
      DO 2000 IV=NSYS, 2, -1
C
C------ eliminate upper VM columns
        VTMP = VDEL(3,1,IV)
        DO 81 KV=IV-1, 1, -1
          VDEL(1,1,KV) = VDEL(1,1,KV) - VM(1,IV,KV)*VTMP
          VDEL(2,1,KV) = VDEL(2,1,KV) - VM(2,IV,KV)*VTMP
          VDEL(3,1,KV) = VDEL(3,1,KV) - VM(3,IV,KV)*VTMP
   81   CONTINUE
C
        VTMP = VDEL(3,2,IV)
        DO 82 KV=IV-1, 1, -1
          VDEL(1,2,KV) = VDEL(1,2,KV) - VM(1,IV,KV)*VTMP
          VDEL(2,2,KV) = VDEL(2,2,KV) - VM(2,IV,KV)*VTMP
          VDEL(3,2,KV) = VDEL(3,2,KV) - VM(3,IV,KV)*VTMP
   82   CONTINUE
C
 2000 CONTINUE
C
      RETURN
      END


      SUBROUTINE UPDATE
C------------------------------------------------------------------
C      Adds on Newton deltas to boundary layer variables.
C      Checks for excessive changes and underrelaxes if necessary.
C      Calculates max and rms changes.
C      Also calculates the change in the global variable "AR".
C        If LALFA=.TRUE. , "AR" is the Reynolds number.
C        If LALFA=.FALSE., "AR" is alpha.
C------------------------------------------------------------------
           implicit real*8 (a-h,o-z) 
      INCLUDE 'XBL.INC'
      INCLUDE 'XFOIL.INC'
      REAL*8 UNEW(IVX,2), U_AR(IVX,2)
      REAL*8 QNEW(IQX),   Q_AR(IQX)
      EQUIVALENCE (VA(1,1,1), UNEW(1,1)) ,
     &            (VB(1,1,1), QNEW(1)  )
      EQUIVALENCE (VA(1,1,IVX), U_AR(1,1)) ,
     &            (VB(1,1,IVX), Q_AR(1)  )
C
C---- calculate new Ue distribution assuming no under-relaxation
C-    also set the sensitivity of Ue wrt to alpha or Re
      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
C
          DUI    = 0.D0
          DUI_AR = 0.D0
          DO 100 JS=1, 2
            DO 1000 JBL=2, NBL(JS)
              J  = IPAN(JBL,JS)
              JV = ISYS(JBL,JS)
              UE_M = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
              DUI    = DUI    + UE_M*(MASS(JBL,JS)+VDEL(3,1,JV))
              DUI_AR = DUI_AR + UE_M*(            -VDEL(3,2,JV))
 1000       CONTINUE
  100     CONTINUE
C
C-------- UINV depends on "AR" only if "AR" is alpha
          IF(LALFA) THEN
           UINV_AR = 0.D0
          ELSE
           UINV_AR = UINV_A(IBL,IS)
          ENDIF
C
          UNEW(IBL,IS) = UINV(IBL,IS) + DUI
          U_AR(IBL,IS) = UINV_AR      + DUI_AR
C
   10   CONTINUE
    1 CONTINUE
C
C---- set new Qtan from new Ue with appropriate sign change
      DO 2 IS=1, 2
        DO 20 IBL=2, IBLTE(IS)
          I = IPAN(IBL,IS)
          QNEW(I) = VTI(IBL,IS)*UNEW(IBL,IS)
          Q_AR(I) = VTI(IBL,IS)*U_AR(IBL,IS)
   20   CONTINUE
    2 CONTINUE
C
C---- calculate new CL from this new Qtan
      SA = DSIN(ARAD)
      CA = DCOS(ARAD)
C
      BETA = DSQRT(1.D0 - MINF**2)
      BFAC = 0.5D0*MINF**2 / (1.D0 + BETA)
C
      CLNEW = 0.D0
      CL_A  = 0.D0
      CL_AR = 0.D0
      QINF = 1.0
C
      I = 1
      CGINC = 1.D0 - (QNEW(I)/QINF)**2
      CPG1 = CGINC / (BETA + BFAC*CGINC)
C
      CPI_Q = -2.D0*QNEW(I)/QINF**2
      CPC_CPI = (1.D0 - BFAC*CPG1)/ (BETA + BFAC*CGINC)
      CPG1_AR = CPC_CPI*CPI_Q*Q_AR(I)
C
      DO 3 I=1, Nvis
        IP = I+1
        IF(I.EQ.Nvis) IP = 1
C
        CGINC = 1.D0 - (QNEW(IP)/QINF)**2
        CPG2 = CGINC / (BETA + BFAC*CGINC)
C
        CPI_Q = -2.D0*QNEW(IP)/QINF**2
        CPC_CPI = (1.D0 - BFAC*CPG1)/ (BETA + BFAC*CGINC)
        CPG2_AR = CPC_CPI*CPI_Q*Q_AR(IP)
C
        DX   =  (X(IP) - X(I))*CA + (Y(IP) - Y(I))*SA
        DX_A = -(X(IP) - X(I))*SA + (Y(IP) - Y(I))*CA
C
        AG    = 0.5D0*(CPG2    + CPG1   )
        AG_AR = 0.5D0*(CPG2_AR + CPG1_AR)
C
        CLNEW = CLNEW + DX  *AG
        CL_A  = CL_A  + DX_A*AG
        CL_AR = CL_AR + DX  *AG_AR
C
        CPG1    = CPG2
        CPG1_AR = CPG2_AR
    3 CONTINUE
C
      IF(LALFA) THEN
C===== alpha is prescribed: AR is Re
C      account for change in Reynolds number due to CL change
C
C----- calculate derivative of CL wrt Reynolds number Re
       CLOLD = CL
       IF(RETYP.EQ.1) THEN
C------ Re = Reinf
        RE_CL = 0.D0
       ELSE IF(RETYP.EQ.2) THEN
C------ Re = Reinf/dsqrt(CL)
        RE_CL = -.5D0*REYBL/CLOLD
       ELSE IF(RETYP.EQ.3) THEN
C------ Re = Reinf/CL
        RE_CL = -REYBL/CLOLD
       ENDIF
C
C----- set change in Re to account for CL changing, since Re = Re(CL)
       DAR = RE_CL * (CLNEW-CLOLD) / (1.D0 - RE_CL*CL_AR)
C
      ELSE
C===== CL is prescribed: AR is alpha
C
C----- set correct sensitivity dCL/dalpha
       CL_AR = CL_AR + CL_A
C
C----- set change in alpha to drive CL to prescribed value
       DAR = (CLSPEC - CLNEW) / CL_AR
C
      ENDIF
ccc  let DAR = 0   to conform to iblint
      DAR=0.D0
C
      RMSBL = 0.D0
      RMXBL = 0.D0
C
      DHI = 1.5D0
      DLO = -.5D0
C
C---- calculate changes in BL variables and under-relaxation factor if necessary
      RLX = 1.D0
      RLX = RLX*BETA_BL
      DO 4 IS=1, 2
        DO 40 IBL=2, NBL(IS)
          IV = ISYS(IBL,IS)
C
C-------- set changes without underrelaxation
          DCTAU = VDEL(1,1,IV) - DAR*VDEL(1,2,IV)
          DTHET = VDEL(2,1,IV) - DAR*VDEL(2,2,IV)
          DMASS = VDEL(3,1,IV) - DAR*VDEL(3,2,IV)
          DUEDG = UNEW(IBL,IS) + DAR*U_AR(IBL,IS)  -  UEDG(IBL,IS)
          DDSTR = (DMASS - DSTR(IBL,IS)*DUEDG)/UEDG(IBL,IS)

C
C-------- normalize changes
          IF(IBL.LT.ITRAN(IS)) DN1 = DCTAU / 10.D0
          IF(IBL.GE.ITRAN(IS)) DN1 = DCTAU / CTAU(IBL,IS)
          DN2 = DTHET / THET(IBL,IS)
          DN3 = DDSTR / DSTR(IBL,IS)
          DN4 = DABS(DUEDG)/0.25D0
C
C-------- accumulate for rms change
          RMSBL = RMSBL + DN1**2 + DN2**2 + DN3**2 + DN4**2
C          
C-------- see if Ctau needs underrelaxation
          RDN1 = RLX*DN1
          IF(DABS(DN1) .GT. DABS(RMXBL)) THEN
           RMXBL = DN1
           IF(IBL.LT.ITRAN(IS)) VMXBL = 'n'
           IF(IBL.GE.ITRAN(IS)) VMXBL = 'C'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN1 .GT. DHI) then
            RLX = DHI/DN1
            RLX = RLX*BETA_BL
          end if
          IF(RDN1 .LT. DLO) then
            RLX = DLO/DN1
            RLX = RLX*BETA_BL
          end if
C
C-------- see if Theta needs underrelaxation
          RDN2 = RLX*DN2
          IF(DABS(DN2) .GT. DABS(RMXBL)) THEN
           RMXBL = DN2
           VMXBL = 'T'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN2 .GT. DHI) then
            RLX = DHI/DN2
            RLX = RLX*BETA_BL
          end if
          IF(RDN2 .LT. DLO) then
            RLX = DLO/DN2
            RLX = RLX*BETA_BL
          end if
C
C-------- see if Dstar needs underrelaxation
          RDN3 = RLX*DN3
          IF(DABS(DN3) .GT. DABS(RMXBL)) THEN
           RMXBL = DN3
           VMXBL = 'D'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN3 .GT. DHI) then
            RLX = DHI/DN3
            RLX = RLX*BETA_BL
          end if
          IF(RDN3 .LT. DLO) then
            RLX = DLO/DN3
            RLX = RLX*BETA_BL
          end if
C
C-------- see if Ue needs underrelaxation
          RDN4 = RLX*DN4
          IF(DABS(DN4) .GT. DABS(RMXBL)) THEN
           RMXBL = DUEDG
           VMXBL = 'U'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN4 .GT. DHI) then
            RLX = DHI/DN4
            RLX = RLX*BETA_BL
          end if
          IF(RDN4 .LT. DLO) then
            RLX = DLO/DN4
            RLX = RLX*BETA_BL
          end if
C
   40   CONTINUE
    4 CONTINUE
C
C---- set true rms change
      RMSBL = DSQRT( RMSBL / (4.D0*DBLE( NBL(1)+NBL(2) )) )
C
C
      IF(LALFA) THEN
C----- set underrelaxed change in Reynolds number from change in lift
       REYBL = REYBL + RLX*DAR
      ELSE
C----- set underrelaxed change in alpha
       ARAD = ARAD + RLX*DAR
       ADEG = ARAD/DTOR
      ENDIF
C
C---- update BL variables with underrelaxed changes
      DO 5 IS=1, 2
        DO 50 IBL=2, NBL(IS)
          IV = ISYS(IBL,IS)
C
          DCTAU = VDEL(1,1,IV) - DAR*VDEL(1,2,IV)
          DTHET = VDEL(2,1,IV) - DAR*VDEL(2,2,IV)
          DMASS = VDEL(3,1,IV) - DAR*VDEL(3,2,IV)
          DUEDG = UNEW(IBL,IS) + DAR*U_AR(IBL,IS)  -  UEDG(IBL,IS)
          DDSTR = (DMASS - DSTR(IBL,IS)*DUEDG)/UEDG(IBL,IS)
C
          CTAU(IBL,IS) = CTAU(IBL,IS) + RLX*DCTAU
          THET(IBL,IS) = THET(IBL,IS) + RLX*DTHET
          DSTR(IBL,IS) = DSTR(IBL,IS) + RLX*DDSTR
          UEDG(IBL,IS) = UEDG(IBL,IS) + RLX*DUEDG
C
C-------- eliminate absurd transients
          IF(IBL.LE.IBLTE(IS))
     &      DSTR(IBL,IS) = DMAX1(DSTR(IBL,IS) , 1.05000D0*THET(IBL,IS))
          IF(IBL.GT.IBLTE(IS))
     &      DSTR(IBL,IS) = DMAX1(DSTR(IBL,IS) , 1.00005D0*THET(IBL,IS))
          IF(IBL.GE.ITRAN(IS))
     &      CTAU(IBL,IS) = DMIN1( CTAU(IBL,IS) , 0.20d0 )
C
C-------- set new mass defect (nonlinear update)
          MASS(IBL,IS) = DSTR(IBL,IS) * UEDG(IBL,IS)
C
   50   CONTINUE
    5 CONTINUE
C
      RETURN
      END !UPDATE


      FUNCTION SEVAL(SS,X,XS,S,N)
      implicit real*8 (a-h,o-z)
      DIMENSION X(N), XS(N), S(N)
C--------------------------------------------------
C     Calculates X(SS)                             |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      ILOW = 1
      I = N
C
   10 IF((I-ILOW) .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      SEVAL = T*X(I) + (1.D0-T)*X(I-1) + (T-T*T)*((1.D0-T)*CX1 - T*CX2)
      RETURN
      END ! SEVAL



      SUBROUTINE SINVRT(SI,XI,X,XS,S,N)
C-------------------------------------------------------
C     Calculates the "inverse" spline function S(X).    |
C     Since S(X) can be multi-valued or not defined,    |
C     this is not a "black-box" routine.  The calling   |
C     program must pass via SI a sufficiently good      |
C     initial guess for S(XI).                          |
C                                                       |
C     XI      specified X value       (input)           |
C     SI      calculated S(XI) value  (input,output)    |
C     X,XS,S  usual spline arrays     (input)           |
C                                                       |
C-------------------------------------------------------
C
      implicit real*8 (a-h,o-z)
      DIMENSION X(N),XS(N),S(N)
      SISAV = SI
C
      DO 10 ITER=1, 10
        RES  = SEVAL(SI,X,XS,S,N) - XI
        RESP = DEVAL(SI,X,XS,S,N)
        DS = -RES/RESP
        SI = SI + DS
        IF(DABS(DS).LT.1.0D-5) RETURN
   10 CONTINUE
      WRITE(*,*)
     &  'SINVRT: spline inversion failed. Input value returned.'
      SI = SISAV
C
      RETURN
      END ! SINVRT







      SUBROUTINE GAUSS(NSIZ,NN,Z,R,NRHS)
C     *******************************************************
C     *                                                     *
C     *   Solves general NxN system in NN unknowns          *
C     *    with arbitrary number (NRHS) of righthand sides. *
C     *   Assumes system is invertible...                   *
C     *    ...if it isn't, a divide by zero will result.    *
C     *                                                     *
C     *   Z is the coefficient matrix...                    *
C     *     ...destroyed during solution process.           *
C     *   R is the righthand side(s)...                     *
C     *     ...replaced by the solution vector(s).          *
C     *                                                     *
C     *                              Mark Drela  1984       *
C     *******************************************************
C
      implicit real*8 (a-h,o-z)
      DIMENSION Z(NSIZ,NSIZ), R(NSIZ,NRHS)
C
      DO 1 NP=1, NN-1
        NP1 = NP+1
C
C------ find max pivot index NX
        NX = NP
        DO 11 N=NP1, NN
          IF(DABS(Z(N,NP))-DABS(Z(NX,NP))) 11,11,111
  111      NX = N
   11   CONTINUE
C
        PIVOT = 1.D0/Z(NX,NP)
C
C------ switch pivots
        Z(NX,NP) = Z(NP,NP)
C
C------ switch rows & normalize pivot row
        DO 12 L=NP1, NN
          TEMP = Z(NX,L)*PIVOT
          Z(NX,L) = Z(NP,L)
          Z(NP,L) = TEMP
   12   CONTINUE
C
        DO 13 L=1, NRHS
          TEMP = R(NX,L)*PIVOT
          R(NX,L) = R(NP,L)
          R(NP,L) = TEMP
   13   CONTINUE
C
C------ forward eliminate everything
        DO 15 K=NP1, NN
          ZTMP = Z(K,NP)
C
C          IF(ZTMP.EQ.0.0) GO TO 15
C
          DO 151 L=NP1, NN
            Z(K,L) = Z(K,L) - ZTMP*Z(NP,L)
  151     CONTINUE
          DO 152 L=1, NRHS
            R(K,L) = R(K,L) - ZTMP*R(NP,L)
  152     CONTINUE
   15   CONTINUE
C
    1 CONTINUE
C
C---- solve for last row
      DO 2 L=1, NRHS
        R(NN,L) = R(NN,L)/Z(NN,NN)
    2 CONTINUE
C
C---- back substitute everything
      DO 3 NP=NN-1, 1, -1
        NP1 = NP+1
        DO 31 L=1, NRHS
          DO 310 K=NP1, NN
            R(NP,L) = R(NP,L) - Z(NP,K)*R(K,L)
  310     CONTINUE
   31   CONTINUE
    3 CONTINUE
C
      RETURN
      END ! GAUSS



      FUNCTION DEVAL(SS,X,XS,S,N)
      implicit real*8 (a-h,o-z)
      DIMENSION X(N), XS(N), S(N)
C--------------------------------------------------
C     Calculates dX/dS(SS)                         |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      ILOW = 1
      I = N
C
   10 IF((I-ILOW) .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      DEVAL = X(I) - X(I-1) + (1.D0-4.D0*T+3.D0*T*T)*CX1 
     &      + T*(3.D0*T-2.D0)*CX2
      DEVAL = DEVAL/DS
      RETURN
      END ! DEVAL
C
C
C    *************   SUBROUTINE UESET  ************ 
      SUBROUTINE UESET
C---------------------------------------------------------
C     Sets Ue from inviscid Ue plus all source influence
C---------------------------------------------------------
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
C
      DO 1 IS=1, 2
c        write(71,*)'is =',is
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
c            write(71,*)'i =',i
C
          DUI = 0.D0
          DO 100 JS=1, 2
            DO 1000 JBL=2, NBL(JS)
              J  = IPAN(JBL,JS)
              UE_M = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
              DUI = DUI + UE_M*MASS(JBL,JS)
c            write(71,*)jbl,UE_M,dui,   mass(jbl,js)
 1000       CONTINUE
  100     CONTINUE
C
          UEDG(IBL,IS) = UINV(IBL,IS) + DUI
C
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END  ! UESET
