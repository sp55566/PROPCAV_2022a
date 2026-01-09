      PROGRAM PROPCAV
************************************************************************
*     PROPCAV MAIN PROGRAM                                             *
*     Unsteady cavitating propeller analysis program                   *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) MASSACHUSETTS INSTITUTE OF TECHNOLOGY              *
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*     MAY 1992                                                         *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*     Version 1.0   November, 1998                                     *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*     Version 1.1   October, 1999                                      *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*     Version 1.2   January, 2001                                      *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*     Version 2.0   December, 2001                                     *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*     Version 2.1   May, 2003                                          *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*     Version 2.1.1 January, 2004                                      *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*     Version 2.3   January, 2007                                      *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*     Version 3.0   May, 2010                                          *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*     Version 3.12  August, 2013                                       *
*----------------------------------------------------------------------*
*     initiated 5-1-92 by Neal Fine using PUF-10 as a base and taking  *
*     cavity bvp subroutines from PSFCAV.                              *
*     Date     Revision or Comment                                     *
*----------------------------------------------------------------------*
*     PROPCAV  v3.3 Release                                            *
*     Unsteady cavitating propeller analysis program                   *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*        October  2016  P.I. Dr. Spyros Kinnas                         *
*----------------------------------------------------------------------*
*     PROPCAV  R2018a Release                                          *
*     Unsteady cavitating propeller analysis program                   *
*----------------------------------------------------------------------*
*     PROPCAV  R2019 Release                                           *
*     Unsteady cavitating propeller analysis program                   *
*----------------------------------------------------------------------*
*     PROPCAV  R2022 Release                                           *
*     Unsteady cavitating propeller analysis program                   *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*         November  2022  P.I. Dr. Spyros Kinnas                       *
*----------------------------------------------------------------------*
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      CHARACTER*29 FNLOG,FNPLT

!     REAL*4 ETIME, ELAPSED,TARRAY(2)!,DTIME
      REAL*4 ELAPSED,TARRAY(2)!,DTIME
      DIMENSION AFM(20,6)
!     EXTERNAL ETIME
!     EXTERNAL DTIME

C-----------------------------------------------------------------------
C     Release version
C-----------------------------------------------------------------------
      WRITE(*,'(A)') '  '
      WRITE(*,'(A)') '  **********************************************'
      WRITE(*,'(A)') '  *                                            *'
      WRITE(*,'(A)') '  *      PROPCAV RELEASE VERSION R2022         *'
      WRITE(*,'(A)') '  *        Released  November, 2022            *'
      WRITE(*,'(A)') '  *                                            *'
      WRITE(*,'(A)') '  **********************************************'
      WRITE(*,'(A)') '  '

C-----------------------------------------------------------------------
C     CONTROL FOR VISCOUS BOUNDARY LAYER ANALYSIS  BY HONG SUN  Jan.,07
C-----------------------------------------------------------------------

      WRITE(*,'(A,$)')
     *  ' PROPCAV> WANT TO INCLUDE VISCOUS RUN? (IVISC=1:YES) '
      WRITE(*,*)
      READ(*,*) IVISC

C-----------------------------------------------------------------------
C     Read and process input data
C-----------------------------------------------------------------------
      WRITE(*,'(A,$)') ' PROPCAV> ENTER DELTAT (in degrees): '
      READ(*,*) DELTAT
      NDLTAT=INT(DELTAT)
      DELTAT=DELTAT*RAD

!s--YE TIAN  for m_param
      NSTEP = 360/NDLTAT
      NNDIM = NSTEP
      NDEL  = NDLTAT
!e--YE TIAN  for m_param

C-----------------------------------------------------------------------
C     What is the cut off for extrapolation?                   JY081600
C-----------------------------------------------------------------------
      WRITE(*,*)
      WRITE(*,'(A,$)') ' PROPCAV> ENTER RADII TO CUT THE TIP (RMRTIP): '
      WRITE(*,*)
      READ(*,*) RMRTIP

C---------This subroutine reads in the propeller geometry---------------

      CALL PROINP

!s--- YE TIAN 07/09/2013----
      NZWSUB=NWZ+NDEL
!e--- YE TIAN 07/09/2013----

!     MODULE m_FOUT
      ALLOCATE(TT(NSTEP),XKT(NSTEP,6),XKTV(NSTEP,6))


C     MODULE GPROD
      ALLOCATE(XB(NBPZ,MBPZ),YB(NBPZ,MBPZ),ZB(NBPZ,MBPZ))
      ALLOCATE(XW(NWPZ,MBPZ),YW(NWPZ,MBPZ),ZW(NWPZ,MBPZ))
      ALLOCATE(XX_TMP(NWPZ,MBPZ),YY_TMP(NWPZ,MBPZ),ZZ_TMP(NWPZ,MBPZ))
c XM YU 12/2011
      allocate(xwvs3d(nwpz,mbpz),ywvs3d(nwpz,mbpz),zwvs3d(nwpz,mbpz))
C/s S.N.KIM | Aug. 2018.
      allocate(xwvs3do(nwpz,mbpz),ywvs3do(nwpz,mbpz),zwvs3do(nwpz,mbpz))
C/e S.N.KIM 
c XM YU 12/2011
      ALLOCATE(XH(NHPZ,MHPZ),YH(NHPZ,MHPZ),ZH(NHPZ,MHPZ))
      ALLOCATE(XW1(NWPZ,MBPZ),YW1(NWPZ,MBPZ),ZW1(NWPZ,MBPZ))
      ALLOCATE(THT(MBPZ),THR(NBPZ),SB(NBHZ))
      ALLOCATE(YTC(NBHZ),YCC(NBHZ),TCC(NBHZ))
      ALLOCATE(XI(NBPZ),ETA(NBPZ))
      ALLOCATE(CONMRX(3,NBPZ,MBPZ),VECMRX(3,NBPZ,MBPZ),
     &         WAKVEC(3,NWPZ,MBPZ))
      ALLOCATE(ARCLNG(NSCZP,MBZ,2))

C     MODULE GPRODW
      ALLOCATE(XWO(NWPZ,MBPZ,KZ),YWO(NWPZ,MBPZ,KZ),ZWO(NWPZ,MBPZ,KZ))
      ALLOCATE(XWW(NWPZ,MBPZ,KZ),YWW(NWPZ,MBPZ,KZ),ZWW(NWPZ,MBPZ,KZ))
C/s S.N.KIM | Aug. 2018.
      ALLOCATE(XWW1(NWPZ,MBPZ,KZ),YWW1(NWPZ,MBPZ,KZ),ZWW1(NWPZ,MBPZ,KZ))
C/e S.N.KIM

C     MODULE GINT
      ALLOCATE(RZ(MBPZ),RZSQ(MBPZ),PITCH(MBPZ),RAKE(MBPZ))
      ALLOCATE(SKEW(MBPZ),CHORD(MBPZ),CAMBR(MBPZ),THICK(MBPZ))

C     MODULE FVEL
      ALLOCATE(VOXW(NSCWZ),VOYW(NSCWZ),VOZW(NSCWZ))
      ALLOCATE(VTOTS(NPANZ),VXIB(NBZ,MBZ),VETAB(NBZ,MBZ))
      ALLOCATE(VINFSB(NBZ,MBZ),VXIH(NHBZ,MHBZ),VETAH(NHBZ,MHBZ))
      ALLOCATE(VINFSH(NHBZ,MHBZ),UXTOT(NBZ,MBZ),UYTOT(NBZ,MBZ))
      ALLOCATE(UZTOT(NBZ,MBZ),UXHTOT(NHBZ,MHBZ),UYHTOT(NHBZ,MHBZ))
      ALLOCATE(UZHTOT(NHBZ,MHBZ))
      ALLOCATE(VOX1(NPANZ),VOY1(NPANZ),VOZ1(NPANZ))

C     MODULE FVEL2
      ALLOCATE(UXTHTOT(NHMX,MCAVM),UYTHTOT(NHMX,MCAVM))
      ALLOCATE(UZTHTOT(NHMX,MCAVM),UXCHTOT(NCAVM,MCAVM))
      ALLOCATE(UYCHTOT(NCAVM,MCAVM),UZCHTOT(NCAVM,MCAVM))

C     MODULE GEOMT
      ALLOCATE(XVP(NPANZ,4),YVP(NPANZ,4),XCT(NPANZ,3),DIR(NPANZ,3,3))
      ALLOCATE(SS(NPANZ,15),SID(NPANZ,4),GSX(NPANZ,4,3))
      ALLOCATE(GSW(NPANZ,4),VEL(NPANZ,6),DELU(NPANZ),DELV(NPANZ))
      ALLOCATE(COSPHI(NPANZ),SINPHI(NPANZ))

C     MODULE GEOMTP
      ALLOCATE(XCTP(NPANZ,3,KZ),DIRCP(NPANZ,3,3,KZ))
      ALLOCATE(XCTPW(NSCWZ,3,KZ),DIRWCP(NSCWZ,3,3,KZ))
      ALLOCATE(XCTPWs(NSCWZ,3,KZ),DIRWsCP(NSCWZ,3,3,KZ))
c XM YU
      allocate(xctpwv(nscwz,3,kz))
      ALLOCATE(XCTPDW(NDWMN,3,KZ),DIRDWCP(NDWMN,3,3,KZ))
      ALLOCATE(XCTPDWs(NSCWZ,3,KZ),DIRDWsCP(NSCWZ,3,3,KZ))
      ALLOCATE(XCPW(NSCWZ,3,KZ))

C     MODULE GEOMTW
      ALLOCATE(XVPW(NSCWZ,4),YVPW(NSCWZ,4),XCTW(NSCWZ,3))
      ALLOCATE(DIRW(NSCWZ,3,3),SSW(NSCWZ,15),SIDW(NSCWZ,4))
      ALLOCATE(GSXW(NSCWZ,4,3),GSWW(NSCWZ,4),VELW(NSCWZ,6))
      ALLOCATE(DELUW(NSCWZ),DELVW(NSCWZ),COSPHIW(NSCWZ))
      ALLOCATE(SINPHIW(NSCWZ))

C     MODULE GEOMTWD
      ALLOCATE(XVPWD(NDWMN,4),YVPWD(NDWMN,4),XCTWD(NDWMN,3))
      ALLOCATE(DIRWD(NDWMN,3,3),SSWD(NDWMN,15),SIDWD(NDWMN,4))
      ALLOCATE(GSXWD(NDWMN,4,3),GSWWD(NDWMN,4),VELWD(NDWMN,6))
      ALLOCATE(DELUWD(NDWMN),DELVWD(NDWMN),COSPHIWD(NDWMN))
      ALLOCATE(SINPHIWD(NDWMN))

C     MODULE GEOMTWO
      ALLOCATE(XVPWO(NSCWZ,4,KZ),YVPWO(NSCWZ,4,KZ),XCTWO(NSCWZ,3,KZ))
      ALLOCATE(DIRWO(NSCWZ,3,3,KZ),SSWO(NSCWZ,15,KZ),SIDWO(NSCWZ,4,KZ))
      ALLOCATE(GSXWO(NSCWZ,4,3,KZ),GSWWO(NSCWZ,4,KZ),VELWO(NSCWZ,6,KZ))
      ALLOCATE(DELUWO(NSCWZ,KZ),DELVWO(NSCWZ,KZ),COSPHIWO(NSCWZ,KZ))
      ALLOCATE(SINPHIWO(NSCWZ,KZ))

C     MODULE TOTG
      ALLOCATE(XG(NPANZ,4,3),CHRLEPS(NPANZ),CHRLEWS(NSCWZ))
      ALLOCATE(XGW(NSCWZ,4,3),chrlewso(NSCWZ,KZ))
      ALLOCATE(XGWD(NDWMN,4,3),CHRLEWSD(NDWMN))

      ALLOCATE(XG1(NPANZ,4,3),XG11(NPANZ,4,3),R111(NPANZ))

C     MODULE SIMQD1
      ALLOCATE(A(NPANZ),B(NPANZ),W(NPANZ,MBZ))
      ALLOCATE(CHDK(NPANZ),WDK(NPANZ),POT(NPANZ),DPDN(NPANZ))
      ALLOCATE(winf(npanz,mbz),CHDKDT(NPANZ),WD(NPANZ,MDMAX))
      ALLOCATE(WKD(NPAWZ,MDMAX))

C     MODULE SIMQD2
      ALLOCATE(AA(NTZ,NTZ),BB(NPANZ,NPANZ),C(NPANZ,NPAWZ))
      ALLOCATE(D(NPAWZ,NPANZ),E(NPAWZ,NPANZ),F(NPAWZ,NPAWZ))
      ALLOCATE(WK(NPAWZ,MBZ),WK2(NPAWZ,MBZ))
      ALLOCATE(WDK1(NPAWZ),CHDK1(NPAWZ),CHDKDT1(NPAWZ))

C     MODULE SIMQ1
      ALLOCATE(NPERB(NBLKMAX))

C     MODULE CORR
      ALLOCATE(DELCP(MBZ),DCPPRE(MBZ))

C     MODULE CIRC
! YE TIAN 04/08/2012
!     ALLOCATE(DELP(MBZ),DELPD(MDMAX))
      ALLOCATE(DELP(MBZ),DELPD(MDMAX),DELPonW(MBZ))
      ALLOCATE(DELPonW2(MBZ))
C/s S.N.KIM | ductwake alignment using FWA | Aug. 2018
      ALLOCATE(DELPonDW(MDMAX))
C/e S.N.KIM | ductwake alignment using FWA | Aug. 2018

C     MODULE PRES
      ALLOCATE(CPB(NBZ,MBZ),CPH(NHBZ,MHBZ),DPDUB(NBZ,MBZ))
      ALLOCATE(DPDVB(NBZ,MBZ),DPDUH(NHBZ,MHBZ),DPDVH(NHBZ,MHBZ))

C     MODULE NORMAL
      ALLOCATE(XCON(NBHZP,MBPZ),YCON(NBHZP,MBPZ),ZCON(NBHZP,MBPZ))

C     MODULE FORCE
      ALLOCATE(BKFX(MBZ),BKFY(MBZ),BKFZ(MBZ),BKMX(MBZ),BKMY(MBZ))
      ALLOCATE(BKMZ(MBZ),HKFX(NHBZ),HKFY(NHBZ),HKFZ(NHBZ))
      ALLOCATE(HKMX(NHBZ),HKMY(NHBZ),HKMZ(NHBZ))

C     MODULE VFOR
      ALLOCATE(FVA(MBZ))
C/s S. Kim | turbine, local friction coefficients
      ALLOCATE(VRR(NBZ,MBZ),RER(NBZ,MBZ))
      ALLOCATE(XCDF1(NBZ,MBZ))
C/e S. Kim | turbine, local friction coefficients

C     MODULE CAVGEOD
      ALLOCATE(DS(NBSCZ,MBZ),DZ(NBSCZ,MBZ))
      ALLOCATE(SZ(NBSCZP,MBZ),SPZ(NBSCZ,MBZ))
      ALLOCATE(RLAM(MBZ,2),CAVL(MBZ,2),CAVLP(MBZ,2),DZW(NWZ,MBZ))
      ALLOCATE(SOP(MBZ,2),CAVLSA(MBZ,NSTEP,2))

C     MODULE CAVGEOI
      ALLOCATE(M0(MBZ,2),JCV(MBZ,2),LCV(MBZ,2))
      ALLOCATE(NOCAV(MBZ,0:NSTEP,2))

C     MODULE CAVST
      ALLOCATE(QC(NSCZP,MBZ,2),PHI0(MBZ,2))
      ALLOCATE(PHI1(NSCZP,MBZ,2),DPHIDS(NSCZP,MBZ,2))
      ALLOCATE(DPHIDV(NBHZ,MBZ,2))

C     MODULE CVSOL2
      ALLOCATE(HT(NSCZP,MBZ,2),DPDNC(NPANZ),DELTA(MBZ,2))
      ALLOCATE(HTW1(MBZ),DELTAP(MBZ,2),SOL(NTZ),SORW(NPAWZ))
      ALLOCATE(POTW(NPAWZ),HTP(NBHZP,MBZ,2),HTWP(NWZ+1,MBZ))

C     MODULE SPLITD
      ALLOCATE(FLP(MBZ,2),FRP(MBZ,2),FLS(MBZ,2),FRS(MBZ,2))
      ALLOCATE(DZL(MBZ,2),DZR(MBZ,2),QSPR(MBZ,2),DLISP(MBZ,2))
      ALLOCATE(DPDVSP(MBZ,2),PHIL(MBZ,2),PHIR(MBZ,2),QSPL(MBZ,2))
      ALLOCATE(QSSR(MBZ))


C     MODULE NWAK
      ALLOCATE(NWC(MBZ,2),NWPAN(MBZ))

C     MODULE UVECB
      ALLOCATE(UL(NPANZ,3),VL(NPANZ,3),WL(NPANZ,3))

C     MODULE UVECBo
      ALLOCATE(ULo(NPANZ,3,KZ),VLo(NPANZ,3,KZ))
      ALLOCATE(WLo(NPANZ,3,KZ))

C     MODULE UVECW
      ALLOCATE(ULW(NSCWZ,3))

C     MODULE XTRAP
      ALLOCATE(CT(MBZ,4,2),DT(MBZ,4,2),QT(MBZ,4,2),QW(MBZ,4,2))

C     MODULE UNSWAKD
      ALLOCATE(WSTINF(NPANZ,MBZ),WSTINFD(NPANZ,MDMAX))

C     MODULE UNSDP
      ALLOCATE(DPHI(MBZ,0:NNDIM))

C     MODULE UNS
      ALLOCATE(POTM(NNDIM,NPANZ),DPDTPRE(NNDIM,NPANZ))
      ALLOCATE(POTWM(NNDIM,NPAWZ),DPDTPREW(NNDIM,NPAWZ))

C     MODULE TMPP1
      ALLOCATE(TEMP1(NPANZ),TEMP2(NPANZ),TEMP3(NPAWZ))
      ALLOCATE(TEMP4(NSCWZ),TEMP5(NPAWZ+MBZ))
      ALLOCATE(STRGTH(NPANZ,KZ),WUSINF(NPANZ,MBZ))
      ALLOCATE(STRGTH1(NSCWZ,KZ),WUSINF1(NSCWZ,MBZ))
      ALLOCATE(TEMPD4(NDWMN),STRGTHD(NPANZ,KZ))

C     MODULE CAVFB
      ALLOCATE(NNWC(MBZ),NWDIR(MBZ))

C     MODULE WAKSUBD
      ALLOCATE(XWS(NZWSUB,MBPZ),YWS(NZWSUB,MBPZ))
      ALLOCATE(ZWS(NZWSUB,MBPZ),WSUBIF(NPANZ,MBZ))
      ALLOCATE(WKFACE(NPAWZ,MBZ,NZWSUB),WSUBIF_IM(NPANZ,MBZ))
C/s S.N.KIM | Aug. 2018.
      ALLOCATE(XWSO(NZWSUB,MBPZ,10),YWSO(NZWSUB,MBPZ,10))
      ALLOCATE(ZWSO(NZWSUB,MBPZ,10))
C/e S.N.KIM

C     MODULE DWAKSUB
      ALLOCATE(XWSD(NZWSUB,MDMAXP),YWSD(NZWSUB,MDMAXP))
      ALLOCATE(ZWSD(NZWSUB,MDMAXP),WSUBIFD(NPANZ,MDMAX))

C     MODULE MEMSOL
      ALLOCATE(ALHS(NTZ,NTZ))

C     MODULE CVRHS
      ALLOCATE(RHS(NTZ))

C     MODULE PRESNO
      ALLOCATE(CPBN(NBZ,MBZ),CPHN(NHBZ,MHBZ))
      ALLOCATE(CPTNN(NTMAXT,NCMAXT),CPDN(NDMAX,MDMAX))
      ALLOCATE(CPTHN(NHMX,MCAVM),CPCN(NCAVM,MCAVM))

C     MODULE PODGEO1,PODGEO3,PODGEO4
      ALLOCATE(HRZ(NBPZ,MBPZ),HRZSQ(NBPZ,MBPZ))
      ALLOCATE(HRZP(NBZ,MBZ),HRZPSQ(NBZ,MBZ))

C     MODULE BKUTTA
      ALLOCATE(POTEMP(NTZ),WWK(NTZ,MBZ))

C     MODULE TOGBFLOW
      ALLOCATE(XCCC(NBHZP,MBPZ),YCCC(NBHZP,MBPZ))
      ALLOCATE(ZCCC(NBHZP,MBPZ),CPBMEAN(NSTEP,NBHZP,MBPZ))

C     MODULE GLOCBL
      ALLOCATE( XIV(NBPZ,MBPZ),ETAV(NBPZ,MBPZ) )

      IF(ISP.EQ.1) THEN

C     MODULE SPINT
        ALLOCATE(ICB(MBZ,2,NSTEP),MSW(MBZ,NSTEP),ICW(MBZ,NSTEP))
        ALLOCATE(IC(2,MBZ,NSTEP),IW(2,MBZ,NSTEP))
        ALLOCATE(ISUBM(NPBZ+NPHZ,NSTEP),ISUBH(NHBZ,NSTEP))
        ALLOCATE(NPERM(MBZ),ISUB1(MBZ,2,NSTEP),ISUWM(NPAWZ,NSTEP))
C     MODULE SPIMG
        ALLOCATE(XCTP_IM(NPBZ+NPHZ,3,KZ),XCPW_IM(NPAWZ,3,KZ))
      ENDIF

      IF(ISC.EQ.1) THEN
C     MODULE FXSR,CAVSR
        ALLOCATE(SSTE(MBZ),VELTE(MBZ,6))
        ALLOCATE(CT2(MBZ,4,2),QCSR(NCSRZP,MBZ,2),PHI0T(MBZ,2))
        ALLOCATE(PSI0T(MBZ,2),PHI2(NCSRZ,MBZ,2))
        ALLOCATE(DPHIDS2(NCSRZ,MBZ,2),DPHIDV2(NZSR2,MBZ,2))
        ALLOCATE(DELTAT1(MBZ),DELTATP(MBZ),CAVLT(MBZ))
        ALLOCATE(CAVLTP(MBZ),CAVLTSA(MBZ,NSTEP))
      ENDIF

      IF(ITUN.EQ.1) THEN
C     MODULE TNGEO
        ALLOCATE(XTUN(NtMAXt+1,NCMAXt+1),YTUN(NtMAXt+1,NCMAXt+1))
        ALLOCATE(ZTUN(NtMAXt+1,NCMAXt+1))
        ALLOCATE(VXITN(NTMAXT,NCMAXT),VETATN(NTMAXT,NCMAXT))
        ALLOCATE(VINFSTN(NTMAXT,NCMAXT))
        ALLOCATE(DPDUTN(NTMAXT,NCMAXT),DPDVTN(NTMAXT,NCMAXT))
        ALLOCATE(UXTNTOT(NTMAXT,NCMAXT),UYTNTOT(NTMAXT,NCMAXT))
        ALLOCATE(UZTNTOT(NTMAXT,NCMAXT),CPTN(NTMAXT,NCMAXT))
      ENDIF

C.....we should read the angle of attack here fore hydrofoil geometry...
C.....so it can be use in gwake.f subroutine. (JY012199)................
C.... Same for ICON=8 as well (JY110100)

      IF((ICON.EQ.5).OR.(ICON.EQ.6).OR.(ICON.EQ.8))THEN
C       WRITE(*,'(A)') '  '
        WRITE(*,*) ' PROPCAV> ENTER ANGLE OF ATTACK (in degrees): '
        READ(*,*) ALPHA
      ENDIF

      IF(ISP .NE. 0) THEN
        WRITE(*,*)
        WRITE(*,'(A,$)') ' PROPCAV> ENTER INCLINED ANGLE (degrees): '
        READ(*,*) SPANGLE
        SPANGLE = SPANGLE * PI / 180.
      ENDIF

C-----------------------------------------------------------------------
C     Check to see if DELK is divisable by DELTAT.              JY111900
C-----------------------------------------------------------------------

      if(360/ndltat*ndltat .ne. 360) then
        write(*,*) ' ================ WARNING ===================='
        write(*,*) ' --> 360 degree is not divided by input DELTAT'
        write(*,*) ' --> Deltat should be one of following numbers'
        write(*,*) '     (ex. 2,3,4,5,6,8,9,10,12,15,18...)       '
        write(*,*) ' --> Solution from the nearest angular        '
        write(*,*) '     will be used!                            '
        write(*,*) ' ============================================='
      endif

C-----------------------------------------------------------------------
C     Generate the geometries of propeller blade (GBLADE), hub (GHUB)
C     and wake (GWAKE)
C-----------------------------------------------------------------------
      CALL CHRLEN(FN,LENCH)

C.... Since M841B is a special geometry, it is hard wired into the code.
      IF(ISC.EQ.1.AND.FN(1:LENCH).EQ.'m841b') THEN
        CALL M841B
        GO TO 1222
      END IF

      IUNS=0

      FNPLT=FN(1:LENCH)//'.plt'
      OPEN(700, FILE=FNPLT, STATUS='UNKNOWN')
      WRITE(700,5112)
 5112 FORMAT(1X,'VARIABLES = "x", "y", "z"')

!s--YE TIAN ---
      NSUB = 5
      IF(IAN.NE.2) THEN
         NWSUB1=4
      ELSE
         NWSUB1=2
      END IF
C    -------------------------------------------------
C           IAN = 2        |       IAN = 1,3,6,7
C    ----------------------|--------------------------
C        VIS   |  INVISCID |     VIS   |  INVISCID
C    ----------|-----------|-----------|--------------
C       5 & 2  |   5 & 2   |    5 & 4  |   5 & 4
C    -------------------------------------------------
C                                  S.N.KIM | Oct. 2018

      IF(IVISC.EQ.1) THEN
        WRITE(*,*)
        WRITE(*,*) ' PROPCAV> INPUT NSUB'
        WRITE(*,*)'          (no. of macro-panels in the wake',
     *            ' to subdivide):'
        READ(*,*) NSUB
        IF(IAN.NE.2) THEN
          WRITE(*,*) ' PROPCAV> INPUT NWSUB1'
          WRITE(*,*) '          (no. of panels each wake macro-panel',
     *               ' to be subdivided into):'
          READ(*,*) NWSUB1
        END IF
      ENDIF
C    -------------------------------------------------
C           IAN = 2        |       IAN = 1,3,6,7
C    ----------------------|--------------------------
C       VIS    | INVISCID  |     VIS       | INVISCID
C    ----------|-----------|---------------|----------
C     NSUB & 2 |  5 & 2    | NSUB & NWSUB1 |  5 & 4
C    -------------------------------------------------
C                                  S.N.KIM | Oct. 2018

      CALL PROGEO

 1222 CONTINUE


C---- Check paramters for dimensioning errors.  Results are printed to
C---- ERR.LOG.

      if(iscav .eq. 0) CALL CHECK_PARAM2

      IF(ICON.EQ.1) THEN
        STOP
      END IF

C-----------------------------------------------------------------------
C     Input unsteady parameters
C-----------------------------------------------------------------------

      CALL CHRLEN(FN,LENCH)

C-----------------------------------------------------------------------
C     How many time steps?
C-----------------------------------------------------------------------
      FNLOG=FN(1:LENCH)//'.log'
      OPEN(26,FILE=FNLOG,STATUS='UNKNOWN')
      WRITE(*,*) ' PROPCAV> ENTER NUMBER OF REVOLUTIONS FOR FULLY',
     *  ' WETTED COMPUTATION'
      WRITE(*,*) '          (ignored if ISP=1 or ISTEADY=0):'
      READ(*,*) NTREVW

C.... If ISP=1, set NTREVW=0
      IF(ISP.EQ.1.OR.ISTEADY.EQ.0) THEN
        IF(ISP.EQ.1) WRITE(*,*) ' Setting NTREVW=0 for ISP=1'
        IF(ISTEADY.EQ.0) WRITE(*,*) ' Setting NTREVW=0 for ISTEADY=0'
        NTREVW=0
      ELSE

C.......If ISC=1, set NTREVW=1
        IF(ISC.EQ.1.AND.ISTEADY.GT.1.AND.NTREVW.NE.1) THEN
          WRITE(*,*) ' Only 1 wetted revolution is needed'
          WRITE(*,*) ' Setting NTREVW=1'
          WRITE(*,*)
          NTREVW=1
        END IF
      END IF

      WRITE(*,*) ' PROPCAV> ENTER NUMBER OF REVOLUTIONS FOR CAVITY',
     *  ' COMPUTATION:'
      READ(*,*) NTREVC

C      WRITE(*,*) NTREVC
C.... If ISC=1, special requirements for NTREVC
      IF(ISP.NE.1.AND.ISC.EQ.1.AND.ISTEADY.GT.1.
     *  AND.NTREVC.LT.6) THEN
        WRITE(*,*) ' At least 6 cavitating revolution is needed'
        WRITE(*,*) ' Setting NTREVC=6'
        WRITE(*,*)
        NTREVC=6
      END IF

      IF(ISC.EQ.1.AND.ISTEADY.LE.1.AND.NTREVC.EQ.0) THEN
        WRITE(*,*) ' You must perform cavitating analysis for ISC=1'
        IF(ISP.NE.1) THEN
          WRITE(*,*) ' Setting NTREVC=6'
          WRITE(*,*)
          NTREVC=6
        ELSE
          WRITE(*,*) ' Setting NTREVC=4'
          WRITE(*,*)
          NTREVC=4
        END IF
      END IF

C     --- Write mean forces with revolution ....

      OPEN(712,FILE='meanforce.plt',STATUS='UNKNOWN')

      IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8) THEN
        WRITE(712,*)
     %    'Title ="Mean total Forces at each revolution" '
        WRITE(712,*)
     %    ' VARIABLES="REV","FX","FY","FZ","MX","MY","MZ"'
      ELSE
        WRITE(712,*)
     %    'Title =" Lift & Drag Coefficients"'
        WRITE(712,*)
     %    ' VARIABLES="REV","CL","CD"'
      ENDIF

      IF(ICON .NE. 5 .AND. NTREVW .NE. 0) THEN
        WRITE(712,*)
     %    ' Zone T="Fully Wetted Forces" '
      ENDIF

      IF(ICON .EQ. 5) THEN
        WRITE(712,*)
     %    ' Zone T="Fully Wetted Forces" '
      ENDIF

      DO I = 1 , 20
        DO J = 1 , 6
          AFM(I,J) = 0.0
        ENDDO
      ENDDO

C     --- Write Circulation distribution on file cir.plt & cir2.plt
C     ---- These outputs replace wet.cir, steady.cir, cav.cir ......
C     -- cir.plt : Non-dim by UR & Vs
C
      OPEN(1700,FILE='cir.plt',STATUS='UNKNOWN')
      OPEN(1701,FILE='cir2.plt',STATUS='UNKNOWN')

      IF(ICON .NE. 5) THEN
        WRITE(1700,*)
     %    'TITLE ="CIRCULATION DISTRIBUTION: Gr=DELP/(2*pi*R*Ur)"'
        WRITE(1700,*) 'Variables ="r/R","100Gr"'
        WRITE(1701,*)
     %    'TITLE ="CIRCULATION DISTRIBUTION: Gs=DELP/(2*pi*R*Vs)"'
        WRITE(1701,*) 'Variables ="r/R","100Gs"'
      ELSE
        WRITE(1700,*)
     %    'TITLE ="CIRCULATION DISTRIBUTION: Ginf=DELP/(2*pi*R*U_inf)"'
        WRITE(1700,*)
     %    'Variables ="r/R","100Ginf"'
        WRITE(1701,*)
     %    'TITLE ="CIRCULATION DISTRIBUTION: Ginf=DELP/(2*pi*R*U_inf)"'
        WRITE(1701,*)
     %    'Variables ="r/R","100Ginf"'
      ENDIF

C-----------------------------------------------------------------------
C     Distance from center of hub to free surface if ISP=1.     JY112299
C-----------------------------------------------------------------------
      IF(ISP.EQ.1) THEN
        WRITE(*,*) '  PROPCAV> ENTER SUBMERGENCE RATIO (h/D)'
        READ(*,*) SRATIO
        YFS=SRATIO*2.-1.
        WRITE(*,*) '  PROPCAV> YFS = ',YFS
      END IF

C-----------------------------------------------------------------------
C     Figure out the number of revolution and the left time steps
C     NTPREV is the number of time steps in one revolution CM
C-----------------------------------------------------------------------
      NTPREV = 360 / NDLTAT !Number of Time Steps in one revolutions (= 360/NDLTAT)
C     WRITE(*,'(A)') ' '
C     WRITE(*,'('' TIME STEPS IN ONE REVOLUTIONS: '',I4)') NTPREV
C     WRITE(*,'(A)') ' '

C-----------------------------------------------------------------------
C     NTREV is the total number of revolutions.
C     NTIME is the total number of timesteps.
C-----------------------------------------------------------------------

      NTREV=NTREVW
      NTIME=NTREV*NTPREV  ! NTPREV = 360 / NDLTAT

      DO 10 N=1,NTPREV
        IF(N.LE.1) THEN
          TT(N)=0.
        ELSE
          TT(N) = REAL(N-1)*REAL(NDLTAT)
        END IF
 10   CONTINUE

C.... Jump directly to the cavity solution if ISP=1. (JY112299)
      IF(ISP.EQ.1) THEN
        NTSTEP=0
        CALL INDPOT(1)
        GO TO 301
      END IF

C-----------------------------------------------------------------------
C---------Computations that only need to be done once or variable-------
C---------assignments------------------------------------------CM100297-
C-----------------------------------------------------------------------

      CAVINI=0.08

      IF((ICON.EQ.4).OR.(ICON.EQ.5).OR.(ICON.EQ.6)
     *  .OR.(ICON.EQ.8))THEN
        WRITE(*,*) '  PROPCAV> ENTER INITIAL CAVITY LENGTH: '
        IF(IFACE.EQ.0.OR.IFACE.EQ.1) THEN
          READ(*,*) CAVINI
        ELSE IF(IFACE.EQ.2) THEN
          READ(*,*) CAVINIB,CAVINIF
        END IF
      ENDIF

      IWET=1

C-----------------------------------------------------------------------
C     READ INPUTS DATA FOR VISCOUS BOUNDARY LAYER ANALYSIS
C     BY HONG SUN   02/26/04
C-----------------------------------------------------------------------

      IF(IVISC.EQ.1) THEN
        WRITE(*,'(A,$)')' PROPCAV> ENTER REYNOLDS NUMBER (=VS*D/v) '
        WRITE(*,*)
        READ(*,*)  REYD
        WRITE(*,*)
     *  ' PROPCAV> ENTER NCRITICAL (9=.09% TURB LEV; 2.623=1%TURB LEV)'
        WRITE(*,*)
        READ(*,*) RVCRIT
C     WRITE(*,'(A,$)')' ENTER FORCED TRANSITION ON THE SUCTION SIDE:'
C     WRITE(*,*)
C     READ(*,*)  XTRANS
C     WRITE(*,'(A,$)')' ENTER FORCED TRANSITION ON THE PRESSURE SIDE:'
C     WRITE(*,*)
C     READ(*,*)  XTRANP
        WRITE(*,*)
     *  ' PROPCAV> ENTER NUMBER OF ITERATIONS FOR VISOUS RUN:'
        WRITE(*,*)
        READ(*,*)  MAXIT
        WRITE(*,*)
     * ' PROPCAV> ENTER THE TOLERANCE FOR VISCOUS RUN (e.g. EPS1=1E-4):'

        WRITE(*,*)
        READ(*,*)  EPS1
        WRITE(*,*)
     * ' PROPCAV> ENTER THE UNDER-RELAXATION FACTOR (e.g. beta=0.5):'
        WRITE(*,*)
        READ(*,*)  BETA_BL_INP
        WRITE(*,*)
C     Read the transition location at both sides for BL analysis
        CALL READXTRAN
      ENDIF
C--------------------------------------------------BY HONG SUN 02/26/04


C-----------------------------------------------------------------------
C     Open files for wetted analysis.
C-----------------------------------------------------------------------
      IF(ISC.NE.1) CALL OPFILEW

C-----------------------------------------------------------------------
C     Open files for hullfpp information by S.H.CHANG 02/25/2010
C-----------------------------------------------------------------------
C      CALL OPFILEH

C-----------------------------------------------------------------------
C     Begin the unsteady process
C     time step  0: steady case
C     time step  1: beginning the unsteady effect, the propeller is
C     at theta=0.0
C-----------------------------------------------------------------------
C
C.....NREV:   current number of revolution
C.....IDXREV: number of time step in current revolution

      NREV=1

      WRITE(*,1004)

 1004 Format('------Start fully-wetted computation---------------')

C     -------------------------------------------------------
C     Open files for tip vortex model
C     -------------------------------------------------------
      IF(IAN.EQ.2) THEN
        CALL OPFILET

C     -- Start run from cavity problem by skipping wetted run.

        if(iscav .ne. 0) then
          rc1 = radini
          if(iscav .eq. 1) call inidetach
          go to 301
        endif
        icavsave = icavmax ! Initial 'icavmax' comes from 'cav.ctr' file.
      END IF

C -- Read effective wake if pc2ns_ueff.dat file exist!
CYiranSu_PC2NS_read_effective wake
      INQUIRE(FILE='pc2ns_ueff.dat',EXIST=IPC2NS)
      IF (IPC2NS.EQ.(.TRUE.)) THEN
        NNN=NTPREV
        IF (ISTEADY.EQ.0) NNN=1
        ALLOCATE(UEFX(NC,MR,NNN),UEFR(NC,MR,NNN),UEFT(NC,MR,NNN))
        OPEN(991,FILE='pc2ns_ueff.dat',STATUS='OLD')
        READ(991,*)
        READ(991,*)
        DO K=1,NNN
          DO J=1,MR
            DO I=1,NC
        READ(991,*) TXX,TYY,TZZ,UEFX(I,J,K),UEFR(I,J,K),UEFT(I,J,K),TPP
            END DO
          END DO
        END DO
        CLOSE(991)
        WRITE(*,*) 'Using effective wake from PC2NS subroutine!'
      ELSE
        WRITE(*,*) 'Using .WAK file!'
      END IF
CYiranSu End

      DO 100 N=0,NTIME  ! NTIME = NTREV * NTPREV

        NTSTEP=N
        IF(N.LE.1) THEN
          TSTEP=0.
          ITSTEP = 0
          IDXREV=N
        ELSE
          TSTEP=TSTEP+DELTAT
          ITSTEP=ITSTEP+NDLTAT
          IF(ITSTEP .EQ. 360) THEN
            NREV=NREV+1
            IDXREV=0
            TSTEP=ZERO
            ITSTEP = 0
          END IF
          IDXREV=IDXREV+1
        END IF

 1010   FORMAT('+ Time Step=',I5,'   Blade angle=',F4.0)
        IF(IDXREV.EQ.0) THEN
          WRITE(*,1010) N,TT(1)
          WRITE(26,1010) N,TT(1)
        ELSE
          WRITE(*,1010) N,TT(IDXREV)
          WRITE(26,1010) N,TT(IDXREV)
        END IF

C     Call wake geometry for the full wake alignment
C     Ippp = 1 : calculate influence coeff. and save into files
C     Ippp = 0 : Do not Calculate Inf. Coeff. and Read from saved files.
C     icavmax : Maximum number of Rev. for the wake alignment (=2 rev)
        if(ian .eq. 2) then
          if(nrev .le. 8 .and. n .ne. 0) then
            call callgeok(1)
          endif

C     --- If (NREV .GT. 8), Do not solve wake align
C     Solve unsteady problem with already aligned unsteady wake.
C     So, Read previous wake geometry corresponding to the blade position.
          if(nrev .gt. 8) call callgeok(2)

C/s S.N.KIM | Tip vortex model is ignored.
c          if(ian.eq.2)then
c            call tcavgeo
c          endif
C/e S.N.KIM | Aug. 2018.

          if(ntstep .eq. 0) then
            icavmax = icavsave
          else
            if(nrev .le. 8) then
              icavmax= 2
              ippp = 2
            else
              icavmax = 1
              ippp = 2
            endif
          endif

C     ------ From here ---------------------------
C     Set for the normal alignment
        else              ! Steady Case ignores input value of icavmax. It sets the value as 1.
          icavmax = 1
          ippp = 1
          if(n .ne. 0) ippp = 0   ! After 0, we only take care of wake related influence coefficient.
        endif

!YE TIAN  04/12/2012
        if (ian.eq.6) then
          icavmax = NALIGN_ITER
          ippp = 1
          nwpanel = WKPAN_INP
C/s S.N.KIM | After n.eq.0, we only take care of wake related inf. coeff. 
C           | Actually, FWA (ian.eq.6) is only for the steady runs (n.eq.0), but
C           | in case it goes unsteady regime, IPPP needs to be set 0 to
C           | avoid inf. coeff. calculation.  
          if(n .ne. 0) then
            icavmax = 1
            ippp = 0
          endif 
C/e S.N.KIM | Aug. 2018.  
        endif

C/s S.N.KIM | Given the 'ktkq-history.plt' file is only for IAN=6,
C           | it is updated only when the FWA is in process in
C           | steady state (ntstep=0). 
        if(ian.eq.6.and.ntstep.eq.0) then
          open(1108,FILE='ktkq-history.plt',STATUS='UNKNOWN')
          write(1108,*) 'TITLE="Force Performance"'
          write(1108,*) 'Variables="Ite","TT","KTP","10KQ","KTPV"
     &,"10KQV"'
          write(1108,*) 'ZONE T="Force Performance History"'
        endif
C/e S.N.KIM | Sep. 2018.

C--   Do iteration for wake alignment up to icavmax rev. -------
        do 1000 icavt = 1 , icavmax  ! icavmax = NALIGN_ITER

!-START(Seung-Nam, Kim: 03/03-2016)-----<REPANELLING PROJECT>----------------------------
          IF(IREPANEL.EQ.1) THEN
            IF(ICAVT .GT. 1) THEN
              CALL DUCTGEO            
              IF(IDFWA.EQ.1) THEN
                CALL DUCTWAKGEO_FIXEDPANEL_FWA ! aligned ductwake is rotated to follow duct.  
              ELSE
                CALL DUCTWAKGEO_FIXEDPANEL ! cylindrical ductwake is rotated to follow duct.  
              ENDIF
              CALL GWAKES_DUCT 
              IF(ICAVT.EQ.ICAVMAX) CALL DUCTPLT 
              WRITE(*,*)'Repaneling on the duct and ductwake is 
     &complete.'
            ENDIF
          ENDIF

C/s S.N.KIM | Considering its huge file size, wake history file is
C           | created only at the 0th time step when FWA is adopted.
          IF (IAN.EQ.6.AND.NTSTEP.EQ.0) THEN
            IF(ICAVT .EQ .1) THEN
              OPEN(3367,FILE='full-wake-his.plt'
     &,STATUS='UNKNOWN')
            WRITE(3367,*)'TITLE="Plot Motion of Wake and Duct Geometry"'
            WRITE(3367,*)'VARIABLES="X","Y","Z"'
            ENDIF
            CALL WAKEPLT_CHECK
            IF (IDUCT.EQ.1) CALL DUCTPLT_CHECK
            IF(ICAVT .EQ. ICAVTMAX) THEN
              CLOSE(3367)
            ENDIF
          ENDIF
C/e S.N.KIM | Sep. 2018.

C/s S.N.KIM | Unsteady Wake Alignment
          IF(IAN .EQ. 2 .AND. IUPLOT .EQ. 1) THEN
            IF(ICAVT .EQ. 1 .AND. NTSTEP .EQ. 0) THEN
              OPEN(8829,FILE='unsteady-wake.plt',STATUS='UNKNOWN')
              WRITE(8829,*) 'TITLE="Unsteady Wake Geometry"'
              WRITE(8829,*) 'VARIABLES = "X","Y","Z","DPHI","InfCoef"'
            ENDIF
          ENDIF
C/e S.N.KIM | Aug. 2018.
!-END(Seung-Nam, Kim: 03/03-2016)-------<REPANELLING PROJECT>----------------------------

          if(ian .eq. 2 .and. ntstep .eq. 0) call savegeow
C-----------------------------------------------------------------------
C     Calculate induced potentials and construct influence
C     coefficient matrices
C-----------------------------------------------------------------------
          if ((ian.eq.6.and.n.eq.0).or.(ian.eq.2.and.ntstep.eq.0)) then 
!Same as the steady runs, blade inf. coeff. only 
!needs to be evaluated at the first ite. for ian=2. | Aug. 2018. S.N.KIM
            if (icavt.eq.1) then
              ippp = 1
            else
              ippp = 2
              if (iduct.eq.1.and.irepanel.eq.1)
     *          ippp = 1
!S.N.Kim 01/03/2016 Changed ippp value from 2 into 1
!since influence coefficients on duct need to be recalculated due to repaneling.
            endif
          endif

          if ((ntstep.eq.0).and.(icavt.eq.1).and.(ian.eq.6)) IAN = 8

          CALL INDPOT(ippp)
          if (n.eq.ntime) call wakeplt
          if (IAN.eq.8) IAN = 6

C-----------------------------------------------------------------------
C     If ICON = 2, only geometries and i.c. matrices are created
C-----------------------------------------------------------------------
          IF(ICON.EQ.2) THEN   ! Geometry and Influence Computations Only
            STOP
          END IF

C-----------------------------------------------------------------------
C     Solve the simultaneous equation and create output files
C-----------------------------------------------------------------------
          CALL SIMEQN

C/s S.N.KIM - Unsteady Wake Alignment
      IF(IAN.EQ.2 .AND. IUPLOT .EQ. 1) THEN
        CALL WAKEPLT_UNSTEADY
        IF(ICAVT.EQ.ICAVMAX.AND.NTSTEP.EQ.NTIME) CLOSE(8829)
      ENDIF
C/e S.N.KIM - Unsteady Wake Alignment     

C-----------------------------------------------------------------------
C     Draw the wetted velocity plot for the last wetted revolution.
C     JY100598
C-----------------------------------------------------------------------
          IF(IDXREV.NE.0.AND.NREV.EQ.NTREV) THEN ! NTREV = NREV = 0.
            IF(IAN .NE. 2) CALL VECPLT
          END IF

          DO 20 M=1,MR
            DPHI(M,IDXREV)=DELP(M)
 20       CONTINUE

!s-YE TIAN evaluate force for each step
            IF(ISC.EQ.1) THEN
              CALL FORCEV_SC
            ELSE
              CALL FORCEV
            END IF

C---------------------------------------------------------------------------------------C
C   S.N. KIM  |^|  This part is the crucial part of the steady full wake alignment in   C
C             |^|  periodic inflow condition. Initial version of this scheme before     C
C  ---------- |^|  PROPCAV V3.3 was initially introduced by Dr. Tian, Y. back in 2012.  C
C  Aug. 2018  |^|  Later, S.N. KIM modified significant part of his scheme to           C
C             |^|  incorporate more general input files with stable outputs and most    C 
C             |^|  importantly corrected the psuedo time step term in the alignment     C
C             |^|  algorithm, which yields better validation with RANS for open and     C
C             |^|  ducted propellers.                                                   C
C             |^|  -------------------------------------------------------------------  C
C             |^|  Reference : Dr. Tian, Y.        Int. J. of Rotating Machinery, 2012  C
C             |^|              KIM, S.      J. of Marine Science and Engineering, 2018  C
C---------------------------------------------------------------------------------------C
            if(ian.eq.6.and.n.eq.0)then
              imrw = mr/ALIGN_DIV
              if(iduct.eq.1) imrwd = mduct/2 !ALIGN_DIV
              call run_alignp(icavt,imrw,imrwd)
              call wrtback_wakgeo   ! After FWA, generate wake coordinates into NC * MR coord.
              if(idfwa.eq.1) call wrtback_ductwakgeo
              call gwakes(NSUB,NWSUB1)
              WRITE(*,*) 'ZONE T="STEADY KT,KQ,KTV,KQV"'
              WRITE(*,7005) TT(1),Nblade*FBXP(1),
     %        Nblade*FBXP(4),Nblade*FBXPV(1),Nblade*FBXPV(4)
              WRITE(1108,7005) icavt,TT(1),Nblade*FBXP(1),
     %        Nblade*FBXP(4),Nblade*FBXPV(1),Nblade*FBXPV(4)
            end if
C---------------------------------------------------------------------------------------C

C-----------------------------------------------------------------------
C     Calculate the forces (wetted)
C     Forcep is potential based forces, forcev is viscous
C-----------------------------------------------------------------------
C-------Modification for full wake alignment.
          if(icavt .eq. icavmax) then

            IF(ISC.EQ.1) THEN
              if(icavt .lt. icavmax) then
              CALL FORCEV_SC
              end if
            ELSE
              if(icavt .lt. icavmax) then
              CALL FORCEV
              end if
            END IF

C-----------------------------------------------------------------------
C     Print steady, fully wetted forces.
C-----------------------------------------------------------------------
            IF(ISC .NE. 1) THEN
              IF(NTSTEP.EQ.0) THEN
                WRITE(73,*) 'ZONE T="STEADY KT,KQ,KTV,KQV"'
                WRITE(73,7005) TT(1),Nblade*FBXP(1),
     %            Nblade*FBXP(4),Nblade*FBXPV(1),Nblade*FBXPV(4)

                IF(IDUCT .NE. 0) THEN
                  WRITE(622,*) 'ZONE T="STEADY KT,KQ,KTV,KQV"'
                  WRITE(622,7005) TT(IDXREV+1),Nblade*FDXP(1),
     %              Nblade*FDXP(4),Nblade*FDXPV(1),
     %              Nblade*FDXPV(4)
                ENDIF

 7005           FORMAT(1X,I2,5(1X,E14.7))

                IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8) THEN
                  DO NN=1,NTPREV
                    DO KK=1,6
                      XKTV(NN,KK)=FBXPV(KK)
                    END DO
                  END DO
                  CALL PUF3HRM(NBLADE,NTPREV,IWET,XKTV,NSTEP)
                END IF

              ELSE
C/s S.N.KIM | wants to see the unsteady forces at every revolution in a seperate zone.
c                IF(NTSTEP.EQ.1)
                IF(ITSTEP.EQ.0) ! Whenever a new revolution starts, ITSTEP sets to be 0, so I use this instead of NESTEP.
     *            WRITE(73,*) 'ZONE T="',NREV,'REVOLUTION"'
c     *            WRITE(73,*) 'ZONE T="UNSTEADY KT,KQ,KTV,KQV"'
C/e S.N.KIM | Aug. 2018.
                IF(NTSTEP.EQ.1 .AND. IDUCT .NE. 0)
     *            WRITE(622,*) 'ZONE T="UNSTEADY KT,KQ,KTV,KQV"'
                IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8) THEN
                  DO IKK=1,6
                    XKT(IDXREV,IKK)=FBXP(IKK)
                    XKTV(IDXREV,IKK)=FBXPV(IKK)
                    AFM(NREV,IKK)=AFM(NREV,IKK)
     %                +FBXPV(IKK)/REAL(NTPREV)
                  ENDDO
                ENDIF
              END IF
            ENDIF

C-----------------------------------------------------------------------
C     Write wetted outputs to files.
C-----------------------------------------------------------------------

            IF(ISC.NE.1) THEN
              IF(NREV.EQ.NTREV.OR.IDXREV.EQ.0) call writefilew ! Force plotting at the last Rev.
!Allen Du 01/11/2018 output the cp_min
              if(opt_cpmin_w<=opt_cpmin_temp) opt_cpmin_w=opt_cpmin_temp
C/s S.N.KIM | Unsteady Force Plotting
              IF(NREV.NE.NTREV.AND.IDXREV.NE.0) THEN
                WRITE(73,7000) TT(IDXREV),FBXP(1),FBXP(4),     ! Force plotting before the last Rev.
     *                         FBXPV(1),FBXPV(4)
 7000      FORMAT(1X,F9.4,2X,F10.6,2X,F10.6,6X,F10.6,2X,F10.6)
              ENDIF
C/e S.N.KIM | Aug. 2018.
            ENDIF

C***********************************************************************
C     COUPLED WITH XFOIL         BY HONG SUN         April, 2006
C***********************************************************************

            IF(IVISC.EQ.1.AND.ISC.NE.1.AND.NTREVC.EQ.0) THEN
              IF(NREV.EQ.NTREV.OR.IDXREV.EQ.0) THEN
                CALL GEO2DW
                IF(IDUCT.EQ.1) CALL GEO2D_D
                CALL INFLOWBL
c XM YU 12/2011-----------------------------------------------------------------
c Iterative scheme for 3D boundary layer
                CALL WAKV 
                CALL BODYV
                umass = 0.0
                do ita=1,3
                   CALL CAVBL2D(ita)
                   residual=0.0
                   do m=2,mr-1
                      do j=1,ntw(m)-1
                       temp=abs(updatem(j,m)-umass(j,m))
                       if (temp.ge.residual) then
                           residual=temp
                           mrecord=m
                           jrecord=j
                       endif
                      enddo
                   enddo
                   write(*,*) 'residual', residual
!Allen 03/21/2019 extension scheme
!                  if (residual.le.1e-3) exit

                   do m=1,mr
                      do j=1,ntw(m)-1
                        umass(j,m)=updatem(j,m)
                      enddo
                   enddo
                enddo
c XM YU 12/2011-----------------------------------------------------------------

C-----------------------------------------------------------------------
C     Recalculate and Print steady and unsteady, fully wetted forces.
C     (with Boundary Layer effects)
C-----------------------------------------------------------------------
                CALL FORCE_BL
                IF(NTSTEP.EQ.0) THEN
                  WRITE(74,*) 'ZONE T="STEADY KT,KQ,KTV,KQV"'
                  WRITE(74,7005) TT(1),Nblade*FBXP(1),Nblade*FBXP(4),
     &              Nblade*FBXPV(1),Nblade*FBXPV(4)

                  IF(IDUCT.NE.0.AND.IDOPT.EQ.1) THEN
                    WRITE(623,*) 'ZONE T="STEADY KT,KQ,KTV,KQV"'
                    WRITE(623,7005) TT(IDXREV+1),Nblade*FDXP(1),
     %                Nblade*FDXP(4),Nblade*FDXPV(1),
     %                Nblade*FDXPV(4)
                  ENDIF

                  IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8) THEN
                    DO NN=1,NTPREV
                      DO KK=1,6
                        XKTV(NN,KK)=FBXPV(KK)
                      END DO
                    END DO
                    CALL PUF3HRM(NBLADE,NTPREV,IWET,XKTV,NSTEP)
                  END IF

                ELSE

                  IF(NTSTEP.EQ.1)
     *              WRITE(74,*) 'ZONE T="UNSTEADY KT,KQ,KTV,KQV"'
                  IF(NTSTEP.EQ.1.AND.IDUCT.NE.0.AND.IDOPT.EQ.1)
     *              WRITE(623,*) 'ZONE T="UNSTEADY KT,KQ,KTV,KQV"'

                  IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8) THEN
                    DO IKK=1,6
                      XKT(IDXREV,IKK)=FBXP(IKK)
                      XKTV(IDXREV,IKK)=FBXPV(IKK)
                    ENDDO
                  ENDIF

              WRITE(74,7005) TT(IDXREV),Nblade*FBXP(1),Nblade*FBXP(4),
     *              Nblade*FBXPV(1),Nblade*FBXPV(4)
                  IF(IDUCT.NE.0.AND.IDOPT.EQ.1)
     *              WRITE(623,7005) TT(IDXREV),Nblade*FBXP(1),
     *              Nblade*FBXP(4),Nblade*FBXPV(1),
     *              Nblade*FBXPV(4)

                END IF

              ENDIF
            ENDIF
C************************************************BY HONG SUN APRIL,06

C-----------------------------------------------------------------------
C     Calculate the initial guess for the detachment line
C     based on the wetted pressure distribution---CM071797
C
C     The DETACH routine should not be called during the steady, fully
C     wetted step (i.e. IDXREV=0).                             JY081498
C----------------------------------------------------------------------
            IF(IDXREV.NE.0) CALL DETACH

          end if

C     -- Begin Tip HSLEE(10/13/99)
C Yiran Su 10/05/2016: Files are not closed at final FWA iteration, so that FWA can be run with cavitation
C Yiran          if(((ian.eq.2).or.(ian.eq.6)) .and. ippp .ne. 0) then
          if (((ian.eq.6).and.(icavt.ne.icavmax)).or.(ian.eq.2)) then
          if (ippp.ne.0) then
C Yiran End 10/05/2016
            close(41)
            close(42)
            close(50)
            close(51)
            close(52)
            close(53)

            do kkk = 1 , nblade+1
              ikkk = 80 + kkk
              close(ikkk)
            enddo

            if(iduct .ne. 0) then
              do kkk = 1, nblade
                 ikkk = 500 + kkk
                 close(ikkk)
              enddo
            endif

            if((icavt .lt. icavmax).and.(ian.eq.2)) then
C---------------------------------------------------------------------------------------C
C   S.N. KIM  |^|  This part is the crucial part of the unsteady wake alignment in      C
C             |^|  periodic inflow condition. Non periodic flow is later covered by     C
C             |^|  Dr. Yiran Su in his Ph.D. diss. as more elaborate case. Initial      C
C  ---------- |^|  version of this scheme before PROPCAV V3.3 was based on tip vortex   C
C     MHL     |^|  model, which was initially intoduced by Dr. H.S. LEE back in 2002.   C
C  Aug. 2018  |^|  Later, S.N. KIM modified significant part of his scheme              C
C             |^|  to incorporate more general input files with stable outputs.         C
C             |^|  Tip vortex model is no longer available with IAN=2,                  C
C             |^|  but remains as future research.                                      C
C             |^|  -------------------------------------------------------------------  C
C             |^|  Reference : Dr. Lee, H.         Ph.D. Diss.,   UT-Austin, 2002.      C
C             |^|              KIM, S.             M.S. Thesis,   UT-Austin, 2017.      C
C             |^|              Dr. Su, Y.          Ph.D. Diss.,   UT-Austin, 2018.      C
C---------------------------------------------------------------------------------------C
              imrw = mr/IUALIGN_DIV ! IUALIGN_DIV comes from 'cav.ctr' file.     
              imrwd = 0             ! Periodic unsteady runs are ONLY available for open propellers as of Aug. 2018.
              if(ntstep.eq.0) then                                  
                call run_alignp(icavt,imrw,imrwd)
              else
                call run_alignp_unsteady(icavt,imrw)
              endif
              call wrtback_wakgeo
              call gwakes(nsub,nwsub1)
C-------------------------------------------------------------------------------------C
            endif
          endif
          endif

          if((ian.eq.2)) then
            if((nrev .le. 5 .or. nrev .eq. ntrev) .and. 
     *          icavt .eq. icavmax .and. ntstep .ne. 0)
     *                        call plot3d2
            if(ntstep .eq. 0) call plot3d2
          end if

 1000   continue

C     ---- Save geometry for the next time step
C     ---- Plot Cp, Vt, and potential on tip vortex at each revolution.
        if(ian .eq. 2) call wtprs
        if(ian .eq. 2 .and. ntstep .eq. 0) call savegeok
        if(ian .eq. 2 .and. nrev .le. 8) call savegeow

 100  CONTINUE

C ======================================
C PC2NS for coupling Fluent with PROPCAV
C                               YIRAN SU
C                             07/07/2016
C ======================================
      IWRT2 = .FALSE.
      INQUIRE (FILE='pc2ns.trigger',EXIST=IWRT2)
      IF (IWRT2) CALL PC2NS

C Yiran Su export wake geometry for IAN=5
        open (1105, file='wakinp.dat', status='unknown')
          do mm= 1 , mr+1
            write(1105, *) nsw(mm)
            do nn= 1, nsw(mm)
              write(1105, '(3F12.6)')xw(nn,mm),yw(nn,mm),zw(nn,mm)
            enddo
          enddo
        close(1105)

C-----------------------------------------------------------------------
C     Calculate the fully wetted blade and shaft harmonics for the last
C     revolution and write the results to harmny.wet.           JY071399
C-----------------------------------------------------------------------
C.... Same for ICON=8 as well (JY110100)

      IF(IDXREV.NE.0.AND.ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8.
     *  AND.ISC.NE.1)
     *  CALL PUF3HRM(NBLADE,NTPREV,IWET,XKTV,NSTEP)

C-----------------------------------------------------------------------
C     caculating the average wetted forces and printing them to
C     *.wktkq.                                                  JY071399
C-----------------------------------------------------------------------
      IF(IDXREV.NE.0.AND.ISC.NE.1) THEN
        AVKT=0
        AVKQ=0
        AVKTV=0
        AVKQV=0
        DO I=1,NTPREV
          AVKT=AVKT+XKT(I,1)
          AVKQ=AVKQ+XKT(I,4)
          AVKTV=AVKTV+XKTV(I,1)
          AVKQV=AVKQV+XKTV(I,4)
        END DO
        AVKT=AVKT/NTPREV
        AVKQ=AVKQ/NTPREV
        AVKTV=AVKTV/NTPREV
        AVKQV=AVKQV/NTPREV
        WRITE(73,*) 'ZONE T="AVERAGED KT,KQ,KTV,KQV"'
        WRITE(73,7005) TT(1),Nblade*AVKT, Nblade*AVKQ,
     &    Nblade*AVKTV, Nblade*AVKQV
        WRITE(73,7005) TT(NTPREV),Nblade*AVKT, Nblade*AVKQ,
     &    Nblade*AVKTV, Nblade*AVKQV

C-----------------------------------------------------------------------
C     Creat new force output file, force.wet, which will contain all
C     six components of the viscous forces at all timesteps in the
C     last revolution.                                          JY011600
C-----------------------------------------------------------------------
        OPEN(711,FILE='force.wet',STATUS='UNKNOWN')
        write(711,*) 'TITLE="Fully Wetted 6-compt. Forces per Blade"'
        WRITE(711,*) 'VARIABLES="ANGLE","FX","FY","FZ","MX","MY","MZ"'

 7111   FORMAT(7(1X,E14.7))
        DO I=1,NTPREV
          WRITE(711,7111) TT(I),(XKTV(I,I1),I1=1,6)
        END DO
        CLOSE(711)


        IF(ICON .NE. 5) THEN
          DO I = 1 , NTREV
            WRITE(712,7112) I, (AFM(I,J)*NBLADE, J = 1, 6)
          ENDDO
        ENDIF

 7112   FORMAT(1x,i5,6(1X,E14.7))
      END IF

C-----------------------------------------------------------------------
C     Close files
C-----------------------------------------------------------------------
      CLOSE(15)
      IF(IHUB.NE.0.AND.IPHUB.EQ.1) CLOSE(17)
      CLOSE(25)
      CLOSE(26)
      CLOSE(73)

      CLOSE(313)
      CLOSE(314)

C-----------------------------------------------------------------------
C     End of fully-wetted time marching
C-----------------------------------------------------------------------

C     -- GBFLOW OUTPUT

C     write(*,*) afm(ntrev,1),afm(ntrev,4)

C     AAKT = AFM(NTREV,1)*NBLADE
C     AAKQ = AFM(NTREV,4)*NBLADE

C     CALL PROP2GBFLOW(AAKT,AAKQ)
      IF(ICON .EQ. 5) THEN
        WRITE(712,*) NREV, CLIFT,CDRAG
      ENDIF

      WRITE(*,1005)
 1005 Format('------End fully-wetted computation-----------------')

C-----------------------------------------------------------------------
C     Time computation stuff
C-----------------------------------------------------------------------

      DELDTI = DTIME ( TARRAY )
      WRITE(*,1001) DELDTI/60
 1001 format(/,'-- Run time for wetted computations is ',F10.2,
     *  ' min. of CPU time.')

 301  CONTINUE

      IF(NTREVC.EQ.0) GOTO 1111

C-----------------------------------------------------------------------
C     Begin cavitating flow time marching
C-----------------------------------------------------------------------
      IWET=0

      IF(NTREVC .NE. 0) THEN
        WRITE(712,*)
     %    ' ZONE T="Cavitating Forces" '
      ENDIF

      DO I = 1 , 20
        DO J = 1 , 6
          AFM(I,J) = 0.0
        ENDDO
      ENDDO

      WRITE(*,1006)
 1006 Format(//,'------Start unsteady cavity computation----------')
C-----------------------------------------------------------------------
C     NTREV is the total number of revolutions.
C     NTIME is the total number of timesteps.
C-----------------------------------------------------------------------

      NTREV=NTREVC
      IF(ISTEADY.EQ.0) THEN
        NCTIME=NTREV
      ELSE
        NCTIME=NTREV*NTPREV
      END IF

C-----------------------------------------------------------------------
C     Zero the stored forces following fully wetted sol'n       JY010200
C-----------------------------------------------------------------------
      DO N=1,NTPREV
        DO L=1,6
          XKT(N,L)=ZERO
          XKTV(N,L)=ZERO
        END DO
      END DO

      NREV=1
      IF(IAN .EQ. 2 .AND. ISCAV .EQ. 2) THEN
        NREV = NNREV
      ENDIF
C     SH--------------------------------------------------------------------
C     SH     interpolate to get the inflow over the rudder control points
C     SH--------------------------------------------------------------------
      if(itergb.eq.1) CALL GBPOST(1)
C     SH-----------------------------rudder----Shreenaath 05/20/03---------------
      
      DO 300 N=1,NCTIME

        IF(ISTEADY.EQ.0) THEN
          NTSTEP=N
          TSTEP=ZERO
          itstep = 0
          IDXREV=1
          NREV=N
        ELSE
          NTSTEP=N
          IF(N.LE.1) THEN
            TSTEP=ZERO
            itstep = 0
            IDXREV=N
          ELSE
            TSTEP=TSTEP+DELTAT
            itstep = itstep + ndltat
            ERRCHK=ABS(TSTEP-TWOPI)
            IF(ERRCHK.LE.0.008) THEN
              NREV=NREV+1
              IDXREV=1
              TSTEP=ZERO
              itstep = 0
            ELSE
              IDXREV=IDXREV+1
            END IF
          END IF
        END IF

        WRITE(*,1010) N,TT(IDXREV)
C-----------------------------------------------------------------------
C     Determine initial detachment here if ISTEADY=0.         JY061300
C     Also plot wetted velocity vectors.                      JY080100
C-----------------------------------------------------------------------

        IF(IAN .EQ. 2 .AND. ISCAV .NE. 0)  GO TO 302

        IF(ISTEADY.EQ.0.AND.NTSTEP.EQ.1.AND.ISP.NE.1) THEN
          CALL DETACH
          IWET=1
          IF(IAN .NE. 2) CALL VECPLT
          IWET=0
        END IF

 302    CONTINUE
C-----------------------------------------------------------------------
C     Cal. the submergence factor of each panel at all time steps.
C     JY112299
C-----------------------------------------------------------------------
        IF(NTSTEP.EQ.1.AND.ISP.EQ.1) CALL SURFACE
C     IF(NTSTEP.EQ.1.AND.ISP.EQ.1) CALL SUBM

C-----------------------------------------------------------------------
C     Full wake alignment algorithms.
C-----------------------------------------------------------------------
C     -- Update geometry of wake and Tip vortex cavity
C     -- For ian.eq.6, xww, yww and zww are already updated in the
C        subroutine, 'wrtback_wakgeo' for all blades. 
C        | S.N.KIM | Aug. 2018.
        if((ian .eq. 2)) then

C     -- Key wake geometry

          call callgeok(2)

C     -- Read Other wakes geometries from saved file

c          call callgeow
c
c          do nn = 1, nwpanel+1
c            do m = 1 , mr + 1
c              xww(nn,m,1) = xw(nn,m)
c              yww(nn,m,1) = yw(nn,m)
c              zww(nn,m,1) = zw(nn,m)
c              do kk = 2 , nblade
c                ik = kk - 1
c                tht1 = -ik * delk
c                yww(nn,m,kk) = ywo(nn,m,ik)*cos(tht1)
c     %            -zwo(nn,m,ik)*sin(tht1)
c                zww(nn,m,kk) = ywo(nn,m,ik)*sin(tht1)
c     %            +zwo(nn,m,ik)*cos(tht1)
c                xww(nn,m,kk) = xwo(nn,m,ik)
c              enddo
c            enddo
c          enddo
          call indpot2
        endif
C-----------------------------------------------------------------------
C     solve for the cavity planform and loading characteristics at
C     this time step
C-----------------------------------------------------------------------
        CALL OPFILEC

!Allen Du 01/15/2018 output the cavity volume
        if(isteady.eq.0) then
            if(n==ntrev) opt_vol_flag=1
        else
            if(n>=(ntrev-1)*ntprev) opt_vol_flag=1
        end if

        IF(IAN .EQ. 2 .AND. NTSTEP .EQ. 1
     %    .AND. ISCAV .EQ. 2) call readforrun

C-----------------------------------------------------------------------
C     CAVITY THE ESSENCE.
C-----------------------------------------------------------------------
        CALL CAVB

!Allen Du 01/11/2018 output the cp_min
        if(opt_cpmin_c<=opt_cpmin_temp) opt_cpmin_c=opt_cpmin_temp

        DO 220 M=1,MR
          DPHI(M,IDXREV)=DELP(M)
 220    CONTINUE

        IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8) THEN
          DO IKK=1,6
            AFM(NREV,IKK) = AFM(NREV,IKK) + FBXPV(IKK)/REAL(NTPREV)
          END DO
        ENDIF

C-----------------------------------------------------------------------
C     Close matrix files for IAN=2.
C-----------------------------------------------------------------------

c        if((ian.eq.2).or.(ian.eq.6)) then
        if((ian.eq.2)) then
          close(41)
          close(42)
          close(50)
          close(51)
          close(52)
          close(53)
          do kkk = 1 , nblade+1
            ikkk = 80 + kkk
            close(ikkk)
          enddo

          do kkk = 1 , nblade
            ik1 = 90 + kkk
            ik2 = 30 + kkk
            close(ik1)
            close(ik2)
          enddo

          close(110)
          close(111)
          close(112)
          close(113)
        endif

        IF(ICON .EQ. 5) THEN
          WRITE(712,*) NREV, CLIFT,CDRAG
        ENDIF

 300  CONTINUE


C-----------------------------------------------------------------------
C     Write cavitating outputs to files.
C-----------------------------------------------------------------------

      CALL WRITEFILEC
      CALL SAVEFORRUN

      IF(ICON .NE. 5) THEN
        DO I = 1 , NTREV
          WRITE(712,7112) I, (AFM(I,J)*NBLADE, J=1,6)
        ENDDO
      ENDIF

C-----------------------------------------------------------------------
C     End of analysis.  Print total computational time
C-----------------------------------------------------------------------
      WRITE(*,1007)
 1007 Format('+------End unsteady cavity computation----------------
     *  ------')
 1008 Format(//,'  ===== TOTAL ELAPSED CPU TIME ',F10.2,' MIN')
 1009 Format(//,'  ===== CPU TIME FOR CAVITATING REVS ',F10.2,' MIN')


C.... Only output total CPU time for ISP=1. (JY112299)
      IF(ISP.EQ.1) THEN
        DELDTI = DTIME ( TARRAY )
        WRITE(*,1009) DELDTI/60.
      ELSE
        ELAPSED = ETIME ( TARRAY )
        DELDTI = DTIME ( TARRAY )

        WRITE(*,1008) ELAPSED/60.
        WRITE(*,1009) DELDTI/60.
      END IF

 1111 CONTINUE


      call write_fpg

! commented by Yiran Su 09/22/2016
!      if (ian.eq.1) call duct_ind

      IF (IAN.EQ.6) CLOSE(1108)

!Allen Du 01/06/2018 add the option of coupling with cavopt3d 
      opt_dpv1=sigma-opt_cpmin_c
      opt_dpv2=sigma-opt_cpmin_w

      open(70401,file="pass2cavopt.dat")
      write(70401,*)"0"
      write(70401,*) opt_vol, "0.0 0.0"
      write(70401,*)"0.0 0.0", opt_dpv1, opt_dpv2
      write(70401,*)opt_kt,opt_kq,"0.0  0.0"
      write(70401,*)"cpmin_w=",opt_cpmin_w
      write(70401,*)"cpmin_c=",opt_cpmin_c
      close(70401)
      END


