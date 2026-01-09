************************************************************************
*                                                                      *
*   PARAMETERS:                                                        *
*     KZ    NUMBER OF BLADES                                           *
*     NBZ   NUMBER OF PANELS OF THE BLADE CHORDWISE                    *
*     MBZ   NUMBER OF PANELS OF THE BLADE THE SPANWASE                 *
*     NWZ   NUMBER OF PANELS IN THE WAKE FOR ONE RADII                 *
*     NHBZ  NUMBER OF PANELS OF THE HUB STREAMWISE                     *
*     MHBZ  NUMBER OF PANELS OF THE HUB CIRCUMFERENTIALLY              *
*     NPANZ TOTAL NUMBER OF PANELS TO BE SOLVED ON THE BLADE AND HUB   *
*     NPAWZ NUMBER OF PANELS IN THE WAKE FOR ONE RADII USED FOR GEOM3D *
*                                                                      *
************************************************************************

C-----------------------------------------------------------------------
C     NEXT COMMON BLOCKS ARE FOR ADMINISTRATION
C-----------------------------------------------------------------------
      USE GPROD
      USE GPRODW
      USE GINT
      USE FVEL
      USE FVEL2
      USE GEOMT
      USE GEOMTP
      USE GEOMTW
      USE GEOMTWD
      USE GEOMTWO
      USE TOTG
      USE SIMQD1
      USE SIMQD2
      USE SIMQ1
      USE CORR
      USE CIRC
      USE PRES
      USE NORMAL
      USE FORCE
      USE VFOR
      USE CAVGEOD
      USE CAVGEOI
      USE CAVST
      USE CVSOL2
      USE SPLITD
      USE SPLITI
      USE NWAK
      USE UVECB
      USE UVECBO
      USE UVECW
      USE XTRAP
      USE SPIMG
      USE SPINT
      USE FXSR
      USE CAVSR
      USE TNGEO
      USE FVELTN1
      USE FVELTN2
      USE DUCTGEOM2
      USE DUCTGEOM3
      USE DUCTGEOM
      USE FVELD1
      USE FVELD2
      USE FVELD3
      USE FVELD4
      USE UNSWAKD
      USE UNS
      USE UNSDP
      USE TMPP1
      USE CAVFB
      USE WAKSUBD
      USE DWAKSUB
      USE MEMSOL
      USE CVRHS
      USE BKUTTA
      USE TOGBFLOW
      USE PODGEO1
      USE PODGEO3
      USE PODGEO4
      USE PRESNO
      USE GLOCBL
      USE FVEL3
      USE FVELW
      USE IVCON
      USE INTGV
      USE PUFBL
      USE VISBL
      USE REBL
      USE ARCLV
      USE INTGVD
      USE PUFBLD
      USE VISBLD
      USE REBLD
      USE ARCLVD
      USE m_GPROI
      USE m_FOUT
      USE m_DISKI
      USE m_DISKD
      USE m_GOUT

      USE m_CAVTY1
      USE m_CAVTY2
      USE m_CAVTY3
      USE m_CAVTY4
      USE m_CAVTY5
      USE m_CAVTY6
      USE m_CAVTY7
      USE m_TIPH1
      USE m_TIPH2
      USE m_TIPH3
      USE m_VELWAKE
      USE m_CAVCENT

      USE SU_INP
      USE PC2NSVEL
      USE UNSPROP
      use m_cavopt  !Allen DU 09.13.2018


      INCLUDE 'PARAM.INC'

      COMMON/CONPAR/ZERO,ONE,TWO,HALF,PI,HALFPI,TWOPI,HUNTPI,RAD
      COMMON/ADMN/LABEL,FN,FNSCR
      COMMON/CONS/ DELK,XUWDK,SIGMA,FROUDE,RMRTIP,RADINI
      COMMON/CTRLI/ICON,wingimag,IHUB, ITERKT,ICSPAC,IRSPAC,ITHK
     *     ,NRECL,NTIME,IWET,ICAM
     *     ,      ISTEADY,ISEARCH,ISC,IAN,IFILE,IPKT,ISP,IMG,IPHUB
     *     ,      IANG,INH,IMOD,NSAPW,IFACE,IFORMAT,ALIGN_DIV
     *     ,      WKPAN_INP, NALIGN_ITER, IPLOT, IDFWA, IREADW
      COMMON/CTRL2/SIGMAB,SIGMAF,XLTRAN,TOLKT

      COMMON /DAT1/ ADVCO,RULT,RHULT,DCD,XULT,RHUB,DTPROP
      COMMON /DAT2/ NX, NBLADE
      COMMON /DAT3/ HPCH1,HPCH2

      COMMON/INTG/ NC,MR,MRTIP,NH,NHP,NHM,NCP,MRP,NPANB
     *     ,      NPANH,NPANW,NPANEL
      COMMON/HUBI/NHBU,NHBDT,NHBX,MHBT
      COMMON/HUBD/XHBU,XHBD,XHBT,XHBTE,XHBFD
      COMMON/RUDDER/ITERGB
      COMMON/RUDDER1/FNRUDDER,IHULL
      COMMON/RUDFLAP/IFLAP,XFLAP,CHFLAP,ANGFLAP

C-----------------------------------------------------------------------
C     NEXT COMMON BLOCKS ARE FOR GEOMETRIES
C-----------------------------------------------------------------------
!     COMMON/DISKI/NHBDK,NUWDK,NSW(MBPZ),NWK
!     COMMON/DISKD/TANBUW,XHDK(MHPZ,2),YHDK(MHPZ,2)
!    *     ,      ZHDK(MHPZ,2),XWDK(MHPZ,2),YWDK(MHPZ,2),ZWDK(MHPZ,2)
!    *     ,      XHDKDT(MHPZ,2),YHDKDT(MHPZ,2),ZHDKDT(MHPZ,2)
      COMMON/GINP/XR(NXMAX),XRSQ(NXMAX),XPI(NXMAX),XRAKE(NXMAX)
     *     ,      XSKEW(NXMAX),XCHD(NXMAX),XCI(NXMAX),XTI(NXMAX)
     *     ,      XVA(NXMAX),XVR(NXMAX),XVT(NXMAX),XUA(NXMAX)
     *     ,      XUAU(NXMAX),XUT(NXMAX),XUTU(NXMAX)
!     COMMON/GOUT/SBP(NBHZ),RZP(MBZ),RZPSQ(MBZ)
      COMMON/GSPL/PICUB(MAXSPL),RKCUB(MAXSPL),SKCUB(MAXSPL)
     *     ,      CHCUB(MAXSPL),CICUB(MAXSPL),TICUB(MAXSPL)
     *     ,      VACUB(MAXSPL),VRCUB(MAXSPL),VTCUB(MAXSPL)
     *     ,      UACUB(MAXSPL),UAUCUB(MAXSPL),UTCUB(MAXSPL)
     *     ,      UTUCUB(MAXSPL)
!     COMMON/GPROI/NDLTAT,NLEP(MBZ,NSTEP,2)

C-----------------------------------------------------------------------
C     NEXT COMMON BLOCKS ARE FOR THE INFLUENCE-COEFFICIENT COMPUTATIONS
C-----------------------------------------------------------------------
      COMMON/GEOM/XV(4),YV(4),S(15),SIDE(4)

C-----------------------------------------------------------------------
C     NEXT COMMON BLOCK IS FOR THE PRESSURE KUTTA CONDITION, FORCE
C     CALCULATIONS AND OUTPUTS
C-----------------------------------------------------------------------
!     COMMON/FOUT/TT(NSTEP),XKT(NSTEP,6),XKTV(NSTEP,6)
!    *    ,       FXP(6),FXV(6),
!    *            FBXP(6),FHXP(6),FDXP(6),FTXP(6),
!    *            FBXV(6),FHXV(6),FDXV(6),
!    *            FBXPV(6),FHXPV(6),FDXPV(6),FTOTAL(6)

CVV -------------------------------------------------------------------
CVV THE COMMON BLOCKS BELOW HAVE BEEN COMMENTED OUT AND REPLACED
CVV WITH MODULES, AND ALLOCATABLE ARRAYS
CVV REVISION : JANUARY 18, 2010 (VIMAL VINAYAN)
CVV -------------------------------------------------------------------

C      COMMON/GINT/RZ(MBPZ),RZSQ(MBPZ),PITCH(MBPZ),RAKE(MBPZ)
C     *     ,      SKEW(MBPZ),CHORD(MBPZ),CAMBR(MBPZ),THICK(MBPZ)
C      COMMON/GPROD/XB(NBPZ,MBPZ),YB(NBPZ,MBPZ),ZB(NBPZ,MBPZ)
C     *     ,      XW(NWPZ,MBPZ),YW(NWPZ,MBPZ),ZW(NWPZ,MBPZ)
C     *     ,      XTW,YTW,ZTW
C     *     ,      XH(NHPZ,MHPZ),YH(NHPZ,MHPZ),ZH(NHPZ,MHPZ)
C     *     ,      ZW1(NWPZ,MBPZ),YW1(NWPZ,MBPZ),XW1(NWPZ,MBPZ)
C     *     ,      THT(MBPZ),THR(NBPZ),SB(NBHZ),YTC(NBHZ),YCC(NBHZ)
C     *     ,      TCC(NBHZ),XI(NBPZ),ETA(NBPZ),CONMRX(3,NBPZ,MBPZ)
C     *     ,     VECMRX(3,NBPZ,MBPZ),WAKVEC(3,NWPZ,MBPZ)
C     *     ,     ARCLNG(NSCZP,MBZ,2),CAVINI
C     *     ,      CAVINIB,CAVINIF
C      COMMON/GPRODW/XWO(NWPZ,MBPZ,KZ),YWO(NWPZ,MBPZ,KZ)
C     *     ,      ZWO(NWPZ,MBPZ,KZ),XWW(NWPZ,MBPZ,KZ)
C     *     ,      YWW(NWPZ,MBPZ,KZ),ZWW(NWPZ,MBPZ,KZ)
CVV -------------------------------------------------------------------

C-----------------------------------------------------------------------
C     NEXT COMMON BLOCK IS FOR VELOCITIES
C-----------------------------------------------------------------------
C      COMMON/FVEL/VOXW(NSCWZ),VOYW(NSCWZ),VOZW(NSCWZ)
C     *     ,      VTOTS(NPANZ),VXIB(NBZ,MBZ),VETAB(NBZ,MBZ)
C     *     ,      VINFSB(NBZ,MBZ),VXIH(NHBZ,MHBZ),VETAH(NHBZ,MHBZ)
C     *     ,      VINFSH(NHBZ,MHBZ),UXTOT(NBZ,MBZ),UYTOT(NBZ,MBZ)
C     *     ,      UZTOT(NBZ,MBZ),UXHTOT(NHBZ,MHBZ),UYHTOT(NHBZ,MHBZ)
C     *     ,      UZHTOT(NHBZ,MHBZ)
C     *     ,      VOX1(NPANZ),VOY1(NPANZ),VOZ1(NPANZ)
C      COMMON/FVEL2/ UXTHTOT(NHMX,MCAVM),UYTHTOT(NHMX,MCAVM)
C     *     ,      UZTHTOT(NHMX,MCAVM), UXCHTOT(NCAVM,MCAVM)
C     *     ,      UYCHTOT(NCAVM,MCAVM),UZCHTOT(NCAVM,MCAVM)

C-----------------------------------------------------------------------
C     NEXT COMMON BLOCKS ARE FOR THE INFLUENCE-COEFFICIENT COMPUTATIONS
C-----------------------------------------------------------------------
C      COMMON/GEOMT/XVP(NPANZ,4),YVP(NPANZ,4),XCT(NPANZ,3),DIR(NPANZ,3,3)
C     *     ,       SS(NPANZ,15),SID(NPANZ,4),GSX(NPANZ,4,3)
C     *     ,       GSW(NPANZ,4),VEL(NPANZ,6),DELU(NPANZ),DELV(NPANZ)
C     *     ,       COSPHI(NPANZ),SINPHI(NPANZ)
C.....HONG ADDED  DIRCP, XCTPW, DIRWCP TO GEOMTP 01/20/2005
C.....HONG ADDED  XCTPWS, DIRWSCP TO GEOMTP 04/06/2005
C      COMMON/GEOMTP/XCTP(NPANZ,3,KZ),DIRCP(NPANZ,3,3,KZ)
C     *     ,        XCTPW(NSCWZ,3,KZ),DIRWCP(NSCWZ,3,3,KZ)
C     *     ,        XCTPWS(NSCWZ,3,KZ),DIRWSCP(NSCWZ,3,3,KZ)
C.....HONG ADDED XCTPDW,DIRDWCP,XCTPDWS,DIRDWSCP     06/04/06
C     *     ,        XCTPDW(NDWMN,3,KZ),DIRDWCP(NDWMN,3,3,KZ)
C     *     ,        XCTPDWS(NSCWZ,3,KZ),DIRDWSCP(NSCWZ,3,3,KZ)
C     *     ,        XCPW(NSCWZ,3,KZ)
C      COMMON/GEOMTW/XVPW(NSCWZ,4),YVPW(NSCWZ,4),XCTW(NSCWZ,3)
C     *     ,        DIRW(NSCWZ,3,3),SSW(NSCWZ,15),SIDW(NSCWZ,4)
C     *     ,        GSXW(NSCWZ,4,3),GSWW(NSCWZ,4),VELW(NSCWZ,6)
C     *     ,        DELUW(NSCWZ),DELVW(NSCWZ),COSPHIW(NSCWZ)
C     *     ,        SINPHIW(NSCWZ)
C      COMMON/GEOMTWD/XVPWD(NDWMN,4),YVPWD(NDWMN,4),XCTWD(NDWMN,3)
C     *     ,        DIRWD(NDWMN,3,3),SSWD(NDWMN,15),SIDWD(NDWMN,4)
C     *     ,        GSXWD(NDWMN,4,3),GSWWD(NDWMN,4),VELWD(NDWMN,6)
C     *     ,        DELUWD(NDWMN),DELVWD(NDWMN),COSPHIWD(NDWMN)
C     *     ,        SINPHIWD(NDWMN)
C      COMMON/GEOMTWO/XVPWO(NSCWZ,4,KZ),YVPWO(NSCWZ,4,KZ)
C     *             ,XCTWO(NSCWZ,3,KZ),DIRWO(NSCWZ,3,3,KZ)
C     *              ,SSWO(NSCWZ,15,KZ),SIDWO(NSCWZ,4,KZ)
C     *              ,GSXWO(NSCWZ,4,3,KZ),GSWWO(NSCWZ,4,KZ)
C     *              ,VELWO(NSCWZ,6,KZ),DELUWO(NSCWZ,KZ)
C     *              ,DELVWO(NSCWZ,KZ),COSPHIWO(NSCWZ,KZ)
C     *              ,SINPHIWO(NSCWZ,KZ)

C      COMMON/TOTG/XG(NPANZ,4,3),CHRLEPS(NPANZ),CHRLEWS(NSCWZ)
C     *     ,      XGW(NSCWZ,4,3),CHRLEWSO(NSCWZ,KZ)
C     *     ,      XGWD(NDWMN,4,3),CHRLEWSD(NDWMN)

C-----------------------------------------------------------------------
C     NEXT COMMON BLOCK IS FOR THE MATRIX SOLUTION
C-----------------------------------------------------------------------
C      COMMON/SIMQD1/A(NPANZ),B(NPANZ),W(NPANZ,MBZ)
C     *     ,      CHDK(NPANZ),WDK(NPANZ),POT(NPANZ),DPDN(NPANZ)
C     *     ,      WINF(NPANZ,MBZ),CHDKDT(NPANZ),WD(NPANZ,MDMAX)
C     *     ,      WKD(NPAWZ,MDMAX)
C      COMMON/SIMQD2/AA(NTZ,NTZ),BB(NPANZ,NPANZ),C(NPANZ,NPAWZ)
C     *     ,        D(NPAWZ,NPANZ),E(NPAWZ,NPANZ),F(NPAWZ,NPAWZ)
C     *     ,        WK(NPAWZ,MBZ),WK2(NPAWZ,MBZ)
C     *     ,        WDK1(NPAWZ),CHDK1(NPAWZ),CHDKDT1(NPAWZ)
C      COMMON/SIMQI/NPERB(NBLKMAX)

C-----------------------------------------------------------------------
C     NEXT COMMON BLOCK IS FOR THE PRESSURE KUTTA CONDITION, FORCE
C     CALCULATIONS AND OUTPUTS
C-----------------------------------------------------------------------
C      COMMON/CORR/DELCP(MBZ),DCPPRE(MBZ)
C      COMMON/CIRC/DELP(MBZ),DELPD(MDMAX)
C      COMMON/PRES/CPB(NBZ,MBZ),CPH(NHBZ,MHBZ),DPDUB(NBZ,MBZ)
C     *     ,      DPDVB(NBZ,MBZ),DPDUH(NHBZ,MHBZ),DPDVH(NHBZ,MHBZ)
C      COMMON/NORMAL/XCON(NBHZP,MBPZ),YCON(NBHZP,MBPZ),ZCON(NBHZP,MBPZ)
C      COMMON/FORCE/BKFX(MBZ),BKFY(MBZ),BKFZ(MBZ),BKMX(MBZ),BKMY(MBZ)
C     *     ,       BKMZ(MBZ),HKFX(NHBZ),HKFY(NHBZ),HKFZ(NHBZ)
C     *     ,       HKMX(NHBZ),HKMY(NHBZ),HKMZ(NHBZ)
C      COMMON/VFOR/FVA(MBZ),XCDF

C-----------------------------------------------------------------------
C    NEXT COMMON BLOCKS ARE FOR THE CAVITY SOLUTION       5-1-92 NF
C-----------------------------------------------------------------------
C      COMMON/CAVGEOD/DS(NBSCZ,MBZ),DZ(NBSCZ,MBZ)
C     *     ,        SZ(NBSCZP,MBZ),SPZ(NBSCZ,MBZ)
C     *     ,        RLAM(MBZ,2),CAVL(MBZ,2),CAVLP(MBZ,2),DZW(NWZ,MBZ)
C     *     ,        SOP(MBZ,2),CAVLSA(MBZ,NSTEP,2)
C      COMMON/CAVGEOI/M0(MBZ,2),ITER,JCV(MBZ,2),LCV(MBZ,2)
C     *     ,        NOCAV(MBZ,0:NSTEP,2),IT2
C      COMMON/CAVST/QC(NSCZP,MBZ,2),PHI0(MBZ,2)
C     *     ,        PHI1(NSCZP,MBZ,2),DPHIDS(NSCZP,MBZ,2)
C     *     ,        DPHIDV(NBHZ,MBZ,2)
C      COMMON/CVSOL2/HT(NSCZP,MBZ,2),DPDNC(NPANZ),DELTA(MBZ,2)
C     *     ,        HTW1(MBZ),DELTAP(MBZ,2),SOL(NTZ),SORW(NPAWZ)
C     *     ,        POTW(NPAWZ),HTP(NBHZP,MBZ,2),HTWP(NWZ+1,MBZ)
C      COMMON/SPLITD/FLP(MBZ,2),FRP(MBZ,2),FLS(MBZ,2),FRS(MBZ,2)
C     *     ,       DZL(MBZ,2),DZR(MBZ,2),QSPR(MBZ,2),DLISP(MBZ,2)
C     *     ,       DPDVSP(MBZ,2),PHIL(MBZ,2),PHIR(MBZ,2),QSPL(MBZ,2)
C     *     ,       QSSR(MBZ)
C      COMMON/SPLITI/NSPP(MBZ,2),NSPS(MBZ,2)
C      COMMON/NWAK/NWC(MBZ,2),NWPAN(MBZ)
C      COMMON/UVECB/UL(NPANZ,3),VL(NPANZ,3),WL(NPANZ,3)
C      COMMON/UVECBO/ULO(NPANZ,3,KZ),VLO(NPANZ,3,KZ)
C     *     ,        WLO(NPANZ,3,KZ)
C      COMMON/UVECW/ULW(NSCWZ,3)
C      COMMON/XTRAP/CT(MBZ,4,2),DT(MBZ,4,2),QT(MBZ,4,2),QW(MBZ,4,2)

      COMMON/NGRN/NBW,NTRA
      COMMON/VOLUME/CVOL,ANGOUT
      COMMON/DELTOL/DTOL
      COMMON/STSL/STPR,SLPR,RLAMDA1,AFJ
      COMMON/CNTRL/ITERMAX
      COMMON/FOIL/ALPHA,CLIFT,CDRAG
      COMMON/VISCOUS/ IVPITCH

C -- BEGIN TIP HSLEE(10/14/99)
C-- DEFINE VARIABLES FOR THE TIP VORTEX CAVITY

      COMMON /MTIP1/  NTHX, NCVX, MCVT,NTHXP,NCVXP,MCVTP
      COMMON /MTIP2/  NPANT, NPANC, NWPANEL, NPWAKE
      COMMON /MTIPCONST/ ICAVT, ICAVMAX, ISCAV, NNREV, IUALIGN_DIV, IUPLOT

      COMMON/DISPL/ DELTAM,RC1,GAMMA
!     COMMON/CAVTY1/ XVC(NCAVMP,MCAVMP),YVC(NCAVMP,MCAVMP),
!    *               ZVC(NCAVMP,MCAVMP)
!     COMMON/CAVTY2/ XIND(NPAWZ),YIND(NPAWZ),ZIND(NPAWZ)
!     COMMON/CAVTY3/ VXMEAN(NCAVM),VYMEAN(NCAVM),VZMEAN(NCAVM)
!    *              ,CMX(NPHC),CMY(NPHC),CMZ(NPHC)
!     COMMON/CAVTY4/ CPC(NCAVM,MCAVM),VTC(NCAVM,MCAVM)
!     COMMON/CAVTY5/ VXIC(NCAVM,MCAVM),VETAC(NCAVM,MCAVM),
!    *               VINFSC(NCAVM,MCAVM)
!     COMMON/CAVTY6/ DPDUC(NCAVM,MCAVM),DPDVC(NCAVM,MCAVM)
!     COMMON/CAVTY7/ UXCTOT(NCAVM,MCAVM),UYCTOT(NCAVM,MCAVM),
!    *               UZCTOT(NCAVM,MCAVM)

!     COMMON /TIPH1/ XCH(NHMXP,MCAVMP),YCH(NHMXP,MCAVMP),
!    *               ZCH(NHMXP,MCAVMP)

!     COMMON /TIPH2/ DPDUTH(NHMX,MCAVM),DPDVTH(NHMX,MCAVM),
!    *               CPTH(NHMX,MCAVM),VTTH(NHMX,MCAVM),
!    *               VINFSTH(NHMX,MCAVM)
!     COMMON /TIPH3/ VXITH(NHMX,MCAVM),VETATH(NHMX,MCAVM)
!     COMMON /VELWAKE/ WINDX(NWMAX,MRMAXP),WINDY(NWMAX,MRMAXP),
!    *                 WINDZ(NWMAX,MRMAXP)
!     COMMON /CAVCENT/ XVCC(NCAVMP),YVCC(NCAVMP),ZVCC(NCAVMP)

C -- END OF TIP (10/14/99)

      CHARACTER*72 LABEL
      CHARACTER*29 FN,FNSCR

C-----------------------------------------------------------------------
C     COMMON BLOCKS FOR ISP=1
C-----------------------------------------------------------------------
      COMMON/SPDAT/YFS,SRATIO,SPANGLE

C      COMMON/SPINT/ICB(MBZ,2,NSTEP),MSW(MBZ,NSTEP),ICW(MBZ,NSTEP),
C     *             IC(2,MBZ,NSTEP),IW(2,MBZ,NSTEP),
C     *             ISUBM(NPBZ+NPHZ,NSTEP),ISUBH(NHBZ,NSTEP),
C     *             NPERM(MBZ),ISUB1(MBZ,2,NSTEP),ISUWM(NPAWZ,NSTEP)
C      COMMON/SPIMG/XCTP_IM(NPBZ+NPHZ,3,KZ),XCPW_IM(NPAWZ,3,KZ)

C-----------------------------------------------------------------------
C     COMMON BLOCKS FOR ISC=1
C-----------------------------------------------------------------------
       COMMON /IXSR/ NSR,NSR2,NTSR,NSRP,NSR2P,NCOLD,NHOLD,N0(2)

C      COMMON /FXSR/ SSTE(MBZ),VELTE(MBZ,6)
!      PARAMETER(NCSRZ=NZSR2+NWZ,NCSRZP=NCSRZ+1)
C      COMMON /CAVSR/ CT2(MBZ,4,2),QCSR(NCSRZP,MBZ,2),PHI0T(MBZ,2)
C     *     ,         PSI0T(MBZ,2),PHI2(NCSRZ,MBZ,2)
C     *     ,         DPHIDS2(NCSRZ,MBZ,2),DPHIDV2(NZSR2,MBZ,2)
C     *     ,         DELTAT1(MBZ),DELTATP(MBZ),CAVLT(MBZ)
C     *     ,         CAVLTP(MBZ),CAVLTSA(MBZ,NSTEP)

C ---------------------------------------------------------------------
C       COMMON BLOCK FOR ITUN = 1
C ---------------------------------------------------------------------
      COMMON /TUNINPUT/ TUNRAD,XFOR,XAFT,XMIDF,XMIDA
      COMMON /ITUNCON1/ ITUNGEO,ITUN
      COMMON /ITUNCON2/ NAXT,MTUNEL,NPANTN
      COMMON /ITUNCON3/ NAX,NAXF,NAXA,NAXM,NSIDE
      COMMON /TUNINFLOW/ AREAR

C       COMMON /TNGEO/ XTUN(NTMAXT+1,NCMAXT+1),YTUN(NTMAXT+1,NCMAXT+1),
C     %                 ZTUN(NTMAXT+1,NCMAXT+1)
C        COMMON /FVELTN1/ VXITN(NTMAXT,NCMAXT),VETATN(NTMAXT,NCMAXT),
C     %                   VINFSTN(NTMAXT,NCMAXT)
C        COMMON /FVELTN2/ DPDUTN(NTMAXT,NCMAXT),DPDVTN(NTMAXT,NCMAXT),
C     %                   UXTNTOT(NTMAXT,NCMAXT),UYTNTOT(NTMAXT,NCMAXT),
C     %                   UZTNTOT(NTMAXT,NCMAXT),CPTN(NTMAXT,NCMAXT)

C ---------------------------------------------------------------------
C       COMMON BLOCK FOR IDUCT = 1
C ---------------------------------------------------------------------

      COMMON /IDUCTCON1/ IDUCT,IDGEO,IDOPT,IDALIGN,IREPANEL
      COMMON /IDUCTCON2/ NDUCT, NDUCTH, NDUCTHP, MDUCT,
     %                   NDUCTP, MDUCTP,NPAND,NWAKEP
      COMMON /IDUCTCON3/ NDFWD, NDAFT, NDBACK, NDDAT
      COMMON /IDUCTCON4/ NDWK, NDWKP,NPWAKED
      COMMON /DUCTGEOM1/ DUCTGAP,DWAKEL,DUCTPT,
C.......HONG ADDED DCHORD,XLED,YLED,DUCTANG TO COMMON BLOCK DUCTGEOM
     *              DCHORD,XLED,YLED,DUCTANG
C       COMMON /DUCTGEOM2/ XXCP(NDMAX),RRCP(NDMAX)
C       COMMON /DUCTGEOM3/ XID(201),ETAD(201)
C       COMMON /DUCTGEOM/ XD(NDMAXP,MDMAXP), YD(NDMAXP,MDMAXP),
C     %                    ZD(NDMAXP,MDMAXP),XDW(NDWMXP,MDMAXP),
C     %                    YDW(NDWMXP,MDMAXP),ZDW(NDWMXP,MDMAXP),
C     %                    XDWDK(MDMAXP,2),YDWDK(MDMAXP,2),
C     %                    ZDWDK(MDMAXP,2),WDKD(NPANZ)
CC    %                    ZDWDK(MDMAXP,2),WDKD(MDMAX)
C        COMMON /FVELD1/ VXID(NDMAX,MDMAX),VETAD(NDMAX,MDMAX),
C     %                   VINFSD(NDMAX,MDMAX)
C        COMMON /FVELD2/ DPDUD(NDMAX,MDMAX),DPDVD(NDMAX,MDMAX),
C     %                   UXDTOT(NDMAX,MDMAX),UYDTOT(NDMAX,MDMAX),
C     %                   UZDTOT(NDMAX,MDMAX),CPD(NDMAX,MDMAX)
C       COMMON /FVELD3/ VVYD(NPANMD),VVZD(NPANMD)
C       COMMON /FVELD4/ DPOTMEAN(NDMAX),DCPMEAN(NDMAX),DPDUMEAN(NDMAX)

C--------------------------------------------------------------
C -- FOR VARYING HUB
C--------------------------------------------------------------

C      COMMON /PODGEO1/ XHIN(200),YHIN(200),XYHINCUB(800)
      COMMON /PODGEO2/ NHIN
C      COMMON /PODGEO3/ HRZ(NBPZ,MBPZ),HRZSQ(NBPZ,MBPZ)
C      COMMON /PODGEO4/ HRZP(NBZ,MBZ),HRZPSQ(NBZ,MBZ)

C ------------------------------------------------------
C       FOR FORCE EXCLUDED STATIC PRESSURE
C ------------------------------------------------------
C      COMMON/PRESNO/CPBN(NBZ,MBZ),CPHN(NHBZ,MHBZ),
C     *           CPTNN(NTMAXT,NCMAXT),CPDN(NDMAX,MDMAX),
C     *           CPTHN(NHMX,MCAVM),CPCN(NCAVM,MCAVM)

C-----------------------------------------------------------------------
C     COMMON BLOCKS FOR IVISC=1
C     HONG ADDED FOR COUPLING WITH XFOIL       APRIL, 2006
C-----------------------------------------------------------------------
      COMMON /IBL/ IVISC, ITER_VI
      COMMON /BLTOL/ REYD,RVCRIT,XTRANS,XTRANP,MAXIT,EPS1
C      COMMON /GLOCBL/ XIV(NBPZ,MBPZ),ETAV(NBPZ,MBPZ)
C     *     ,          XIVW(NWPZ,MBPZ),ETAVW(NWPZ,MBPZ)
C     *     ,          XIVWS(NWZ+NDEL,MBPZ),ETAVWS(NWZ+NDEL,MBPZ)
C     *     ,          XIVIS(NZVIS,MBZ),ETAVIS(NZVIS,MBZ)
C.....FOR DUCT CASE
C     *     ,          XIVD(NDMAXP,MDMAXP),ETAVD(NDMAXP,MDMAXP)
C     *     ,          XIVDW(NDMAXP,MDMAXP),ETAVDW(NDMAXP,MDMAXP)
C     *     ,          XIVDWS(NWZ+NDEL,MDMAXP),ETAVDWS(NWZ+NDEL,MDMAXP)
C     *     ,          XIVISD(NZVISD,MDMAX),ETAVISD(NZVISD,MDMAX)

C.... FVEL3 BLOCK FOR WAKE PERTURBATION POTENTIAL  04/13/2005
C      COMMON/FVEL3/PPOTW(NSCWZ),PPOTWS(NSCWZ)
C     *     ,       VTOTS_V(NPANZ)
C     *     ,       UXTOT_V(NBZ,MBZ),UYTOT_V(NBZ,MBZ),UZTOT_V(NBZ,MBZ)
C..... FOR DUCT CASE
C     *     ,       PPOTDW(NSCWZ),PPOTDWS(NSCWZ)
C     *     ,       UXDTOT_V(NDMAX,MDMAX),UYDTOT_V(NDMAX,MDMAX)
C     *     ,       UZDTOT_V(NDMAX,MDMAX)

C      COMMON/FVELW/VXIW(NWZ,MBZ),VETAW(NWZ,MBZ),VINFSW(NWZ,MBZ)
C     *     ,       VXIWS(NWZ+NDEL,MBZ),VETAWS(NWZ+NDEL,MBZ)
C     *     ,       VINFSWS(NWZ+NDEL,MBZ)
C     *     ,       CPBL(NBZ,MBZ)
C.....FOR DUCT CASE
C     *     ,       VXIDW(NWZ,MDMAX),VETADW(NWZ,MDMAX)
C     *     ,       VXIDWS(NWZ+NDEL,MDMAX),VETADWS(NWZ+NDEL,MDMAX)
C     *     ,       VINFSDWS(NWZ+NDEL,MDMAX),VINFSDW(NWZ,MDMAX)
C     *     ,       CPBLD(NDMAX,MDMAX)

C      COMMON/IVCON/LVFLAG(MBPZ),LVFLAGD(MDMAX)

C<<<<<<<<<<<<<<<<<<<<<<<<<END OF PUFCAV.INC>>>>>>>>>>>>>>>>>>>>>>>>>>>>>







