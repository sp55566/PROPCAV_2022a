************************************************************************
*     PROPCAV  R2019 Release                                           *
*     Unsteady cavitating propeller analysis program                   *
*----------------------------------------------------------------------*
*     PROPCAV  R2022 Release                                           *
*     Unsteady cavitating propeller analysis program                   *
*----------------------------------------------------------------------*
*     COPYRIGHT (C) THE UNIVERSITY OF TEXAS AT AUSTIN                  *
*        November  2022  P.I. Dr. Spyros A. Kinnas                     *
*            MARINE HYDRODYNAMICS LABORATORY                           *
*----------------------------------------------------------------------*
      MODULE m_param
      implicit none
      integer KZ
      integer NSTEP, NNDIM
      integer NDEL
      integer NBZ,MBZ,NBPZ,MBPZ,NBHZ,NBHZP
      integer NSCZ,NSCZP,NBSCZ,NBSCZP
      integer NPBZ,NPANZ
      integer NTZ
      integer NBLKMAX
      integer NPHZ
      integer nphc, NPHT
      integer NPAWZ,NPAWZ2
      integer NSCWZ
      integer mrmax , nwmax
      integer ncavm , ncavmp
      integer mcavm, mcavmp
      integer nhmx, nhmxp
      integer mrmaxp, nwmaxp
      integer NWZVIS, NZVIS
      integer NDWZVIS,NZVISD
      integer NHUZ,NHDZ,NHBZ,NHPZ
      integer MHBZ,MHPZ
      integer NZWSUB
      integer NWZ,NWPZ
      integer NCSRZ,NCSRZP
      integer NtMAXt, NCMAXT
      integer NPANMt
      integer NDMAX, MDMAX,NDWMX
      integer NDMAXP, MDMAXP, NDWMXP
      integer NPANMD, NDWMN
      END MODULE m_param

      MODULE su_inp
      implicit none
      integer IDRT
      integer INRPAN
      END MODULE su_inp

!Allen Du 01/06/2018 add the option of coupling with cavopt3d 
      MODULE m_cavopt
          implicit none
          real :: opt_kt=0.0
          real :: opt_kq=0.0
          real :: opt_cpmin_w=0.0
          real :: opt_cpmin_c=0.0
          real :: opt_cpmin_temp=0.0
          real :: opt_dpv=0.0
          real :: opt_vol=0.0
          integer:: opt_vol_flag=0
      END MODULE m_cavopt

      MODULE m_gproi
      implicit none
      integer NDLTAT
      integer,allocatable :: NLEP(:,:,:)
      END MODULE m_gproi

      MODULE m_FOUT
        implicit none
!       real TT(NSTEP),XKT(NSTEP,6),XKTV(NSTEP,6)
        real,allocatable :: TT(:),XKT(:,:),XKTV(:,:)
        real FXP(6),FXV(6)
        real FBXP(6),FHXP(6),FDXP(6),FTXP(6)
        real FBXV(6),FHXV(6),FDXV(6)
        real FBXPV(6),FHXPV(6),FDXPV(6),FTOTAL(6)

      END MODULE m_FOUT

      MODULE m_P3DAT1
        implicit none
        integer,parameter :: N4=45 ! 30 S.N.KIM | Aug. 2018
        real,allocatable  :: F(:,:),C(:,:),S(:,:),FT(:,:)
        real              :: AO(6),A(N4,6),B(N4,6)
        real AMP(N4,6),PH(N4,6),HAMP(N4),HPH(N4),VAMP(N4)
        real VPH(N4)
      END MODULE m_P3DAT1

      MODULE m_PSXYZ
        implicit none
        real,allocatable :: xps(:,:,:),yps(:,:,:),zps(:,:,:)
        real,allocatable :: ipsn(:,:),nel(:)
      END MODULE m_PSXYZ

      module m_DISKI
        implicit none
        integer NHBDK,NUWDK,NWK
        integer,allocatable :: NSW(:)
      end module m_DISKI

      module m_DISKD
        implicit none
        real TANBUW
        real,allocatable :: XHDK(:,:),YHDK(:,:)
        real,allocatable :: ZHDK(:,:),XWDK(:,:)
        real,allocatable :: YWDK(:,:),ZWDK(:,:)
        real,allocatable :: XHDKDT(:,:),YHDKDT(:,:),ZHDKDT(:,:)
      end module m_DISKD

      module m_GOUT
        implicit none
        real,allocatable::SBP(:),RZP(:),RZPSQ(:)
      end module m_GOUT

      module m_CAVTY1
        implicit none
        real,allocatable ::XVC(:,:),YVC(:,:),ZVC(:,:)
      end module m_CAVTY1

      module m_CAVTY2
        implicit none
        real,allocatable ::XIND(:),YIND(:),ZIND(:)
      end module m_CAVTY2

      module m_CAVTY3
        implicit none
        real,allocatable ::VXMEAN(:),VYMEAN(:),VZMEAN(:)
        real,allocatable ::CMX(:),CMY(:),CMZ(:)
      end module m_CAVTY3

      module m_CAVTY4
        implicit none
        real,allocatable :: CPC(:,:),VTC(:,:)
      end module

      module m_CAVTY5
        implicit none
        real,allocatable ::VXIC(:,:),VETAC(:,:), VINFSC(:,:)
      end module m_CAVTY5

      module m_CAVTY6
        implicit none
        real,allocatable :: DPDUC(:,:),DPDVC(:,:)
      end module m_CAVTY6

      module m_CAVTY7
        implicit none
        real,allocatable :: UXCTOT(:,:),UYCTOT(:,:),UZCTOT(:,:)
      end module

      module  m_TIPH1
        implicit none
        real,allocatable ::XCH(:,:),YCH(:,:),  ZCH(:,:)
      end module

      module  m_TIPH2
        implicit none
        real,allocatable ::DPDUTH(:,:),DPDVTH(:,:)
        real,allocatable ::CPTH(:,:),VTTH(:,:), VINFSTH(:,:)
      end module

      module  m_TIPH3
        implicit none
        real,allocatable :: VXITH(:,:),VETATH(:,:)
      end module

      module  m_VELWAKE
        implicit none
        real,allocatable :: WINDX(:,:),WINDY(:,:),WINDZ(:,:)
      end module

      module  m_CAVCENT
        implicit none
        real,allocatable :: XVCC(:),YVCC(:),ZVCC(:)
      end module

      module m_WAKGEO
        implicit none
        real,allocatable :: RTBL(:,:),XTBL(:),RW(:),VAR(:),VTR(:),
     *     UASTAR(:),UTSTAR(:),UAUSTR(:),UTUSTR(:)
      end module m_WAKGEO

      module m_INPGEO2
       implicit none
       real,allocatable :: CUBTMP(:),XCS1(:,:),XTS1(:,:),XCS0(:),XTS0(:)
      end module m_INPGEO2

      module m_BLADEM
        implicit none
        real,allocatable ::XBM(:,:),YBM(:,:),ZBM(:,:)
      end module m_BLADEM

      module m_NACAIO
        implicit none
        real,allocatable :: XKEY(:),YKEY(:),ZKEY(:)
     #     ,XNEAR(:),YNEAR(:),ZNEAR(:),CROSS(:,:)
      end module m_NACAIO

      module m_WKNP
        implicit none
        integer,allocatable ::NWIDX(:),NWSEC(:)
      end module m_WKNP

      MODULE GPROD
      ALLOCATABLE :: XB(:,:),YB(:,:),ZB(:,:)
      ALLOCATABLE :: XW(:,:),YW(:,:),ZW(:,:)
      ALLOCATABLE :: XX_TMP(:,:),YY_TMP(:,:),ZZ_TMP(:,:)
c XM YU 12/2011
c add for wake geometry including subpanels
      allocatable :: xwvs3d(:,:),ywvs3d(:,:),zwvs3d(:,:)
C/s S.N.KIM | Aug. 2018.
      allocatable :: xwvs3do(:,:),ywvs3do(:,:),zwvs3do(:,:)
C/e S.N.KIM 
c XM YU 12/2011
      ALLOCATABLE :: XH(:,:),YH(:,:),ZH(:,:)
      ALLOCATABLE :: XW1(:,:),YW1(:,:),ZW1(:,:)
      ALLOCATABLE :: THT(:),THR(:),SB(:),YTC(:),YCC(:),TCC(:)
      ALLOCATABLE :: XI(:),ETA(:)
      ALLOCATABLE :: CONMRX(:,:,:),VECMRX(:,:,:),WAKVEC(:,:,:)
      ALLOCATABLE :: ARCLNG(:,:,:)
      REAL CAVINI,CAVINIB,CAVINIF
      END MODULE GPROD

      MODULE GPRODW
      ALLOCATABLE XWO(:,:,:),YWO(:,:,:),ZWO(:,:,:)
      ALLOCATABLE XWW(:,:,:),YWW(:,:,:),ZWW(:,:,:)
C/s S.N.KIM | Aug. 2018.
      ALLOCATABLE XWW1(:,:,:),YWW1(:,:,:),ZWW1(:,:,:)
C/e S.N.KIM
      END MODULE GPRODW

      MODULE GINT
      ALLOCATABLE RZ(:),RZSQ(:),PITCH(:),RAKE(:)
      ALLOCATABLE SKEW(:),CHORD(:),CAMBR(:),THICK(:)
      END MODULE GINT

      MODULE FVEL
      ALLOCATABLE VOXW(:),VOYW(:),VOZW(:),VTOTS(:)
      ALLOCATABLE VXIB(:,:),VETAB(:,:)
      ALLOCATABLE VINFSB(:,:),VXIH(:,:),VETAH(:,:)
      ALLOCATABLE VINFSH(:,:)
      ALLOCATABLE UXTOT(:,:),UYTOT(:,:),UZTOT(:,:)
      ALLOCATABLE UXHTOT(:,:),UYHTOT(:,:),UZHTOT(:,:)
      ALLOCATABLE VOX1(:),VOY1(:),VOZ1(:)
      END MODULE FVEL

      MODULE FVEL2
      ALLOCATABLE UXTHTOT(:,:),UYTHTOT(:,:)
      ALLOCATABLE UZTHTOT(:,:),UXCHTOT(:,:)
      ALLOCATABLE UYCHTOT(:,:),UZCHTOT(:,:)
      END MODULE FVEL2

      MODULE GEOMT
      ALLOCATABLE XVP(:,:),YVP(:,:),XCT(:,:),DIR(:,:,:)
      ALLOCATABLE SS(:,:),SID(:,:),GSX(:,:,:)
      ALLOCATABLE GSW(:,:),VEL(:,:),DELU(:),DELV(:)
      ALLOCATABLE COSPHI(:),SINPHI(:)
      END MODULE GEOMT

      MODULE GEOMTP
      ALLOCATABLE XCTP(:,:,:),DIRCP(:,:,:,:)
      ALLOCATABLE XCTPW(:,:,:),DIRWCP(:,:,:,:)
      ALLOCATABLE XCTPWs(:,:,:),DIRWsCP(:,:,:,:)
c XM YU
      allocatable xctpwv(:,:,:)
      ALLOCATABLE XCTPDW(:,:,:),DIRDWCP(:,:,:,:)
      ALLOCATABLE XCTPDWs(:,:,:),DIRDWsCP(:,:,:,:)
      ALLOCATABLE XCPW(:,:,:)
      END MODULE GEOMTP

      MODULE GEOMTW
      ALLOCATABLE XVPW(:,:),YVPW(:,:),XCTW(:,:)
      ALLOCATABLE DIRW(:,:,:),SSW(:,:),SIDW(:,:)
      ALLOCATABLE GSXW(:,:,:),GSWW(:,:),VELW(:,:)
      ALLOCATABLE DELUW(:),DELVW(:),COSPHIW(:)
      ALLOCATABLE SINPHIW(:)
      END MODULE GEOMTW

      MODULE GEOMTWD
      ALLOCATABLE XVPWD(:,:),YVPWD(:,:),XCTWD(:,:)
      ALLOCATABLE DIRWD(:,:,:),SSWD(:,:),SIDWD(:,:)
      ALLOCATABLE GSXWD(:,:,:),GSWWD(:,:),VELWD(:,:)
      ALLOCATABLE DELUWD(:),DELVWD(:),COSPHIWD(:)
      ALLOCATABLE SINPHIWD(:)
      END MODULE GEOMTWD

      MODULE GEOMTWO
      ALLOCATABLE XVPWO(:,:,:),YVPWO(:,:,:),XCTWO(:,:,:)
      ALLOCATABLE DIRWO(:,:,:,:),SSWO(:,:,:),SIDWO(:,:,:)
      ALLOCATABLE GSXWO(:,:,:,:),GSWWO(:,:,:),VELWO(:,:,:)
      ALLOCATABLE DELUWO(:,:),DELVWO(:,:),COSPHIWO(:,:)
      ALLOCATABLE SINPHIWO(:,:)
      END MODULE GEOMTWO

      MODULE TOTG
      ALLOCATABLE XG(:,:,:),CHRLEPS(:),CHRLEWS(:)
      ALLOCATABLE XGW(:,:,:),chrlewso(:,:)
      ALLOCATABLE XGWD(:,:,:),CHRLEWSD(:)

      ALLOCATABLE XG1(:,:,:),XG11(:,:,:),R111(:)
      END MODULE TOTG

      MODULE SIMQD1
      ALLOCATABLE A(:),B(:),W(:,:)
      ALLOCATABLE CHDK(:),WDK(:),POT(:),DPDN(:)
      ALLOCATABLE winf(:,:),CHDKDT(:),WD(:,:)
      ALLOCATABLE WKD(:,:)
! YE TIAN 04/11/2012 add following 2 lines
      allocatable gamonb_s(:,:) ! dimension(1:nc,1:mr) spanwise vorticity
      allocatable gamonb_c(:,:) ! dimension(1:nc,1:mr+1) chordwise vorticity
C/s S.N.KIM | includes HUB effect | 15-08-2016
      allocatable gamonh_s(:,:)
      allocatable gamonh_c(:,:)
C/e S.N.KIM | includes HUB effect | 15-08-2016
      allocatable gamond_s(:,:)
      allocatable gamond_c(:,:)
      allocatable gamonwd_s(:,:)
      allocatable gamonwd_c(:,:)
C/s S.N.KIM | ductwake alignment using FWA | Aug. 2018
      allocatable gamonw_s(:,:)
      allocatable gamonw_c(:,:)
C/e S.N.KIM | ductwake alignment using FWA | Aug. 2018
      END MODULE SIMQD1

      MODULE SIMQD2
      ALLOCATABLE AA(:,:),BB(:,:),C(:,:)
      ALLOCATABLE D(:,:),E(:,:),F(:,:)
      ALLOCATABLE WK(:,:),WK2(:,:)
      ALLOCATABLE WDK1(:),CHDK1(:),CHDKDT1(:)
      END MODULE SIMQD2

      MODULE SIMQ1
      ALLOCATABLE NPERB(:)
      END MODULE SIMQ1

      MODULE CORR
      ALLOCATABLE DELCP(:),DCPPRE(:)
      END MODULE CORR

      MODULE CIRC
! YE TIAN 04/08/2012
      ALLOCATABLE DELPonW(:)
C/s S.N.KIM | unsteady wake strength | Aug. 2019
      ALLOCATABLE DELPonW2(:)
C/e S.N.KIM | unsteady wake strength | Aug. 2019
C/s S.N.KIM | ductwake alignment using FWA | Aug. 2018
      ALLOCATABLE DELPonDW(:)
C/e S.N.KIM | ductwake alignment using FWA | Aug. 2018
      ALLOCATABLE DELP(:),DELPD(:)
      END MODULE CIRC

      MODULE PRES
      ALLOCATABLE CPB(:,:),CPH(:,:),DPDUB(:,:)
      ALLOCATABLE DPDVB(:,:),DPDUH(:,:),DPDVH(:,:)
      END MODULE PRES

      MODULE NORMAL
      ALLOCATABLE XCON(:,:),YCON(:,:),ZCON(:,:)
      END MODULE NORMAL

      MODULE FORCE
      ALLOCATABLE BKFX(:),BKFY(:),BKFZ(:),BKMX(:),BKMY(:)
      ALLOCATABLE BKMZ(:),HKFX(:),HKFY(:),HKFZ(:)
      ALLOCATABLE HKMX(:),HKMY(:),HKMZ(:)
      END MODULE FORCE

      MODULE VFOR
      ALLOCATABLE FVA(:)
      REAL XCDF
C/s S. Kim | turbine, local friction coefficients
      ALLOCATABLE VRR(:,:),RER(:,:)
      ALLOCATABLE XCDF1(:,:)
      REAL REI
C/e S. Kim | turbine, local friction coefficients
      END MODULE VFOR

      MODULE CAVGEOD
      ALLOCATABLE DS(:,:),DZ(:,:)
      ALLOCATABLE SZ(:,:),SPZ(:,:)
      ALLOCATABLE RLAM(:,:),CAVL(:,:),CAVLP(:,:),DZW(:,:)
      ALLOCATABLE SOP(:,:),CAVLSA(:,:,:)
      END MODULE CAVGEOD

      MODULE CAVGEOI
      ALLOCATABLE M0(:,:)
      INTEGER ITER
      ALLOCATABLE JCV(:,:),LCV(:,:)
      ALLOCATABLE NOCAV(:,:,:)
      INTEGER IT2
      END MODULE CAVGEOI

      MODULE CAVST
      ALLOCATABLE QC(:,:,:),PHI0(:,:)
      ALLOCATABLE PHI1(:,:,:),DPHIDS(:,:,:)
      ALLOCATABLE DPHIDV(:,:,:)
      END MODULE CAVST

      MODULE CVSOL2
      ALLOCATABLE HT(:,:,:),DPDNC(:),DELTA(:,:)
      ALLOCATABLE HTW1(:),DELTAP(:,:),SOL(:),SORW(:)
      ALLOCATABLE POTW(:),HTP(:,:,:),HTWP(:,:)
      END MODULE CVSOL2

      MODULE SPLITD
      ALLOCATABLE FLP(:,:),FRP(:,:),FLS(:,:),FRS(:,:)
      ALLOCATABLE DZL(:,:),DZR(:,:),QSPR(:,:),DLISP(:,:)
      ALLOCATABLE DPDVSP(:,:),PHIL(:,:),PHIR(:,:),QSPL(:,:)
      ALLOCATABLE QSSR(:)
      END MODULE SPLITD

      MODULE SPLITI
      ALLOCATABLE NSPP(:,:),NSPS(:,:)
      END MODULE SPLITI

      MODULE NWAK
      ALLOCATABLE NWC(:,:),NWPAN(:)
      END MODULE NWAK

      MODULE UVECB
      ALLOCATABLE UL(:,:),VL(:,:),WL(:,:)
      END MODULE UVECB

      MODULE UVECBo
      ALLOCATABLE ULo(:,:,:),VLo(:,:,:)
      ALLOCATABLE WLo(:,:,:)
      END MODULE UVECBo

      MODULE UVECW
      ALLOCATABLE ULW(:,:)
      END MODULE UVECW

      MODULE UVECWo
      ALLOCATABLE ULWo(:,:,:)
      END MODULE UVECWo

      MODULE XTRAP
      ALLOCATABLE CT(:,:,:),DT(:,:,:),QT(:,:,:),QW(:,:,:)
      END MODULE XTRAP

      MODULE SPINT
      ALLOCATABLE ICB(:,:,:),MSW(:,:),ICW(:,:)
      ALLOCATABLE IC(:,:,:),IW(:,:,:)
      ALLOCATABLE ISUBM(:,:),ISUBH(:,:)
      ALLOCATABLE NPERM(:),ISUB1(:,:,:),ISUWM(:,:)
      END MODULE SPINT

      MODULE SPIMG
      ALLOCATABLE XCTP_IM(:,:,:),XCPW_IM(:,:,:)
      END MODULE SPIMG

      MODULE FXSR
      ALLOCATABLE SSTE(:),VELTE(:,:)
      END MODULE FXSR

      MODULE CAVSR
      ALLOCATABLE CT2(:,:,:),QCSR(:,:,:),PHI0T(:,:)
      ALLOCATABLE PSI0T(:,:),PHI2(:,:,:)
      ALLOCATABLE DPHIDS2(:,:,:),DPHIDV2(:,:,:)
      ALLOCATABLE DELTAT1(:),DELTATP(:),CAVLT(:)
      ALLOCATABLE CAVLTP(:),CAVLTSA(:,:)
      END MODULE CAVSR

      MODULE TNGEO
      ALLOCATABLE XTUN(:,:),YTUN(:,:)
      ALLOCATABLE ZTUN(:,:)
      END MODULE TNGEO

      MODULE FVELTN1
      ALLOCATABLE VXITN(:,:),VETATN(:,:)
      ALLOCATABLE VINFSTN(:,:)
      END MODULE FVELTN1

      MODULE FVELTN2
      ALLOCATABLE DPDUTN(:,:),DPDVTN(:,:)
      ALLOCATABLE UXTNTOT(:,:),UYTNTOT(:,:)
      ALLOCATABLE UZTNTOT(:,:),CPTN(:,:)
      END MODULE FVELTN2

      MODULE DUCTGEOM2
      ALLOCATABLE XXCP(:),RRCP(:)
      END MODULE DUCTGEOM2

      MODULE DUCTGEOM3
      ALLOCATABLE XID(:),ETAD(:)
      END MODULE DUCTGEOM3

      MODULE DUCTGEOM
      ALLOCATABLE XD(:,:), YD(:,:)
      ALLOCATABLE ZD(:,:),XDW(:,:)
      ALLOCATABLE YDW(:,:),ZDW(:,:)
      ALLOCATABLE XDWDK(:,:),YDWDK(:,:)
      ALLOCATABLE ZDWDK(:,:),WDKD(:)
      END MODULE DUCTGEOM

      MODULE FVELD1
      ALLOCATABLE VXID(:,:),VETAD(:,:)
      ALLOCATABLE VINFSD(:,:)
      END MODULE FVELD1

      MODULE FVELD2
      ALLOCATABLE DPDUD(:,:),DPDVD(:,:)
      ALLOCATABLE UXDTOT(:,:),UYDTOT(:,:)
      ALLOCATABLE UZDTOT(:,:),CPD(:,:)
      END MODULE FVELD2

      MODULE FVELD3
      ALLOCATABLE VVYD(:),VVZD(:)
      END MODULE FVELD3

      MODULE FVELD4
      ALLOCATABLE DPOTMEAN(:),DCPMEAN(:),DPDUMEAN(:),DPDVMEAN(:)
      END MODULE FVELD4

      MODULE UNSWAKD
      ALLOCATABLE WSTINF(:,:),WSTINFD(:,:)
      END MODULE UNSWAKD

      MODULE UNSDP
      ALLOCATABLE DPHI(:,:)
      END MODULE UNSDP

      MODULE UNS
      ALLOCATABLE POTM(:,:),DPDTPRE(:,:)
      ALLOCATABLE POTWM(:,:),DPDTPREW(:,:)
      END MODULE UNS

      MODULE TMPP1
      ALLOCATABLE TEMP1(:),TEMP2(:),TEMP3(:)
      ALLOCATABLE TEMP4(:),TEMP5(:)
      ALLOCATABLE STRGTH(:,:),WUSINF(:,:)
      ALLOCATABLE STRGTH1(:,:),WUSINF1(:,:)
      ALLOCATABLE TEMPD4(:),STRGTHD(:,:)
      END MODULE TMPP1

      MODULE CAVFB
      ALLOCATABLE NNWC(:),NWDIR(:)
      END MODULE CAVFB

      MODULE WAKSUBD
      ALLOCATABLE XWS(:,:),YWS(:,:)
      ALLOCATABLE ZWS(:,:),WSUBIF(:,:)
      ALLOCATABLE WKFACE(:,:,:),WSUBIF_IM(:,:)
C/s S.N.KIM | Aug. 2018.
      ALLOCATABLE XWSO(:,:,:),YWSO(:,:,:),ZWSO(:,:,:)
C/e S.N.KIM
      END MODULE WAKSUBD

      MODULE DWAKSUB
      ALLOCATABLE XWSD(:,:),YWSD(:,:)
      ALLOCATABLE ZWSD(:,:),WSUBIFD(:,:)
      END MODULE DWAKSUB

      MODULE MEMSOL
      ALLOCATABLE ALHS(:,:)
      END MODULE MEMSOL

      MODULE CVRHS
      ALLOCATABLE RHS(:)
      END MODULE CVRHS

      MODULE BKUTTA
      ALLOCATABLE POTEMP(:),WWK(:,:)
      END MODULE BKUTTA

      MODULE TOGBFLOW
      ALLOCATABLE XCCC(:,:),YCCC(:,:),ZCCC(:,:),CPBMEAN(:,:,:)
      END MODULE TOGBFLOW

      MODULE PODGEO1
      ALLOCATABLE XHIN(:),YHIN(:),XYHINCUB(:)
      END MODULE PODGEO1

      MODULE PODGEO3
      ALLOCATABLE HRZ(:,:),HRZSQ(:,:)
      END MODULE PODGEO3

      MODULE PODGEO4
      ALLOCATABLE HRZP(:,:),HRZPSQ(:,:)
      END MODULE PODGEO4

      MODULE PRESNO
      ALLOCATABLE CPBN(:,:),CPHN(:,:)
      ALLOCATABLE CPTNN(:,:),CPDN(:,:)
      ALLOCATABLE CPTHN(:,:),CPCN(:,:)
      END MODULE PRESNO

      MODULE GLOCBL
      ALLOCATABLE XIV(:,:),ETAV(:,:)
      ALLOCATABLE XIVW(:,:),ETAVW(:,:)
      ALLOCATABLE XIVWS(:,:),ETAVWS(:,:)
      ALLOCATABLE XIVIS(:,:),ETAVIS(:,:)
C     FOR DUCT CASE
      ALLOCATABLE XIVD(:,:),ETAVD(:,:)
      ALLOCATABLE XIVDW(:,:),ETAVDW(:,:)
      ALLOCATABLE XIVDWS(:,:),ETAVDWS(:,:)
      ALLOCATABLE XIVISD(:,:),ETAVISD(:,:)
      END MODULE GLOCBL

      MODULE FVEL3
      ALLOCATABLE PPOTW(:),PPOTWS(:)
      ALLOCATABLE VTOTS_V(:)
      ALLOCATABLE UXTOT_V(:,:),UYTOT_V(:,:),UZTOT_V(:,:)
C     FOR DUCT CASE
      ALLOCATABLE PPOTDW(:),PPOTDWS(:)
      ALLOCATABLE UXDTOT_V(:,:),UYDTOT_V(:,:)
      ALLOCATABLE UZDTOT_V(:,:)
      END MODULE FVEL3

      MODULE FVELW
      ALLOCATABLE VXIW(:,:),VETAW(:,:),VINFSW(:,:)
      ALLOCATABLE VXIWS(:,:),VETAWS(:,:)
      ALLOCATABLE VINFSWS(:,:)
      ALLOCATABLE CPBL(:,:)
C     FOR DUCT CASE
      ALLOCATABLE VXIDW(:,:),VETADW(:,:)
      ALLOCATABLE VXIDWS(:,:),VETADWS(:,:)
      ALLOCATABLE VINFSDWS(:,:),VINFSDW(:,:)
      ALLOCATABLE CPBLD(:,:)
      END MODULE FVELW

      MODULE IVCON
      ALLOCATABLE LVFLAG(:)
C     FOR DUCT CASE
      ALLOCATABLE LVFLAGD(:)
      END MODULE IVCON

C     ALL THE MODULES LISTED BELOW ARE PART OF
C     PUFBL.INC
      MODULE INTGV
      ALLOCATABLE NTW(:)
      END MODULE INTGV

      MODULE PUFBL
      ALLOCATABLE AV(:,:),AVINV(:,:)
      ALLOCATABLE BETA(:,:),CSIGMA(:,:)
      ALLOCATABLE AWV(:,:),CA(:,:)
C XM YU 12/2011
c add for influence coefficients and mass defect of viscous run
      allocatable AVL(:,:,:),BV(:,:,:)
      allocatable csig(:,:,:),awv3(:,:,:)
      allocatable bsbb(:,:),bswb(:,:)
      allocatable duvl(:,:),umass(:,:)
      allocatable bsbw(:,:),bsww(:,:)
      allocatable updatem(:,:),savemass(:,:)
C XM YU 12/2011
      ALLOCATABLE CB(:,:),DH1(:,:)
      ALLOCATABLE DH(:,:),DD(:,:)
      ALLOCATABLE REC(:),CS(:,:)
      ALLOCATABLE XTRANLP(:),XTRANLS(:)
      END MODULE PUFBL

      MODULE VISBL
      ALLOCATABLE UVL(:),SVL(:),XNVIS(:)
      ALLOCATABLE YNVIS(:),XCVIS(:),YCVIS(:)
      ALLOCATABLE DVL(:,:)
      END MODULE VISBL

      MODULE REBL
      ALLOCATABLE BSRCB(:,:),BSRCWS(:,:),BSRCW(:,:)
      END MODULE REBL

      MODULE ARCLV
      ALLOCATABLE DELTZ(:,:),DELTS(:,:)
      ALLOCATABLE SC(:,:),SP(:,:),CDV(:)
c XM YU 03/2012
c add for local skin friction coefficients
      allocatable cfskin(:,:)
c XM YU 03/2012

      END MODULE ARCLV

C     ALL THE MODULES LISTED BELOW ARE PART OF
C     PUFBL.INC : DUCT CASE
      MODULE INTGVD
      ALLOCATABLE NTWD(:)
      END MODULE INTGVD

      MODULE PUFBLD
      ALLOCATABLE AVD(:,:),AVDINV(:,:)
      ALLOCATABLE BETAD(:,:),CSIGMAD(:,:)
      ALLOCATABLE AWVD(:,:),CAD(:,:)
      ALLOCATABLE CBD(:,:),DH1D(:,:)
      ALLOCATABLE DHD(:,:),DD_D(:,:)
      ALLOCATABLE RECD(:),CSD(:,:)
      ALLOCATABLE XTRANDP(:),XTRANDS(:)
      END MODULE PUFBLD

      MODULE VISBLD
      ALLOCATABLE UVLD(:),SVLD(:),XNVISD(:)
      ALLOCATABLE YNVISD(:),XCVISD(:),YCVISD(:)
      ALLOCATABLE DVLD(:,:)
      END MODULE VISBLD

      MODULE REBLD
      ALLOCATABLE BSRCBD(:,:),BSRCWSD(:,:),BSRCWD(:,:)
      END MODULE REBLD

      MODULE ARCLVD
      ALLOCATABLE DELTZD(:,:),DELTSD(:,:)
      ALLOCATABLE SCD(:,:),SPD(:,:),CDVD(:)
      END MODULE ARCLVD

      MODULE PC2NSVEL
      LOGICAL IPC2NS,IWRT2
      REAL,ALLOCATABLE:: UEFX(:,:,:),UEFR(:,:,:),UEFT(:,:,:)
      END MODULE PC2NSVEL

      MODULE UNSPROP
      INTEGER IUNS
      INTEGER NTSTEP_KEY
      INTEGER IDXREV_KEY
      END MODULE
