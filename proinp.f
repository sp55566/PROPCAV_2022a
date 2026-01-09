      SUBROUTINE PROINP
************************************************************************
*     PROINP: PROpeller geometry INPut                                 *
*     --- Read input data from two input files                        *
*     --- Calculate SPLINE cubic coefficients for geometries          *
*     *
*     Date of last revision      Revision                                 *
*     ---------------------      --------                                 *
*     10/18/2001                 Undergone major clean up for v2.0        *
*     *
************************************************************************

      use m_INPGEO2
      use m_BLADEM
      use m_NACAIO
      use m_WKNP
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'INPGEO.INC'
      COMMON /MODGEO1/ DUMXC(MXINP),PC(16),
     % YP(MXINP,NXMAX),YS(MXINP,NXMAX)
      COMMON /HUB5TMP/RHULT1
      COMMON/HIMAGE2/RRHUB

      COMMON /IFORT2/ NPOINTS

      DIMENSION RLVLRX(NXMAX),XR1(NXMAX),VLRCUB(MAXSPL)
      DIMENSION CUBIC(4*(MXINP-1)),DUMCAM(MXINP),DUMTHK(MXINP)
      DIMENSION YC1(MXINP),TC1(MXINP),YTMP(16),XRSQ1(1)
      DATA PC/0.01,0.025,0.05,0.1,0.2,0.3,
     * 0.4,0.5,0.6,0.7,0.8,0.9,0.95,
     * 0.975,0.99,1.00/
      CHARACTER*29 FNADM,FNGEO,FNOUT
      LOGICAL IDSTP

C-----------------------------------------------------------------------
C     IPHUB=1   Print hub pressures
C     IPHUB=0   Do not print hub pressures
C-----------------------------------------------------------------------
      IPHUB=1

C-----------------------------------------------------------------------
C     Read a universal filename and open input and general output files
C-----------------------------------------------------------------------
C     write(*,'(a)') fn,fnscr
      WRITE(*,'(A,$)') ' PROPCAV> ENTER UNIVERSAL FILENAME: '
      READ(*,'(A)') FN
      CALL CHRLEN(FN,LENCH)
      FNADM=FN(1:LENCH)//'.adm'
      FNGEO=FN(1:LENCH)//'.geo'
      FNOUT=FN(1:LENCH)//'.out'
      OPEN(1,FILE=FNADM,STATUS='OLD',FORM='FORMATTED')
      OPEN(3,FILE=FNGEO,STATUS='OLD',FORM='FORMATTED')
      OPEN(2,FILE=FNOUT,STATUS='UNKNOWN',FORM='FORMATTED')

C-----------------------------------------------------------------------
C     Read run-control options
C-----------------------------------------------------------------------
C     -- IAN : 0 = no align , 1 = align, 2= unsteady align
C     -- ISP : 0 = fully submerged, 1 = partially submerged
C     -- IMG : 0 = no image, 1 = image (used only if ISP=1)

C     write(*,*)
C     write(*,*) ' READ FROM ADM FILE '
C     write(*,*)

!     READ(1,*) ICON,wingimag,IHUB,ISTEADY,IAN,IFILE,ISP,IMG
      READ(1,*) ICON,IHUB,ISTEADY,IAN,IFILE,ISP,IMG
      READ(1,*) ITERKT,TOLKT,IPK,IDRT,INRPAN
      READ(1,*) ICSPAC,IRSPAC,ITHK,ICAM,ISEARCH,IFACE
      READ(1,*) IVPITCH,TFACTOR

C/s S.N.KIM | modifies IAN to incorporate duct wake alignment with FWA | Aug. 2018
      IF (IAN.EQ.7) THEN
        IAN = 6
        IDFWA = 1
      ENDIF 
C/e S.N.KIM | modifies IAN to incorporate duct wake alignment with FWA | Aug. 2018

C Modify parameters for unsteady case -- Yiran Su 09/19/2017 --
      IF (IUNS.EQ.1) THEN
        ICON = 0
        ISTEADY = 2
        IAN = 1
        ISP = 0
        IMG = 0
        INRPAN = 0
        ISC = 0
        IAN = 1
        IF (ITUN .NE.0) THEN
          WRITE(*,*) 'Tunnel not supported in unsteady mode.'
          STOP
        END IF
        IF (IDUCT.NE.0) THEN
          WRITE(*,*) 'Duct not supported in unsteady mode.'
          STOP
        END IF
      END IF


C---- Yiran panel input method 10/05/2016
      IF ((IDRT.EQ.1).AND.((IAN.EQ.0).OR.(IAN.EQ.1)))  THEN
        WRITE(*,*) "Warning: IAN=5 or 6 is recommended for IDRT=1!"
      END IF

C---- Yiran input wake geometry option 10/05/2016
      INQUIRE(FILE='wakinp.dat',EXIST=IDSTP)
      IF (IAN.EQ.5) THEN
        IAN=1
        IREADW=1 ! S.N.KIM | for cavity runs in case wakes are read in | Aug. 2018.
        WRITE(*,*) 'IAN=5: use wakinp.dat to determine wake geometry!'
        IF (IDSTP.EQ.(.FALSE.)) THEN
          WRITE(*,*) 'File wakinp.dat does not exist! PROPCAV stopped!!'
        ENDIF
      ELSEIF (IDSTP.EQ.(.TRUE.)) THEN
        OPEN(UNIT=721,FILE='wakinp.dat',STATUS='UNKNOWN')
        CLOSE(UNIT=721,STATUS='DELETE')
      ENDIF

!s--- YE TIAN ---  09/02/2013 --
      if (IVISC .eq. 1) then
        IVPITCH = 0
      end if
!s--- YE TIAN ---  09/02/2013 --
      wingimag = 0
C     sh
C
C     temporary change in angle of attack
C
C     sh
C-----------------------------------------------------------------------
C     The next three options are only read if ITHK=99.          JY111600
C
C     IFORMAT: 0 = input thickness distribution is in T/C (PROPCAV)
C     format.
C     1 = input thickness distribution is in T/Tmax (MPUF3A)
C     format.
C     2 = input chord-wise locations x/C and pressure side
C     y_p/C + suction side y_s/C coordinates instead of
C     the camber and thickness distributions.
C
C     IMOD: 0 = use input camber and thickness as given in *.geo.
C     1 = modify suction side geometry so that the T.E. thickness
C     is zero.
C-----------------------------------------------------------------------
      READ(1,*) IFORMAT,IMOD
!s--- YE TIAN --- 07/04/2013----
!     New parameter for T.E. closure
!     x_trunc: the location of truncation after which the section is
!     modified in order to close the T.E.
      if (IMOD .eq. 2) read(1,*) x_trunc
!e--- YE TIAN --- 07/04/2013----



C-----------------------------------------------------------------------
C     New parameter for supercavitating sections.               JY060901
C
C     ISC: 0 = not supercavitating section
C     1 = supercavitating section (Note: IMOD should be 0)
C-----------------------------------------------------------------------
      READ(1,*) ISC
      READ(1,*) IFLAP,IHULL,ITUN, IDUCT
      READ(1,*) ITERGB

C/s S.N.KIM | duct wake alignment with FWA | Aug. 2018.
      IF (IDFWA.EQ.1) THEN
        IF (IDUCT.EQ.0) THEN
          WRITE(*,*) 'Duct should be included (IDUCT=1) if 
     &duct wake alignment is incldued (IDFWA=1).'
          STOP
        ENDIF
      ENDIF
C/e S.N.KIM | duct wake alignment with FWA | Aug. 2018.

C-----------------------------------------------------------------------
C     Read the head of the input file for geometries
C-----------------------------------------------------------------------
C     write(*,*) ' READ FROM GEO FILE '
C     write(*,*)

      READ(3,'(A)') LABEL
      READ(3,*) NX,NBLADE,NC,MR
      READ(3,*) RHUB, NHBU,MHBT, XHBU,XHBD,XHBT
      READ(3,*) ADVCO, RULT,RHULT,DCD,XULT,XUWDK
      READ(3,*) FROUDE,SIGMA,DTOL,RLAMDA1,ITERMAX

C-----Finished building the blade geometry
C-----If IDRT=1, use node points location on panel surface
C-----Yiran panel input
      IF (IDRT.EQ.1) THEN
         OPEN(197,FILE='panel_inp.dat',STATUS='OLD')
         READ(197,*) NC,MR
         NC=NC-1
         MR=MR-1
         CLOSE(197)
      ENDIF
C-----YIRAN FINISHED

C.... If it's a SC propeller, then the following variables need to be
C.... set.
      IF(ISC.EQ.1) THEN
       ITERKT=0
       IF(ISP.NE.1)ISEARCH=1
       IFACE=2
       IMOD=0
       WRITE(*,*) '--------------------------------------------------'
       WRITE(*,*) 'For ISC=1 (supercavitating propeller sections)'
       WRITE(*,*) '--------------------------------------------------'
       WRITE(*,*) ' The following input variables are used:'
       WRITE(*,*) ' ITERKT=0,ISEARCH=1,IFACE=2,IMOD=0'
       WRITE(*,*)
       IF(IAN.EQ.2) THEN
        WRITE(*,*) ' IAN=2 is not allowed.  Setting IAN=1.'
        IAN=1
       END IF
      END IF

      READ(1,*) XCDF, DAMP
      READ(1,'(A)') FNSCR

C---- New variables for special fan-like SP propeller (JY090100)
      IF(ICON.EQ.7) READ(1,*) IANG,INH

C---- New variables for special foil with different cavitation number
C---- on the face and on the back (JY110100)
      IF(ICON.EQ.8) READ(1,*) SIGMAB,SIGMAF,XLTRAN

      IF(IHUB .EQ. 6 .OR. IHUB .EQ. 7) THEN
C     HPCH1,HPCH2 : Initial Pitch angle of hub grid at LE & TE (45deg. recommended)
C     If Pitch is too low => Increase angles
C
       READ(1,*) HPCH1,HPCH2
       HPCH1 = HPCH1 * PI / 180.
       HPCH2 = HPCH2 * PI / 180.
      ENDIF

!s--- YE TIAN ----
      if (IAN.eq.6) then
        READ(1,*) ALIGN_DIV, WKPAN_INP, NALIGN_ITER, IPLOT
      end if
!e--- YE TIAN ----
      CLOSE (1)

!-----------------------------------------------------------------------
!     Read the input file for geometries
!-----------------------------------------------------------------------
!     write(*,*) ' READ FROM GEO FILE '
!     write(*,*)

!     READ(3,'(A)') LABEL
!     READ(3,*) NX,NBLADE,NC,MR
!     READ(3,*) RHUB, NHBU,MHBT, XHBU,XHBD,XHBT
!     READ(3,*) ADVCO, RULT,RHULT,DCD,XULT,XUWDK
!     READ(3,*) FROUDE,SIGMA,DTOL,RLAMDA1,ITERMAX

!s-- YE TIAN 07/05/2013---- assign values of the parameters
      KZ = NBLADE
      if (KZ .le. 3) KZ = 3
      NBZ= NC+NZSR
      MBZ= MR
      if (IAN .eq. 2) MBZ = MBZ + 1
      NBPZ=NBZ+1
      MBPZ=MBZ+1
      NBHZ=NBZ/2
      NBHZP=NBHZ+1
!e-- YE TIAN 07/05/2013---- assign values of the parameters

      NWZ=360
      NWPZ = NWZ + 1

      NSCZ=NBHZ+NWZ
      NSCZP=NSCZ+1
      NBSCZ=NBZ+NWZ
      NBSCZP=NBSCZ+1

      NHUZ=60;
      if (IAN.eq.2) then
        NHDZ = NWZ
      else
        NHDZ=200;
      end if
      NHBZ=NHUZ+NHDZ+NBHZ;NHPZ=NHBZ+1
      MHBZ=40;MHPZ=MHBZ+1

      mrmax = mbz; nwmax=nwz
      ncavm = nwz; ncavmp = ncavm+1
      mcavm = 10 ; mcavmp=mcavm + 1
      nhmx=nbz/2; nhmxp = nhmx + 1
      mrmaxp = mrmax + 1; nwmaxp=nwmax + 1

      NPHZ=NHBZ*MHBZ
      nphc=mcavm*ncavm

      if (ITUN .eq. 0) then
        NtMAXt=3; NCMAXT=3
      else
        NtMAXt=200; NCMAXT=20
      end if
      NPANMT=NTMAXT*NCMAXT

      if (IDUCT .eq. 0) then
        NDMAX=3;MDMAX=3;NDWMX=3
      else
        NDMAX=200; MDMAX=80;NDWMX=1000
      end if
      NDMAXP= NDMAX+1; MDMAXP = MDMAX + 1; NDWMXP=NDWMX+1
      NPANMD=NDMAX*MDMAX; NDWMN=NDWMX*MDMAX



      NPHT = MCAVM*NHMX
      NPAWZ=MBZ*(NWZ+NDEL)
      NPAWZ2=NPAWZ
      NSCWZ=NPAWZ

      NPBZ=NBZ*MBZ
      NPANZ=NPBZ+NPHZ+NPANMT+NPHC+NPHT+NPANMD
      NTZ=NPANZ+NPAWZ2
      NBLKMAX=MBZ+NHBZ+NTMAXT+NCAVM+NHMX+MDMAX








      NTZ=NPANZ+NPAWZ2


      NWZVIS=2*NWPZ+NDEL;NZVIS=NBPZ+NWZVIS
      NDWZVIS=NDWMX+NWZ+NDEL;NZVISD=NDMAXP+NDWZVIS


      NCSRZ=NZSR2+NWZ;NCSRZP=NCSRZ+1



!e-- YE TIAN 07/05/2013---- assign values of the parameters

!     module m_INPGEO2
      allocate(CUBTMP(MAXSPL),
     *		XCS1(15,MBPZ),XTS1(16,MBPZ),XCS0(MBPZ),XTS0(MBPZ))

!     module m_BLADEM
      allocate(XBM(NBHZP,MBPZ),YBM(NBHZP,MBPZ),ZBM(NBHZP,MBPZ))

!     module m_NACAIO
      allocate(XKEY(NBPZ),YKEY(NBPZ),ZKEY(NBPZ)
     &     ,XNEAR(NBPZ),YNEAR(NBPZ),ZNEAR(NBPZ),CROSS(3,NBPZ))

!     module m_WKNP
      allocate(NWIDX(MBZ),NWSEC(MBZ))

      IF(IVISC.EQ.1) THEN

C     MODULE GLOCBL
        ALLOCATE( XIVW(NWPZ,MBPZ),ETAVW(NWPZ,MBPZ) )
        ALLOCATE( XIVWS(NWZ+NDEL,MBPZ),ETAVWS(NWZ+NDEL,MBPZ) )
        ALLOCATE( XIVIS(NZVIS,MBZ),ETAVIS(NZVIS,MBZ) )

C     MODULE FVEL3
        ALLOCATE (PPOTW(NSCWZ),PPOTWS(NSCWZ))
        ALLOCATE (VTOTS_V(NPANZ))
        ALLOCATE (UXTOT_V(NBZ,MBZ),UYTOT_V(NBZ,MBZ),UZTOT_V(NBZ,MBZ))

C     MODULE FVELW
        ALLOCATE (VXIW(NWZ,MBZ),VETAW(NWZ,MBZ),VINFSW(NWZ,MBZ))
        ALLOCATE (VXIWS(NWZ+NDEL,MBZ),VETAWS(NWZ+NDEL,MBZ))
        ALLOCATE (VINFSWS(NWZ+NDEL,MBZ))
        ALLOCATE (CPBL(NBZ,MBZ))

C     MODULE IVCON
        ALLOCATE (LVFLAG(MBPZ))

C     MODULE INTGV
        ALLOCATE (NTW(MBZ))

C     MODULE PUFBL
        ALLOCATE (AV(NZVIS,NZVIS),AVINV(NZVIS,NZVIS))
        ALLOCATE (BETA(NZVIS,NZVIS),CSIGMA(NZVIS,NZVIS))
        ALLOCATE (AWV(NZVIS,NZVIS),CA(NZVIS,NZVIS))
c XM YU 12/2011
        allocate (AVL(nzvis,nzvis,nzvis),BV(nzvis,nzvis,nzvis))
        allocate (CSIG(nzvis,nzvis,nzvis),awv3(nzvis,nzvis,nzvis))
        allocate (bsbb(npanz,npanz),bswb(nscwz,npanz))
        allocate (bsbw(npanz,nscwz),bsww(nscwz,nscwz))
        allocate (duvl(nzvis,mbpz),umass(nzvis,mbpz))
        allocate (updatem(nzvis,mbpz),savemass(nzvis,mbpz))
c XM YU 12/2011
        ALLOCATE (CB(NZVIS,NZVIS),DH1(NZVIS,NZVIS))
        ALLOCATE (DH(NZVIS,NZVIS),DD(NZVIS,NZVIS))
        ALLOCATE (REC(MBZ),CS(NZVIS,NZVIS))
        ALLOCATE (XTRANLP(MBZ),XTRANLS(MBZ))

C     MODULE VISBL
        ALLOCATE (UVL(NZVIS),SVL(NZVIS),XNVIS(NZVIS))
        ALLOCATE (YNVIS(NZVIS),XCVIS(NZVIS),YCVIS(NZVIS))
        ALLOCATE (DVL(NZVIS,NZVIS))

C     MODULE REBL
        ALLOCATE (BSRCB(NPANZ,KZ),BSRCWS(NSCWZ,KZ),BSRCW(NSCWZ,KZ))

C     MODULE ARCLV
        ALLOCATE (DELTZ(NZVIS,MBZ),DELTS(NZVIS,MBZ))
        ALLOCATE (SC(NZVIS,MBZ),SP(NZVIS,MBZ),CDV(MBZ))
c XM YU 03/2012
        allocate (cfskin(nzvis,mbz))
c XM YU 03/2012
      ENDIF

      IF(IHUB.EQ.6) THEN
        ALLOCATE(XHIN(200),YHIN(200),XYHINCUB(800))
      END IF

      IF(IDUCT.EQ.1) THEN

C     MODULE DUCTGEOM2,DUCTGEOM3,DUCTGEOM,FVELD1,FVELD2,FVELD3,FVELD4
        ALLOCATE(XXCP(NDMAX),RRCP(NDMAX))
        ALLOCATE(XID(201),ETAD(201))
        ALLOCATE(XD(NDMAXP,MDMAXP), YD(NDMAXP,MDMAXP))
        ALLOCATE(ZD(NDMAXP,MDMAXP),XDW(NDWMXP,MDMAXP))
        ALLOCATE(YDW(NDWMXP,MDMAXP),ZDW(NDWMXP,MDMAXP))
        ALLOCATE(XDWDK(MDMAXP,2),YDWDK(MDMAXP,2))
        ALLOCATE(ZDWDK(MDMAXP,2),WDKD(NPANZ))
        ALLOCATE(VXID(NDMAX,MDMAX),VETAD(NDMAX,MDMAX))
        ALLOCATE(VINFSD(NDMAX,MDMAX))
        ALLOCATE(DPDUD(NDMAX,MDMAX),DPDVD(NDMAX,MDMAX))
        ALLOCATE(UXDTOT(NDMAX,MDMAX),UYDTOT(NDMAX,MDMAX))
        ALLOCATE(UZDTOT(NDMAX,MDMAX),CPD(NDMAX,MDMAX))
        ALLOCATE(VVYD(NPANMD),VVZD(NPANMD))
        ALLOCATE(DPOTMEAN(NDMAX),DCPMEAN(NDMAX),DPDUMEAN(NDMAX),
     &          DPDVMEAN(NDMAX))

        IF(IVISC.EQ.1) THEN

C     MODULE GLOCBL : FOR DUCT CASE
          ALLOCATE( XIVD(NDMAXP,MDMAXP),ETAVD(NDMAXP,MDMAXP) )
          ALLOCATE( XIVDW(NDMAXP,MDMAXP),ETAVDW(NDMAXP,MDMAXP) )
          ALLOCATE( XIVDWS(NWZ+NDEL,MDMAXP),ETAVDWS(NWZ+NDEL,MDMAXP) )
          ALLOCATE( XIVISD(NZVISD,MDMAX),ETAVISD(NZVISD,MDMAX) )

C     MODULE FVEL3 : FOR DUCT CASE
          ALLOCATE (PPOTDW(NSCWZ),PPOTDWS(NSCWZ))
          ALLOCATE (UXDTOT_V(NDMAX,MDMAX),UYDTOT_V(NDMAX,MDMAX))
          ALLOCATE (UZDTOT_V(NDMAX,MDMAX))

C     MODULE FVELW : FOR DUCT CASE
          ALLOCATE (VXIDW(NWZ,MDMAX),VETADW(NWZ,MDMAX))
          ALLOCATE (VXIDWS(NWZ+NDEL,MDMAX),VETADWS(NWZ+NDEL,MDMAX))
          ALLOCATE (VINFSDWS(NWZ+NDEL,MDMAX),VINFSDW(NWZ,MDMAX))
          ALLOCATE (CPBLD(NDMAX,MDMAX))

C     MODULE IVCON : FOR DUCT CASE
          ALLOCATE (LVFLAGD(MDMAX))
C     MODULE INTGVD
          ALLOCATE (NTWD(MDMAX))

C     MODULE PUFBLD
          ALLOCATE(AVD(NZVISD,NZVISD),AVDINV(NZVISD,NZVISD))
          ALLOCATE(BETAD(NZVISD,NZVISD),CSIGMAD(NZVISD,NZVISD))
          ALLOCATE(AWVD(NZVISD,NZVISD),CAD(NZVISD,NZVISD))
          ALLOCATE(CBD(NZVISD,NZVISD),DH1D(NZVISD,NZVISD))
          ALLOCATE(DHD(NZVISD,NZVISD),DD_D(NZVISD,NZVISD))
          ALLOCATE(RECD(MDMAX),CSD(NZVISD,NZVISD))
          ALLOCATE(XTRANDP(MBZ),XTRANDS(MBZ))

C      MODULE VISBLD
          ALLOCATE(UVLD(NZVISD),SVLD(NZVISD),XNVISD(NZVISD))
          ALLOCATE(YNVISD(NZVISD),XCVISD(NZVISD),YCVISD(NZVISD))
          ALLOCATE(DVLD(NZVISD,NZVISD))

C      MODULE REBLD
          ALLOCATE(BSRCBD(NPANZ,KZ),BSRCWSD(NSCWZ,KZ))
          ALLOCATE(BSRCWD(NSCWZ,KZ))

C      MODULE ARCLVD
          ALLOCATE(DELTZD(NZVISD,MDMAX),DELTSD(NZVISD,MDMAX))
          ALLOCATE(SCD(NZVISD,MDMAX),SPD(NZVISD,MDMAX),CDVD(MDMAX))
        ENDIF
      ENDIF

!     module m_CAVTY1
      allocate(XVC(NCAVMP,MCAVMP),YVC(NCAVMP,MCAVMP),
     *               ZVC(NCAVMP,MCAVMP))
!     module m_CAVTY2
      allocate(XIND(NPAWZ),YIND(NPAWZ),ZIND(NPAWZ))

!     module m_CAVTY3
      allocate(VXMEAN(NCAVM),VYMEAN(NCAVM),VZMEAN(NCAVM)
     *              ,CMX(NPHC),CMY(NPHC),CMZ(NPHC))

!     module m_CAVTY4
      allocate(CPC(NCAVM,MCAVM),VTC(NCAVM,MCAVM))

!     module m_CAVTY5
      allocate(VXIC(NCAVM,MCAVM),VETAC(NCAVM,MCAVM),
     *               VINFSC(NCAVM,MCAVM))

!     module m_CAVTY6
      allocate(DPDUC(NCAVM,MCAVM),DPDVC(NCAVM,MCAVM))

!     module m_CAVTY7
      allocate(UXCTOT(NCAVM,MCAVM),UYCTOT(NCAVM,MCAVM),
     *               UZCTOT(NCAVM,MCAVM))

!     module  m_TIPH1
      allocate(XCH(NHMXP,MCAVMP),YCH(NHMXP,MCAVMP),
     *               ZCH(NHMXP,MCAVMP))

!     module  m_TIPH2
      allocate(DPDUTH(NHMX,MCAVM),DPDVTH(NHMX,MCAVM),
     *               CPTH(NHMX,MCAVM),VTTH(NHMX,MCAVM),
     *               VINFSTH(NHMX,MCAVM))

!     module  m_TIPH3
      allocate(VXITH(NHMX,MCAVM),VETATH(NHMX,MCAVM))

!     module  m_VELWAKE
      allocate(WINDX(NWMAX,MRMAXP),WINDY(NWMAX,MRMAXP),
     *                 WINDZ(NWMAX,MRMAXP))

!     module  m_CAVCENT
      allocate(XVCC(NCAVMP),YVCC(NCAVMP),ZVCC(NCAVMP))





      MRP = MR + 1

      RRHUB = RHUB

!s--- YE TIAN 08/13/2013-----
!     module m_DISKI
      allocate(NSW(MBPZ))

!     module m_DISKD
      allocate(XHDK(MHPZ,2),YHDK(MHPZ,2)
     *     ,  ZHDK(MHPZ,2),XWDK(MHPZ,2),YWDK(MHPZ,2),ZWDK(MHPZ,2)
     *     ,  XHDKDT(MHPZ,2),YHDKDT(MHPZ,2),ZHDKDT(MHPZ,2))

!     module m_GOUT
      allocate(SBP(NBHZ),RZP(MBZ),RZPSQ(MBZ))

!     MODULE m_gproi
      if (.NOT.allocated(NLEP)) then
        ALLOCATE(NLEP(MBZ,NSTEP,2))
      end if

!     MODULE SPLITI
      if (.NOT.allocated(NSPP)) then
        ALLOCATE(NSPP(MBZ,2),NSPS(MBZ,2))
      end if

!s--- YE TIAN remarks 08/13/2013----
C.... Since M841B is a special geometry, it is hard wired into the code.
      IF(ISC.EQ.1.AND.FN(1:LENCH).EQ.'m841b') RETURN

!     This is ugly...... there is a version with more general treatment
!     but it is not included in this release.
!e--- YE TIAN remarks 08/13/2013----


C     IF(ISC.EQ.1.AND.FN(1:LENCH).EQ.'m841b') THEN
C     CALL M841B
C     GO TO 2000
C     END IF

      READ(3,*) (XR(N),N=1,NX)
      READ(3,*) (XPI(N),N=1,NX)
      READ(3,*) (XRAKE(N),N=1,NX)
      READ(3,*) (XSKEW(N),N=1,NX)
      READ(3,*) (XCHD(N),N=1,NX)
      READ(3,*) (XCI(N),N=1,NX)
      READ(3,*) (XTI(N),N=1,NX)

c YIRAN panel input
      IF (IDRT.EQ.1) THEN
        DO N=1,NX
          XR(N)=RHUB+(REAL(N-1))/(REAL(NX-1))*(XR(NX)-RHUB)
        END DO
      END IF
c Yiran Finished

C---- Check paramters for dimensioning errors.  Results are printed to
C---- ERR.LOG.

      IF(IDUCT .NE. 0) CALL DUCTINPUT

C     -- Check PARAM.INC

      CALL CHECK_PARAM0

C     -- Check input parameters in geometry file

      CALL CHECK_PARAM


!s--- YE TIAN ---- 07/05/2013---- allocation



C     -------------------------------------

C     T--- 082600 Tempo. Change (same as MPUF3a)
C     T--- This should not be applied for ISP=1 (JY091400)
C.....Same for ICON=8 as well (JY110100)
      IF(ICON.NE.5.AND.ICON.NE.6.AND.ISP.NE.1.AND.ICON.NE.8) THEN
       BUG=FILLIN(0.95,XR,XCHD,NX)
       BUG1=0.1*BUG
       IF(XCHD(NX).LT.BUG1) XCHD(NX)=BUG1
      END IF
C     T----------------------------------------------

      DO N=1,NX
       XUA(N)=ZERO
       XUAU(N)=ZERO
       XUT(N)=ZERO
       XUTU(N)=ZERO
      END DO

      IF(IRSPAC.EQ.3) THEN
       READ(3,*) (RZ(M),M=1,MRP)
      END IF

C.....For the VLR section .input the LE radii
      IF(ITHK.EQ.7) THEN
       READ(3,*) (RLVLRX(N),N=1,NX)
      END IF

C-----------------------------------------------------------------------
C     For ICAM=99, read camber distribution from .geo file.    JY112798
C-----------------------------------------------------------------------
c      IF(ICAM.EQ.99.AND.IFORMAT.NE.2)
c     * READ(3,*) ((XCS(N,J),N=1,NX),J=1,15)
      IF(ICAM.EQ.99.AND.IFORMAT.NE.2) THEN
         READ(3,*) (XCS(N,1),N=1,NX)
         READ(3,*) (XCS(N,2),N=1,NX)
         READ(3,*) (XCS(N,3),N=1,NX)
         READ(3,*) (XCS(N,4),N=1,NX)
         READ(3,*) (XCS(N,5),N=1,NX)
         READ(3,*) (XCS(N,6),N=1,NX)
         READ(3,*) (XCS(N,7),N=1,NX)
         READ(3,*) (XCS(N,8),N=1,NX)
         READ(3,*) (XCS(N,9),N=1,NX)
         READ(3,*) (XCS(N,10),N=1,NX)
         READ(3,*) (XCS(N,11),N=1,NX)
         READ(3,*) (XCS(N,12),N=1,NX)
         READ(3,*) (XCS(N,13),N=1,NX)
         READ(3,*) (XCS(N,14),N=1,NX)
         READ(3,*) (XCS(N,15),N=1,NX)
      ENDIF

C-----------------------------------------------------------------------
C     For ITHK=99, read thickness distribution from .geo file. JY112798
C-----------------------------------------------------------------------
      IF(ITHK.EQ.99.AND.IFORMAT.NE.2) THEN
c        READ(3,*) ((XTS(N,J),N=1,NX),J=1,16)
        READ(3,*) (XTS(N,1),N=1,NX)
        READ(3,*) (XTS(N,2),N=1,NX)
        READ(3,*) (XTS(N,3),N=1,NX)
        READ(3,*) (XTS(N,4),N=1,NX)
        READ(3,*) (XTS(N,5),N=1,NX)
        READ(3,*) (XTS(N,6),N=1,NX)
        READ(3,*) (XTS(N,7),N=1,NX)
        READ(3,*) (XTS(N,8),N=1,NX)
        READ(3,*) (XTS(N,9),N=1,NX)
        READ(3,*) (XTS(N,10),N=1,NX)
        READ(3,*) (XTS(N,11),N=1,NX)
        READ(3,*) (XTS(N,12),N=1,NX)
        READ(3,*) (XTS(N,13),N=1,NX)
        READ(3,*) (XTS(N,14),N=1,NX)
        READ(3,*) (XTS(N,15),N=1,NX)
        READ(3,*) (XTS(N,16),N=1,NX)

!s--- YE TIAN --- 07/04/2013----
        IF (IMOD .eq. 2) then
          call close_te(nx, nxmax, xts, x_trunc)
          IMOD = 0
        end if
!e--- YE TIAN --- 07/04/2013----
      END IF

C-----------------------------------------------------------------------
C     The next if statement is used to convert the thickness
C     distribution input from MPUF3A (T/Tmax) to PROPCAV (T/C) format.
C     JY111600
C-----------------------------------------------------------------------
      IF(ITHK.EQ.99.AND.ICON.NE.5.AND.IFORMAT.EQ.1) THEN
       DO J=1,16
        DO N=1,NX
         IF(XCHD(N).EQ.0.0) THEN
          XTS(N,J)=XTS(N,J)*XTI(N)/0.1
         ELSE
          XTS(N,J)=XTS(N,J)*XTI(N)/XCHD(N)
         END IF
        END DO
       END DO
      END IF

C-----------------------------------------------------------------------
C     For IFORMAT=2, read input chord-wise locations and corresponding
C     pressure side and suction side coordinates.  Then convert them
C     to XCS and XTS if IMOD=0.                                 JY111600
C     THis routine is wrong!
C     However, I'll keep this routine for the moment to use for PSF2ALIGN.F
C     (Hope this doesnot affect much on calculation of induced velocity),
C     and I'll directly calculate blade geometry in <gblade.f>
C
C-----------------------------------------------------------------------
      IF(IFORMAT.EQ.2) THEN
        READ(3,*) NPOINTS
!s--- Allen Du 06/07/2018 add the following line to check the array size
        if(npoints.ge.mxinp) then
          write(*,*)"Increase MXINP in PARAM.INC file."
          write(*,*)"Will stop here."
          stop
        endif
!e--- Allen Du 06/07/2018 add the following line to check the array size
        READ(3,*) (DUMXC(N),N=1,NPOINTS)
        DO I=1,NX
         READ(3,*) (YP(N,I),N=1,NPOINTS)
         READ(3,*) (YS(N,I),N=1,NPOINTS)
        END DO

!s--- YE TIAN --- 07/04/2013----
        IF (IMOD .eq. 2) then
          call  close_te2(nx, NPOINTS, MXINP, DUMXC, YP, x_trunc)
          call  close_te2(nx, NPOINTS, MXINP, DUMXC, YS, x_trunc)
          IMOD = 0
        end if
!e--- YE TIAN --- 07/04/2013----
        IF(IMOD.EQ.0) THEN
          DO I=1,NX
            DO N=1,NPOINTS
             DUMCAM(N)=(YP(N,I)+YS(N,I))/2.
             DUMTHK(N)=YS(N,I)-YP(N,I)
            ENDDO

            CALL UGLYDK(NPOINTS,1,1,DUMXC,DUMCAM,ZERO,ZERO,CUBIC)
            CALL EVALDK(NPOINTS,15,DUMXC,PC,YTMP,CUBIC)
            DO N=1,15
             XCS(I,N)=YTMP(N)
            END DO

            CALL UGLYDK(NPOINTS,1,1,DUMXC,DUMTHK,ZERO,ZERO,CUBIC)
            CALL EVALDK(NPOINTS,16,DUMXC,PC,YTMP,CUBIC)
            DO N=1,16
             XTS(I,N)=YTMP(N)
            END DO
          END DO
        END IF
      END IF

C-----------------------------------------------------------------------
C     Close data file.
C-----------------------------------------------------------------------
      CLOSE(3)

C.... Change the input data so all the variables are given starting
C.... from RHUB.  This is needed when IAN=1; otherwise, RHUB will
C.... be change to XR(1) inside psf2align.f (HSL,JY030200)

C Yiran 20170526 Add IAN.eq.6 because psf2 alignment is also called when ian=6
C .............. Not sure about IAN.eq.2 case. (test ian.eq.2 and xr(1).ne.rhub)

      IF(XR(1).NE.RHUB.AND.(IAN.EQ.1 .OR. IAN.EQ.6)) THEN

       DO N=1,NX
        XRSQ(N)=1.0-SQRT(ABS(1.0-XR(N)))
       ENDDO
       CALL UGLYDK(NX,1,1,XR,XPI,0.0,0.0,PICUB)
       CALL UGLYDK(NX,1,1,XR,XRAKE,0.0,0.0,RKCUB)
       CALL UGLYDK(NX,1,1,XR,XSKEW,0.0,0.0,SKCUB)
       CALL UGLYDK(NX,1,1,XRSQ,XCHD,0.0,0.0,CHCUB)
       CALL UGLYDK(NX,1,1,XR,XCI,0.0,0.0,CICUB)
       CALL UGLYDK(NX,1,1,XR,XTI,0.0,0.0,TICUB)
       IF(ITHK.EQ.7) THEN
        CALL UGLYDK(NX,1,1,XRSQ,RLVLRX,0.0,0.0,VLRCUB)
       END IF
       DO I=1,NX
        XR1(I)=XR(I)
       ENDDO

       XR(1) = RHUB
       XRSQ1(1)=1.0-SQRT(ABS(1.0-RHUB))

       CALL EVALDK(NX,1,XR1,XR(1),XPI(1),PICUB)
       CALL EVALDK(NX,1,XR1,XR(1),XRAKE(1),RKCUB)
       CALL EVALDK(NX,1,XR1,XR(1),XSKEW(1),SKCUB)
       CALL EVALDK(NX,1,XRSQ,XRSQ1,XCHD(1),CHCUB)
       CALL EVALDK(NX,1,XR1,XR(1),XCI(1),CICUB)
       CALL EVALDK(NX,1,XR1,XR(1),XTI(1),TICUB)
       IF(ITHK .EQ. 7) THEN
        CALL EVALDK(NX,1,XRSQ,XRSQ1,RLVLRX(1),VLRCUB)
       ENDIF

       IF(ICAM.EQ.99) THEN
        DO JJ=1,15
         CALL UGLYDK(NX,1,1,XR1,XCS(1,JJ),0.0,0.0,CUBTMP)
         CALL EVALDK(NX,1,XR1,XR(1),XCS(1,JJ),CUBTMP)
        ENDDO
       END IF

       IF(ITHK.EQ.99) THEN
        DO JJ=1,16
         CALL UGLYDK(NX,1,1,XR1,XTS(1,JJ),0.0,0.0,CUBTMP)
         CALL EVALDK(NX,1,XR1,XR(1),XTS(1,JJ),CUBTMP)
        ENDDO
       END IF

      ENDIF
C.....End of changes (HSL,JY030200)

C     --120999 (HSLEE) include to change thickness/chord ratio for the ITHK=99 case
C     write(*,*) 'Type the factor you want to increase t/c ratio !!!'
C     wtite(*,*) ' (nochange=1.0 , 50% decrease : 0.5 .... ) '
C     include viscous pitch correction
C     ivpitch = 0 : no  other: yes

      if(ITHK .eq. 99) then
       do n = 1 , nx
        do j = 1 , 16
         xts(n,j) = xts(n,j) * tfactor
        enddo
       enddo
      else
       do n = 1 , nx
        xti(n) = xti(n) * tfactor
       enddo
      endif

C--   120999 (HSLEE)


C--   120999 (HSLEE)
C     -- Read HUB data if IHUB .eq. 6

      IF(IHUB .EQ. 6) THEN
       OPEN(903,FILE="hub.dat",STATUS='UNKNOWN')

       READ(903,*) NHIN
       DO I = 1,NHIN
          READ(903,*) XHIN(I),YHIN(I)
       END DO
       CLOSE(903)
       CALL UGLYDK(NHIN,1,1,XHIN,YHIN,0.0,0.0,XYHINCUB)
      ENDIF

C-----------------------------------------------------------------------
C     Read tip vortex cavity data
C-----------------------------------------------------------------------
      IF(IAN.EQ.2) THEN
       CALL CTR
      ELSE
       ISCAV = 0
      ENDIF

C-----------------------------------------------------------------------
C     Establish the NLEP(MR,NTPREV,2) matrix as the "initial guess"
C     of detachment at NLE on each strip ...CM 032097
C     What is done here is that NLEP(MR,0) is established as
C     a constant leading edge detachment line at the 0th timestep
C-----------------------------------------------------------------------
      DO 7 I=1,2
       DO 5 M=1, MR
        NLEP(M,1,I)=0
 5     CONTINUE
 7    CONTINUE

C-----------------------------------------------------------------------
C     --- CASE FOR ICON = 5  & 8
C     Next 5 lines added by HL for hydrofoil case             CM010598
C-----------------------------------------------------------------------
      IF(rhub.eq.0.0.or. ICON.EQ.5.OR.ICON.EQ.8) THEN
       DO N = 1 , NX
        XPI(N) = XPI(N) * PI / 180.
       ENDDO
      ENDIF

C-----------------------------------------------------------------------
C     Read name of wake file (this part is moved from propcav.f), then
C     read and interpolate data from wake file by calling subroutine
C     readwak.f.
C     JY090799
C-----------------------------------------------------------------------
      WRITE(*,*) ' '
      WRITE(*,*) ' PROPCAV> ENTER INFLOW WAKE FILENAME: '
      READ(*,*) WKFILE
C     SH---------------------ICON = 5------RUDDER----SHREENAATH--04/02/2003--
      IF(ICON.EQ.5)THEN
       WRITE(*,*) ' '
       WRITE(*,*)'RUDDER> ENTER FROUDE NUMBER OF RUDDER: '
       READ(*,*) FNRUDDER

       IF(IFLAP.EQ.1) THEN
        WRITE(*,*)''
        WRITE(*,*)'RUD FLAP> ENTER FLAP PIVOT RATIO TO CHORD, FLAP ANG:'
        READ(*,*) CHFLAP, ANGFLAP
       ENDIF
      ENDIF
C     SH--------------------------------------------------------------------
      CALL READWAK

C-----End of modification (JY090799)------------------------------------

C-----------------------------------------------------------------------
C     Use un-aligned helix wake for fan-line SP-propeller.      JY090300
C-----------------------------------------------------------------------
      IF(((IAN.EQ.1).or.(IAN.eq.6).or.(IAN.eq.2)) .AND. ICON.NE.7 .AND.
C      IF(((IAN.EQ.1).or.(IAN.eq.2).or.(ISC.EQ.1)) .AND. ICON.NE.7 .AND.
     & ISP.NE.1) THEN
       RHULT1=RHULT
       IHTMP = IHUB
       IDTMP = IDUCT
       IDUCT = 0
       CALL PSF2ALIGN(ITHK,ICAM,DAMP,IHUB,IDUCT)

       IHUB = IHTMP
       IDUCT = IDTMP
      ENDIF

C-----------------------------------------------------------------------
C     Modify suction side geometry to give zero T.E. thickness if
C     IMOD=1.                                                   JY111600
C
C     This part of the subroutine used to be called before PSF2ALIGN,
C     I'm moving it here so that the wake geometry remains the same
C     regardless of IMOD.                                       JY071601
C-----------------------------------------------------------------------
      IF(IMOD.EQ.1) THEN

       IF(IFORMAT.NE.2) THEN
        NPOINTS=16
        DO N=1,NPOINTS
         DUMXC(N)=PC(N)
        END DO
        DO J=1,NX
         IF(ICAM.NE.99) THEN
          CRT=XCI(J)
          IF(ICAM.EQ.0) THEN
           CALL A8ML(NPOINTS,CRT,DUMXC,YC1,TC1)
          ELSE IF(ICAM.EQ.1) THEN
           CALL CBPB(NPOINTS,CRT,DUMXC,YC1,TC1)
          END IF
          DO N=1,NPOINTS-1
           XCS(J,N)=YC1(N)
          END DO
         END IF

         DO N=1,NPOINTS
          IF(N.LT.NPOINTS) THEN
           YP(N,J)=XCS(J,N)-XTS(J,N)/2.
           YS(N,J)=XCS(J,N)+XTS(J,N)/2.
          ELSE
           YP(N,J)=-XTS(J,N)/2.
           YS(N,J)=XTS(J,N)/2.
          END IF
         END DO
        END DO
       END IF

       CALL MODGEO(IFORMAT,ICAM,ICON,NX,NPOINTS)

      END IF

C-----------------------------------------------------------------------
C     Calculate SPLINE cubic coefficients for later interpolations
C-----------------------------------------------------------------------
      DO 10 N=1,NX
       XRSQ(N)=1.0-SQRT(ABS(1.0-XR(N)))
 10   CONTINUE
      CALL UGLYDK(NX,1,1,XR,XPI,0.0,0.0,PICUB)
      CALL UGLYDK(NX,1,1,XR,XRAKE,0.0,0.0,RKCUB)
      CALL UGLYDK(NX,1,1,XR,XSKEW,0.0,0.0,SKCUB)
      CALL UGLYDK(NX,1,1,XRSQ,XCHD,0.0,0.0,CHCUB)
      CALL UGLYDK(NX,1,1,XR,XCI,0.0,0.0,CICUB)
      CALL UGLYDK(NX,1,1,XR,XTI,0.0,0.0,TICUB)
      CALL UGLYDK(NX,1,1,XR,XVA,0.0,0.0,VACUB)
      CALL UGLYDK(NX,1,1,XR,XVR,0.0,0.0,VRCUB)
      CALL UGLYDK(NX,1,1,XR,XVT,0.0,0.0,VTCUB)
      CALL UGLYDK(NX,1,1,XR,XUA,0.0,0.0,UACUB)
      CALL UGLYDK(NX,1,1,XR,XUAU,0.0,0.0,UAUCUB)
      CALL UGLYDK(NX,1,1,XR,XUT,0.0,0.0,UTCUB)
      CALL UGLYDK(NX,1,1,XR,XUTU,0.0,0.0,UTUCUB)
      IF(ITHK.EQ.7) THEN
       CALL UGLYDK(NX,1,1,XRSQ,RLVLRX,0.0,0.0,VLRCUB)
      END IF


C-----------------------------------------------------------------------
C     Write input data to the output file
C-----------------------------------------------------------------------
      WRITE(2,900)

      WRITE(2,902)

      WRITE(2,910) ICON,IHUB,ISTEADY,IAN,IFILE,ISP,IMG
     % ,  ITERKT,IPK,ICSPAC,IRSPAC,ITHK
     % ,  ICAM,ISEARCH,IFACE,IVPITCH
     % ,  IFORMAT,IMOD,ISC

      WRITE(2,930)

      WRITE(2,935) NX
      WRITE(2,940) NBLADE,NC,MR

      WRITE(2,950) RHUB, NHBU,MHBT, XHBU,XHBD,XHBT
      WRITE(2,960) ADVCO, RULT,RHULT,DCD,XULT,XUWDK
      WRITE(2,965) FROUDE,SIGMA,DTOL,RLAMDA1,ITERMAX

      WRITE(2,970) (XR(N),XPI(N),XRAKE(N),XSKEW(N),XCHD(N),
     * XCI(N),XTI(N),N=1,NX)

      WRITE(2,980) (XR(N),XVA(N),XVR(N),XVT(N),XUA(N),XUT(N),
     * XUAU(N),XUTU(N),N=1,NX)

C-----------------------------------------------------------------------
C------------------------------FORMAT ---------------------------------
C-----------------------------------------------------------------------
 900  FORMAT(
     % '=========================================================',/
     % '               PROPCAV-R2022 : PANEL METHOD              ',/
     % '                                                         ',/
     % '            The University of Texas at Austin            ',/
     % '          Computational Hydrodynamic Laboratory          ',/
     % '                                                         ',/
     % '========================================================='//)

 902  FORMAT(
     % '---------------- RUN CONTROL OPTIONS -----------------------')

 910  FORMAT(/
     * '     ICON    =  ',I3,'               IHUB    =  ',I3,/
     * '     ISTEADY =  ',I3,'               IAN     =  ',I3,/
     * '     IFILE   =  ',I3,'               ISP     =  ',I3,/
     * '     IMG     =  ',I3,'               ITERKT  =  ',I3,/
     * '     IPK     =  ',I3,'               ICSPAC  =  ',I3,/
     * '     IRSPAC  =  ',I3,'               ITHK    =  ',I3,/
     * '     ICAM    =  ',I3,'               ISEARCH =  ',I3,/
     * '     IFACE   =  ',I3,'               IVPITCH =  ',I3,/
     * '     IFORMAT =  ',I3,'               IMOD    =  ',I3,/
     * '     ISC     =  ',I3//)

 930  FORMAT(
     % '---------------- GEOMETRY INPUT DATA -----------------------')

 935  FORMAT(
     1 /5X,'NX=    ',I2,8X,'  NUMBER OF INPUT RADII'/)

 940  FORMAT(5X,'BLADE PANEL:'
     2 /5X,'NBLADE=',I3,8X,'  N0 OF BLADES',
     3 /5X,'NC    =',I3,8X,'  CHORDWISE PANELS ON THE BLADE'
     4 /5X,'MR    =',I3,8X,'  SPANWISE PANELS ON THE BLADE'/)

 950  FORMAT(5X,'HUB PANEL:'
     1 /5X,'RHUB = ',F6.3,4X,'  HUB RADIUS'
     2 /5X,'NHBU = ',I2,8X, '  NO. OF UPSTREAM PANELS OF HUB'
     3 /5X,'MHBT = ',I2,8X, '  NO. OF CIRCUMFERENTIAL PANELS OF HUB'
     4 /5X,'XHBU = ',F6.3,4X,'  UPSTREAM HUB LENGTH/RADIUS'
     5 /5X,'XHBD = ',F6.3,4X,'  DOWNSTREAM HUB LENGTH/RADIUS'
     6 /5X,'XHBT = ',F6.3,4X,'  HUB TAIL LENGTH/RADIUS')

 960  FORMAT(//
     % 5X,'-----ADVANCE COEFFICIENT AND SLIPSTREAM INFORMATION-----'//
     % 5X,'ADVCO=',F6.3,'   ADVANCE COEFFICIENT,J, BASED ON SHIP SPEED'/
     % 5X,'RULT =',F6.3,'   ULTIMATE TIP VORTEX RADIUS'/
     % 5X,'RHULT=',F6.3,'   ULTIMATE HUB VORTEX RADIUS'/
     % 5X,'DCD  =',F6.3,'   PROP ROTATION INCREMENT FOR WAKE CONVECTION'/
     % 5X,'XULT =',F6.2,'   AXIAL POSITION OF ULTIMATE WAKE'/
     % 5X,'XUWDK=',F6.2,'   AXIAL POSITION OF ULT. WAKE SINK DISK')

 965  FORMAT(/5X,'-------CAVITY FLOW INFORMATION-------'//5X,
     * 'FROUDE =',F8.3,'  FROUDE NUMBER BASED ON ROTATION'/5X,
     * 'SIGMA  =',F8.3,'  CAVITATION NUMBER BASED ON ROTATION'/5X,
     * 'DTOL   =',F8.3,'  TOLERANCE FOR CAVITY OPENNESS'/5X,
     * 'RLAMDA1=',F8.3,'  LAMDA/CAVL PRESSURE MODEL ZONE LENGTH'/5X,
     * 'ITERMAX=',I4,'      MAX NUMBER OF PLANFORM ITERATIONS FOR GIVEN
     * CROSSFLOW'/5X)

 970  FORMAT(/5X,
     % '-------------------INPUT BLADE GEOMETRY------------------------'/
     % /5X,' R/RO      P/D      XS/D     SKEW   C/D      FO/C     TO/D'/
     % (3X,F7.3,F10.4,F9.4,F9.3,F8.4,2F9.4))
 980  FORMAT(1H1,/5X,'-----INFLOW VELOCITIES-----    TRANS WAKE INDUCED
     * VELOCITIES'/5X,' R/RO   VX/VS  VR/VS  VT/VS    UANW   UTNW      UA
     * UW   UTUW'/(F11.4,3F7.3,2X,2F7.3,3X,2F7.3))

 2000 CONTINUE

      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine PROINP>>>>>>>>>>>>>>>>>>>>>>>>>
      END


!s--- YE TIAN 07/04/2013-----
      subroutine  close_te(nx, nxmax, y, x_trunc)
      implicit none
      integer nx, nxmax
      real y(nxmax,16), x_trunc
      real x(16)
      data x/0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,
     &       0.8 ,0.9,  0.95,0.975, 0.99, 1.0/
      real y_trunc,dx,dy,ydiff

      integer i, j, k

      do j = 1, 16
        if (x(j).ge.x_trunc) then
          k = j
          exit
        end if
      end do
      dx = 1.0 - x_trunc
      do i = 1, nx
        dy = y(i,16)
        do j = k, 16
          ydiff = ((x(j)-x_trunc)/dx)**2*dy
          y(i,j) = y(i,j) - ydiff
        end do
      end do
      end subroutine close_te

      subroutine  close_te2(nx, npt, mxinp, x, y, x_trunc)
      implicit none
      integer nx, mxinp, npt
      real x(npt), y(mxinp,nx), x_trunc
      real dx,dy,ydiff

      integer i, j, k

      do j = 1,npt
        if (x(j).ge.x_trunc) then
          k = j
          exit
        end if
      end do
      dx = 1.0 - x_trunc
      do i = 1, nx
        dy = y(npt,i)
        do j = k, npt
          ydiff = ((x(j)-x_trunc)/dx)**2*dy
          y(j,i) = y(j,i) - ydiff
        end do
      end do
      end subroutine close_te2

!e--- YE TIAN 07/04/2013-----
