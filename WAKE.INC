C**************************** MIT PSF10 ********************************
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
C***********************************************************************

C-----------------------------------------------------------------------
C     Next common blocks are for parameters
C-----------------------------------------------------------------------
      PARAMETER (MAXB=12, NXMAX=17,MCMAX=21,NXMAX2=NXMAX**2,NGPMAX=161,
     *           NGMAX2=NGPMAX*NGPMAX,NG2MAX=3*NGPMAX,MAX=15000)

      COMMON/PARAM/ NPAN(MAXB),MPAN(MAXB),NGRID(MAXB),MGRID(MAXB), 
     *              NGPSUM(MAXB),NCPSUM(MAXB),NCVSUM(MAXB),
     *		    NTWSUM(MAXB),NSW(MAXB,MCMAX),KSPACE,
     *              KTHICK,NCASE,IRUN,IPASS,IDAT,JDAT,KDAT,NCPTOT

C-----------------------------------------------------------------------
C     Next common blocks are for administration
C-----------------------------------------------------------------------
      COMMON/ADMNWAK1/ VOLMNV,ADVCOM,DIAM,RADIUS,RPM,RPS,OMEGA,VSKTS,
     *         VSFPS,RHO,CDRAG,EAR,KSEPTV,DISPN,DISPR,SFC,UWLAM,
     *         VA(MCMAX),VT(MCMAX),VR(MCMAX),XVSTAR(NXMAX),VSTAR(MCMAX),
     *         XCDRAG(NXMAX),XSFC(NXMAX),XUAN(NXMAX),XUTN(NXMAX)
      COMMON/ADMNWAK2/ NPUW, MLTYPE, MTHICK,ICDRAG,ISFC

C-----------------------------------------------------------------------
C     Next common blocks are for geometries
C----------------------------------------------------------------------
      COMMON/GEOINP/ XRFAC(NXMAX2),IPVXR(NXMAX),XRSFAC(NXMAX2),
     *               IPVXRS(NXMAX),XCAM(NXMAX),XTHICK(NXMAX),
     *               PER(17),PSQ(17),XTW
      COMMON /INPGEO1/ XCS(NXMAX,15),XTS(NXMAX,16)

      COMMON/GEOWAK/ X(1200),Y(1200),Z(1200),XP(144),YP(144),ZP(144),
     *               XON(144),YON(144),ZON(144),XCOLL,YCOLL,ZCOLL,
     *               XTIP(18,19),YTIP(18,19),ZTIP(18,19),XW(8640),
     *               YW(8640),ZW(8640),XU(61,12),YU(61,12),ZU(61,12),
     *               X1HV,Y1HV,Z1HV,X2HV,Y2HV,Z2HV,RZ(21),RZSQ(21),
     *               R(20),RFAC(400),IPVR(20),RSQ(20),PITCH(MCMAX),
     *               RAKE(MCMAX),SKEW(MCMAX),CAMBR(MCMAX),THICK(MCMAX),
     *               CHORZ(MCMAX),CHORCP(20),XV(MCMAX),XC(MCMAX),
     *               HC(MCMAX),HT(MCMAX),FC(MCMAX),CPDC(MCMAX)

C-----------------------------------------------------------------------
C     Next common block is for the interation
C-----------------------------------------------------------------------
      COMMON/INTERP1/ AS(20,10),AC(20,10),DIST(40)
      COMMON/INTERP2/ LKB(40),INDS(10),INDC(10)

C-----------------------------------------------------------------------
C     Next common blocks are for the influence-coefficient computations
C-----------------------------------------------------------------------
      COMMON/QVEL/ UBN(20),UTN(30,2),SOURCN(144)
      COMMON/VORSOR/ GB(768),GC(864),SB(768),GAMTW(12,21),
     *               SFACT(20),GAMULT

C-----------------------------------------------------------------------
C     Next common block is for the matrix solution
C-----------------------------------------------------------------------
      COMMON/MATRIX/ AA(20736),AB(144,144),AOBC(144,72),RHS(144)
      COMMON/MATRIX2/ IPIVA(144)
      COMMON/VISCOUS/ IVPITCH
      COMMON/CTICHK/ CT
      COMMON/CHARS/  ADATE,ATIME,PSFGEO,PSFADM,PSFF14,TDATE

      CHARACTER*9    ADATE
      CHARACTER*11   ATIME
      CHARACTER*23   PSFGEO,PSFADM,PSFF14,TDATE


