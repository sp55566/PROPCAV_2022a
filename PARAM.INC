************************************************************************
*   Parameters:                                                        *
*     NSTEP number of blade angles (timesteps) per revolution          *
*     KZ    number of blades                                           *
*     NBZ   number of panels of the blade chordwise                    *
*     MBZ   number of panels of the blade the spanwase                 *
*     NWZ   number of panels in the wake for one radii                 *
*     NHDZ  number of streamwise panels in hub downstream              *
*     NHUZ  number of streamwise panels in hub upstream                *
*     MHBZ  number of panels of the hub circumferentially              *
*     NZSR  number of panels in the SR region                          *
*     MCAVM number of circumferential panels on the tip vortex         *
*                                                                      *    
*   IMPORTANT notes to minimize meory requirement:                     *
*   (Please refer to Table 1 in page 6 of PROPCAV v2.0 manual for      *
*    additional recommendations.)                                      *
*   ------------------------------------------------------------------ *
*     If ISC=0, then set NZSR=0 (i.e. use C1A instead of C1B).         *
*                                                                      *
*     If IHUB=0 & IAN=0 or 1, then set NHDZ=0 & NHUZ=0 (i.e. use C2A   *
*                             instead of C2C).                         *
*                                                                      *
*     If IHUB>0 & IAN=2, then set NHDZ=NWZ & NHUZ=10 (i.e. use C2B     *
*                        instead of C2C).                              *
*                                                                      *
*     If IAN=0 or 1, then set NPHC=0 and NPHT=0 (i.e. use C3A instead  *
*                    of C3B).                                          *
*                                                                      *
*   VERY IMPORTANT:                                                    *
*   ------------------------------------------------------------------ *
*     If you make any change to this file, you need to RECOMPILE all   *
*     the source codes for PROPCAV.                                    *
*                                                                      *    
************************************************************************


!s--- YE TIAN ----- 07/05/2013-----
!     Try to make those parameters determined dynamically.

      use m_param
!e--- YE TIAN ----- 07/05/2013-----


!     PARAMETER(KZ=5)

C.....dimension variable for timestep size (added, JY052600)
C     PARAMETER(NSTEP=90,NDEL=360/NSTEP)
!     PARAMETER(NSTEP=60,NDEL=360/NSTEP)
!     PARAMETER(NDEL=6)
      
C.....Maximum Number of Input radii
      PARAMETER (NXMAX=17,MAXSPL=4*(NXMAX-1) )

C.....Maximum Number of Input data for Modified geometry
!s--- Allen Du 06/07/2018 change the parameter from 20 to 43 so some cases will
!run with out a problem.
!     PARAMETER (MXINP=20)
      PARAMETER (MXINP=43)
!e--- Allen Du 06/07/2018

C.....dimension variable for SR panels 
C1A   PARAMETER(NZSR=0,NZSR2=NZSR/2,NZSRP=NZSR+1,NZSR2P=NZSR2+1)	
      PARAMETER(NZSR=10,NZSR2=NZSR/2,NZSRP=NZSR+1,NZSR2P=NZSR2+1)
C     PARAMETER(NZSR=2,NZSR2=NZSR/2,NZSRP=NZSR+1,NZSR2P=NZSR2+1)

C.....dimension variables for the blade
!     PARAMETER(NBZ=80+NZSR,MBZ=110,NBPZ=NBZ+1,MBPZ=MBZ+1,NBHZ=NBZ/2,
!    *          NBHZP=NBHZ+1)

C.....dimension variables for the wake
C     PARAMETER(NWZ=300,NWPZ=NWZ+1)
C     PARAMETER(NWZ=240,NWPZ=NWZ+1)
!     PARAMETER(NWZ=180,NWPZ=NWZ+1)
C     PARAMETER(NWZ=150,NWPZ=NWZ+1)
C     PARAMETER(NWZ=120,NWPZ=NWZ+1)
C     PARAMETER(NWZ=90,NWPZ=NWZ+1)

C.....dimension variables for the supercavity stuff
!     PARAMETER(NSCZ=NBHZ+NWZ,NSCZP=NSCZ+1,NBSCZ=NBZ+NWZ,NBSCZP=NBSCZ+1)

C.....dimension variables for the hub
C2A   PARAMETER(NHUZ=0,NHDZ=0,NHBZ=NHUZ+NHDZ+NBHZ,NHPZ=NHBZ+1,
C2B   PARAMETER(NHUZ=10,NHDZ=NWZ,NHBZ=NHUZ+NHDZ+NBHZ,NHPZ=NHBZ+1,
C2C   PARAMETER(NHUZ=10,NHDZ=40,NHBZ=NHUZ+NHDZ+NBHZ,NHPZ=NHBZ+1,
C     PARAMETER(NHUZ=10,NHDZ=30,NHBZ=NHUZ+NHDZ+NBHZ,NHPZ=NHBZ+1,

C     PARAMETER(NHUZ=20,NHDZ=NWZ,NHBZ=NHUZ+NHDZ+NBHZ,NHPZ=NHBZ+1,
C     PARAMETER(NHUZ=10,NHDZ=NWZ,NHBZ=NHUZ+NHDZ+NBHZ,NHPZ=NHBZ+1,
C     PARAMETER(NHUZ=20,NHDZ=20,NHBZ=NHUZ+NHDZ+NBHZ,NHPZ=NHBZ+1,
!     PARAMETER(NHUZ=20,NHDZ=100,NHBZ=NHUZ+NHDZ+NBHZ,NHPZ=NHBZ+1,
C     PARAMETER(NHUZ=20,NHDZ=110,NHBZ=NHUZ+NHDZ+NBHZ,NHPZ=NHBZ+1,
C     PARAMETER(NHUZ=10,NHDZ=20,NHBZ=NHUZ+NHDZ+NBHZ,NHPZ=NHBZ+1,
C     PARAMETER(NHUZ=5,NHDZ=30,NHBZ=1,NHPZ=NHBZ+1,
!    *          MHBZ=20,MHPZ=MHBZ+1)

C.....for tunnel
!     PARAMETER(NtMAXt=3, NCMAXT=3)
C     PARAMETER(NtMAXt=200, NCMAXT=20)
!     PARAMETER(NPANMt=NTMAXt*nCmaxt)

C.....for Duct
C     PARAMETER (NDMAX=160, MDMAX=20,NDWMX=100)
C     PARAMETER (NDMAX=200, MDMAX=20,NDWMX=100)
!     PARAMETER (NDMAX=3, MDMAX=3,NDWMX=3)
!     PARAMETER (NDMAXP= NDMAX+1, MDMAXP = MDMAX + 1, NDWMXP=NDWMX+1)
!     PARAMETER (NPANMD=NDMAX*MDMAX, NDWMN=NDWMX*MDMAX)

C.....dimension variables for the tip vortex cavity
C      parameter (mrmax = 3, nwmax=3)	
C      parameter (mcavm = 3 , mcavmp=3)
C      parameter (ncavm = 3, ncavmp = 3)
C      parameter (nhmx=3, nhmxp = 3)
!      parameter (mrmax = mbz, nwmax=nwz)
!      parameter (ncavm = nwz, ncavmp = ncavm+1)
!      parameter (mcavm = 10 , mcavmp=mcavm + 1)
!      parameter (nhmx=nbz/2, nhmxp = nhmx + 1)
!      parameter (mrmaxp = mrmax + 1, nwmaxp=nwmax + 1)

	
C.....dimension variables for solution
C     PARAMETER(NPHZ=3)
!     PARAMETER(NPHZ=NHBZ*MHBZ)
!     parameter(nphc=mcavm*ncavm, NPHT = MCAVM*NHMX)
!     PARAMETER(NPAWZ=MBZ*(NWZ+NDEL),NPAWZ2=NPAWZ)
!     PARAMETER(NSCWZ=NPAWZ)

!     PARAMETER(NPBZ=NBZ*MBZ,NPANZ=NPBZ+NPHZ+NPANMT+NPHC+NPHT+NPANMD)
!     PARAMETER(NTZ=NPANZ+NPAWZ2)
!     PARAMETER (NBLKMAX=MBZ+NHBZ+NTMAXT+NCAVM+NHMX+MDMAX)

C.....dimension variables for BL code   by Hong Sun
!     PARAMETER(NWZVIS=2*NWPZ+NDEL,NZVIS=NBPZ+NWZVIS)
c-----for duct case
!     PARAMETER(NDWZVIS=NDWMX+NWZ+NDEL,NZVISD=NDMAXP+NDWZVIS)

