C ====================================
      SUBROUTINE OPFILET
C ====================================
      
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      
c      open(21,file='wprsc.plt',status='unknown')
      open(910,file='3dgeom-1rev.dat',status='unknown')
      open(911,file='3dgeom-2rev.dat',status='unknown')
C/s S.N.KIM | increased geometry files up to 5th and the last rev. 
      open(913,file='3dgeom-3rev.dat',status='unknown')
      open(914,file='3dgeom-4rev.dat',status='unknown')
      open(915,file='3dgeom-5rev.dat',status='unknown')
      open(916,file='3dgeom-Lrev.dat',status='unknown')
C/e S.N.KIM | Aug. 2018.

c      write(21,210)
      write(910,*) 'Variables ="x","y","z"'
      write(911,*) 'Variables ="x","y","z"'
C/s S.N.KIM | increased geometry files up to 5th and the last rev.
      write(913,*) 'Variables ="x","y","z"'
      write(914,*) 'Variables ="x","y","z"'
      write(915,*) 'Variables ="x","y","z"'
      write(916,*) 'Variables ="x","y","z"'
C/e S.N.KIM | Aug. 2018.

c 210  FORMAT('TITLE = "CP and VT on cavity" '/
c     %     'VARIABLES="X","-Cp_n","Vt","POT"')

C***********************************************************************
C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
C***********************************************************************
c      ncvx = nwk-1
c      ncvxp = ncvx + 1
c      nthx = nc/2
c      nthxp = nthx + 1
c      mcvtp = mcvt + 1
c      npanc = ncvx * mcvt
c      npant = nthx * mcvt
c      npanel = npanb + npanh + npanc + npant
c      nsapw = nwk*mrp
      
      ncvx = 0 !nwk-1
      ncvxp = 0 !ncvx + 1
      nthx = 0 !nc/2
      nthxp = 0 !nthx + 1
      mcvtp = 0 !mcvt + 1
      npanc = 0 !ncvx * mcvt
      npant = 0 !nthx * mcvt
      npanel = npanb + npanh !+ npanc + npant
      nsapw = nwk*mrp

      iunsplot = 1 ! parameter for plotting Unsteady wake history   
C***********************************************************************
C/e S.N.KIM | Aug. 2018.
C***********************************************************************
      
      nrecl2 = nsapw*8
      n2 = nrecl2/1024
      nleft2=nrecl2 - n2*1024
      if(nleft2 .eq. 0) then
         nreclw = 1024 * n2
      else
         nreclw = 1024 * (n2+1)
      endif
      
C-- Open Files to save geometry at each time step

C -- Files for the influence of Blade & Hub

C      OPEN(67,FILE='AA_BH.DAT',STATUS='UNKNOWN'
C     %     ,FORM='UNFORMATTED')
C      OPEN(68,FILE='BB_BH.DAT',STATUS='UNKNOWN'
C     %     ,FORM='UNFORMATTED')

C ---------------------------------------------------
      
      OPEN(143,FILE='GEOWX.DAT',STATUS='UNKNOWN',
     *     FORM='UNFORMATTED',
     *     ACCESS='DIRECT',RECL=NRECLW)
      OPEN(144,FILE='GEOWY.DAT',STATUS='UNKNOWN',
     *     FORM='UNFORMATTED',
     *     ACCESS='DIRECT',RECL=NRECLW)
      OPEN(145,FILE='GEOWZ.DAT',STATUS='UNKNOWN',
     *     FORM='UNFORMATTED',
     *     ACCESS='DIRECT',RECL=NRECLW)
      OPEN(146,FILE='GEOWK.DAT',STATUS='UNKNOWN',
     *     FORM='UNFORMATTED',
     *     ACCESS='DIRECT',RECL=NRECLW)
      
      RETURN
      END
