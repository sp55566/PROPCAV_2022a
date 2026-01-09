C ----------------------------------------------
      subroutine ctr
C ----------------------------------------------
      include 'PUFCAV.INC'
      
      open(3,file='cav.ctr',status='old')

C--------------------------------------------------------------------------
C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
C           | I would not delete the old version of input file, rather
C           | comment out for later use in any case.
C--------------------------------------------------------------------------
cC --- mc : number of circumferential panel of tip vortex cavity
cC     icavmax : Maximum iteration of tip vortex cavity shape
cC     deltam  : distance from the C.P. to misplaced point
cC     nnrev : starting revolution for cavity run
c
c      read(3,*) mcvt, nwk
c
c      IF(NWK.GT.NWPZ) THEN
c         WRITE(*,*) 'Need to increase NWPZ TO NWK=',NWK
c         STOP
c      END IF
c
cC/s S.N.KIM | includes a new parameter 'iualign_div' to use in FWA at
cC           | the steady revolution of unsteady runs.
cC           | Note: if IAN=2, Steady Rev. -> Unsteady Rev. -> Fully Unsteady Rev.
cC           | Aug. 2018.
c      read(3,*) icavmax, iscav, iualign_div
c      read(3,*) deltam
c      read(3,*) radini
c      if(iscav .eq. 2) read(3,*) nnrev
C-------------------------------------------------------------------------
C/e S.N.KIM | Aug. 2018.
C-------------------------------------------------------------------------


C --- mc : number of circumferential panel of tip vortex cavity
C     icavmax : Maximum iteration of tip vortex cavity shape
C     deltam  : distance from the C.P. to misplaced point
C     nnrev : starting revolution for cavity run
      read(3,*) nwk, iuplot
      IF(NWK.GT.NWPZ) THEN
         WRITE(*,*) 'Need to increase NWPZ TO NWK=',NWK
         STOP
      END IF
      mcvt = 0
      deltam = 0.0
      radini = 0.0
C/s S.N.KIM | includes a new parameter 'iualign_div' to use in FWA at
C           | the steady revolution of unsteady runs.
C           | Note: if IAN=2, Steady Rev. -> Unsteady Rev. -> Fully Unsteady Rev.
C           | Aug. 2018.
      read(3,*) icavmax
      read(3,*) iualign_div
      read(3,*) iscav
      if(iscav .eq. 2) read(3,*) nnrev
      close(3)

      return
      end
