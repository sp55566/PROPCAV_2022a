C ------------------------------------------
      subroutine callgeok(ik)
C ------------------------------------------
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'
      include 'PUFCAVC.INC'

!     parameter (nmz=mbpz*nwpz)
!     dimension geox(nmz),geoy(nmz),geoz(nmz)

!s--- YE TIAN 07/07/2013-------
      integer nmz
      real,allocatable :: geox(:),geoy(:),geoz(:)

      nmz = mbpz*nwpz

      if (.NOT. allocated(geox)) then
        allocate(geox(nmz),geoy(nmz),geoz(nmz))
      end if
!e--- YE TIAN 07/07/2013-------

C ---------------------------------------------------------
C  Call wake geometry data From files
C       of kk-th blade angle      
C ---------------------------------------------------------

      if(ik .eq. 1) then
c        call read2(146,1,geox,nsapw)
c        call read2(146,2,geoy,nsapw)
c        call read2(146,3,geoz,nsapw)
c/s S.N.KIM | Key wake was read from the steady results previously, 
C             making the key wake same for all time step locations. I modified it such that 
C             the key wake at certain time step is read from 
C             the preceding time step in unsteady case.
        icc = 0
        do n = 1, nwk
          do m = 1, mrp
            icc = icc + 1
            geox(icc) = xw(n,m)
            geoy(icc) = yw(n,m)
            geoz(icc) = zw(n,m)
          enddo
        enddo
C/e S.N.KIM | Unsteady Wake Alignment. | Aug. 2018.
      elseif(ik .eq. 2) then
        irec = ntpos(1)
          call read2(143,irec,geox,nsapw)
          call read2(144,irec,geoy,nsapw)
          call read2(145,irec,geoz,nsapw)
      endif

      icc = 0      
      do n = 1 , nwk
        do m = 1 , mrp
         icc = icc + 1
         xw(n,m) = geox(icc)
         yw(n,m) = geoy(icc)
         zw(n,m) = geoz(icc)
        enddo
      enddo

C/s S.N.KIM | Once the key wake is updated, subpanels need to be regenerated accordingly.
      if(ik.eq.2) call gwakes(nsub,nwsub1)
C/e S.N.KIM | Unsteady Wake Alignment. | Aug. 2018.        

      return
      end
 

C ------------------------------------------
      subroutine savegeok
C ------------------------------------------
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'
!     parameter (nmz=mbpz*nwpz)
!     dimension geox(nmz),geoy(nmz),geoz(nmz)
!s--- YE TIAN 07/07/2013-------
      integer nmz
      real,allocatable :: geox(:),geoy(:),geoz(:)

      nmz = mbpz*nwpz

      if (.NOT. allocated(geox)) then
        allocate(geox(nmz),geoy(nmz),geoz(nmz))
      end if
!e--- YE TIAN 07/07/2013-------

      icc = 0
      do n = 1 , nwk
        do m = 1 , mrp
           icc = icc + 1
           geox(icc) = xw(n,m) 
           geoy(icc) = yw(n,m) 
           geoz(icc) = zw(n,m) 
        enddo
      enddo

      call write2(146,1,geox,nsapw)
      call write2(146,2,geoy,nsapw)
      call write2(146,3,geoz,nsapw)
 
      return
      end
 
