C ------------------------------------------
      subroutine callgeow
C ------------------------------------------
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'
!       parameter (nmz=mbpz*nwpz)
!       dimension geox(nmz),geoy(nmz),geoz(nmz)
!s--- YE TIAN 07/07/2013-------
      integer nmz
      real,allocatable :: geox(:),geoy(:),geoz(:)

      nmz = mbpz*nwpz

      if (.NOT. allocated(geox)) then
        allocate(geox(nmz),geoy(nmz),geoz(nmz))
      end if
!e--- YE TIAN 07/07/2013-------
C
C
C  Call Tip hub, Tip Vortex and wake geometry data From files
C       of kk-th blade angle      
C ----------------------------------------------------------

      do kk = 2 , nblade

        ik = kk - 1

        irec = ntpos(kk)

        call read2(143,irec,geox,nsapw)
        call read2(144,irec,geoy,nsapw)
        call read2(145,irec,geoz,nsapw)
      
        icc = 0
        do n = 1 , nwk
          do m = 1 , mrp
           icc = icc + 1
           xwo(n,m,ik) = geox(icc)
           ywo(n,m,ik) = geoy(icc)
           zwo(n,m,ik) = geoz(icc)
          enddo
        enddo

      enddo

100      format(2i5,3(1x,f12.5))

      return
      end
 









