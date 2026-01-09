C ------------------------------------------
      subroutine savegeow
C ------------------------------------------
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'
!     parameter (nmz=mbpz*nwpz)
      dimension geox(mbpz*nwpz),geoy(mbpz*nwpz),geoz(mbpz*nwpz)
      integer nmz
      nmz=mbpz*nwpz
C
C
C  Save Blade, Tip hub, Tip Vortex and wake geometry data into files
C      for each time steps
C ----------------------------------------------------------

      nrec = ntprev
C      write(*,*) ' nrec = ', nrec
      icc = 0
      do n = 1 , nwk
        do m = 1 , mrp
           icc = icc + 1
           geox(icc) = xw(n,m) 
           geoy(icc) = yw(n,m) 
           geoz(icc) = zw(n,m) 
        enddo
      enddo

C ----- Save Steady Results ------------------
C
      if(ntstep .eq. 0) then
C         write(*,*) ' save geo : NTSTEP = 0 ' 
C
C ---------------------------------------------
            
        do kk = 1 , nrec 
          call write2(143,kk,geox,nsapw)
          call write2(144,kk,geoy,nsapw)
          call write2(145,kk,geoz,nsapw)
        enddo

C --------Update Unsteady geometry -------
C
      Else
C
C ----------------------------------------

         irec = itstep / ndltat + 1
                   
C         write(*,*) ' save geo irec = ', irec,nsapw

         call write2(143,irec,geox,nsapw)
         call write2(144,irec,geoy,nsapw)
         call write2(145,irec,geoz,nsapw)

C -----------------------------------------------
      Endif
C -----------------------------------------------
 
      return
      end
 
