      module uns_output
        real,allocatable :: unsforce(:,:)
      end module

      subroutine writefileuns(n_step,n_bld0)
************************************************************************
*                                                                      *
*     writefileuns                                                     *
*                                                                      *
*     --- Outputs for unsteady version of PROPCAV                      *
*                                                                      *
*     By Yiran Su 09/19/2017                                           *
*                                                                      *
************************************************************************
      use uns_output
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'
      include 'PUFCAVC.INC'
      character*20 uns_fn

      n_bld = mod(n_bld0, nblade)
      if (n_bld.eq.0) n_bld = nblade

************************************************************************
*    OPEN file at the beginning (once)
************************************************************************
      if (n_step.eq.0) then

        !! ===== #1: Unsteady forces all blades =====
        do i = 1,nblade
          write(uns_fn, "(A9,I1,A4)") 'UNS_ktkq_', i, '.plt'
          open(1080+i, file = uns_fn, status = 'unknown')
          write(1080+i, *) 'Variables = Tstep, kt, kq, ktv, kqv'
          write(1080+i, *) 'zone T = "Unsteady forces for blade',i,'"'
        end do

        !! ===== #2: Blade Circulation =====
        do i = 1,nblade
          write(uns_fn, "(A8,I1,A4)") 'UNS_cir_', i, '.plt'
          open(1180+i, file = uns_fn, status = 'unknown')
          write(1180+i,*) 'Title="CIRCULATION:Gs=100*DELP/(2*pi*R*Vs)"'
          write(1180+i,*) 'Variables="r/R","100*DELP/(2*pi*R*Vs)"'
        end do

        !! == end subroutine
        return

      end if

************************************************************************
*    CLOSE file and exit at final step (once)
************************************************************************
      if (n_step.lt.0) then

        !! ===== #1: Unsteady forces all blades =====
        do i = 1,nblade
          close(1080+i)
        end do

        !! ===== #2: Blade Circulation =====
        do i = 1,nblade
          close(1180+i)
        end do

        !! == end subroutine
        return

      end if

************************************************************************
*    Main part of the iteration (every unsteady inner iteration)
************************************************************************

      !! ===== #1: Unsteady forces all blades =====
      write(1080+n_bld, *) n_step,fbxp(1),fbxp(4),fbxpv(1),fbxpv(4)

      !! ===== #2: Blade Circulation =====
      write(1180+n_bld, *) 'zone T = "blade ',n_bld,'"'
      do i = 1,mr
        write(1180+n_bld, *) hrzp(1,i), delp(i)*100.0/2.0/3.14159265359
      end do

      end subroutine


      subroutine writepressure(n_step)
        include 'PUFCAV.INC'
        include 'PUFCAVB.INC'
        include 'PUFCAVC.INC'
        character*20 uns_fn

        write(uns_fn, "(A9,I4.4,A4)") 'UNS_wprs_', n_step, '.plt'
        open(1021, file = uns_fn, status = 'unknown')
        write(1021, *) 'TITLE="Wetted Pressures on Blade"'
        do m = 1,mr,3
          j  = indexb(nc/2, m)
          rr = sqrt(xct(j,2)**2 + xct(j,3)**2)
          write(1021, "(A12,F4.2,A1)") 'Zone T="r/R=', rr, '"'
          do n = 1, nc
            j = indexb(n, m)
            write(1021, *) xct(j, 1), -cpb(n, m)*advco**2
          end do
        end do
        close(1021)
      end subroutine
