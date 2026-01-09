! ----------------------------------------------------------------------------- !
!                       CREATED BY YIRAN SU 10/04/2017                          !
!                    Store, calculation, and Write of body Force!               !
! ----------------------------------------------------------------------------- !

      module M_UNS_PC2NS_BF
      use M_UNS_GLOBAL_PAR

      public  :: UNS_BF_INIT, UNS_BF_CALC, UNS_BF_EXTRA, UNS_WRT_BF
      private :: bfcalc_2sid, test_plot

      ! Third index "it" is PROPCAV blade angle index ! Not the theta-index for RANS mesh.
      double precision,dimension(:,:,:),allocatable::pbfx,pbfr,pbft,psrc

      contains

! Initialization
        subroutine UNS_BF_INIT
          implicit none
          allocate(pbfx(pc_nh,pc_mr,pc_ntpr))
          allocate(pbfr(pc_nh,pc_mr,pc_ntpr))
          allocate(pbft(pc_nh,pc_mr,pc_ntpr))
          allocate(psrc(pc_nh,pc_mr,pc_ntpr))
          pbfx = 0.0
          pbfr = 0.0
          pbft = 0.0
          psrc = 0.0
        end subroutine

! Calculate Body Force -- wrapper --
        subroutine UNS_BF_CALC(idxrev)
          implicit none
          integer idxrev
          double precision,dimension(pc_nh, pc_mr)::bfx,bfr,bft,src
          integer i, j, k

          call bfcalc_2sid(bfx, bfr, bft, src)

          if (idxrev.eq.0) then
            do k=1,pc_ntpr
              do i=1,pc_nh
                do j=1,pc_mr
                  pbfx(i,j,k) = bfx(i,j)
                  pbfr(i,j,k) = bfr(i,j)
                  pbft(i,j,k) = bft(i,j)
                  psrc(i,j,k) = src(i,j)
                end do
              end do
            end do
          else
            do j=1,pc_mr
              do i=1,pc_nh
                pbfx(i,j,idxrev) =  bfx(i,j)
                pbfr(i,j,idxrev) =  bfr(i,j)
                pbft(i,j,idxrev) =  bft(i,j)
                psrc(i,j,idxrev) =  src(i,j)
              end do
            end do
          end if

          ! Write for Debugging Purpose
          ! write(*,*) 'idxrev = ', idxrev

        end subroutine

! Extrapolate from the previous two time steps to get current body force field
        subroutine UNS_BF_EXTRA
          use M_UNS_THETA, only : it_prop
          implicit none
          integer ib, it, it1, it2
          do ib = 1, pc_nbld
            it = it_prop(ib)
            it1 = it - 1
            it2 = it - 2
            if (it1 < 1) it1 = it1 + pc_ntpr
            if (it2 < 1) it2 = it2 + pc_ntpr
            pbfx(:,:,it) = pbfx(:,:,it1) * 2.0 - pbfx(:,:,it2)
            pbfr(:,:,it) = pbfr(:,:,it1) * 2.0 - pbfr(:,:,it2)
            pbft(:,:,it) = pbft(:,:,it1) * 2.0 - pbft(:,:,it2)
            psrc(:,:,it) = psrc(:,:,it1) * 2.0 - psrc(:,:,it2)
          end do
        end subroutine

! Calculate body force -- core --
        subroutine bfcalc_2sid(bfx,bfr,bft,src)
        use GEOCAMB, only : dirc, indexCamber
        include 'PUFCAV.INC'
        include 'PUFCAVB.INC'
        double precision,dimension(nh,mr)::bfx,bfr,bft,src
        double precision :: tc_fx,tc_fy,tc_fz,tc_sc,tc_tt,tc_rr,tc_ftot
        double precision :: totfx, bf_tot
        bfx = 0.0
        bfr = 0.0
        bft = 0.0
        src = 0.0
        ! Generate the forces
        do m = 1,mr
          do n = 1,nc
            j = indexb(n,m)
            n1 = nhp - n
            if (n .gt. nh) n1 = n - nh
            tc_rr=sqrt(xct(j,2)**2+xct(j,3)**2)
            tc_tt=atan2(xct(j,3),xct(j,2))
            ! force on a bem panel (-1 means propeller to flow)
            ! << cpb is p*2/(rho*vs^2) >> need to /2.0e0
            tc_fx = -1.0*cpb(n,m)*vel(j,1)*ss(j,1)/2.0
            tc_fy = -1.0*cpb(n,m)*vel(j,2)*ss(j,1)/2.0
            tc_fz = -1.0*cpb(n,m)*vel(j,3)*ss(j,1)/2.0
            tc_sc = -1.0*dpdn(j)*ss(j,1) ! unit kg/s

            ! Test_Normal: Use only the normal component
            ! j = indexCamber(n1, m)
            ! bf_tot = tc_fx * dirc(j,3,1) + tc_fy * dirc(j,3,2)
            ! &                                   + tc_fz * dirc(j,3,3)
            ! tc_fx = bf_tot * dirc(j,3,1)
            ! tc_fy = bf_tot * dirc(j,3,2)
            ! tc_fz = bf_tot * dirc(j,3,3)

            ! force and mass source on a panel! Not density
            bfx(n1,m) = bfx(n1,m) + (tc_fx)
            bfr(n1,m) = bfr(n1,m) + (tc_fy*cos(tc_tt)+tc_fz*sin(tc_tt))
            bft(n1,m) = bft(n1,m) + (tc_fz*cos(tc_tt)-tc_fy*sin(tc_tt))
            src(n1,m) = src(n1,m) + (tc_sc)
          end do
        end do

        ! Skin friction
        do m = 1,mr
          do n = 1,nc
            j = indexb(n,m)
            n1 = nhp - n
            if (n .gt. nh) n1 = n - nh
            tc_tt=atan2(xct(j,3),xct(j,2))

            tc_ftot = xcdf * 0.5 * vtots(j) * ss(j,1)
            tc_fx = -1.0 * uxtot(n,m)/sqrt(vtots(j)) * tc_ftot
            tc_fy = -1.0 * uytot(n,m)/sqrt(vtots(j)) * tc_ftot
            tc_fz = -1.0 * uztot(n,m)/sqrt(vtots(j)) * tc_ftot

            ! force and mass source on a panel! Not density
            bfx(n1,m) = bfx(n1,m) + (tc_fx)
            bfr(n1,m) = bfr(n1,m) + (tc_fy*cos(tc_tt)+tc_fz*sin(tc_tt))
            bft(n1,m) = bft(n1,m) + (tc_fz*cos(tc_tt)-tc_fy*sin(tc_tt))
          end do
        end do

        ! Debug purpose: plot total axial forces
*        totfx = 0.0
*        do m = 1,mr
*          do n = 1,nh
*            j = indexb(n,m)
*            totfx = totfx + bfx(n,m)
*          end do
*        end do
*        write(*,*) 'Total PROPCAV Fx =', totfx

        return
        end subroutine

! Write body force
        subroutine UNS_WRT_BF
          use M_LNK_LST_P2M
          use M_UNS_F2M2P_CORR
          use M_UNS_THETA
          implicit none
          type(COR_ENTRY),pointer :: cor_p2m, temp
          integer :: ixp, irp, ixm, irm
          integer :: ix, ir, ib, it, it1
          double precision,dimension(mesh_nx, mesh_nr, pc_nbld) ::
     &                                          mbfx,mbfr,mbft,msrc
          double precision :: weight, tt, bfx, bfy, bfz, bfr, bft, src

          ! Write it info
          open(285, file='bfout.dat', status='unknown')
          do ib = 1, pc_nbld
            write(285,*) it_rans(ib) - 1   ! -1 means from Fortran index to C index  ! Future time step location
          end do

          mbfx = 0.0; mbfr = 0.0; mbft = 0.0; msrc = 0.0

          cor_p2m => getCorr()
          do while(associated(cor_p2m))
            ixp = cor_p2m%ixp
            irp = cor_p2m%irp
            ixm = cor_p2m%ixm
            irm = cor_p2m%irm
            weight = cor_p2m%weight
            temp => cor_p2m%next
            cor_p2m => temp
            do ib = 1, pc_nbld
              it  = it_prop(ib)      ! Current time step
              ! it1 = it - 1           ! Previous time step
              ! if (it1 <= 0) it1 = it1 + pc_ntpr
              bfx = pbfx(ixp,irp,it) ! *0.5 + pbfx(ixp,irp,it1)*0.5 ! 2nd order accuracy if add extra term
              bfr = pbfr(ixp,irp,it) ! *0.5 + pbfr(ixp,irp,it1)*0.5
              bft = pbft(ixp,irp,it) ! *0.5 + pbft(ixp,irp,it1)*0.5
              src = psrc(ixp,irp,it) ! *0.5 + psrc(ixp,irp,it1)*0.5
              mbfx(ixm,irm,ib) = mbfx(ixm,irm,ib) + bfx * weight
              mbfr(ixm,irm,ib) = mbfr(ixm,irm,ib) + bfr * weight
              mbft(ixm,irm,ib) = mbft(ixm,irm,ib) + bft * weight
              msrc(ixm,irm,ib) = msrc(ixm,irm,ib) + src * weight
            end do
          end do

          ! Write out
          do ib = 1, pc_nbld
          do ir = 1, mesh_nr
          do ix = 1, mesh_nx
              it = it_rans(ib) - (ix - 1)
              if (it < 1) it = it + mesh_nt
              tt = atan2( sum(zNode(ix:ix+1, ir:ir+1, it)),
     &                    sum(yNode(ix:ix+1, ir:ir+1, it))  )  ! Nodes with the same index are at t+1/2 step
              tt = tt - rot_rans
              bfy = mbfr(ix,ir,ib)*cos(tt) - mbft(ix,ir,ib)*sin(tt)
              bfz = mbfr(ix,ir,ib)*sin(tt) + mbft(ix,ir,ib)*cos(tt)
              write(285,*) mbfx(ix,ir,ib), bfy, bfz, msrc(ix,ir,ib)
          end do
          end do
          end do
          close(285)

          ! Debugging plot Db3
          call test_plot

        contains

          subroutine test_plot
            use UNSPROP, only : NTSTEP_KEY
            character*20 uns_fn
            integer ix, ir, it
            double precision :: totbfx,xCent,yCent,zCent,bfx,bfr,bft,src
            write(uns_fn, "(A8,I4.4,A4)") 'Test_bf_',NTSTEP_KEY,'.plt'
            open(285, file = uns_fn, status = 'unknown')
            write(285, *) 'Variables = x,y,z,bfx,bfr,bft,src'
            write(285, *) 'zone T = "blade1",I=',mesh_nx,',J=',mesh_nr
            totbfx= 0.0
            do ir = 1, mesh_nr
            do ix = 1, mesh_nx
              it = mesh_nt - (ix - 1)
              xCent = sum(xNode(ix:ix+1, ir:ir+1, it:it+1))/8.0
              yCent = sum(yNode(ix:ix+1, ir:ir+1, it:it+1))/8.0
              zCent = sum(zNode(ix:ix+1, ir:ir+1, it:it+1))/8.0
              bfx = mbfx(ix,ir,1)
              bfr = mbfr(ix,ir,1)
              bft = mbft(ix,ir,1)
              src = msrc(ix,ir,1)
              totbfx = totbfx + meshVol(ix,ir) * bfx
              write(285, *) xCent, yCent, zCent, bfx, bfr, bft, src
            end do
            end do
            close(285)

            ! write(*,*) 'Integrated body force of the 1st blade:',totbfx

          end subroutine

        end subroutine


      end module
