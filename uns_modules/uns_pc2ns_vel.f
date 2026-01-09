! ----------------------------------------------------------------------------- !
!                       CREATED BY YIRAN SU 10/04/2017                          !
!                    Store and calculation of effective wake!                   !
! ----------------------------------------------------------------------------- !

      module M_UNS_PC2NS_VEL
        use M_UNS_GLOBAL_PAR
        use PC2NSVEL

        public :: UNS_VEL_INIT, UNS_UE_CALC, UNS_UE_EXTRA, UNS_UI_CALC,
     &             UNS_UT_CALC, WRT_VEL_FIELD
        private:: wrt_2d_vector, indv_inner1, indv_inner2

        ! Induced velocity and total velocity at the normal offset locations
        double precision,allocatable,dimension(:,:,:,:):: tUx,tUr,tUt ! NC*MR*NTPREV*NSEED
        double precision,allocatable,dimension(:,:,:,:):: iUx,iUr,iUt ! NC*MR*NTPREV*NSEED


      contains

! Initialized
        subroutine UNS_VEL_INIT
          implicit none

          ! Final data
          allocate(uefx(pc_nc,pc_mr,pc_ntpr))
          allocate(uefr(pc_nc,pc_mr,pc_ntpr))
          allocate(ueft(pc_nc,pc_mr,pc_ntpr))
          uefx = 1.0d0
          uefr = 0.0d0
          ueft = 0.0d0
          ! This module init
          allocate(tUx(pc_nc,pc_mr,pc_ntpr,nSeed))
          allocate(tUr(pc_nc,pc_mr,pc_ntpr,nSeed))
          allocate(tUt(pc_nc,pc_mr,pc_ntpr,nSeed))
          allocate(iUx(pc_nc,pc_mr,pc_ntpr,nSeed))
          allocate(iUr(pc_nc,pc_mr,pc_ntpr,nSeed))
          allocate(iUt(pc_nc,pc_mr,pc_ntpr,nSeed))
          tUx = 1.0
          tUr = 0.0
          tUt = 0.0
          iUx = 0.0
          iUr = 0.0
          iUt = 0.0
        end subroutine

! Calculate effective wake
        subroutine UNS_UE_CALC
          use M_UNS_THETA, only : it_prop
          use M_UNS_F2M2P_CORR, only : mesh2pc
          implicit none
          integer :: ib, it, ix
          double precision :: af, bt

          do ib = 1, pc_nbld
            it = it_prop(ib)

            uefx(:,:,it) = uefx(:,:,it) * (1.0 - ueReF) +
     &            (tUx(:,:,it,useLayer) - iUx(:,:,it,useLayer)) * ueReF
            uefr(:,:,it) = uefr(:,:,it) * (1.0 - ueReF) +
     &            (tUr(:,:,it,useLayer) - iUr(:,:,it,useLayer)) * ueReF
            ueft(:,:,it) = ueft(:,:,it) * (1.0 - ueReF) +
     &            (tUt(:,:,it,useLayer) - iUt(:,:,it,useLayer)) * ueReF

            ! Not relaxed and faster version
            ! uefx(:,:,it) = tUx(:,:,it,useLayer) - iUx(:,:,it,useLayer)
            ! uefr(:,:,it) = tUr(:,:,it,useLayer) - iUr(:,:,it,useLayer)
            ! ueft(:,:,it) = tUt(:,:,it,useLayer) - iUt(:,:,it,useLayer)
          end do

          ! Extrapolate the hub station because of the hub panel saw tooth effect
          do ib = 1, pc_nbld
          it = it_prop(ib)
          do ix = 1, pc_nc
            af = (mesh2pc(ix,3,useLayer)%R - mesh2pc(ix,1,useLayer)%R) /
     &           (mesh2pc(ix,3,useLayer)%R - mesh2pc(ix,2,useLayer)%R)
            bt = af - 1.0
! Version 1
!            uefx(ix,1,it) = uefx(ix,2,it)*af - uefx(ix,3,it)*bt
!            uefr(ix,1,it) = uefr(ix,2,it)*af - uefr(ix,3,it)*bt
!            ueft(ix,1,it) = ueft(ix,2,it)*af - ueft(ix,3,it)*bt
! Version 2 D3P
            uefx(ix,1,it) = uefx(ix,3,it)
            uefr(ix,1,it) = uefr(ix,3,it)
            ueft(ix,1,it) = ueft(ix,3,it)
            uefx(ix,2,it) = uefx(ix,3,it)
            uefr(ix,2,it) = uefr(ix,3,it)
            ueft(ix,2,it) = ueft(ix,3,it)
          end do
          end do

        end subroutine

! Extrapolate from the previous two time steps to get current body force field
        subroutine UNS_UE_EXTRA
          use M_UNS_THETA, only : it_prop
          implicit none
          integer :: ib, it, it1, it2
          do ib = 1, pc_nbld
            it = it_prop(ib)
            it1 = it - 1
            it2 = it - 2
            if (it1 < 1) it1 = it1 + pc_ntpr
            if (it2 < 1) it2 = it2 + pc_ntpr
            uefx(:,:,it) = uefx(:,:,it1) * 2.0 - uefx(:,:,it2)
            uefr(:,:,it) = uefr(:,:,it1) * 2.0 - uefr(:,:,it2)
            ueft(:,:,it) = ueft(:,:,it1) * 2.0 - ueft(:,:,it2)
          end do
        end subroutine

! Calculate Ind-Vel inside the fluid domain for all the blades
        subroutine UNS_UI_CALC
          use M_UNS_THETA, only : it_prop
          use M_UNS_F2M2P_CORR, only : mesh2pc
          implicit none
          integer :: ix, ir, it, isd, ib, ct
          double precision,dimension(3,pc_nc*pc_mr*nSeed)         ::loc  ! X  Y  Z
          double precision,dimension(3,pc_nc*pc_mr*nSeed,pc_nbld) ::indv ! Ux Ur Ut
          double precision :: af, bt

          ! Gather all points
          ct = 0
          do isd = 1, nSeed
            do ir = 1, pc_mr
              do ix = 1, pc_nc
                ct = ct + 1
                loc(1,ct) = mesh2pc(ix,ir,isd)%X
                loc(2,ct) = mesh2pc(ix,ir,isd)%Y
                loc(3,ct) = mesh2pc(ix,ir,isd)%Z
              end do
            end do
          end do

          ! Calculate
          indv = 0.0d0
          call indv_inner2(loc, indv, pc_nc*pc_mr*nSeed)
          call indv_inner1(loc, indv, pc_nc*pc_mr*nSeed)
          ! indv = 0.1

          ! call test_Bstrength
          ! call test_Hstrength
          ! call test_Wstrength

          ! Write back
          do ib = 1, pc_nbld
            it = it_prop(ib)
            ct = 0
            do isd = 1, nSeed
              do ir = 1, pc_mr
                do ix = 1, pc_nc
                  ct = ct + 1
                  iUx(ix,ir,it,isd) = indv(1,ct,ib)
                  iUr(ix,ir,it,isd) = indv(2,ct,ib)
                  iUt(ix,ir,it,isd) = indv(3,ct,ib)
                end do
              end do
            end do
          end do

          ! Extrapolate the hub station because of the hub panel saw tooth effect
          do ib = 1, pc_nbld
          it = it_prop(ib)
          do isd = 1, nSeed
          do ix = 1, pc_nc
            af = (mesh2pc(ix,3,isd)%R - mesh2pc(ix,1,isd)%R) /
     &           (mesh2pc(ix,3,isd)%R - mesh2pc(ix,2,isd)%R)
            bt = af - 1.0
            iUx(ix,1,it,isd) = iUx(ix,2,it,isd)*af - iUx(ix,3,it,isd)*bt
            iUr(ix,1,it,isd) = iUx(ix,2,it,isd)*af - iUx(ix,3,it,isd)*bt
            iUt(ix,1,it,isd) = iUx(ix,2,it,isd)*af - iUx(ix,3,it,isd)*bt
          end do
          end do
          end do

        end subroutine

! Read Total Flow for all the blade
        subroutine  UNS_UT_CALC
          use M_UNS_THETA, only : it_prop, it_rans, rot_rans, isRotate
          use M_UNS_F2M2P_CORR, only : yNode, zNode, mesh2pc
          implicit none
          integer :: nlayer1, nlayer2, ntotal
          integer :: ib, ilay, ix, ir, it, it1, it2, px, pr, isd
          double precision :: uxtemp, uytemp, uztemp, urtemp, uttemp, tt
          double precision, allocatable, dimension(:,:,:) :: uxx,urr,utt
          double precision, dimension(3)::v1, v2, v3, v4, vt1, vt2, vtar
          double precision :: xper, rper, tper

          open(285, file = 'total_vel.dat', status = 'old')
          read(285, *) nlayer1, nlayer2, ntotal

          allocate(uxx(mesh_nx, mesh_nr, mesh_nt),
     &             urr(mesh_nx, mesh_nr, mesh_nt),
     &             utt(mesh_nx, mesh_nr, mesh_nt) )
          uxx = -999999.0
          urr = -999999.0
          utt = -999999.0

          do ib = 1, pc_nbld

            ! Read U_total on all layers of a blades
            ! Convert Uy and Uz to Ur and Ut
            ! Store them on a mesh-type index matrix
            ! Store everything on the key-blade's location
            do ilay = nlayer1, nlayer2
            do ir = 1, mesh_nr
            do ix = 1, mesh_nx
              read(285, *) uxtemp, uytemp, uztemp

                                  ! Adjust for different Ix     - (ix-1)
                                  ! Adjust for different layer  + ilay
                                  ! Adjust for rotation         - rot_rans
              it = it_rans(ib) - (ix-1) + ilay
              if (it < 1) it = it + mesh_nt
              if (it > mesh_nt) it = it - mesh_nt
              tt = atan2( sum(zNode(ix:ix+1, ir:ir+1, it:it+1)),
     &                    sum(yNode(ix:ix+1, ir:ir+1, it:it+1))  )
              tt = tt - rot_rans
              urtemp = uytemp * cos(tt) + uztemp * sin(tt)
              uttemp = uztemp * cos(tt) - uytemp * sin(tt)

              !! All interpolation points are defined around the key blade at t = 0
              !! Therefore, we put the n-th blade's total velocity to it_rans(1)
              it = it_rans(1) - (ix-1) + ilay
              if (it < 1) it = it + mesh_nt
              if (it > mesh_nt) it = it - mesh_nt
              uxx(ix,ir,it) = uxtemp
              urr(ix,ir,it) = urtemp
              utt(ix,ir,it) = uttemp
            end do
            end do
            end do

            ! Interpolate and set uttx/r/t
            do isd = 1,nSeed
            do pr = 1,pc_mr
            do px = 1,pc_nc
                ix   = mesh2pc(px, pr, isd)%ix
                ir   = mesh2pc(px, pr, isd)%ir
                it   = mesh2pc(px, pr, isd)%it
                xper = mesh2pc(px, pr, isd)%xper
                rper = mesh2pc(px, pr, isd)%rper
                tper = mesh2pc(px, pr, isd)%tper

                v1(1) = uxx(ix  , ir  , it);
                v1(2) = urr(ix  , ir  , it);
                v1(3) = utt(ix  , ir  , it);
                v2(1) = uxx(ix+1, ir  , it);
                v2(2) = urr(ix+1, ir  , it);
                v2(3) = utt(ix+1, ir  , it);
                v3(1) = uxx(ix  , ir+1, it);
                v3(2) = urr(ix  , ir+1, it);
                v3(3) = utt(ix  , ir+1, it);
                v4(1) = uxx(ix+1, ir+1, it);
                v4(2) = urr(ix+1, ir+1, it);
                v4(3) = utt(ix+1, ir+1, it);
                vt1   = v4*xper*rper + v1*(1.0-xper)*(1.0-rper)
     &                + v3*(1.0-xper)*rper + v2*xper*(1.0-rper)

                it = it + 1
                if (it > mesh_nt) it = it - mesh_nt
                v1(1) = uxx(ix  , ir  , it);
                v1(2) = urr(ix  , ir  , it);
                v1(3) = utt(ix  , ir  , it);
                v2(1) = uxx(ix+1, ir  , it);
                v2(2) = urr(ix+1, ir  , it);
                v2(3) = utt(ix+1, ir  , it);
                v3(1) = uxx(ix  , ir+1, it);
                v3(2) = urr(ix  , ir+1, it);
                v3(3) = utt(ix  , ir+1, it);
                v4(1) = uxx(ix+1, ir+1, it);
                v4(2) = urr(ix+1, ir+1, it);
                v4(3) = utt(ix+1, ir+1, it);
                vt2   = v4*xper*rper + v1*(1.0-xper)*(1.0-rper)
     &                + v3*(1.0-xper)*rper + v2*xper*(1.0-rper)

                vtar = vt1*(1.0-tper) + vt2*tper

                ! Set values
                it1 = it_prop(ib)
                tUx(px,pr,it1,isd) = vtar(1)
                tUr(px,pr,it1,isd) = vtar(2)
                tUt(px,pr,it1,isd) = vtar(3)
            end do
            end do
            end do
          end do
          close(285)
          deallocate(uxx, urr, utt)
        end subroutine

! Write two files
        subroutine wrt_2d_vector
        use UNSPROP, only : NTSTEP_KEY
        use M_UNS_F2M2P_CORR, only : mesh2pc
        use M_UNS_THETA, only : it_prop
        implicit none

        integer :: ix, ir, isd, it
        double precision :: xx, yy, zz, rr, tt, ux, uy, uz, ur, ut
        character(len = 100) :: name

        write(name, "(A11,I4.4,A4)") 'Test_V_vec_',NTSTEP_KEY,'.plt'
        open(287,file = name,status = 'unknown')
        write(287,*) 'Variables = X, Y, Z, Ux, Uy, Uz, iLay'
        ir = pc_mr/2
        it = it_prop(1)
        ! Zone 1: Foil surface
        write(287,*) 'Zone T = "1_Cylindrical_Cut"'
        do ix = 1, pc_nc
          write(287,*) mesh2pc(ix,ir,1)%X, mesh2pc(ix,ir,1)%Y, mesh2pc(ix,ir,1)%Z
          write(287,*) 0.0, 0.0, 0.0, 0
        end do
        ! Zone 2: Foil surface
        do ix = 1, pc_nc
          write(287,*) 'Zone T = "2_Utotal_RANS"'
          do isd = 1, nSeed
            xx = mesh2pc(ix,ir,isd)%X
            yy = mesh2pc(ix,ir,isd)%Y
            zz = mesh2pc(ix,ir,isd)%Z
            rr = sqrt(yy**2 + zz**2)
            tt = atan2(zz, yy)
            ux = tUx(ix,ir,it,isd)
            ur = tUr(ix,ir,it,isd)
            ut = tUt(ix,ir,it,isd)
            uy = ur*cos(tt) - ut*sin(tt)
            uz = ur*sin(tt) + ut*cos(tt)
            write(287,*) xx, yy, zz, ux, uy, uz, isd
          end do
        end do
        ! Zone 3: Foil surface
        do ix = 1, pc_nc
          write(287,*) 'Zone T = "3_Uind_BEM"'
          do isd = 1, nSeed
            xx = mesh2pc(ix,ir,isd)%X
            yy = mesh2pc(ix,ir,isd)%Y
            zz = mesh2pc(ix,ir,isd)%Z
            rr = sqrt(yy**2 + zz**2)
            tt = atan2(zz, yy)
            ux = iUx(ix,ir,it,isd)
            ur = iUr(ix,ir,it,isd)
            ut = iUt(ix,ir,it,isd)
            uy = ur*cos(tt) - ut*sin(tt)
            uz = ur*sin(tt) + ut*cos(tt)
            write(287,*) xx, yy, zz, ux, uy, uz, isd
          end do
        end do
        close(287)
        end subroutine

! Write two files
        subroutine WRT_VEL_FIELD
        use UNSPROP, only : NTSTEP_KEY
        use M_UNS_F2M2P_CORR, only : mesh2pc
        use M_UNS_THETA, only : it_prop
        implicit none
        integer :: ix, ir, isd, it
        double precision :: xx, yy, zz, rr, tt, ux, uy, uz, ur, ut
        character(len = 100) :: name
        ! Vector plot
        ! call wrt_2d_vector

        it = it_prop(1)
        write(name, "(A10,I4.4,A4)") 'Test_V_3d_',NTSTEP_KEY,'.plt'
        open(287,file = name,status = 'unknown')
        write(287,*) 'Variables = X,Y,Z,UTx,UTy,UTz,UIx,UIy,UIz'
        write(287,*) 'Zone T="3d vel",I=',pc_nc,',J=',pc_mr,',K=',nSeed

        do isd = 1, nSeed
          do ir = 1, pc_mr
            do ix = 1, pc_nc
              xx = mesh2pc(ix,ir,isd)%X
              yy = mesh2pc(ix,ir,isd)%Y
              zz = mesh2pc(ix,ir,isd)%Z
              rr = sqrt(yy**2 + zz**2)
              tt = atan2(zz, yy)
              ux = tUx(ix,ir,it,isd)
              ur = tUr(ix,ir,it,isd)
              ut = tUt(ix,ir,it,isd)
              uy = ur*cos(tt) - ut*sin(tt)
              uz = ur*sin(tt) + ut*cos(tt)
              write(287,*) xx, yy, zz, ux, uy, uz
              ux = iUx(ix,ir,it,isd)
              ur = iUr(ix,ir,it,isd)
              ut = iUt(ix,ir,it,isd)
              uy = ur*cos(tt) - ut*sin(tt)
              uz = ur*sin(tt) + ut*cos(tt)
              write(287,*) ux, uy, uz
            end do
          end do
        end do
        close(287)
        end subroutine

      subroutine test_Bstrength
      use M_UNS_THETA, only : it_prop
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'

      integer :: i,j,k,ib,ib1
      character*20 uns_fn
      integer ix, ir, it

      real :: tb_pot(nblade, npanel),tb_dpdn(nblade, npanel)
      real :: tb_temp5(nblade, nwmin*mr)

! - Read source/dipole strength
      nread = nwmin * mr
      do ib = 1, nblade
        irec = it_prop(ib)
        call read2(45, irec, pot,   npanel)
        call read2(47, irec, dpdn,  npanel)
        call read2(46, irec, temp5, nread )
        tb_pot  (ib, 1:npanel) =   pot(1:npanel)
        tb_dpdn (ib, 1:npanel) =  dpdn(1:npanel)
        tb_temp5(ib, 1:nread ) = temp5(1:nread )
      end do

! - Write output
      ! double precision :: totbfx,xCent,yCent,zCent,bfx,bfr,bft,src
      write(uns_fn, "(A12,I4.4,A4)") 'Test_Bstren_',NTSTEP_KEY,'.plt'
      open(299, file = uns_fn, status = 'unknown')
      write(299, *) 'variables = x,y,z,pot,dpdn'
      write(299,*)'Zone T="1",i=',ncp,',j=',mrp,',datapacking=block,varlocation=(4=cellcentered,5=cellcentered)'
      write(299,*) ((xb(n,m), n=1,ncp),m=1,mrp)
      write(299,*) ((yb(n,m), n=1,ncp),m=1,mrp)
      write(299,*) ((zb(n,m), n=1,ncp),m=1,mrp)
      do m = 1,mr
        do n = 1,nc
          l1 = indexb(n,m)
          write(299,*) tb_pot(1,l1)
        end do
      end do
      do m = 1,mr
        do n = 1,nc
          l1 = indexb(n,m)
          write(299,*) tb_dpdn(1,l1)
        end do
      end do
      write(299,*)'zone t="2",i=',ncp,',j=',mrp,',datapacking=block,varlocation=(4=cellcentered,5=cellcentered)'
      write(299,*) ((xb(n,m), n=1,ncp),m=1,mrp)
      write(299,*) ((yb(n,m), n=1,ncp),m=1,mrp)
      write(299,*) ((zb(n,m), n=1,ncp),m=1,mrp)
      do m = 1,mr
        do n = 1,nc
          l1 = indexb(n,m)
          write(299,*) tb_pot(2,l1)
        end do
      end do
      do m = 1,mr
        do n = 1,nc
          l1 = indexb(n,m)
          write(299,*) tb_dpdn(2,l1)
        end do
      end do
      write(299,*)'zone t="3",i=',ncp,',j=',mrp,',datapacking=block,varlocation=(4=cellcentered,5=cellcentered)'
      write(299,*) ((xb(n,m), n=1,ncp),m=1,mrp)
      write(299,*) ((yb(n,m), n=1,ncp),m=1,mrp)
      write(299,*) ((zb(n,m), n=1,ncp),m=1,mrp)
      do m = 1,mr
        do n = 1,nc
          l1 = indexb(n,m)
          write(299,*) tb_pot(3,l1)
        end do
      end do
      do m = 1,mr
        do n = 1,nc
          l1 = indexb(n,m)
          write(299,*) tb_dpdn(3,l1)
        end do
      end do
      write(299,*)'zone t="4",i=',ncp,',j=',mrp,',datapacking=block,varlocation=(4=cellcentered,5=cellcentered)'
      write(299,*) ((xb(n,m), n=1,ncp),m=1,mrp)
      write(299,*) ((yb(n,m), n=1,ncp),m=1,mrp)
      write(299,*) ((zb(n,m), n=1,ncp),m=1,mrp)
      do m = 1,mr
        do n = 1,nc
          l1 = indexb(n,m)
          write(299,*) tb_pot(4,l1)
        end do
      end do
      do m = 1,mr
        do n = 1,nc
          l1 = indexb(n,m)
          write(299,*) tb_dpdn(4,l1)
        end do
      end do
      close(299)
      end subroutine
      subroutine test_Hstrength
      use M_UNS_THETA, only : it_prop
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'

      integer :: i,j,k,ib,ib1
      character*20 uns_fn
      integer ix, ir, it

      real :: tb_pot(nblade, npanel),tb_dpdn(nblade, npanel)
      real :: tb_temp5(nblade, nwmin*mr)

! - Read source/dipole strength
      nread = nwmin * mr
      do ib = 1, nblade
        irec = it_prop(ib)
        call read2(45, irec, pot,   npanel)
        call read2(47, irec, dpdn,  npanel)
        call read2(46, irec, temp5, nread )
        tb_pot  (ib, 1:npanel) =   pot(1:npanel)
        tb_dpdn (ib, 1:npanel) =  dpdn(1:npanel)
        tb_temp5(ib, 1:nread ) = temp5(1:nread )
      end do

! - Write output
      ! double precision :: totbfx,xCent,yCent,zCent,bfx,bfr,bft,src
      write(uns_fn, "(A12,I4.4,A4)") 'Test_Hstren_',NTSTEP_KEY,'.plt'
      open(299, file = uns_fn, status = 'unknown')
      write(299, *) 'variables = x,y,z,pot,dpdn'
      write(299,*)'Zone T="1",i=',nhbx,',j=',mhbt+1,',datapacking=block,varlocation=(4=cellcentered,5=cellcentered)'
      write(299,*) ((xh(n,m), n=1,nhbx),m=1,mhbt+1)
      write(299,*) ((yh(n,m), n=1,nhbx),m=1,mhbt+1)
      write(299,*) ((zh(n,m), n=1,nhbx),m=1,mhbt+1)
      do m = 1,mhbt
        do n = 1,nhbx-1
          l1 = indexh(n,m)
          write(299,*) tb_pot(1,l1)
        end do
      end do
      do m = 1,mhbt
        do n = 1,nhbx-1
          l1 = indexh(n,m)
          write(299,*) tb_dpdn(1,l1)
        end do
      end do
      write(299,*)'Zone T="2",i=',nhbx,',j=',mhbt+1,',datapacking=block,varlocation=(4=cellcentered,5=cellcentered)'
      write(299,*) ((xh(n,m), n=1,nhbx),m=1,mhbt+1)
      write(299,*) ((yh(n,m), n=1,nhbx),m=1,mhbt+1)
      write(299,*) ((zh(n,m), n=1,nhbx),m=1,mhbt+1)
      do m = 1,mhbt
        do n = 1,nhbx-1
          l1 = indexh(n,m)
          write(299,*) tb_pot(2,l1)
        end do
      end do
      do m = 1,mhbt
        do n = 1,nhbx-1
          l1 = indexh(n,m)
          write(299,*) tb_dpdn(2,l1)
        end do
      end do
      write(299,*)'Zone T="3",i=',nhbx,',j=',mhbt+1,',datapacking=block,varlocation=(4=cellcentered,5=cellcentered)'
      write(299,*) ((xh(n,m), n=1,nhbx),m=1,mhbt+1)
      write(299,*) ((yh(n,m), n=1,nhbx),m=1,mhbt+1)
      write(299,*) ((zh(n,m), n=1,nhbx),m=1,mhbt+1)
      do m = 1,mhbt
        do n = 1,nhbx-1
          l1 = indexh(n,m)
          write(299,*) tb_pot(3,l1)
        end do
      end do
      do m = 1,mhbt
        do n = 1,nhbx-1
          l1 = indexh(n,m)
          write(299,*) tb_dpdn(3,l1)
        end do
      end do
      write(299,*)'Zone T="4",i=',nhbx,',j=',mhbt+1,',datapacking=block,varlocation=(4=cellcentered,5=cellcentered)'
      write(299,*) ((xh(n,m), n=1,nhbx),m=1,mhbt+1)
      write(299,*) ((yh(n,m), n=1,nhbx),m=1,mhbt+1)
      write(299,*) ((zh(n,m), n=1,nhbx),m=1,mhbt+1)
      do m = 1,mhbt
        do n = 1,nhbx-1
          l1 = indexh(n,m)
          write(299,*) tb_pot(4,l1)
        end do
      end do
      do m = 1,mhbt
        do n = 1,nhbx-1
          l1 = indexh(n,m)
          write(299,*) tb_dpdn(4,l1)
        end do
      end do
      close(299)
      end subroutine
      subroutine test_Wstrength
      use M_UNS_THETA, only : it_prop
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'

      integer :: i,j,k,ib,ib1
      character*20 uns_fn
      integer ix, ir, it

      real :: tb_temp5(nblade, nwmin*mr)

! - Read source/dipole strength
      nread = nwmin * mr
      do ib = 1, nblade
        irec = it_prop(ib)
        call read2(46, irec, temp5, nread )
        tb_temp5(ib, 1:nread ) = temp5(1:nread )
      end do

! - Write output
      ! double precision :: totbfx,xCent,yCent,zCent,bfx,bfr,bft,src
      write(uns_fn, "(A12,I4.4,A4)") 'Test_Wstren_',NTSTEP_KEY,'.plt'
      open(299, file = uns_fn, status = 'unknown')
      write(299, *) 'variables = x,y,z,delp'
      write(299,*)'Zone T="1",i=',nwmin+1,',j=',mrp,',datapacking=block,varlocation=(4=cellcentered)'
      write(299,*) ((xw(n,m), n=1,nwmin+1),m=1,mrp)
      write(299,*) ((yw(n,m), n=1,nwmin+1),m=1,mrp)
      write(299,*) ((zw(n,m), n=1,nwmin+1),m=1,mrp)
      do m = 1,mr
        do n = 1,nwmin
          l1 = indexw2(n,m)
          write(299,*) tb_temp5(1,l1)
        end do
      end do

      write(299,*)'Zone T="2",i=',nwmin+1,',j=',mrp,',datapacking=block,varlocation=(4=cellcentered)'
      write(299,*) ((xw(n,m), n=1,nwmin+1),m=1,mrp)
      write(299,*) ((yw(n,m), n=1,nwmin+1),m=1,mrp)
      write(299,*) ((zw(n,m), n=1,nwmin+1),m=1,mrp)
      do m = 1,mr
        do n = 1,nwmin
          l1 = indexw2(n,m)
          write(299,*) tb_temp5(2,l1)
        end do
      end do

      write(299,*)'Zone T="3",i=',nwmin+1,',j=',mrp,',datapacking=block,varlocation=(4=cellcentered)'
      write(299,*) ((xw(n,m), n=1,nwmin+1),m=1,mrp)
      write(299,*) ((yw(n,m), n=1,nwmin+1),m=1,mrp)
      write(299,*) ((zw(n,m), n=1,nwmin+1),m=1,mrp)
      do m = 1,mr
        do n = 1,nwmin
          l1 = indexw2(n,m)
          write(299,*) tb_temp5(3,l1)
        end do
      end do

      write(299,*)'Zone T="4",i=',nwmin+1,',j=',mrp,',datapacking=block,varlocation=(4=cellcentered)'
      write(299,*) ((xw(n,m), n=1,nwmin+1),m=1,mrp)
      write(299,*) ((yw(n,m), n=1,nwmin+1),m=1,mrp)
      write(299,*) ((zw(n,m), n=1,nwmin+1),m=1,mrp)
      do m = 1,mr
        do n = 1,nwmin
          l1 = indexw2(n,m)
          write(299,*) tb_temp5(4,l1)
        end do
      end do
      close(299)
      end subroutine

      subroutine indv_inner1(p_loc_d, v_ind, n_pt)
! ===========================================================
! This subroutine calculates the induced velocities.
! Currently working on!!
! ===========================================================
      use omp_lib
      use M_UNS_THETA, only : it_prop
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'

      integer :: n_pt
      double precision, dimension(3, n_pt) :: p_loc_d
      real, dimension(3, n_pt) :: p_loc
      double precision, dimension(3, n_pt, nblade) :: v_ind, v_tmp

      real :: ssx,ssy,ssz,ssr,sst,st0,p_chl
      real :: pv_xv(4),pv_yv(4),pv_side(4),pv_s(15)
      integer :: i,j,k,ib,ib1

      real :: tb_pot(nblade, npanel),tb_dpdn(nblade, npanel)
      real :: tb_temp5(nblade, nwmin*mr)

      integer :: om_ct1, om_ct2, om_ctr
      real :: om_dt
      call system_clock(om_ct1,om_ctr)

! - Induced velocity set to zero
      p_loc  = sngl(p_loc_d)
      v_tmp  = 0.0d0

! - Read source/dipole strength
      nread = nwmin * mr
      do ib = 1, nblade
        irec = it_prop(ib)
        call read2(45, irec, pot,   npanel)
        call read2(47, irec, dpdn,  npanel)
        call read2(46, irec, temp5, nread )
        tb_pot  (ib, 1:npanel) =   pot(1:npanel)
        tb_dpdn (ib, 1:npanel) =  dpdn(1:npanel)
        tb_temp5(ib, 1:nread ) = temp5(1:nread )
      end do

! =============================================================================
!                    begin wall panel induced velocity!
! =============================================================================
!      write(*,*) 'begin calculating blade/hub/casing induced velocity! '

!$omp parallel private(i,j,k,ib,ib1,
!$omp&                 pv_xv,pv_yv,pv_side,pv_s,ssx,ssy,ssz,ssr,sst,st0,
!$omp&                 xloc,yloc,zloc,fs,fd,fsx,fsy,fdx,fdy,fdz,p_chl,
!$omp&                 vdx,vdy,vdz,vsx,vsy,vsz,vx0,vy0,vz0,vr0,vt0)
!$omp&         shared (npanb,npanel,n_pt,nblade,pi,p_loc,
!$omp&                 xvp,yvp,sid,ss,dir,xct,chrleps,
!$omp&                 tb_pot,tb_dpdn)
!$omp& reduction (+:v_tmp)
!$omp do collapse(2) schedule(dynamic,500)
      do i = (npanb+1), npanel
      do j = 1, n_pt
        pv_xv(1:4)   = xvp(i, 1:4)
        pv_yv(1:4)   = yvp(i, 1:4)
        pv_side(1:4) = sid(i, 1:4)
        pv_s(1:15)   = ss(i, 1:15)
        p_chl = chrleps(i)
        ssx = p_loc(1, j)
        ssr = sqrt (p_loc(2,j)**2 + p_loc(3,j)**2)
        st0 = atan2(p_loc(3,j), p_loc(2,j))

        do ib = 1,nblade
          sst = st0 + 2.0 * pi / real(nblade) * real(ib - 1)
          ssy = ssr * cos(sst)
          ssz = ssr * sin(sst)

          xloc = (ssx-xct(i,1))*dir(i,1,1) + (ssy-xct(i,2))*dir(i,1,2)
     &                                     + (ssz-xct(i,3))*dir(i,1,3)
          yloc = (ssx-xct(i,1))*dir(i,2,1) + (ssy-xct(i,2))*dir(i,2,2)
     &                                     + (ssz-xct(i,3))*dir(i,2,3)
          zloc = (ssx-xct(i,1))*dir(i,3,1) + (ssy-xct(i,2))*dir(i,3,2)
     &                                     + (ssz-xct(i,3))*dir(i,3,3)

          call rpan_pv(xloc,yloc,zloc,p_chl,fs,fd,fsx,fsy,
     &                 fdx,fdy,fdz,1,1,pv_xv,pv_yv,pv_s,pv_side)

          vdx = fdx * dir(i,1,1) + fdy * dir(i,2,1) + fdz * dir(i,3,1)
          vdy = fdx * dir(i,1,2) + fdy * dir(i,2,2) + fdz * dir(i,3,2)
          vdz = fdx * dir(i,1,3) + fdy * dir(i,2,3) + fdz * dir(i,3,3)
          vsx = fsx * dir(i,1,1) + fsy * dir(i,2,1) - fd  * dir(i,3,1)
          vsy = fsx * dir(i,1,2) + fsy * dir(i,2,2) - fd  * dir(i,3,2)
          vsz = fsx * dir(i,1,3) + fsy * dir(i,2,3) - fd  * dir(i,3,3)

          do ib1 = 1, nblade
            k = ib + ib1 - 1
            if (k .gt. nblade) k = k - nblade
            vx0 = -1.0* (tb_pot(k,i)* vdx - tb_dpdn(k,i)* vsx) /4.0 /pi
            vy0 = -1.0* (tb_pot(k,i)* vdy - tb_dpdn(k,i)* vsy) /4.0 /pi
            vz0 = -1.0* (tb_pot(k,i)* vdz - tb_dpdn(k,i)* vsz) /4.0 /pi
            vr0 = vz0 * sin(sst) + vy0 * cos(sst)
            vt0 = vz0 * cos(sst) - vy0 * sin(sst)

            v_tmp(1,j,ib1) = v_tmp(1,j,ib1) + dble(vx0)
            v_tmp(2,j,ib1) = v_tmp(2,j,ib1) + dble(vr0)
            v_tmp(3,j,ib1) = v_tmp(3,j,ib1) + dble(vt0)
          end do

        end do
      end do
      end do
!$omp end do
!$omp end parallel
      v_ind  = v_ind + v_tmp
      v_tmp = 0.0d0

! =============================================================================
!                    begin wake induced velocity!
! =============================================================================
!      write(*,*) 'begin calculating wake induced velocity! '

!$omp parallel private(j,m,n,ib,i,i2,ib1,k,
!$omp&                 pv_xv,pv_yv,pv_s,pv_side,ssx,ssr,sst,st0,ssy,ssz,
!$omp&                 xloc,yloc,zloc,fs,fd,fsx,fsy,fdx,fdy,fdz,p_chl,
!$omp&                 vdx,vdy,vdz,vx0,vy0,vz0,vr0,vt0,tb_t5)
!$omp&          shared(n_pt,mr,nwmin,nblade,nread,pi,p_loc,
!$omp&                 xvpw,yvpw,sidw,ssw,
!$omp&                 dirw,xctw,chrlews,tb_temp5)
!$omp& reduction (+:v_tmp)
!$omp do collapse(3) schedule(dynamic,50)
      do j = 1, n_pt
      do m = 1, mr
      do n = 1, nwmin
        do ib = 1,nblade
          i  = idxwak(n,m)
          i2 = indexw2(n,m)
          pv_xv(1:4)   = xvpw(i, 1:4)
          pv_yv(1:4)   = yvpw(i, 1:4)
          pv_side(1:4) = sidw(i, 1:4)
          pv_s(1:15)   = ssw (i, 1:15)
          p_chl = chrlews(i)
          ssx = p_loc(1, j)
          ssr = sqrt (p_loc(2,j)**2 + p_loc(3,j)**2)
          st0 = atan2(p_loc(3,j), p_loc(2,j))
          sst = st0 + 2.0 * pi / real(nblade) * real(ib - 1)
          ssy = ssr * cos(sst)
          ssz = ssr * sin(sst)

        xloc = (ssx-xctw(i,1))*dirw(i,1,1) + (ssy-xctw(i,2))*dirw(i,1,2)
     *                                     + (ssz-xctw(i,3))*dirw(i,1,3)
        yloc = (ssx-xctw(i,1))*dirw(i,2,1) + (ssy-xctw(i,2))*dirw(i,2,2)
     *                                     + (ssz-xctw(i,3))*dirw(i,2,3)
        zloc = (ssx-xctw(i,1))*dirw(i,3,1) + (ssy-xctw(i,2))*dirw(i,3,2)
     *                                     + (ssz-xctw(i,3))*dirw(i,3,3)

          call rpan_pv(xloc,yloc,zloc,p_chl,fs,fd,fsx,fsy,
     *                 fdx,fdy,fdz,1,1,pv_xv,pv_yv,pv_s,pv_side)

          vdx = fdx*dirw(i,1,1) + fdy*dirw(i,2,1) + fdz*dirw(i,3,1)
          vdy = fdx*dirw(i,1,2) + fdy*dirw(i,2,2) + fdz*dirw(i,3,2)
          vdz = fdx*dirw(i,1,3) + fdy*dirw(i,2,3) + fdz*dirw(i,3,3)

          do ib1 = 1, nblade
            k = ib + ib1 - 1
            if (k .gt. nblade) k = k - nblade
            tb_t5 = tb_temp5(k,i2)

            vx0 = -1.0* (tb_t5 * vdx) /4.0 /pi
            vy0 = -1.0* (tb_t5 * vdy) /4.0 /pi
            vz0 = -1.0* (tb_t5 * vdz) /4.0 /pi
            vr0 = vz0 * sin(sst) + vy0 * cos(sst)
            vt0 = vz0 * cos(sst) - vy0 * sin(sst)

            v_tmp(1,j,ib1) = v_tmp(1,j,ib1) + dble(vx0)
            v_tmp(2,j,ib1) = v_tmp(2,j,ib1) + dble(vr0)
            v_tmp(3,j,ib1) = v_tmp(3,j,ib1) + dble(vt0)
          end do

        end do
      end do
      end do
      end do
!$omp end do
!$omp end parallel

      v_ind  = v_ind + v_tmp

      call system_clock(om_ct2)
      om_dt=(real(om_ct2)-real(om_ct1))/real(om_ctr)
      write(*, "(A1,F4.1)") " ", om_dt


      return
      end subroutine

      subroutine indv_inner2(p_loc, v_ind, n_pt)
      ! Blade perturbation velocity induced by the camber surface
      use omp_lib
      use M_UNS_THETA, only: it_prop
      use GEOCAMB, only: indexCamber,xvpc,yvpc,ssc,sidc,chlenc,xctc,dirc
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'

      integer :: n_pt, n_cb
      double precision, dimension(3, n_pt) :: p_loc
      double precision, dimension(3, n_pt, nblade) :: v_ind, v_tmp

      double precision :: ssx,ssy,ssz,ssr,sst,st0,p_chl,xloc,yloc,zloc
      double precision :: fs,fd,fsx,fsy,fdx,fdy,fdz
      double precision :: pv_xv(4),pv_yv(4),pv_side(4),pv_s(15)
      double precision :: vdx,vdy,vdz,vsx,vsy,vsz,vx0,vy0,vz0,vr0,vt0
      integer :: i,j,k,ib,ib1

      double precision :: tb_pot(nblade, npanel),tb_dpdn(nblade, npanel)
      double precision :: tc_pot(nblade, nh*mr),tc_dpdn(nblade, nh*mr)
      double precision :: ndot1, ndot2

      integer :: om_ct1, om_ct2, om_ctr
      real :: om_dt
      call system_clock(om_ct1,om_ctr)

      v_tmp  = 0.0d0
      n_cb = nh * mr
! - Read source/dipole strength
      do ib = 1, nblade
        irec = it_prop(ib)
        call read2(45, irec, pot,   npanel)
        call read2(47, irec, dpdn,  npanel)
        tb_pot  (ib, 1:npanel) =   pot(1:npanel)
        tb_dpdn (ib, 1:npanel) =  dpdn(1:npanel)
      end do

! - Calculate camber source/dipole strength
      ! Index ix, ir : camber panel index in chord-wise and span-wise direction
      ! Index ix1, ix2 : blade panel chord-wise index (ix1 suction side, ix2 pres side)
      ! Index i0 : camber panel index
      ! Index i1, i2 : corresponding blade panel index
      do ib = 1, nblade
        do ix = 1, nh
          do ir = 1, mr
            i0  = indexCamber(ix, ir)
            ix1 = nh  + ix   ! Suction side
            ix2 = nhp - ix   ! Pressure side, opposite normal vector relative to camber panel
            i1  = indexb(ix1, ir)
            i2  = indexb(ix2, ir)
            tc_pot(ib,i0)  = (tb_pot(ib,i1) *ss(i1,1) -
     &                        tb_pot(ib,i2) *ss(i2,1)   )/ssc(i0,1)
            tc_dpdn(ib,i0) = (tb_dpdn(ib,i1)*ss(i1,1) +
     &                        tb_dpdn(ib,i2)*ss(i2,1)   )/ssc(i0,1)

            ! Test_Normal: Use only the normal component
            ! ndot1 = dirc(i0,3,1)*dir(i1,3,1) + dirc(i0,3,2)*dir(i1,3,2)
            ! &                                       + dirc(i0,3,3)*dir(i1,3,3)
            ! ndot2 = dirc(i0,3,1)*dir(i2,3,1) + dirc(i0,3,2)*dir(i2,3,2)
            ! &                                       + dirc(i0,3,3)*dir(i2,3,3)
            ! tc_pot(ib,i0) = (tb_pot(ib,i1)*ss(i1,1)*ndot1 +
            ! &                       tb_pot(ib,i2)*ss(i2,1)*ndot2  )/ssc(i0,1)

          end do
        end do
      end do

!$omp parallel private(i,j,k,ib,ib1,
!$omp&                 pv_xv,pv_yv,pv_side,pv_s,p_chl,ssx,ssy,ssz,ssr,sst,st0,
!$omp&                 xloc,yloc,zloc,fs,fd,fsx,fsy,fdx,fdy,fdz,
!$omp&                 vdx,vdy,vdz,vsx,vsy,vsz,vx0,vy0,vz0,vr0,vt0)
!$omp&         shared (n_pt,n_cb,nblade,p_loc,pi,
!$omp&                 xvpc,yvpc,sidc,ssc,chlenc,dirc,xctc,
!$omp&                 tc_pot,tc_dpdn)
!$omp& reduction (+:v_tmp)
!$omp do collapse(2) schedule(dynamic,100)
      do i = 1, n_cb
      do j = 1, n_pt
        pv_xv(1:4)   = xvpc(i, 1:4)
        pv_yv(1:4)   = yvpc(i, 1:4)
        pv_side(1:4) = sidc(i, 1:4)
        pv_s(1:15)   = ssc(i, 1:15)
        p_chl = chlenc(i)
        ssx = p_loc(1,j)
        ssr = sqrt (p_loc(2,j)**2 + p_loc(3,j)**2)
        st0 = atan2(p_loc(3,j), p_loc(2,j))

        do ib = 1,nblade
          sst = st0 + 2.0 * pi / dble(nblade) * dble(ib - 1)
          ssy = ssr * cos(sst)
          ssz = ssr * sin(sst)

          xloc = (ssx-xctc(i,1))*dirc(i,1,1)+(ssy-xctc(i,2))*dirc(i,1,2)
     &                                      +(ssz-xctc(i,3))*dirc(i,1,3)
          yloc = (ssx-xctc(i,1))*dirc(i,2,1)+(ssy-xctc(i,2))*dirc(i,2,2)
     &                                      +(ssz-xctc(i,3))*dirc(i,2,3)
          zloc = (ssx-xctc(i,1))*dirc(i,3,1)+(ssy-xctc(i,2))*dirc(i,3,2)
     &                                      +(ssz-xctc(i,3))*dirc(i,3,3)

          call rpan_dbl(xloc,yloc,zloc,p_chl,fs,fd,fsx,fsy,
     &                 fdx,fdy,fdz,1,1,pv_xv,pv_yv,pv_s,pv_side)

          vdx = fdx*dirc(i,1,1) + fdy*dirc(i,2,1) + fdz*dirc(i,3,1)
          vdy = fdx*dirc(i,1,2) + fdy*dirc(i,2,2) + fdz*dirc(i,3,2)
          vdz = fdx*dirc(i,1,3) + fdy*dirc(i,2,3) + fdz*dirc(i,3,3)
          ! Temporarily turn off the source effect
          vsx = 0.0 !fsx*dirc(i,1,1) + fsy*dirc(i,2,1) - fd *dirc(i,3,1)
          vsy = 0.0 !fsx*dirc(i,1,2) + fsy*dirc(i,2,2) - fd *dirc(i,3,2)
          vsz = 0.0 !fsx*dirc(i,1,3) + fsy*dirc(i,2,3) - fd *dirc(i,3,3)

          do ib1 = 1, nblade
            k = ib + ib1 - 1
            if (k .gt. nblade) k = k - nblade
            vx0 = -1.0* (tc_pot(k,i)* vdx - tc_dpdn(k,i)* vsx) /4.0 /pi
            vy0 = -1.0* (tc_pot(k,i)* vdy - tc_dpdn(k,i)* vsy) /4.0 /pi
            vz0 = -1.0* (tc_pot(k,i)* vdz - tc_dpdn(k,i)* vsz) /4.0 /pi
            vr0 = vz0 * sin(sst) + vy0 * cos(sst)
            vt0 = vz0 * cos(sst) - vy0 * sin(sst)

            v_tmp(1,j,ib1) = v_tmp(1,j,ib1) + vx0
            v_tmp(2,j,ib1) = v_tmp(2,j,ib1) + vr0
            v_tmp(3,j,ib1) = v_tmp(3,j,ib1) + vt0
          end do

        end do
      end do
      end do
!$omp end do
!$omp end parallel
      v_ind  = v_ind + v_tmp

      call system_clock(om_ct2)
      om_dt=(real(om_ct2)-real(om_ct1))/real(om_ctr)
      write(*, "(A1,F4.1)", advance = "no") " ", om_dt


      return
      end subroutine

      end module
