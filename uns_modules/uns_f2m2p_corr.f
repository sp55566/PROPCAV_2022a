! ---------------------------------------------------------------------------- !
!                       Created by Yiran Su 10/xx/2017                         !
!                 Build Correlation between Fluent Cell/Mesh/PROPCAV           !
! ---------------------------------------------------------------------------- !
      module M_UNS_F2M2P_CORR
        use M_UNS_GLOBAL_PAR

        public  :: UNS_READ_RANS, UNS_READ_MESH,UNS_BLD_R2M,UNS_WRT_R2M,
     &             UNS_READ_PROPCAV, UNS_BLD_P2M, UNS_BLD_M2P

        private :: find_lower, find_higher, quad_area, debugging_plt,
     &             calc_eta, iisd

        ! RANS Cells + its correlation to Mesh
        type :: CELL
          integer           :: pid, cid       ! Partition ID and cell ID
          double precision  :: xc,  yc,  zc  ! Centroid coordinates
          integer           :: idx, idr, idt  ! Corresponding Ix Ir It in Mesh
        end type CELL
        type(CELL),allocatable :: ransMesh(:)
        integer totCells

        ! Structured Mesh
        double precision, dimension(:,:,:),allocatable ::
     &                                xNode,yNode,zNode,rNode ! Node coordinate
        double precision, dimension(:,:), allocatable :: meshVol  ! mesh volume

        ! PROPCAV panels
        double precision,dimension(:,:),allocatable:: bxx, brr, byy, bzz

        ! Data structure for interpolating RANS total flow to PROPCAV mesh (with offset)
        type PC_MSH
          integer :: ix, ir, it
          double precision :: X, Y, Z, R, T
          double precision :: xper, rper, tper
        end type
        type(PC_MSH),allocatable :: mesh2pc(:,:,:)  ! pc_nc * pc_mr * nSeed

      contains

! ----- Read RANS cell information
        subroutine UNS_READ_RANS
          implicit none
          integer pid, numCells
          integer i, j, ierr

          ! get number of cells in total
          open(285, file = 'bf_zone_info.dat', status = 'old')
          totCells = 0
          do
            read(285, *, iostat = ierr) pid, numCells
            if (ierr /= 0) exit
            totCells = totCells + numCells
            do i = 1, numCells
              read(285, *)
            end do
          end do
          close(285)
          write(*,*) 'Total number of cells:',totCells

          ! read in the cell info from RANS
          allocate(ransMesh(totCells));
          open(285, file = 'bf_zone_info.dat', status = 'old')
          j = 0
          do
            read(285, *, iostat = ierr) pid, numCells
            if (ierr /= 0) exit
            do i = 1, numCells
              j = j + 1;
              ransMesh(j)%pid = pid;
              read(285, *) ransMesh(j)%cid, ransMesh(j)%xc,
     &                     ransMesh(j)%yc,  ransMesh(j)%zc
            end do
          end do
          close(285)
        end subroutine

! ----- Read structured mesh information
        subroutine UNS_READ_MESH
          implicit none
          integer :: nId(8), ix, ir, it, i, j
          double precision xxx, yyy, zzz, rr, rmin, rmax

          allocate (xNode(mesh_nx+1,mesh_nr+1,mesh_nt+1))
          allocate (yNode(mesh_nx+1,mesh_nr+1,mesh_nt+1))
          allocate (zNode(mesh_nx+1,mesh_nr+1,mesh_nt+1))
          allocate (rNode(mesh_nx+1,mesh_nr+1,mesh_nt+1))

          open(285, file = 'mesh3d.plt', status = 'old')
          read(285, *)
          read(285, *)
          read(285, *)

          ! Read node information
          do it = 1, mesh_nt+1
          do ir = 1, mesh_nr+1
          do ix = 1, mesh_nx+1
            read(285,*) xNode(ix,ir,it),yNode(ix,ir,it),zNode(ix,ir,it)
          end do
          end do
          end do
          close(285)

          rNode = sqrt(yNode**2 + zNode**2)

        end subroutine

! ----- Read propcav panel info to this module
        subroutine UNS_READ_PROPCAV
          use M_BLADEM, only : xbm, ybm, zbm
          use GEOCAMB, only : xctc, dirc, indexCamber
          implicit none
          integer :: i, j, i1, i2, i0, j0, m, k
          allocate (bxx(pc_nhp,pc_mrp), byy(pc_nhp,pc_mrp))
          allocate (bzz(pc_nhp,pc_mrp), brr(pc_nhp,pc_mrp))
          do j = 1, pc_mrp
            do i = 1, pc_nhp
              bxx(i,j) = xbm(i,j)
              byy(i,j) = ybm(i,j)
              bzz(i,j) = zbm(i,j)
              brr(i,j) = sqrt(byy(i,j)**2 + bzz(i,j)**2)
            end do
          end do

          allocate(mesh2pc(pc_nc, pc_mr, nSeed))
          do j = 1, pc_mr
            do i = 1, pc_nh
              i1 = pc_nh + i   ! Suction side
              i2 = pc_nhp - i  ! Pressure side, opposite normal vector

              ! The normal vectors at the leading edge and the tip are not good
              i0 = i
              j0 = j
              if (i0.le.2) i0 = 3
              if (j0.ge.(pc_mr - 1)) j0 = pc_mr - 2
              m = indexCamber(i0, j0)

              do k = 1, nSeed
                mesh2pc(i1,j,k)%X = xctc(m,1) - dirc(m,3,1) * hSeed * k
                mesh2pc(i1,j,k)%Y = xctc(m,2) - dirc(m,3,2) * hSeed * k
                mesh2pc(i1,j,k)%Z = xctc(m,3) - dirc(m,3,3) * hSeed * k
                mesh2pc(i1,j,k)%R = sqrt ( mesh2pc(i1,j,k)%Y**2 +
     &                                     mesh2pc(i1,j,k)%Z**2   )
                mesh2pc(i1,j,k)%T = atan2( mesh2pc(i1,j,k)%Z,
     &                                     mesh2pc(i1,j,k)%Y      )

                mesh2pc(i2,j,k)%X = xctc(m,1) + dirc(m,3,1) * hSeed * k
                mesh2pc(i2,j,k)%Y = xctc(m,2) + dirc(m,3,2) * hSeed * k
                mesh2pc(i2,j,k)%Z = xctc(m,3) + dirc(m,3,3) * hSeed * k
                mesh2pc(i2,j,k)%R = sqrt ( mesh2pc(i2,j,k)%Y**2 +
     &                                     mesh2pc(i2,j,k)%Z**2   )
                mesh2pc(i2,j,k)%T = atan2( mesh2pc(i2,j,k)%Z,
     &                                     mesh2pc(i2,j,k)%Y      )
              end do
            end do
          end do
        end subroutine

! ----- Build RANS to mesh correlation
        subroutine UNS_BLD_R2M
          implicit none
          integer ic, ix, ir, it, itt
          integer xmin, xmax, rmin, rmax, tmin, tmax
          double precision dist, dist0, distp, rc, tc
          double precision, dimension(mesh_nr) :: blkRmin, blkRmax
          double precision                     :: blkX(mesh_nx, mesh_nr)
          double precision, dimension(mesh_nx,mesh_nr,mesh_nt) ::
     &                                           xCent,yCent,zCent,dCent
          dCent = 1.0

          ! Calculate centroid
          do it = 1, mesh_nt
          do ir = 1, mesh_nr
          do ix = 1, mesh_nx
            xCent(ix,ir,it) = sum(xNode(ix:ix+1, ir:ir+1, it:it+1))/8.0
            yCent(ix,ir,it) = sum(yNode(ix:ix+1, ir:ir+1, it:it+1))/8.0
            zCent(ix,ir,it) = sum(zNode(ix:ix+1, ir:ir+1, it:it+1))/8.0
          end do
          end do
          end do

          ! helper to make the scheme more efficient
          do ir = 1, mesh_nr
            blkRmin(ir) = minval(rNode(:,ir  ,1))
            blkRmax(ir) = maxval(rNode(:,ir+1,1))
            do ix = 1, mesh_nx
              blkX(ix,ir) = xCent(ix,ir,1)
            end do
          end do

          ! Begin searching
          do ic = 1, totCells
            if (mod(ic,30000)==0) write(*,*) ic,'/',totCells,'completed'
            dist0 = 99999.
            rc = sqrt(ransMesh(ic)%yc**2 + ransMesh(ic)%zc**2)
            tc = atan2(ransMesh(ic)%zc, ransMesh(ic)%yc)

            rmin = find_lower (blkRmin, mesh_nr, rc)
            rmax = find_higher(blkRmax, mesh_nr, rc)

            do ir = rmin, rmax

              xmin = find_lower(blkX(:,ir), mesh_nx, ransMesh(ic)%xc)
              xmax = xmin + 1
              if (xmax > mesh_nx) xmax = xmax - 1

              do ix = xmin, xmax

                do it = 1, mesh_nt
                  dist = sqrt( (xCent(ix,ir,it)-ransMesh(ic)%xc)**2 +
     &                         (yCent(ix,ir,it)-ransMesh(ic)%yc)**2 +
     &                         (zCent(ix,ir,it)-ransMesh(ic)%zc)**2   )
                  if (dist < dist0) then
                    dist0 = dist
                    ransMesh(ic)%idx = ix
                    ransMesh(ic)%idr = ir
                    ransMesh(ic)%idt = it
                  end if
                end do

              end do
            end do
      dCent(ransMesh(ic)%idx,ransMesh(ic)%idr,ransMesh(ic)%idt)=dist0
          end do

*          ! test plot the cell to mesh correlation
*           open(285, file='Test_mesh_centroid.plt', status='unknown')
*           write(285,*) 'Variables = x,y,z,dist'
*           write(285,*) 'Zone T="mesh", I=',mesh_nx,',J=',mesh_nr,',K=',mesh_nt+1
*           write(*,*) 'Writing test files for debugging...'
*           do it = 1, mesh_nt + 1
*             itt = it
*             if (itt == mesh_nt + 1) itt = itt - mesh_nt
*             do ir = 1, mesh_nr
*               do ix = 1, mesh_nx
*                 write(285,*) xCent(ix,ir,itt),yCent(ix,ir,itt),
*     &                        zCent(ix,ir,itt),dCent(ix,ir,itt)
*               end do
*             end do
*           end do
*           close(285)

        end subroutine

! ----- Find in list x the element that is closest to t! == Binary search
        integer function find_lower(x, nx, t)
          implicit none
          integer :: nx
          double precision, dimension(nx) :: x
          double precision :: t
          integer :: nstart, nend, nmid
          nstart = 1
          nend = nx
          do while (nend - nstart >= 2)
            nmid = (nstart + nend) / 2
            if (t >= x(nmid)) then
              nstart = nmid
            else
              nend = nmid
            end if
          end do
          find_lower = nstart
        end function
        integer function find_higher(x, nx, t)
          implicit none
          integer :: nx
          double precision, dimension(nx) :: x
          double precision :: t
          integer :: nstart, nend, nmid
          nstart = 1
          nend = nx
          do while (nend - nstart >= 2)
            nmid = (nstart + nend) / 2
            if (t >= x(nmid)) then
              nstart = nmid
            else
              nend = nmid
            end if
          end do
          find_higher = nend
        end function

! ----- Write correlation info to Fluent
        subroutine UNS_WRT_R2M
          implicit none
          integer i
          write(*,*) 'Writing cell to mesh correlation for Fluent ...'
          open(285, file='corr.dat', status='unknown')
          write(285,*) mesh_nx, mesh_nr, mesh_nt, pc_nbld, totCells
          do i = 1, totCells
            write(285,*) ransMesh(i)%pid, ransMesh(i)%cid,
     &           ransMesh(i)%idx-1, ransMesh(i)%idr-1, ransMesh(i)%idt-1
          end do
          close(285)
        end subroutine

! ----- Build PROPCAV to mesh correlation (for body force conserved interpolation)
        subroutine UNS_BLD_P2M
          use M_LNK_LST_P2M
          use M_UNS_2D_INTERSECT, only : UNS_XR_CORR
          implicit none
          integer px, pr, ix, ir, icount
          double precision :: rCent, deltT, area, weight
          double precision,allocatable :: cor(:,:,:,:), propVol(:,:)
          allocate(cor(mesh_nx,mesh_nr,pc_nh,pc_mr)) ! intersecting volume from nc.mr panel to nx.nr cell
          allocate(meshVol(mesh_nx,mesh_nr), propVol(pc_nh,pc_mr))
          cor = 0.0
          propVol = 0.0

          ! Calculate
          call UNS_XR_CORR(bxx, brr, xNode(:,:,1), rNode(:,:,1), cor)

          ! Add propVol together
          do pr = 1,pc_mr
          do px = 1,pc_nh
            do ir = 1,mesh_nr
            do ix = 1,mesh_nx
              propVol(px,pr) = propVol(px,pr) + cor(ix,ir,px,pr)
            end do
            end do
          end do
          end do

          ! Calculate mesh cell volume meshVol(ix,ir)
          do ir = 1,mesh_nr
            do ix = 1,mesh_nx
              rCent = sum(rNode(ix:ix+1, ir:ir+1, 1))/4.0
              deltT = 2.0d0 * uns_pi / mesh_nt
              area  = quad_area( xNode(ix:ix+1, ir:ir+1, 1),
     &                           rNode(ix:ix+1, ir:ir+1, 1)  )
              meshVol(ix,ir) = area * deltT * rCent
            end do
          end do

          ! Final output
          icount = 0
          do pr = 1,pc_mr
          do px = 1,pc_nh
          do ir = 1,mesh_nr
          do ix = 1,mesh_nx
            if (cor(ix,ir,px,pr).gt.0.0) then
              icount = icount + 1
              weight = cor(ix,ir,px,pr)/propVol(px,pr)/meshVol(ix,ir)
              call insertCor(px,pr,ix,ir,weight)
            end if
          end do
          end do
          end do
          end do
          write(*,*) 'Data compressed ratio:', 1.0*icount/pc_mr/pc_nh/mesh_nr/mesh_nx

          ! Plot for debugging purpose, Plot Db2 - see if propVol is correct - see if summation of weight is 1
          call debugging_plt(propVol)

          ! No longer needed
          deallocate(cor, propVol)

        end subroutine

! ----- Calculate Quad area
        function quad_area(x, y)
          implicit none
          integer i, j
          double precision, dimension(2,2) :: x,y
          double precision :: quad_area
          quad_area =   (x(2,1)-x(1,1)) * (y(2,2)-y(1,1))
     &                - (x(2,2)-x(1,1)) * (y(2,1)-y(1,1))
     &                + (x(2,2)-x(1,1)) * (y(1,2)-y(1,1))
     &                - (x(1,2)-x(1,1)) * (y(2,2)-y(1,1))
          quad_area = abs(quad_area) * 0.5
        end function

! ----- Build mesh to PROPCAV correlation (for total flow interpolation)
        subroutine UNS_BLD_M2P
        implicit none
        integer :: px, pr, ps, ix, ir, it, rmin, rmax, xmin, xmax, ct
        logical :: isFound
        double precision :: xprop, rprop, tprop, etaX, etaR, etaT
        double precision, dimension(mesh_nr)       :: blkRmin, blkRmax
        double precision, dimension(mesh_nx, mesh_nr) :: minX, maxX
        double precision, dimension(2):: p1,p2,p3,p4,pp

          ! helper to make the scheme more efficient
          do ir = 1, mesh_nr
            blkRmin(ir) = minval(rNode(:,ir  ,1)) + 0.00001
            blkRmax(ir) = maxval(rNode(:,ir+1,1)) - 0.00001
          end do
          do ir = 1, mesh_nr
            do ix = 1, mesh_nx
              minX(ix,ir) = minval(xNode(ix:ix+1, ir:ir+1, 1)) - 0.00001
              maxX(ix,ir) = maxval(xNode(ix:ix+1, ir:ir+1, 1)) + 0.00001
            end do
          end do

          ! Begin searching
          do ps = 1,nSeed
          do pr = 1,pc_mr
          do px = 1,pc_nc
            ! PROPCAV point location
            xprop= mesh2pc(px,pr,ps)%X
            rprop= mesh2pc(px,pr,ps)%R
            tprop= mesh2pc(px,pr,ps)%T

            isFound = .false.
            ct = 0
911         continue

            ! Search ir and ix
            rmin = find_lower (blkRmin, mesh_nr, rprop)
            rmax = find_higher(blkRmax, mesh_nr, rprop)
            do ir = rmin, rmax
            do ix = 1, mesh_nx
            if(xprop.ge.minX(ix,ir) .and. xprop.le.maxX(ix,ir)) then
              p1 = (/xNode(ix  ,ir  ,1), rNode(ix  ,ir  ,1)/)
              p2 = (/xNode(ix+1,ir  ,1), rNode(ix+1,ir  ,1)/)
              p3 = (/xNode(ix+1,ir+1,1), rNode(ix+1,ir+1,1)/)
              p4 = (/xNode(ix  ,ir+1,1), rNode(ix  ,ir+1,1)/)
              pp = (/xprop,              rprop             /)
              !Found
              if (iisd(p1,p2,p3,p4,pp) == .true.) then
                call calc_eta(p1, p2, p3, p4, pp, etaX, etaR)
                if (etaX > 0.5) then
                    mesh2pc(px,pr,ps)%ix = ix
                    mesh2pc(px,pr,ps)%xper = etaX - 0.5
                else
                    mesh2pc(px,pr,ps)%ix = ix - 1
                    mesh2pc(px,pr,ps)%xper = etaX + 0.5
                end if
                if (etaR > 0.5 .or. ir == 1) then
                    mesh2pc(px,pr,ps)%ir = ir
                    mesh2pc(px,pr,ps)%rper = etaR - 0.5
                else
                    mesh2pc(px,pr,ps)%ir = ir - 1
                    mesh2pc(px,pr,ps)%rper = etaR + 0.5
                end if
                isFound = .true.
                ! Search it
                p1 = (/yNode(ix  ,ir  ,1), zNode(ix  ,ir  ,1)/)
                p2 = (/yNode(ix+1,ir  ,1), zNode(ix+1,ir  ,1)/)
                p3 = (/yNode(ix+1,ir+1,1), zNode(ix+1,ir+1,1)/)
                p4 = (/yNode(ix  ,ir+1,1), zNode(ix  ,ir+1,1)/)
                pp =   p3*etaX*etaR + p1*(1.0-etaX)*(1.0-etaR)
     &               + p4*(1.0-etaX)*etaR + p2*etaX*(1.0-etaR)
                etaT = atan2(pp(2), pp(1))
                etaT = etaT - tprop
                etaT = etaT/(uns_2pi/mesh_nt) - 0.5
                if (etaT < 0) etaT = etaT + mesh_nt
                it = 1 + floor(etaT)
                if (it > mesh_nt) it = it - mesh_nt
                etaT = etaT - floor(etaT)
                mesh2pc(px,pr,ps)%it = it
                mesh2pc(px,pr,ps)%tper = etaT

                exit
              end if
            end if
            end do
            end do

            if (isFound .eq. .false. .and. ct .le. 2) then
              ct = ct + 1
              write(*,*) 'Warning in UNS_BLD_M2P', ct
              write(*,*) xprop, rprop
              write(*,*) rmin, rmax
              xprop = xprop + 1d-4 * ct
              rprop = rprop + 1d-4 * ct
              write(*,*) xprop, rprop, '!!'
              goto 911
            end if

          end do
          end do
          end do
        end subroutine
! ----- Find the local coordinate of a point in a cell (0-1)->(0-1)
        subroutine calc_eta(p1, p2, p3, p4, pp, ebs, eta)
          !----------------------------------------
          !   p4----p7---- p3  eta (0-1)
          !   |  .        |   |
          !   |   pp      |   |
          !   p8    p9    p6  |
          !   |           |   |
          !   |           |   |
          !   p1----p5----p2   -------> ebs (0-1)
          !----------------------------------------
          implicit none
          double precision,dimension(2)::p1,p2,p3,p4,pp
          double precision,dimension(2)::p5,p6,p7,p8,p9
          double precision,dimension(2)::lp1,lp2,lp3,lp4
          double precision,dimension(2)::lp5,lp6,lp7,lp8,lp9
          double precision :: ebs,eta
          integer:: nitr,i,j
          nitr=20  !number of subdivision times
          if (iisd(p1,p2,p3,p4,pp).eq.(.false.)) then
            return
          end if
          lp1 = (/0.0, 0.0/)
          lp2 = (/1.0, 0.0/)
          lp3 = (/1.0, 1.0/)
          lp4 = (/0.0, 1.0/)
          do i = 1,nitr
            p5  = (p1  + p2 ) / 2.0d0
            p6  = (p2  + p3 ) / 2.0d0
            p7  = (p3  + p4 ) / 2.0d0
            p8  = (p4  + p1 ) / 2.0d0
            p9  = (p8  + p6 ) / 2.0d0
            lp5 = (lp1 + lp2) / 2.0d0
            lp6 = (lp2 + lp3) / 2.0d0
            lp7 = (lp3 + lp4) / 2.0d0
            lp8 = (lp4 + lp1) / 2.0d0
            lp9 = (lp8 + lp6) / 2.0d0
            if (iisd(p1,p5,p9,p8,pp).eq.(.true.)) then
              p1  = p1
              p2  = p5
              p3  = p9
              p4  = p8
              lp1 = lp1
              lp2 = lp5
              lp3 = lp9
              lp4 = lp8
              cycle
            else if(iisd(p5,p2,p6,p9,pp).eq.(.true.)) then
              p1  = p5
              p2  = p2
              p3  = p6
              p4  = p9
              lp1 = lp5
              lp2 = lp2
              lp3 = lp6
              lp4 = lp9
              cycle
            else if(iisd(p9,p6,p3,p7,pp).eq.(.true.)) then
              p1  = p9
              p2  = p6
              p3  = p3
              p4  = p7
              lp1 = lp9
              lp2 = lp6
              lp3 = lp3
              lp4 = lp7
              cycle
            else if(iisd(p8,p9,p7,p4,pp).eq.(.true.)) then
              p1  = p8
              p2  = p9
              p3  = p7
              p4  = p4
              lp1 = lp8
              lp2 = lp9
              lp3 = lp7
              lp4 = lp4
              cycle
            end if
          end do

          ebs = (lp1(1)+lp2(1)+lp3(1)+lp4(1)) / 4.0d0
          eta = (lp1(2)+lp2(2)+lp3(2)+lp4(2)) / 4.0d0
          return
        end subroutine

! ----- Find the local coordinate of a point in a cell (0-1)->(0-1)
        logical function iisd(p1,p2,p3,p4,pp)
          !----------------------------------------
          !               t3
          !           p4------ p3
          !           |  .     |
          !        t4 |   pp   | t2
          !           |        |
          !           p1------ p2
          !               t1
          !----------------------------------------
          implicit none
          double precision,dimension(2) :: p1,p2,p3,p4,pp,ps,pe
          double precision :: t1,t2,t3,t4,a,b,c

          iisd=(.false.)
          ! - equation
          ! (y-ys)(xe-xs)=(x-xs)(ye-ys)
          !  ax+by+c=0  a=ye-ys  b=-1(xe-xs) c=ys(xe-xs)-xs(ye-ys)
          ps=p1
          pe=p2
          a=pe(2)-ps(2)
          b=ps(1)-pe(1)
          c=ps(2)*(pe(1)-ps(1))-ps(1)*(pe(2)-ps(2))
          t1=a*pp(1)+b*pp(2)+c

          ps=p2
          pe=p3
          a=pe(2)-ps(2)
          b=ps(1)-pe(1)
          c=ps(2)*(pe(1)-ps(1))-ps(1)*(pe(2)-ps(2))
          t2=a*pp(1)+b*pp(2)+c

          ps=p3
          pe=p4
          a=pe(2)-ps(2)
          b=ps(1)-pe(1)
          c=ps(2)*(pe(1)-ps(1))-ps(1)*(pe(2)-ps(2))
          t3=a*pp(1)+b*pp(2)+c

          ps=p4
          pe=p1
          a=pe(2)-ps(2)
          b=ps(1)-pe(1)
          c=ps(2)*(pe(1)-ps(1))-ps(1)*(pe(2)-ps(2))
          t4=a*pp(1)+b*pp(2)+c

          if (((t1*t3).ge.0).and.((t2*t4).ge.0)) iisd=(.true.)

          return
        end function

! ----- Plot for debugging purpose
        subroutine debugging_plt(propVol)
          include 'PUFCAV.INC'
          include 'PUFCAVB.INC'
          double precision :: propVol(pc_nh, pc_mr)
          open(285, file='Test_prop_vol.plt', status='unknown')
          write(285,*) 'Variables = x,y,z,vol_sum'
          write(285,*) 'Zone T="mesh", I=',pc_nh,',J=',pc_mr
          write(*,*) 'Writing test file Db2 for debugging...'
          do ir = 1, pc_mr
          do ix = 1, pc_nh
            write(285,*) bxx(ix,ir),byy(ix,ir),bzz(ix,ir),propVol(ix,ir)
          end do
          end do
          close(285)
        end subroutine

      end module

