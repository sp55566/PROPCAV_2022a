      module M_UNS_2D_INTERSECT
        use M_UNS_GLOBAL_PAR
        implicit none

        public  :: UNS_XR_CORR
        private :: set_point, intArea, its_pnt, polySeq,
     &             polyArea, find_lower, find_higher, iisd

        type POINT
          double precision :: x, y
        end type
                      !! ===========================
                      !!
                      !!      P4  ------- P3   o R (y)
                      !!        /       /     /|\
                      !!       /  QUAD /       |
                      !!      /       /        |
                      !!  P1  ------- P2       o---> X
                      !!
        type QUAD     !! ===========================
          type(POINT) :: p(4)
          double precision :: xmin, xmax, rmin, rmax
        end type

      contains

! ----- Add intersecting volumes to the matrix
        subroutine UNS_XR_CORR(xProp, rProp, xNode, rNode, cor)
          implicit none
          double precision :: xProp(pc_nhp,pc_mrp), rProp(pc_nhp,pc_mrp)
          double precision :: xNode(mesh_nx+1,mesh_nr+1)
          double precision :: rNode(mesh_nx+1,mesh_nr+1)
          double precision :: cor(mesh_nx,mesh_nr,pc_nh,pc_mr)

          integer px, pr, ix, ir, rmin, rmax
          double precision :: pminR, pmaxR, pminX, pmaxX
          double precision :: delt, rmid, areaTotal, areaSum
          double precision, dimension(mesh_nr) :: blkRmin, blkRmax
          type(QUAD), dimension(mesh_nx, mesh_nr) :: mesh
          type(QUAD) :: panel

          delt = uns_pi * 2.0 / mesh_nt

          ! Generate RANS mesh data structure
          do ir = 1, mesh_nr
            do ix = 1, mesh_nx
              ! Set all 4 points
      call set_point(mesh(ix,ir)%p(1),xNode(ix  ,ir  ),rNode(ix  ,ir  ))
      call set_point(mesh(ix,ir)%p(2),xNode(ix+1,ir  ),rNode(ix+1,ir  ))
      call set_point(mesh(ix,ir)%p(3),xNode(ix+1,ir+1),rNode(ix+1,ir+1))
      call set_point(mesh(ix,ir)%p(4),xNode(ix  ,ir+1),rNode(ix  ,ir+1))
              mesh(ix,ir)%xmin = minval(xNode(ix:ix+1, ir:ir+1))
              mesh(ix,ir)%xmax = maxval(xNode(ix:ix+1, ir:ir+1))
              mesh(ix,ir)%rmin = minval(rNode(ix:ix+1, ir:ir+1))
              mesh(ix,ir)%rmax = maxval(rNode(ix:ix+1, ir:ir+1))
            end do
          end do

          ! helper to make the scheme more efficient
          do ir = 1, mesh_nr
            blkRmin(ir) = minval(mesh(:,ir)%rmin)
            blkRmax(ir) = maxval(mesh(:,ir)%rmax)
          end do

          ! Begin searching
          ! Begin searching
          do pr = 1,pc_mr
          do px = 1,pc_nh
            ! PROPCAV panel setup
            call set_point(panel%p(1),xProp(px  ,pr  ),rProp(px  ,pr  ))
            call set_point(panel%p(2),xProp(px+1,pr  ),rProp(px+1,pr  ))
            call set_point(panel%p(3),xProp(px+1,pr+1),rProp(px+1,pr+1))
            call set_point(panel%p(4),xProp(px  ,pr+1),rProp(px  ,pr+1))
            panel%xmin = minval(xProp(px:px+1, pr:pr+1))
            panel%xmax = maxval(xProp(px:px+1, pr:pr+1))
            panel%rmin = minval(rProp(px:px+1, pr:pr+1))
            panel%rmax = maxval(rProp(px:px+1, pr:pr+1))
            rmid = sum(rProp(px:px+1, pr:pr+1))/4.0

            areaTotal = polyArea(panel%p,4)             ! Debugging
            areaSum = 0.0d0                             ! Debugging

            ! Search ir
            rmin = find_lower (blkRmin, mesh_nr, panel%rmin)
            rmax = find_higher(blkRmax, mesh_nr, panel%rmax)
            do ir = rmin,rmax
            do ix = 1, mesh_nx
              if ( panel%xmin < mesh(ix,ir)%xmax .and.
     &             panel%xmax > mesh(ix,ir)%xmin .and.
     &             panel%rmin < mesh(ix,ir)%rmax .and.
     &             panel%rmax > mesh(ix,ir)%rmin       ) then
                cor(ix,ir,px,pr) = delt*rmid*intArea(panel,mesh(ix,ir))
                areaSum = areaSum + intArea(panel,mesh(ix,ir))   ! Debugging
              end if
            end do
            end do

            if (abs(areaTotal-areaSum) .gt. 0.0003) then           ! Debugging
              write (*,*) areaTotal, areaSum, px, pr               ! Debugging
              call debug_intersection(panel)                       ! Debugging
            end if                                                 ! Debugging
            if (pr == 2 .and. px == pc_nh/3-1) then
              call debug_intersection(panel)
            end if
          end do
          end do

        contains

          subroutine debug_intersection(panel)
            type(QUAD) :: panel
            open(284, file='Test_int_2D.plt', status='unknown')
            write(284,*) 'Variables = x,r'
            ! Print panel
            call print_polygon(284, panel%p, 4, 'PANEL')
            ! Print sub-panel
            do ir = rmin,rmax
            do ix = 1, mesh_nx
              if ( panel%xmin < mesh(ix,ir)%xmax .and.
     &             panel%xmax > mesh(ix,ir)%xmin .and.
     &             panel%rmin < mesh(ix,ir)%rmax .and.
     &             panel%rmax > mesh(ix,ir)%rmin       ) then
                call print_inter_poly(284, panel, mesh(ix,ir))
              end if
            end do
            end do
            close(284)
          end subroutine

          subroutine print_polygon(fh, pint, siz, tname)
            integer siz, fh, i
            character(len=5) :: tname
            type(POINT) :: pint(siz)
            write(fh,*) 'zone T = "', tname, '", I = ', siz+1
            do i = 1, siz
              write(fh, *) pint(i)%x, pint(i)%y
            end do
            write(fh, *) pint(1)%x, pint(1)%y
          end subroutine

          subroutine print_inter_poly(fh, quad1, quad2)
            integer :: fh
            type(QUAD) :: quad1, quad2
            type(POINT):: points(24), ptemp
            logical :: ltemp
            integer i, j, npts, i1, j1
            ! Collect all inner points
            npts = 0
            do i = 1,4
              if (iisd(quad1%p(i),quad2)) then
                npts = npts + 1
                points(npts) = quad1%p(i)
              end if
            end do
            do i = 1,4
              if (iisd(quad2%p(i),quad1)) then
                npts = npts + 1
                points(npts) = quad2%p(i)
              end if
            end do
            ! Collect intersection between edges
            do i = 1,4
              i1 = i + 1
              if (i1 == 5) i1 = 1
              do j = 1,4
                j1 = j + 1
                if (j1 == 5) j1 = 1
      ltemp=its_pnt(quad1%p(i),quad1%p(i1),quad2%p(j),quad2%p(j1),ptemp)
                if (ltemp) then
                  npts = npts + 1
                  points(npts) = ptemp
                end if
              end do
            end do
            if (npts >= 3) then
              ! Rearrange the sequence
              call polySeq(points(1:npts),npts)
              call print_polygon(fh, points, npts, 'parts')
            end if
          end subroutine

        end subroutine

! ----- Calculate intersecting area
        function intArea(quad1, quad2)
          implicit none
          double precision :: intArea
          type(QUAD) quad1, quad2
          type(POINT):: points(24), ptemp
          logical :: ltemp
          integer i, j, npts, i1, j1

          ! Collect all inner points
          npts = 0
          do i = 1,4
            if (iisd(quad1%p(i),quad2)) then
              npts = npts + 1
              points(npts) = quad1%p(i)
            end if
          end do
          do i = 1,4
            if (iisd(quad2%p(i),quad1)) then
              npts = npts + 1
              points(npts) = quad2%p(i)
            end if
          end do

          ! Collect intersection between edges
          do i = 1,4
            i1 = i + 1
            if (i1 == 5) i1 = 1

            do j = 1,4
              j1 = j + 1
              if (j1 == 5) j1 = 1

      ltemp=its_pnt(quad1%p(i),quad1%p(i1),quad2%p(j),quad2%p(j1),ptemp)

              if (ltemp) then
                npts = npts + 1
                points(npts) = ptemp
              end if

            end do
          end do

          if (npts >= 3) then
            ! Rearrange the sequence
            call polySeq(points(1:npts),npts)
            ! Calculate polygon area
            intArea = polyArea(points(1:npts),npts)
          else
            intArea = 0.0d0
          end if
        end function

! ----- Get the intersection point between two line segments
        function its_pnt(p1, p2, p3, p4, pint)
          ! Check this website for equations
          ! http://www-cs.ccny.cuny.edu/~wolberg/capstone/intersection/Intersection%20point%20of%20two%20lines.html
          logical its_pnt
          type(POINT) :: p1, p2, p3, p4, pint
          double precision :: ua, ub, below
          double precision, parameter :: tol = 1.0d-8

          below=(p4%y - p3%y)*(p2%x - p1%x)-(p4%x - p3%x)*(p2%y - p1%y)
          if (dabs(below) < tol) then
            its_pnt = .false.
            return
          end if
          ua = (p4%x - p3%x)*(p1%y - p3%y)-(p4%y - p3%y)*(p1%x - p3%x)
          ub = (p2%x - p1%x)*(p1%y - p3%y)-(p2%y - p1%y)*(p1%x - p3%x)
          ua = ua/below
          ub = ub/below
          if (ua<0.0 .or. ua>1.0 .or. ub<0.0 .or. ub>1.0) then
            its_pnt = .false.
            return
          end if

          its_pnt = .true.
          pint%x = p1%x + ua * (p2%x - p1%x)
          pint%y = p1%y + ua * (p2%y - p1%y)

        end function

! ----- Arrange the sequence of several points in counter clockwise direction
        subroutine polySeq(pin, n)
          implicit none
          integer n, i, jmax, imax
          type(POINT) :: pin(n), ptemp
          double precision :: ttt(n), xmid, tmid, ttemp

          xmid = 0.0
          tmid = 0.0
          do i = 1,n
            xmid = xmid + pin(i)%x
            tmid = tmid + pin(i)%y
          end do
          xmid = xmid / n
          tmid = tmid / n

          do i = 1,n
            ttt(i) = atan2(pin(i)%y - tmid, pin(i)%x - xmid)
          end do

          imax = n
          do i = 1,(n-1)
            jmax = maxloc(ttt(1:imax),1)
            ttemp = ttt(imax); ttt(imax) = ttt(jmax); ttt(jmax) = ttemp
            ptemp = pin(imax); pin(imax) = pin(jmax); pin(jmax) = ptemp
            imax = imax - 1
          end do

        end subroutine

! ----- Calculate Polygon area
        function polyArea(pin, n)
          implicit none
          integer n, i, j
          type(POINT) :: p(n), pin(n)
          double precision :: polyArea

          do i = 1,n
            p(i)%x  = pin(i)%x  - pin(1)%x
            p(i)%y = pin(i)%y - pin(1)%y
          end do

          polyArea = 0.0
          do i = 1,n
            j = i + 1
            if (j.gt.n) j = 1
            polyArea = polyArea + p(i)%x * p(j)%y - p(j)%x * p(i)%y
          end do
          polyArea = abs(polyArea) * 0.5
        end function

! ----- Set value for
        subroutine set_point(p, x, y)
          implicit none
          type(POINT) :: p
          double precision :: x, y
          p%x = x
          p%y = y
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

! ----- Find the local coordinate of a point in a cell (0-1)->(0-1)
        logical function iisd(pp, panel)
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
          type(POINT) :: ps, pe, pp
          type(QUAD) :: panel
          double precision :: t1, t2, t3, t4, a, b, c

          iisd=(.false.)
          ! - equation
          ! (y-ys)(xe-xs)=(x-xs)(ye-ys)
          !  ax+by+c=0  a=ye-ys  b=-1(xe-xs) c=ys(xe-xs)-xs(ye-ys)
          ps = panel%p(1)
          pe = panel%p(2)
          a  = pe%y - ps%y
          b  = ps%x - pe%x
          c  = ps%y * (pe%x - ps%x) - ps%x * (pe%y - ps%y)
          t1 = a * pp%x + b * pp%y + c

          ps = panel%p(2)
          pe = panel%p(3)
          a  = pe%y - ps%y
          b  = ps%x - pe%x
          c  = ps%y * (pe%x - ps%x) - ps%x * (pe%y - ps%y)
          t2 = a * pp%x + b * pp%y + c

          ps = panel%p(3)
          pe = panel%p(4)
          a  = pe%y - ps%y
          b  = ps%x - pe%x
          c  = ps%y * (pe%x - ps%x) - ps%x * (pe%y - ps%y)
          t3 = a * pp%x + b * pp%y + c

          ps = panel%p(4)
          pe = panel%p(1)
          a  = pe%y - ps%y
          b  = ps%x - pe%x
          c  = ps%y * (pe%x - ps%x) - ps%x * (pe%y - ps%y)
          t4 = a * pp%x + b * pp%y + c

          if (((t1*t3).ge.0).and.((t2*t4).ge.0)) iisd = .true.
*          if (dabs(t1)<1e-8 .or. dabs(t1)<1e-8) iisd = .true.
*          if (dabs(t3)<1e-8 .or. dabs(t4)<1e-8) iisd = .true.

          return
        end function



      end module
