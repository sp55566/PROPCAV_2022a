! the alignment methods is the same as that used in Propcav

      module m_wake
      implicit none
      integer  mrw, ncw, nbld
      integer  mrpw, ncpw !mrpw = mrw+1, ncpw = ncw+1
      integer ndte
      real delta_w
      real omega_p  !the angular velocity in rad/s of the propeller
                    ! = 2*PI*Vs/J/D = PI/J
      real dth
      real theta_db
      real, allocatable:: wpt(:,:,:,:) !3, ncpw,mrpw, nbld
      real, allocatable:: wr(:,:,:)   ! ncpw,mrpw, nbld   wake radius
      real, allocatable:: wth(:,:,:)   ! ncpw,mrpw, nbld  wake angle
      real, allocatable:: wv(:,:,:,:)  !3, ncpw,mrpw, nbld
      real, allocatable:: wvi(:,:,:) 
      real, allocatable:: wqt(:,:,:)   ! ncpw,mrpw, nbld  wake angle
      real, allocatable:: wqr(:,:,:)   ! ncpw,mrpw, nbld  wake angle
      real, allocatable:: womg(:,:)    ! ncpw,mrpw (for angular velocity from inflow)
!     real, allocatable:: phiw(:,:,:)  !ncw,mrw, nbld

      real, allocatable:: wpx_dumy(:,:) !ncpw,mrpw
      real, allocatable:: wpy_dumy(:,:) !ncpw,mrpw
      real, allocatable:: wpz_dumy(:,:) !ncpw,mrpw
      real, allocatable:: wvxb(:,:) ! ncpw,mrpw,
      real, allocatable:: wvrb(:,:)
      real, allocatable:: wvtb(:,:)
      integer, allocatable:: ivbmiss(:,:)

      real, allocatable:: phiw0(:)     !mrw
      real, allocatable:: gamw0(:)     !mrpw
      real, allocatable:: gamw1(:,:,:) !nbld,nwpanel,mrpw 
      real, allocatable:: gamw2(:,:,:) !nbld,nwpanel+1,mrw
      real, allocatable:: err(:)       !error
      real err_m
      end module m_wake

      module m_dwake
      implicit none
      integer  mrdw, ncdw
      integer  mrpdw, ncpdw !mrpw = mrw+1, ncpw = ncw+1

      real, allocatable:: dwpt(:,:,:,:) !3, ncpdw,mrpdw, nbld
      real, allocatable:: dwr(:,:,:)   ! ncpdw,mrpdw, nbld   wake radius
      real, allocatable:: dwth(:,:,:)   ! ncpdw,mrpdw, nbld  wake angle
      real, allocatable:: dwv(:,:,:,:)  !3, ncpdw,mrpdw, nbld
      real, allocatable:: dwvi(:,:,:)
      real, allocatable:: dwqt(:,:,:)   ! ncpdw,mrpdw, nbld  wake angle
      real, allocatable:: dwqr(:,:,:)   ! ncpdw,mrpdw, nbld  wake angle
!      real, allocatable:: phiw(:,:,:)  !ncw,mrw, nbld
c
      real, allocatable:: dwpx_dumy(:,:) !ncpdw,mrpdw
      real, allocatable:: dwpy_dumy(:,:) !ncpdw,mrpdw
      real, allocatable:: dwpz_dumy(:,:) !ncpdw,mrpdw
      real, allocatable:: wvxd(:,:) ! ncpdw,mrpdw
      real, allocatable:: wvrd(:,:) ! ncpdw,mrpdw
      real, allocatable:: wvtd(:,:) ! ncpdw,mrpdw
      integer, allocatable:: ivdmiss(:,:) ! ncpdw,mrpdw
c
      real, allocatable:: phidw0(:)     !mrdw
      real, allocatable:: phidw1(:)     !mrdw
      real, allocatable:: gamdw0(:)     !mrpdw
      real, allocatable:: derr(:)       !error
      real derr_m
      end module m_dwake

      subroutine wk_interp1(x1,y1,z1,phi1,n1,
     &                     x2,y2,z2,phi2,n2)
      implicit none
      integer :: n1, n2
      real,dimension(n1) :: x1,y1,z1,s1,idx1
      real,dimension(n2) :: x2,y2,z2,s2,idx2
      real sp1(n1-1), sp2(n2-1)
      real phi1(n1-1), phi2(n2-1)

      real,dimension(4*(n1-1)) :: cubs,cubx,cuby,cubz,cubphi
      integer :: i
      real ds, vtmp(3)

      do i = 1, n1
        idx1(i) = real(i-1)/real(n1-1)
      end do
      do i = 1, n2
        idx2(i) = real(i-1)/real(n2-1)
      end do

      s1(1) = 0.0d0
      do i = 2,n1
        vtmp(1) = x1(i)-x1(i-1)
        vtmp(2) = y1(i)-y1(i-1)
        vtmp(3) = z1(i)-z1(i-1)
        ds = sqrt(dot_product(vtmp,vtmp))
        s1(i) = s1(i-1) + ds
      end do

      call uglydk(n1,0,0,idx1,s1,0.0,0.0,cubs)
      call evaldk(n1,n2,idx1,idx2,s2,cubs)

      sp1(1:n1-1) = 0.5*(s1(1:n1-1)+s1(2:n1))
      sp2(1:n2-1) = 0.5*(s2(1:n2-1)+s2(2:n2))

      call uglydk(n1,1,1,s1,x1,0.0,0.0,cubx)
      call uglydk(n1,1,1,s1,y1,0.0,0.0,cuby)
      call uglydk(n1,1,1,s1,z1,0.0,0.0,cubz)
      call evaldk(n1,n2,s1,s2,x2,cubx)
      call evaldk(n1,n2,s1,s2,y2,cuby)
      call evaldk(n1,n2,s1,s2,z2,cubz)

      call uglydk(n1-1,1,1,sp1,phi1,0.0,0.0,cubphi)
      call evaldk(n1-1,n2-1,sp1,sp2,phi2,cubphi)
      return
      end subroutine wk_interp1

      subroutine wk_interp2(x1,y1,z1,n1,
     &                      x2,y2,z2,n2)
      implicit none
      integer :: n1, n2
      real,dimension(n1) :: x1,y1,z1,s1,idx1  ! ,rr1,th1
      real,dimension(n2) :: x2,y2,z2,s2,idx2  ! ,rr2,th2

      real,dimension(4*(n1-1)) :: cubs,cubx,cuby,cubz,cubphi  ! ,cubr,cubth
      integer :: i
      real ds, vtmp(3)

      do i = 1, n1
        idx1(i) = real(i-1)/real(n1-1)
      end do
      do i = 1, n2
        idx2(i) = real(i-1)/real(n2-1)
      end do

      s1(1) = 0.0d0

      do i = 2,n1
        vtmp(1) = x1(i)-x1(i-1)
c        vtmp(1) = 0.0d0
        vtmp(2) = y1(i)-y1(i-1)
        vtmp(3) = z1(i)-z1(i-1)
        ds = sqrt(dot_product(vtmp,vtmp))
        s1(i) = s1(i-1) + ds

      end do

      call uglydk(n1,0,0,idx1,s1,0.0,0.0,cubs)
      call evaldk(n1,n2,idx1,idx2,s2,cubs)

      call uglydk(n1,0,0,s1,x1,0.0,0.0,cubx)
      call uglydk(n1,0,0,s1,y1,0.0,0.0,cuby)
      call uglydk(n1,0,0,s1,z1,0.0,0.0,cubz)
      call evaldk(n1,n2,s1,s2,x2,cubx)
      call evaldk(n1,n2,s1,s2,y2,cuby)
      call evaldk(n1,n2,s1,s2,z2,cubz)

      return
      end subroutine wk_interp2


      subroutine alloc_wk_mem(imrw,incw, inbld)
      use m_wake
      implicit none
      integer imrw, incw, inbld
      mrw = imrw
      ncw = incw
      nbld = inbld
      mrpw = mrw + 1
      ncpw = ncw + 1
      if (.NOT. allocated(wpt)) then
        allocate(wpt(3,ncpw,mrpw,nbld))
        allocate(wv(3,ncpw,mrpw,nbld))
        allocate(wvi(3,ncpw,mrpw))
        allocate(wr(ncpw,mrpw,nbld))
        allocate(wth(ncpw,mrpw,nbld))
        allocate(wqt(ncpw,mrpw,nbld))
        allocate(wqr(ncpw,mrpw,nbld))
        allocate(womg(ncpw,mrpw))

        allocate(wpx_dumy(ncpw,mrpw))
        allocate(wpy_dumy(ncpw,mrpw))
        allocate(wpz_dumy(ncpw,mrpw))
        allocate(wvxb(ncpw,mrpw))
        allocate(wvrb(ncpw,mrpw))
        allocate(wvtb(ncpw,mrpw))
        allocate(ivbmiss(ncpw,mrpw))

        allocate(phiw0(mrw))
        allocate(gamw0(mrpw))
        allocate(gamw1(nbld,ncpw,mrpw))
        allocate(gamw2(nbld,ncpw,mrw))
        allocate(err(mrpw))
      end if
      return
      end subroutine alloc_wk_mem

      subroutine alloc_dwk_mem(imrwd,incdw,inbld)
      use m_dwake
      implicit none
      integer imrwd, incdw, inbld
      mrdw = imrwd
      ncdw = incdw
      mrpdw = mrdw + 1
      ncpdw = ncdw + 1
      if (.NOT. allocated(dwpt)) then
        allocate(dwpt(3,ncpdw,mrpdw,inbld))
        allocate(dwv(3,ncpdw,mrpdw,inbld))
        allocate(dwvi(3,ncpdw,mrpdw))
        allocate(dwr(ncpdw,mrpdw,inbld))
        allocate(dwth(ncpdw,mrpdw,inbld))
        allocate(dwqt(ncpdw,mrpdw,inbld))
        allocate(dwqr(ncpdw,mrpdw,inbld))

        allocate(dwpx_dumy(ncpdw,mrpdw))
        allocate(dwpy_dumy(ncpdw,mrpdw))
        allocate(dwpz_dumy(ncpdw,mrpdw))
        allocate(wvxd(ncpdw,mrpdw))
        allocate(wvrd(ncpdw,mrpdw))
        allocate(wvtd(ncpdw,mrpdw))
        allocate(ivdmiss(ncpdw,mrpdw))

        allocate(phidw0(mrdw))
        allocate(phidw1(mrdw))
        allocate(gamdw0(mrpdw))
        allocate(derr(mrpdw))
      end if
      return
      end subroutine alloc_dwk_mem

      subroutine symm_wkp
      use m_wake
      use m_constant
      implicit none
      integer k
      do k = 2, nbld
!$OMP parallel shared(wpt,wr,wth,k,theta_db)
!$OMP workshare
        wpt(1,:,:,k)=wpt(1,:,:,1);
        wr(:,:,k)=wr(:,:,1);
        wth(:,:,k)=wth(:,:,1)-real(k-1)*theta_db;
!$OMP end workshare
!$OMP end parallel
      end do
      return
      end subroutine symm_wkp

      subroutine symm_dwkp
      use m_wake
      use m_dwake
      use m_constant
      implicit none
      integer k
      do k = 2, nbld
!$OMP parallel shared(wpt,wr,wth,k,theta_db)
!$OMP workshare
        dwpt(1,:,:,k)=dwpt(1,:,:,1);
        dwr(:,:,k)=dwr(:,:,1);
        dwth(:,:,k)=dwth(:,:,1)-real(k-1)*theta_db;
!$OMP end workshare
!$OMP end parallel
      end do
      return
      end subroutine symm_dwkp

      subroutine cyl2cart
      use m_wake
      use m_constant
      implicit none
      integer i, j, k
      do k = 1, nbld
!$omp parallel
!$omp&shared(mrpw, ncpw, wpt,wr, wth,k) private(i,j)
!$omp do
        do i = 1,mrpw
          do j = 1,ncpw
            wpt(2,j,i,k) = wr(j,i,k)*cos(wth(j,i,k))
            wpt(3,j,i,k) = wr(j,i,k)*sin(wth(j,i,k))
          end do
        end do
!$omp end do
!$omp end parallel
      end do
      return
      end subroutine cyl2cart

      subroutine cyl2cart_duct
      use m_wake
      use m_dwake
      use m_constant
      implicit none
      integer i, j, k
      do k = 1, nbld
!$omp parallel
!$omp&shared(mrpdw, ncpdw, dwpt,dwr, dwth,k) private(i,j)
!$omp do
        do i = 1,mrpdw
          do j = 1,ncpdw
            dwpt(2,j,i,k) = dwr(j,i,k)*cos(dwth(j,i,k))
            dwpt(3,j,i,k) = dwr(j,i,k)*sin(dwth(j,i,k))
          end do
        end do
!$omp end do
!$omp end parallel
      end do
      return
      end subroutine cyl2cart_duct

      subroutine trans_delponw(delponw,imrw)
      use m_wake
      real delponw(imrw)
      do i = 1, imrw
       phiw0(i) = (- DELPonW(i))*0.6+phiw0(i)*0.4
C       phiw0(i) = (- DELPonW(i))
      end do

      gamw0(1) = phiw0(1)
      do i = 2, mrw
        gamw0(i) = phiw0(i)-phiw0(i-1)
      end do
      gamw0(mrpw) = - phiw0(mrw)

      return
      end subroutine

      subroutine trans_delponw_uns(delponw,imrw)
      use m_wake
      real delponw(imrw)
      do i = 1, imrw
       phiw0(i) = - DELPonW(i)
      end do

      gamw0(1) = phiw0(1)
      do i = 2, mrw
        gamw0(i) = phiw0(i)-phiw0(i-1)
      end do
      gamw0(mrpw) = - phiw0(mrw)
      return
      end subroutine

      subroutine trans_delpondw(delpondw,imrwd)
      use m_wake
      use m_dwake
      real delpondw(imrwd)
      do i = 1, imrwd
       phidw0(i) = (-DELPonDW(i))*0.6+phidw0(i)*0.4
       phidw1(i) = DELPonDW(i)
      end do

      gamdw0(1) = phidw0(1)
      do i = 2, mrdw
        gamdw0(i) = phidw0(i)-phidw0(i-1)
      end do
      gamdw0(mrpdw) = - phidw0(mrdw)

      return
      end subroutine

      subroutine init_wkp2s_propcav(x,y,z,delponw,imrw,JS,DELTAT,ihub)! initialize the wake with circulation
      !cosine sps
      use m_wake
      use m_constant
      !implicit none
      !INCLUDE 'PUFCAVB.INC'
      implicit none
      real dt,dr,r,cir_max,th,JS,dl
      real DELTAT
      integer i,imrw,ihub
      real,dimension(imrw+1,2)::x,y,z
      real delponw(imrw)

      delta_w = 0.075
!     delta_w = 0.05
!     delta_w = 0.03
      omega_p = PI/JS
      dth = DELTAT
      dt = dth/omega_p

      dt = dt*2.0
      theta_db = 2*PI/real(nbld)

      do i = 1, imrw+1
        wpt(1,1,i,1)=x(i,1)
        wpt(2,1,i,1)=y(i,1)
        wpt(3,1,i,1)=z(i,1)

        r = sqrt(y(i,1)**2+z(i,1)**2)
        dl=sqrt((r*dth)**2+dt**2)

        wpt(1,2,i,1)=x(i,1)+x(i,2)*dl!/2.0 !Here we get half wake length for several strips, which helps the convergence of cavity runs.   
        wpt(2,2,i,1)=y(i,1)+y(i,2)*dl!/2.0
        wpt(3,2,i,1)=z(i,1)+z(i,2)*dl!/2.0
      end do

      do i = 1, imrw+1
        wth(1,i,1)=atan2(z(i,1),y(i,1))
         wr(1,i,1)=sqrt(z(i,1)*z(i,1)+y(i,1)*y(i,1))
        wth(2,i,1)=atan2(wpt(3,2,i,1),wpt(2,2,i,1))
         wr(2:ncpw,i,1)=sqrt(wpt(3,2,i,1)**2+wpt(2,2,i,1)**2)
      end do
      do i = 2, ncw
        if (i+1.lt.6) then 
          wth(i+1,:,1) = wth(i,:,1) + dth!/2.0
        else 
          wth(i+1,:,1) = wth(i,:,1) + dth
        endif
      end do
      do i = 2, ncw
        if (i+1.lt.6) then 
          wpt(1,i+1,:,1) = wpt(1,i,:,1) + dt!/2.0
        else
          wpt(1,i+1,:,1) = wpt(1,i,:,1) + dt
        endif
      end do

      do i = 1, imrw
        phiw0(i) = - DELPonW(i)
      end do

      gamw0(1) = phiw0(1)
      do i = 2, mrw
        gamw0(i) = phiw0(i)-phiw0(i-1)
      end do
      gamw0(mrpw) = - phiw0(mrw)
      call symm_wkp
      call cyl2cart ! After this, all the wpt(:,:,:,:) are set.
      if (ihub.eq.6) call podgwake_fwa
      end subroutine init_wkp2s_propcav

      subroutine init_dwkp2s_propcav(x2,y2,z2,delpondw,imrwd,JS,DELTAT)! initialize the wake with circulation
      !cosine sps
      use m_wake
      use m_dwake
      use m_constant
      !implicit none
      !INCLUDE 'PUFCAVB.INC'
      implicit none
      real dt,dr,r,cir_max,th,JS,dl
      real DELTAT
      integer i,imrwd
      real,dimension(imrwd+1,2)::x2,y2,z2
      real delpondw(imrwd)

      delta_w = 0.075
!     delta_w = 0.05
!     delta_w = 0.03
      omega_p = PI/JS
      dth = DELTAT
      dt = dth/omega_p

      dt = dt*2.0
      theta_db = 2*PI/real(nbld)

      do i = 1, imrwd+1
        dwpt(1,1,i,1)=x2(i,1)
        dwpt(2,1,i,1)=y2(i,1)
        dwpt(3,1,i,1)=z2(i,1)

        r = sqrt(y2(i,1)**2+z2(i,1)**2)
        dl=sqrt((r*dth)**2+dt**2)

        dwpt(1,2,i,1)=x2(i,1)+x2(i,2)*dl!/2.0
        dwpt(2,2,i,1)=y2(i,1)+y2(i,2)*dl!/2.0
        dwpt(3,2,i,1)=z2(i,1)+z2(i,2)*dl!/2.0
      end do

      do i = 1, imrwd+1
        dwth(1,i,1)=atan2(z2(i,1),y2(i,1))
         dwr(1,i,1)=sqrt(z2(i,1)*z2(i,1)+y2(i,1)*y2(i,1))
        dwth(2,i,1)=atan2(dwpt(3,2,i,1),dwpt(2,2,i,1))
         dwr(2:ncpdw,i,1)=sqrt(dwpt(3,2,i,1)**2+dwpt(2,2,i,1)**2)
      end do

      do i = 2, ncdw
        if (i+1.lt.6) then ! 3, 4, 5
          dwth(i+1,:,1) = dwth(i,:,1) + dth!/2.0
        else ! 6, 7, 8...
          dwth(i+1,:,1) = dwth(i,:,1) + dth
        endif
      end do
      do i = 2, ncdw
        if (i+1.lt.6) then ! 3, 4, 5
          dwpt(1,i+1,:,1) = dwpt(1,i,:,1) + dt!/2.0
        else ! 6, 7, 8...
          dwpt(1,i+1,:,1) = dwpt(1,i,:,1) + dt
        endif
      end do

      do i = 1, imrwd
        phidw0(i) = - DELPonDW(i)
        phidw1(i) = DELPonDW(i)
      end do

      gamdw0(1) = phidw0(1)
      do i = 2, mrdw
        gamdw0(i) = phidw0(i)-phidw0(i-1)
      end do
      gamdw0(mrpdw) = - phidw0(mrdw)
      call symm_dwkp
      call cyl2cart_duct ! After this, all the wpt(:,:,:,:) are set.

      end subroutine init_dwkp2s_propcav

! =======================================================================
      subroutine vseg_dfast(xp, x1, x2, v, inf, delta)
      use m_constant
      implicit none
      real , intent(in)   :: xp(dm), x1(dm), x2(dm)
      integer, intent(in) :: inf
      real,    intent(out):: v(dm)
      real,    intent(in) :: delta
      real,  parameter    :: radius = 5.0
      real,  parameter  :: tol =    1.0e-8

      real :: r0(dm), r1(dm), r2(dm), r1xr2(dm), fvec(dm), h(dm)
      real :: lr0,l2r0, lr1, lr2, l2r1xr2, tmp, f

      r0 = x2 - x1
      l2r0= dot_product(r0,r0)
      lr0= sqrt(l2r0)
      fvec =HALF * (x2+x1)-xp
      f=sqrt(dot_product(fvec,fvec))

      if((f/lr0.ge.radius).and.(inf .ne. 0)) then
        tmp=HALF/f**3
        h(1) =r0(3)*fvec(2)-r0(2)*fvec(3)
        h(2) =r0(1)*fvec(3)-r0(3)*fvec(1)
        h(3) =r0(2)*fvec(1)-r0(1)*fvec(2)
        v=tmp*h
        return
      endif


      r1 = xp - x1
      lr1= sqrt(dot_product(r1,r1))
      r2 = xp - x2
      lr2= sqrt(dot_product(r2,r2))
      r1xr2(1) =r2(3)*r1(2)-r2(2)*r1(3)
      r1xr2(2) =r2(1)*r1(3)-r2(3)*r1(1)
      r1xr2(3) =r2(2)*r1(1)-r2(1)*r1(2)


      l2r1xr2= dot_product(r1xr2,r1xr2)
     &         + delta*delta *l2r0

      if (lr1 .lt. tol) then
        tmp = ZERO
      else
        tmp = dot_product(r1,r0)/lr1
      endif
      if (inf .eq. 1) then
        tmp = tmp + lr0
      else
        if (lr2 .gt. tol) then
          tmp= tmp - dot_product(r2,r0)/lr2
        endif
      endif
      v   = HALF* tmp / l2r1xr2 * r1xr2

      return
      end subroutine vseg_dfast

      subroutine indvel3p
      use m_constant
      use m_wake
      use omp_lib
      implicit none
      integer i, j, k, l, inf, n1, n2, m
      real v(3),delta_w2,x_tmp,vy,vz
      wv = 0.0
      inf = 0
      do n2 = 1, 1
!$omp parallel
!$omp&shared(n2,wpt,wv,gamw0,delta_w,inf)
!$omp& private(k,n1,i,l,j,v,x_tmp,delta_w2)
!$omp do
        do k = 1, mrpw
          do l = 2, ncpw
            do n1 = 1, nbld
              do i = 1, mrpw
                do j = 1, ncw
!                 if ((i.eq.k).and.((j.eq.l).or.(j+1.eq.l)))cycle
                  if ((n1.eq.1).and.(i.eq.k)
     &            .and.((j.eq.l).or.(j+1.eq.l)))cycle
                  x_tmp = 0.5*(wpt(1,j,i,n1)+wpt(1,j+1,i,n1));
                  if (x_tmp .gt. 0.5) then
                    delta_w2 = delta_w + (x_tmp-0.5)*0.05
                  else
                    delta_w2 = delta_w
                  endif
                  call vseg_dfast(wpt(:,l,k,n2),
     &                            wpt(:,j,i,n1),
     &                            wpt(:,j+1,i,n1),
     &                            v, inf,
     &                     delta_w2)
!    &                     delta_w+(0.1*real(j-1)/real(ncw-1)))
                  wv(:,l,k,n2) = wv(:,l,k,n2)+v*gamw0(i)/TWOPI
                end do
              end do
            end do
          end do
        end do
!$omp end do
!$omp end parallel
      end do

      wv(1,:,:,1) = wv(1,:,:,1) + 1.0
!$omp parallel
!$omp&shared(n2,wpt,wv,gamw0,delta_w,inf) private(k,n1,i,l,j,vy,vz)
!$omp do
      do l = 1, mrpw
        do k = 1, ncw
!         vy = 0.5*(wv(2,k,l,1)+wv(2,k+1,l,1))
!         vz = 0.5*(wv(3,k,l,1)+wv(3,k+1,l,1))
!         wqr(k,l,1) = vy*cos(wth(k,l,1)+0.5*dth) +
!    &                 vz*sin(wth(k,l,1)+0.5*dth)
!         wqt(k,l,1) = (-vy*sin(wth(k,l,1)+0.5*dth) +
!    &                 vz*cos(wth(k,l,1)+0.5*dth) ) /
!    &                 wr(k,l,1) + omega_p
          wqr(k,l,1) = wv(2,k,l,1)*cos(wth(k,l,1)) +
     &                 wv(3,k,l,1)*sin(wth(k,l,1))
          wqt(k,l,1) = (-wv(2,k,l,1)*sin(wth(k,l,1)) +
     &                 wv(3,k,l,1)*cos(wth(k,l,1)) ) /
     &                 wr(k,l,1) + omega_p
!         wv(2,k,l,1) = wv(2,k,l,1)-omega_p*wpt(3,k,l,1)
!         wv(3,k,l,1) = wv(3,k,l,1)+omega_p*wpt(2,k,l,1)
        end do
      end do
!$omp end do
!$omp end parallel
      return
      end subroutine indvel3p

      subroutine align_wkps(beta,beta2) !the symmetrical version of align_wkp
      use m_wake
      implicit none
      integer i, j, nal,nend
      real vx,vr,vt,v(3),v1(3),vall,pt(3), dx, dl,dt,dt0
      real xtp, beta, beta2
      dt0 = dth/omega_p
!     nend = min(nend,ncpw)
!$omp parallel
!$omp&shared(wpt,wv,dt0,dth,err,beta)
!$omp&private(i,j,vx,vr,vt,dl,v,dt,xtp)
!$omp do
      do i = 1, mrpw
!       do j = ncpw,2,-1
        do j = 2,ncpw
!       do j = 2,nend
!         v = 0.6*wv(:,j-1,i,1)
!    &        +0.4*wv(:,j,i,1)
          vx = beta2*wv(1,j-1,i,1)
     &        +(1.0-beta2)*wv(1,j,i,1)
          vt = beta2*wqt(j-1,i,1)
     &        +(1.0-beta2)*wqt(j,i,1)
          vr = beta2*wqr(j-1,i,1)
     &        +(1.0-beta2)*wqr(j,i,1)
!         vx =  wv(1,j-1,i,1)
!         vt = wqt(j-1,i,1)
!         vr = wqr(j-1,i,1)
          dt = dth/vt
!         dt = dt0/vx
          vx = vx*dt
!         vt = vt*dt
          vr = vr*dt
          if(wpt(1,j-1,i,1).gt. 1.0) vr = 0.0
          if (j.eq.ncpw) err(i)=wpt(1,j,i,1)
          wpt(1,j,i,1) = (wpt(1,j-1,i,1)+vx)*beta+
     &                   (1.0-beta)* wpt(1,j,i,1)
!         wth(j,i,1)    = (wth(j-1,i,1)+vt)*0.95 +
!    &                   0.05* wth(j,i,1)
          wr(j,i,1)    = (wr(j-1,i,1)+vr)*beta+
     &                   (1.0-beta)* wr(j,i,1)
          if (j.eq.ncpw) err(i)=err(i)-wpt(1,j,i,1)
        end do
      end do
!$omp end do
!$omp end parallel

!      if (nend .lt. ncpw) then
!!$omp parallel
!!$omp&shared(wpt,wr,nend,mrpw,dt0) private(i,j)
!!$omp do
!        do i = 1, mrpw
!          do j = nend+1, ncpw
!            wpt(1,j,i,1)=wpt(1,j-1,i,1)+dt0
!            wr(j,i,1) = wr(nend,i,1)
!          end do
!        end do
!!$omp end do
!!$omp end parallel
!      endif
      err = abs(err)
      err_m = sum(err)/real(mrpw)
      call symm_wkp
      call cyl2cart
      return
      end subroutine align_wkps

      subroutine indvel4p
      use m_constant_s
      use m_wake
      use omp_lib
      implicit none
      integer i, j, k, l, inf, n1, n2, m
      real v(3),delta_w2,x_tmp

      wv = 0.0
      inf = 0
*      do n2 = 1, 1
!$omp parallel
!$omp&shared(n2,wpt,wv,gamw0,delta_w,inf)
!$omp& private(k,n1,i,l,j,v,x_tmp,delta_w2)

!$omp do
        do k = 1, mrpw
          do l = 1, ncpw
            do n1 = 1, nbld
              do i = 1, mrpw
                do j = 1, ncw
!                 if ((i.eq.k).and.((j.eq.l).or.(j+1.eq.l)))cycle
                  if ((n1.eq.1).and.(i.eq.k)
     &            .and.((j.eq.l).or.(j+1.eq.l)))cycle
                  x_tmp = 0.5*(wpt(1,j,i,n1)+wpt(1,j+1,i,n1));
                  if (x_tmp .gt. 0.5) then
                    delta_w2 = delta_w + (x_tmp-0.5)*0.05
                  else
                    delta_w2 = delta_w    ! delta_w = 0.075
                  endif
                  call vseg_dfast(wpt(:,l,k,1),
     &                            wpt(:,j,i,n1),
     &                            wpt(:,j+1,i,n1),
     &                            v, 0,
     &                     delta_w2)
!    &                     delta_w+(0.1*real(j-1)/real(ncw-1)))
                  wv(:,l,k,1) = wv(:,l,k,1)+v*gamw0(i)/TWOPI
                end do
              end do
            end do
          end do
        end do
!$omp end do
!$omp end parallel
*      end do

C      wv(1,:,:,1) = wv(1,:,:,1) + 1.0 ! 1.0 means inflow velocity
      return
      end subroutine indvel4p

      subroutine indvel4pd
      use m_constant_s
      use m_wake
      use m_dwake
      use omp_lib
      implicit none
      integer i, j, k, l, inf, n1, n2, m
      real v(3),delta_w2,x_tmp
      dwv = 0.0
      inf = 0
!      do n2 = 1, 1
!$omp parallel
!$omp&shared(n2,dwpt,dwv,gamdw0,delta_w,inf)
!$omp& private(k,n1,i,l,j,v,x_tmp,delta_w2)
!$omp do
        do k = 1, mrpdw
          do l = 1, ncpdw
            do n1 = 1, nbld
              do i = 1, mrpdw
                do j = 1, ncdw
!                 if ((i.eq.k).and.((j.eq.l).or.(j+1.eq.l)))cycle
                  if ((n1.eq.1).and.(i.eq.k)
     &            .and.((j.eq.l).or.(j+1.eq.l)))cycle
                  x_tmp = 0.5*(dwpt(1,j,i,n1)+dwpt(1,j+1,i,n1));
                  if (x_tmp .gt. 0.5) then
                    delta_w2 = delta_w + (x_tmp-0.5)*0.05
                  else
                    delta_w2 = delta_w    ! delta_w = 0.075
                  endif
                  call vseg_dfast(dwpt(:,l,k,1),
     &                            dwpt(:,j,i,n1),
     &                            dwpt(:,j+1,i,n1),
     &                            v, 0,
     &                     delta_w2)
!    &                     delta_w+(0.1*real(j-1)/real(ncw-1)))
                  dwv(:,l,k,1) = dwv(:,l,k,1)-v*gamdw0(i)/TWOPI
                end do
              end do
            end do
          end do
        end do
!$omp end do
!$omp end parallel
!      end do

C      dwv(1,:,:,1) = dwv(1,:,:,1) !+ 1.0 ! 1.0 means inflow velocity

      return
      end subroutine indvel4pd

      subroutine indvel4pU(nwpanel,mr)
      use m_constant_s
      use m_wake
      use omp_lib
      use GPRODW
      implicit none
      integer nwpanel, mr
      integer i, j, k, l, inf, n1, n, m
      real v(3),delta_w2,x_tmp
      real xww_tmp(mr+1),yww_tmp(mr+1),zww_tmp(mr+1)
      real xww_tmp2(mrpw),yww_tmp2(mrpw),zww_tmp2(mrpw)

!$OMP parallel
!$OMP&shared(xww,yww,zww,wpt)
!$OMP&private(n1,n,m,xww_tmp,yww_tmp,zww_tmp,xww_tmp2,yww_tmp2,zww_tmp2)
!$OMP DO
      do n1 = 1, nbld
        do n = 1, nwpanel + 1
          do m = 1, mr + 1
            xww_tmp(m) = xww(n,m,n1)
            yww_tmp(m) = yww(n,m,n1)
            zww_tmp(m) = zww(n,m,n1)
          enddo
          call wk_interp2(xww_tmp,yww_tmp,zww_tmp,mr+1,
     *                    xww_tmp2,yww_tmp2,zww_tmp2,mrw+1)
          do m = 1, mrw+1
            wpt(1,n,m,n1) = xww_tmp2(m)
            wpt(2,n,m,n1) = yww_tmp2(m)
            wpt(3,n,m,n1) = zww_tmp2(m)
          enddo  
        enddo
      enddo
!$OMP END DO
!$OMP END parallel

C---------------------------------------
C/s S.N.KIM - for FWA unsteady alignment
!$OMP parallel
!$OMP&shared(mrpw,ncpw,wpt,wr,wth) private(n,m)
!$OMP DO
      do m = 1, mrpw
        do n = 1, ncpw
          wr(n,m,1) = sqrt(wpt(2,n,m,1)**2 + wpt(3,n,m,1)**2)
          wth(n,m,1) = atan2(wpt(3,n,m,1),wpt(2,n,m,1))
        enddo
      enddo
!$OMP END DO
!$OMP END parallel
C/s S.N.KIM - for FWA unsteady alignment
C---------------------------------------

      wv = 0.0
      inf = 0
C SPANWISE DIRECTION 
!$omp parallel
!$omp&shared(wpt,wv,gamw1,delta_w,inf)
!$omp& private(k,n1,i,l,j,v,x_tmp,delta_w2)
!$omp do
        do k = 1, mrpw
          do l = 1, nwpanel+1
            do n1 = 1, nbld
              do i = 1, mrpw
                do j = 1, nwpanel
!                 if ((i.eq.k).and.((j.eq.l).or.(j+1.eq.l)))cycle
                  if ((n1.eq.1).and.(i.eq.k)
     &            .and.((j.eq.l).or.(j+1.eq.l)))cycle
                  x_tmp = 0.5*(wpt(1,j,i,n1)+wpt(1,j+1,i,n1));
                  if (x_tmp .gt. 0.5) then
                    delta_w2 = delta_w + (x_tmp-0.5)*0.05
                  else
                    delta_w2 = delta_w    ! delta_w = 0.075
                  endif
                  call vseg_dfast(wpt(:,l,k,1),
     &                            wpt(:,j,i,n1),
     &                            wpt(:,j+1,i,n1),
     &                            v, 0,
     &                     delta_w2)
!    &                     delta_w+(0.1*real(j-1)/real(ncw-1)))
                  wv(:,l,k,1) = wv(:,l,k,1)+v*gamw1(n1,j,i)/TWOPI
                end do
              end do
            end do
          end do
        end do
!$omp end do
!$omp end parallel

C STREAMWISE DIRECTION (WE DO NOT NEED THIS IN STEADY CASE THOUGH)
!$omp parallel
!$omp&shared(wpt,wv,gamw2,delta_w,inf)
!$omp& private(k,n1,i,l,j,v,x_tmp,delta_w2)
!$omp do
        do k = 1, mrpw
          do l = 1, nwpanel+1
            do n1 = 1, nbld
              do i = 1, mrw
                do j = 2, nwpanel
                  if ((n1.eq.1).and.(i.eq.k)
     &            .and.((j.eq.l).or.(j+1.eq.l)))cycle
                  x_tmp = 0.5*(wpt(1,j,i,n1)+wpt(1,j+1,i,n1));
                  if (x_tmp .gt. 0.5) then
                    delta_w2 = delta_w + (x_tmp-0.5)*0.05
                  else
                    delta_w2 = delta_w    ! delta_w = 0.075
                  endif
                  call vseg_dfast(wpt(:,l,k,1),
     &                            wpt(:,j,i+1,n1),
     &                            wpt(:,j,i,n1),
     &                            v, 0,
     &                     delta_w2)
!    &                     delta_w+(0.1*real(j-1)/real(ncw-1)))
                  wv(:,l,k,1) = wv(:,l,k,1)+v*gamw2(n1,j,i)/TWOPI
                end do
              end do
            end do
          end do
        end do
!$omp end do
!$omp end parallel

      return
      end subroutine indvel4pU

      subroutine align_wkps2(beta,beta2) !the symmetrical version of align_wkp
      use m_wake
      implicit none
      integer i, j, nal
      real vx,vr,vt,v(3),v1(3),vall,pt(3), dx, dl,dt,dt0
      real beta, beta2
      dt0 = dth/omega_p
!$omp parallel
!$omp&shared(wpt,wv,dt0,beta) private(i,j,vx,vr,vt,dl,v,v1,vall,pt,dt)
!$omp do
      do i = 1, mrpw
!       do j = ncpw,2,-1
        do j = 2,ncpw
!         v = wv(:,j-1,i,1)
          v = beta2*wv(:,j-1,i,1)
     &        +(1.0-beta2)*wv(:,j,i,1)
!         vx = 0.6*wv(1,j-1,i,1)
!    &        +0.4*wv(1,j,i,1)
!         vt = 0.6*wqt(j-1,i,1)
!    &        +0.4*wqt(j,i,1)
!         vr = 0.6*wqr(j-1,i,1)
!    &        +0.4*wqr(j,i,1)
!         dl = sqrt(dot_product(v,v))
!         dt =dth/dl
          v1=v
          dl = sqrt((wr(j-1,i,1)*dth)**2+(dt0*1.0)**2)
          v1(2) = v1(2) - omega_p*wpt(3,j-1,i,1)
          v1(3) = v1(3) + omega_p*wpt(2,j-1,i,1)
          vall = sqrt(dot_product(v1,v1))
          dt = dl/vall
          v  = v*dt
          pt(1) = wpt(1,j-1,i,1)
          pt(2) = wr(j-1,i,1)*cos(wth(j-1,i,1)+omega_p*dt)
          pt(3) = wr(j-1,i,1)*sin(wth(j-1,i,1)+omega_p*dt)
          pt = pt + v
!         v(2) = v(2) - wpt(3,j-1,i,1)*dth
!         v(3) = v(3) + wpt(2,j-1,i,1)*dth
!         vx = vx*dt
!         vt = vt*dt
!         vr = vr*dt
!         dt = dth/vt
!         vx = vx*dt
!         vr = vr*dt
          if (j.eq.ncpw) err(i)=wpt(1,j,i,1)
          wpt(:,j,i,1) = pt*beta +
     &                   (1.0-beta)* wpt(:,j,i,1)
          if (j.eq.ncpw) err(i)=err(i)-wpt(1,j,i,1)
        end do
      end do
!$omp end do
!$omp end parallel
      err = abs(err)
      err_m = sum(err)/real(mrpw)

      call symm_wkp2
      call cyl2cart
      return
      end subroutine align_wkps2

      subroutine align_wkps3(ihub,idfwa) !the symmetrical version of align_wkp
      use m_wake
      use m_dwake
      implicit none
      integer i, j, ihub, idfwa
      real dtstar, vall
      real v(3),v1(3),v2(3),pt(3),dl,dt,dt0
      real beta3, epsil

      dt0 = dth/omega_p  ! dth is in 'm_wake'
                         ! dth = DELTAT
      if(idfwa.eq.1) then
        epsil = 0.25 !Under-relaxation factor for ducted propellers.
      else
        epsil = 0.50 !Under-relaxation factor for open propellers.
      endif

      call inflow2fwa

C-------------------------------------------------------------------------------------
C   S.N. KIM  |  This part is the essence of the full wake alignment. Initial
C  ---------- |  version of this scheme before PROPCAV V3.3 was based on the several 
C   Aug-2018  |  simplifications using constant values for some parameters, 
C             |  therefore it was a lot different from what 
C             |  was meant to be by the equations. I modified most of this part 
C             |  to fully incorporate the equations derived by Dr. Y. Tian, with 
C             |  the corrected dtstar which was in typo in his journal.
C             |  ---------------------------------------------------------------------
C             |  Reference : Tian, Y.; Kinnas, S. A.         Int. J. Rotat. 2012  
C             |              Kim, S.; Kinnas, S. A.; Du, W.  JMSE. 2018
C-------------------------------------------------------------------------------------
!$omp parallel
!$omp&shared(wpt,wv,wvi,dt0)
!$omp&private(i,j,dl,v,v1,v2,vall,pt,dt,beta3,dtstar)
!$omp do
      do i = 1, mrpw
        do j = ncpw, 2, -1

          dt = dt0

          v = (0.5*wv(:,j-1,i,1) + 0.5*wv(:,j,i,1))
          v1 = (0.5*wvi(:,j-1,i) + 0.5*wvi(:,j,i))*dt

          vall = sqrt(dot_product(v1,v1))
          v1 = v1/vall

          dl = dot_product(v,v1)
          v2 = v - v1*dl

          pt(1) = wpt(1,j-1,i,1) + v1(1)*vall
          pt(2) = wpt(2,j-1,i,1) + v1(2)*vall
          pt(3) = wpt(3,j-1,i,1) + v1(3)*vall

          dtstar = vall/(vall/dt + dl)
          beta3 = dt/dtstar

          if (j.eq.ncpw) err(i) = wpt(1,j,i,1)
          wpt(:,j,i,1) = (pt * beta3 + v2 * dt +
     &             (1.0 - beta3) * wpt(:,j,i,1)) * epsil +  
     &             wpt(:,j,i,1) * (1.0 - epsil) ! 0.25 above means under relaxation factor.
          if (j.eq.ncpw) err(i) = err(i) - wpt(1,j,i,1)
        end do
      end do
!$omp end do
!$omp end parallel
C-------------------------------------------------------------------------------------

      err = abs(err)
      err_m = sum(err)/real(mrpw)

      if (ihub.eq.6) call bldwk2hub_pen
      call symm_wkp2
      call cyl2cart
      return
      end subroutine align_wkps3

      subroutine align_wkps3_duct !the symmetrical version of align_wkp
      use m_wake
      use m_dwake
      implicit none
      integer i, j
      real dtstar,vall
      real v(3),v1(3),v2(3),pt(3),dl,dt,dt0
      real beta3

      dt0 = dth/omega_p  ! dth is in 'm_wake'
                         ! dth = DELTAT

      CALL inflow2fwa_duct

C-------------------------------------------------------------------------------------
C   S.N. KIM  |  This part is written for the ductwake alignment using FWA.  
C   --------  |  --------------------------------------------------------------------- 
C   Aug-2018  |  Reference : Tian, Y.; Kinnas, S. A.         Int. J. Rotat. 2012  
C             |              Kim, S.; Kinnas, S. A.; Du, W.  JMSE. 2018
C-------------------------------------------------------------------------------------
!$omp parallel
!$omp&shared(dwpt,dwv,dt0,dwvi)
!$omp&private(i,j,dl,v,v1,v2,vall,pt,dt,beta3,dtstar)
!$omp do
      do i = 1, mrpdw
        do j = ncpdw, 2, -1

          dt = dt0

          v = 0.5*dwv(:,j-1,i,1) + 0.5*dwv(:,j,i,1)
          v1 = (0.5*dwvi(:,j-1,i) + 0.5*dwvi(:,j,i))*dt

          vall = sqrt(dot_product(v1,v1))
          v1 = v1/vall

          dl = dot_product(v,v1)
          v2 = v - v1*dl

          pt(1) = dwpt(1,j-1,i,1) + v1(1)*vall
          pt(2) = dwpt(2,j-1,i,1) + v1(2)*vall
          pt(3) = dwpt(3,j-1,i,1) + v1(3)*vall

          dtstar = vall/(vall/dt + dl)
          beta3 = dt/dtstar

          if (j.eq.ncpdw) derr(i) = dwpt(1,j,i,1)
          dwpt(:,j,i,1) = (v2 * dt + pt * beta3 +
     &                    (1.0 - beta3) * dwpt(:,j,i,1)) * 0.25 +
     &                    dwpt(:,j,i,1) * 0.75
          if (j.eq.ncpdw) derr(i) = derr(i) - dwpt(1,j,i,1)
        enddo
      enddo
!$omp end do
!$omp end parallel

      derr = abs(derr)
      derr_m = sum(derr)/real(mrpdw)

      call symm_dwkp2
      call cyl2cart_duct
      return
      end subroutine align_wkps3_duct

      subroutine symm_wkp2
      use m_wake
      use m_constant_s
      implicit none
      integer k,i,j
!$OMP parallel private(j,i) shared(wpt,wr,wth)
!$OMP DO
      do j = 1, mrpw
        do i = 1, ncpw
          wr(i,j,1)=
     &    sqrt(
     &    wpt(2,i,j,1)**2+wpt(3,i,j,1)**2
     &    )
          wth(i,j,1) =
     &    atan2(wpt(3,i,j,1),wpt(2,i,j,1))
        end do
      end do
!$OMP END DO
!$OMP end parallel
      do k = 2, nbld
!$OMP parallel shared(wpt,wr)
!$OMP workshare
        wpt(1,:,:,k)=wpt(1,:,:,1);
        wr(:,:,k)=wr(:,:,1);
        wth(:,:,k)=wth(:,:,1)-real(k-1)*theta_db;
!$OMP end workshare
!$OMP end parallel
      end do
      return
      end subroutine symm_wkp2

      subroutine symm_dwkp2
      use m_wake
      use m_dwake
      use m_constant_s
      implicit none
      integer k,i,j
!$OMP parallel private(j,i) shared(dwpt,dwr,dwth)
!$OMP DO
      do j = 1, mrpdw
        do i = 1, ncpdw
          dwr(i,j,1)=
     &    sqrt(
     &    dwpt(2,i,j,1)**2+dwpt(3,i,j,1)**2
     &    )
          dwth(i,j,1) =
     &    atan2(dwpt(3,i,j,1),dwpt(2,i,j,1))
        end do
      end do
!$OMP END DO
!$OMP end parallel
      do k = 2, nbld
!$OMP parallel shared(dwpt,dwr)
!$OMP workshare
        dwpt(1,:,:,k)=dwpt(1,:,:,1);
        dwr(:,:,k)=dwr(:,:,1);
        dwth(:,:,k)=dwth(:,:,1)-real(k-1)*theta_db;
!$OMP end workshare
!$OMP end parallel
      end do
      return
      end subroutine symm_dwkp2

      subroutine out_wk
      use m_wake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
!     implicit none
      integer i, j, k,n
      open(100, file='wake.dat',status='unknown')
      write(100,*) 'variables = "x","y","z", "u", "v", "w"'
!     write(100,*) 'variables = "x","y","z"'
      do n = 1,nbld
      write(100,*) 'zone datapacking=block i=', nwpanel+1, ',j=', mrpw
      do k = 1, 3
        do i = 1, mrpw
          do j = 1, nwpanel+1
            write(100,*) wpt(k,j,i,n)
          end do
        end do
      end do
      do k = 1, 3
        do i = 1, mrpw
          do j = 1, nwpanel+1
            write(100,*) wv(k,j,i,n)
          end do
        end do
      end do
      end do
      close(100)
      return
      end subroutine out_wk

      subroutine out_bld
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      integer i, j, k,n
      open(100, file='bld.dat',status='unknown')
      write(100,*) 'variables = "x","y","z"'
!     write(100,*) 'variables = "x","y","z"'
      write(100,*) 'zone datapacking=block i=', nc+1, ',j=', mr+1
      do i = 1,mr+1
        do j = 1, nc+1
          write(100,*) xb(j,i)
        end do
      end do
      do i = 1,mr+1
        do j = 1, nc+1
          write(100,*) yb(j,i)
        end do
      end do
      do i = 1,mr+1
        do j = 1, nc+1
          write(100,*) zb(j,i)
        end do
      end do
      close(100)
      return
      end subroutine out_bld

      subroutine get_err(oerr1,oerr2)
      use m_wake
      use m_dwake
      implicit none
      real oerr1, oerr2
      oerr1 = err_m
      oerr2 = derr_m
      return
      end subroutine get_err

      subroutine run_alignp(icav,imrw,imrwd)
      use m_wake
      use m_dwake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
C/s S.N.KIM | FWA on blade wake
      real,allocatable::x(:,:),y(:,:),z(:,:)
      real,allocatable:: xbtmp_o1(:),ybtmp_o1(:),zbtmp_o1(:)
      real,allocatable:: xbtmp_o2(:),ybtmp_o2(:),zbtmp_o2(:)
      real xbtmp1(imrw+1),ybtmp1(imrw+1),zbtmp1(imrw+1)
      real xbtmp2(imrw+1),ybtmp2(imrw+1),zbtmp2(imrw+1)
      real phiwtmp(imrw)
      real vtmp(3), dl
      real beta1,err1,err2
      integer i,nipt,nmiss,icav, imrw, ii, iii
      character*29 :: wkcheck2,filename2
C/e S.N.KIM | FWA on blade wake
C/s S.N.KIM | FWA on duct wake
      real,allocatable::x2(:,:),y2(:,:),z2(:,:)
      real,allocatable::xdtmp_o1(:),ydtmp_o1(:),zdtmp_o1(:)
      real,allocatable::xdtmp_o2(:),ydtmp_o2(:),zdtmp_o2(:)
      real xdtmp1(imrwd+1),ydtmp1(imrwd+1),zdtmp1(imrwd+1)
      real xdtmp2(imrwd+1),ydtmp2(imrwd+1),zdtmp2(imrwd+1)
      real phidwtmp(imrwd)
      integer imrwd,niptd,nmissd
C/e S.N.KIM | FWA on duct wake
C-------------------------------------------------------------
C/s S.N.KIM | FWA on blade wake
      write(*,*) 'icav=',icav,'out of',icavmax,'FWA iterations'
      if(.NOT.allocated(x)) then
        allocate(x(imrw+1,2),y(imrw+1,2),z(imrw+1,2))
        allocate(xbtmp_o1(mr+1),ybtmp_o1(mr+1),zbtmp_o1(mr+1))
        allocate(xbtmp_o2(mr+1),ybtmp_o2(mr+1),zbtmp_o2(mr+1))
      end if
C/e S.N.KIM | FWA on blade wake
C/s S.N.KIM | FWA on duct wake
      if(iduct.eq.1) then
        if(.NOT.allocated(x2)) then
         allocate(x2(imrwd+1,2),y2(imrwd+1,2),z2(imrwd+1,2))
         allocate(xdtmp_o1(mduct+1),ydtmp_o1(mduct+1),zdtmp_o1(mduct+1))
         allocate(xdtmp_o2(mduct+1),ydtmp_o2(mduct+1),zdtmp_o2(mduct+1))
        end if
      endif
C/e S.N.KIM | FWA on duct wake

      do i = 1, mr+1
        xbtmp_o1(i) = xb(1,i)
        ybtmp_o1(i) = yb(1,i)
        zbtmp_o1(i) = zb(1,i)
        xbtmp_o2(i) = xb(nc,i)
        ybtmp_o2(i) = yb(nc,i)
        zbtmp_o2(i) = zb(nc,i)
      end do

      if(iduct.eq.1) then
        do i = 1, mduct + 1
          xdtmp_o1(i) = xd(1,i)
          ydtmp_o1(i) = yd(1,i)
          zdtmp_o1(i) = zd(1,i)
          xdtmp_o2(i) = xd(nduct,i)
          ydtmp_o2(i) = yd(nduct,i)
          zdtmp_o2(i) = zd(nduct,i)
        enddo
      endif

      IF(ICAV.EQ.1) THEN
C----------------------------------------------------------------------
C/s S.N.KIM - FWA on the blade wake.
        if(ian.eq.2) then
          call alloc_wk_mem(imrw,nwpanel,nblade)    !For unsteady, we use exact wake length.
        else
          call alloc_wk_mem(imrw,nwpanel+20,nblade) !For steady, we use wake 20 panels longer to alleviate wake expanding around downstream.
        endif
        call wk_interp1(xbtmp_o1,ybtmp_o1,zbtmp_o1,delponw,mrp,
     &                  xbtmp1,ybtmp1,zbtmp1,phiwtmp,imrw+1)
        call wk_interp2(xbtmp_o2,ybtmp_o2,zbtmp_o2,mrp,
     &                  xbtmp2,ybtmp2,zbtmp2,imrw+1)
        do i = 1, imrw+1
          x(i,1) = xbtmp1(i)
          y(i,1) = ybtmp1(i)
          z(i,1) = zbtmp1(i)

          vtmp(1)=(xbtmp1(i)-xbtmp2(i))
          vtmp(2)=(ybtmp1(i)-ybtmp2(i))
          vtmp(3)=(zbtmp1(i)-zbtmp2(i))
          dl = sqrt(dot_product(vtmp,vtmp))

          if (dl .gt. 1e-5) then
            vtmp = vtmp / dl
            x(i,2) = vtmp(1)
            y(i,2) = vtmp(2)
            z(i,2) = vtmp(3)
          else
            x(i,2) = x(i-1,2)
            y(i,2) = y(i-1,2)
            z(i,2) = z(i-1,2)
          end if
        end do
        call init_wkp2s_propcav(x,y,z,phiwtmp,imrw,ADVCO,DELTAT,ihub) ! initialize blade wake with the blade circulation
C/e S.N.KIM - FWA on the blade wake.
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C/s S.N.KIM - FWA om the duct wake.
        if(iduct.eq.1) then
          call alloc_dwk_mem(imrwd,ndwk+20,nblade)
          call wk_interp1(xdtmp_o1,ydtmp_o1,zdtmp_o1,delpondw,mduct+1,
     &                    xdtmp1,ydtmp1,zdtmp1,phidwtmp,imrwd+1)
          call wk_interp2(xdtmp_o2,ydtmp_o2,zdtmp_o2,mduct+1,
     &                    xdtmp2,ydtmp2,zdtmp2,imrwd+1)
          do i = 1, imrwd+1
            x2(i,1) = xdtmp1(i)
            y2(i,1) = ydtmp1(i)
            z2(i,1) = zdtmp1(i)

            vtmp(1) = (xdtmp1(i) - xdtmp2(i))
            vtmp(2) = (ydtmp1(i) - ydtmp2(i))
            vtmp(3) = (zdtmp1(i) - zdtmp2(i))
            dl = sqrt(dot_product(vtmp,vtmp))

            if (dl .gt. 1e-5) then
              vtmp = vtmp / dl
              x2(i,2) = vtmp(1)
              y2(i,2) = vtmp(2)
              z2(i,2) = vtmp(3)
            else
              x2(i,2) = x2(i-1,2)
              y2(i,2) = y2(i-1,2)
              z2(i,2) = z2(i-1,2)
            endif
          enddo
          call init_dwkp2s_propcav(x2,y2,z2,phidwtmp,imrwd,ADVCO,DELTAT)! initialize duct wake with duct circulation.
        endif
C/e S.N.KIM - FWA on the duct wake.
C----------------------------------------------------------------------

      ELSE

C----------------------------------------------------------------------
C/s S.N.KIM - FWA on the blade wake.
        call wk_interp1(xbtmp_o1,ybtmp_o1,zbtmp_o1,delponw,mrp,
     &                  xbtmp1,ybtmp1,zbtmp1,phiwtmp,imrw+1)
        call trans_delponw(phiwtmp,imrw)
C/e S.N.KIM - FWA on the blade wake.
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C/s S.N.KIM - FWA on the duct wake.
        if(iduct.eq.1) then
          call wk_interp1(xdtmp_o1,ydtmp_o1,zdtmp_o1,delpondw,mduct+1,
     &                    xdtmp1,ydtmp1,zdtmp1,phidwtmp,imrwd+1)
          call trans_delpondw(phidwtmp,imrwd)
        endif
C/e S.N.KIM - FWA on the duct wake.
C----------------------------------------------------------------------
      ENDIF

      nipt = mrpw*ncpw                   ! = (mrw+1)*(ncw+1),   mrw = MR/ALIGN_DIV,     ncw = nwpanel+20.
      if(iduct.eq.1) niptd = mrpdw*ncpdw ! = (mrdw+1)*(ncdw+1), mrdw = MDUCT/ALIGN_DIV, ncdw = ndwk+20.  

      call cal_bld_vort

C**********************************************************************
C              INNER ITERATION OF THE FULL WAKE ALIGNMENT  
C**********************************************************************
      if(iduct.eq.1.and.icavt.eq.1.and.idfwa.ne.1) 
     &  write(*,*) 'Ductwake is assumed to be cylindrical.'
      do i = 1, 40
        call indvel4p 
        if (idfwa.eq.1) call indvel4pd
        call dump_wakgeo
        if (idfwa.eq.1) call dump_dwakgeo
        call velfrombld(nipt,wpx_dumy,wpy_dumy,wpz_dumy
     &                 ,wvxb,wvrb,wvtb,nmiss,ivbmiss,delta_w
     &                 ,iduct,ihub)
        call merge_vel2
        call align_wkps3(ihub,idfwa)
        if (idfwa.eq.1) then
          call veltoduct(niptd,dwpx_dumy,dwpy_dumy,dwpz_dumy
     &                  ,wvxd,wvrd,wvtd,nmissd,ivdmiss,delta_w,idopt)
          call merge_vel2_duct
          call align_wkps3_duct
        endif
        if(iplot.eq.1) call wpt_check_plt(i) ! S.N.KIM - includes plotting option to reduce computing time. 
        call get_err(err1,err2)
        if(idfwa.eq.1) then
          write(*,7005) i,err1,err2
        else
          write(*,7006) i,err1
        endif
        if (err1 .lt. 1e-4 .and. err2 .lt. 1e-4) exit
 7005   FORMAT(1x,'Residuals at ite.=',I2,'  Blade Wake=',E13.6,
     &'  Duct Wake=',E13.6)
 7006   FORMAT(1x,'Residuals at ite.=',I2,'  Blade Wake=',E13.6)
      end do
      call out_wk  ! 'wake.dat' file is generated.
      call out_bld ! 'bld.dat' file is generated.
C**********************************************************************
C/e S.N.KIM - NOTE: 'outer' loops are determined by 'icavt', Aug. 2018.
C**********************************************************************

      end subroutine run_alignp

      subroutine merge_vel
      use m_wake
      implicit none
      integer i, j
      do i = 1, mrpw
        do j = 1, ncpw
         wv(1,j,i,1) = wv(1,j,i,1) + wvxb(j,i)
         wqr(j,i,1)  = wqr(j,i,1)  + wvrb(j,i)
         wqt(j,i,1)  = wqt(j,i,1)  + wvtb(j,i)/wr(j,i,1)
        end do
      end do
      end subroutine merge_vel

      subroutine merge_vel2
      use m_wake
      implicit none
      integer i, j
      real wr_tmp,wt_tmp

      do i = 1, mrpw
        do j = 1, ncpw
         wv(1,j,i,1) = wv(1,j,i,1) + wvxb(j,i)
         wv(2,j,i,1) = wv(2,j,i,1) + wvrb(j,i)*cos(wth(j,i,1))
     &                             - wvtb(j,i)*sin(wth(j,i,1))
         wv(3,j,i,1) = wv(3,j,i,1) + wvrb(j,i)*sin(wth(j,i,1))
     &                             + wvtb(j,i)*cos(wth(j,i,1))
        enddo
      enddo
      end subroutine merge_vel2

      subroutine merge_vel2_duct
      use m_wake
      use m_dwake
      implicit none
      integer i, j
      real wr_tmp,wt_tmp

      do i = 1, mrpdw
        do j = 1, ncpdw
         dwv(1,j,i,1) = dwv(1,j,i,1) + wvxd(j,i) 
         dwv(2,j,i,1) = dwv(2,j,i,1) + wvrd(j,i)*cos(dwth(j,i,1))
     &                - wvtd(j,i)*sin(dwth(j,i,1)) 
         dwv(3,j,i,1) = dwv(3,j,i,1) + wvrd(j,i)*sin(dwth(j,i,1))
     &                + wvtd(j,i)*cos(dwth(j,i,1)) 
        enddo
      enddo
      end subroutine merge_vel2_duct

      subroutine merge_vel2_uns
      use m_wake
      implicit none
      integer i, j
      real wr_tmp,wt_tmp

      call inflowk3_no_wrovs 

      do i = 1, mrpw
        do j = 1, ncpw
         wv(1,j,i,1) = wv(1,j,i,1) + wvxb(j,i) + wvi(1,j,i)
         wv(2,j,i,1) = wv(2,j,i,1) + wvrb(j,i)*cos(wth(j,i,1))
     &                             + wvi(2,j,i)
         wv(3,j,i,1) = wv(3,j,i,1) + wvrb(j,i)*sin(wth(j,i,1))
     &                             + wvi(3,j,i)
        enddo
      enddo
      end subroutine 

      subroutine dump_wakgeo
      use m_wake
      implicit none
      integer i, j
      real beta
      beta=0.3
      do i = 1, mrpw
        do j = 1, ncpw
          wpx_dumy(j,i) = wpt(1,j,i,1)
          wpy_dumy(j,i) = wpt(2,j,i,1)
          wpz_dumy(j,i) = wpt(3,j,i,1)
        end do
      end do
      do i = 1, mrpw
        wpx_dumy(1,i)=(1-beta)*wpx_dumy(1,i)+beta*wpx_dumy(2,i)
        wpy_dumy(1,i)=(1-beta)*wpy_dumy(1,i)+beta*wpy_dumy(2,i)
        wpz_dumy(1,i)=(1-beta)*wpz_dumy(1,i)+beta*wpz_dumy(2,i)
      end do

      return
      end subroutine dump_wakgeo

      subroutine dump_dwakgeo
      use m_wake
      use m_dwake
      implicit none
      integer i, j
      real beta
      beta=0.3
      do i = 1, mrpdw
        do j = 1, ncpdw
          dwpx_dumy(j,i) = dwpt(1,j,i,1)
          dwpy_dumy(j,i) = dwpt(2,j,i,1)
          dwpz_dumy(j,i) = dwpt(3,j,i,1)
        end do
      end do
      do i = 1, mrpdw
        dwpx_dumy(1,i)=(1-beta)*dwpx_dumy(1,i)+beta*dwpx_dumy(2,i)
        dwpy_dumy(1,i)=(1-beta)*dwpy_dumy(1,i)+beta*dwpy_dumy(2,i)
        dwpz_dumy(1,i)=(1-beta)*dwpz_dumy(1,i)+beta*dwpz_dumy(2,i)
      end do
      return
      end subroutine dump_dwakgeo

      subroutine wrtback_wakgeo
      use m_wake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      integer n,m,k
      real ttm1,ttm2
      real, allocatable:: xx(:),yy(:),zz(:) ,RRW(:),ThetaW(:)
      real, allocatable:: xxo(:),yyo(:),zzo(:)
      if (.NOT.allocated(xx)) then
        allocate(xx(mrp),yy(mrp),zz(mrp) ,RRW(mrp),ThetaW(mrp))
        allocate(xxo(mrpw),yyo(mrpw),zzo(mrpw))
      end if

      do n = 2, nwpanel + 1
        do m = 1, mrpw
          xxo(m) = wpt(1,n,m,1)
          yyo(m) = wpt(2,n,m,1)
          zzo(m) = wpt(3,n,m,1)
        end do

        call wk_interp2(xxo,yyo,zzo,mrpw,xx,yy,zz,mrp)

        do m = 1, mrp
           xw(n,m) = xx(m)
           yw(n,m) = yy(m)
           zw(n,m) = zz(m)
        end do
      end do

      if (iduct.eq.1) then
        CALL bldwk_pen757 ! Subroutine that controls the penetration of wake on duct.
      endif

      if (ian.ne.2) then ! for ian=2, wake geo is read-in at 'infwak2.f'.
        do k = 1, nblade
          do m = 1 , mrp
            do n = 1 , nwpanel + 1
              rr = sqrt(yw(n,m)**2+zw(n,m)**2)
              th = atan2(zw(n,m),yw(n,m))-real(k-1)*theta_db
              xww(n,m,k) = xw(n,m)
              yww(n,m,k) = rr*cos(th)
              zww(n,m,k) = rr*sin(th)
            enddo
          enddo
        enddo
      endif

      do m = 1 , mr+1
        NSW(m) = nwpanel + 1
      enddo

      return
      end subroutine wrtback_wakgeo

      subroutine wrtback_ductwakgeo
      use m_wake
      use m_dwake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      integer n,m,k
      real angtmp
      real, allocatable:: xx2(:),yy2(:),zz2(:),rrr(:)
      real, allocatable:: xxo2(:),yyo2(:),zzo2(:)
      if (.NOT.allocated(xx2)) then
        allocate(xx2(mduct+1),yy2(mduct+1),zz2(mduct+1))
        allocate(xxo2(mrpdw),yyo2(mrpdw),zzo2(mrpdw))
        allocate(rrr(mduct+1))
      end if

      do n = 2, ndwk + 1
        do m = 1, mrpdw
          xxo2(m) = dwpt(1,n,m,1)
          yyo2(m) = dwpt(2,n,m,1)
          zzo2(m) = dwpt(3,n,m,1)
        end do

        call wk_interp2(xxo2,yyo2,zzo2,mrpdw,xx2,yy2,zz2,mduct+1)

        do m = 1, mduct+1
           xdw(n,m) = xx2(m)
           ydw(n,m) = yy2(m)
           zdw(n,m) = zz2(m)
        end do
      end do

      return
      end subroutine wrtback_ductwakgeo

      subroutine wk_interp3(n1,x1,y1,z1,r2,x2)
      implicit none
      integer n1
      real x1(n1),y1(n1),z1(n1),r1(n1)
      real x2,r2
      real wcx(4*n1-4)
      integer i

      do i=1,n1
        r1(i)=sqrt(y1(i)**2+z1(i)**2)
      enddo

      call uglydk(n1,1,1,r1,x1,0.0,0.0,wcx)
      call evaldks(n1,1,r1,r2,x2,wcx)

      return
      end subroutine wk_interp3

      subroutine wk_interp4(n1,ip,x1,y1,z1,r2,s2)
      implicit none
      integer n1,ip
      real x1(n1),y1(n1),z1(n1)
      real xtmp(ip),ytmp(ip),ztmp(ip),rtmp(ip),stmp(ip)
      real x2,y2,z2,r2,s2
      real wcs(4*ip-4)
      real ds,vtmp(3)
      integer i

      stmp(1)=0.0D0
      do i=1,ip
        xtmp(i)=x1(i)
        ytmp(i)=y1(i)
        ztmp(i)=z1(i)
      enddo

      do i=2,ip
        vtmp(1) = xtmp(i)-xtmp(i-1)
        vtmp(2) = ytmp(i)-ytmp(i-1)
        vtmp(3) = ztmp(i)-ztmp(i-1)
        ds = sqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
        stmp(i) = stmp(i-1) + ds
      enddo

      do i=1,ip
        rtmp(i)=sqrt(ytmp(i)**2+ztmp(i)**2)
      enddo

      call uglydk(ip,1,1,rtmp,stmp,0.0,0.0,wcs)
      call evaldks(ip,1,rtmp,r2,s2,wcs)

      return
      end subroutine wk_interp4


      subroutine wk_interp5(n1,x1,y1,z1,ss2,
     &                        x3,y3,z3)
      implicit none
      integer n1
      real x1(n1),y1(n1),z1(n1),s1(n1),ss1
      real ss2,s2(n1)
      real x3(n1),y3(n1),z3(n1)
      real wcs(4*n1-4),wcx(4*n1-4),wcy(4*n1-4),wcz(4*n1-4)
      integer i
      real ds,vtmp(3)

      s1(1) = 0.0d0
      do i = 2,n1
        vtmp(1) = x1(i)-x1(i-1)
        vtmp(2) = y1(i)-y1(i-1)
        vtmp(3) = z1(i)-z1(i-1)
        ds = sqrt(dot_product(vtmp,vtmp))
        s1(i) = s1(i-1) + ds
      enddo
      ss1=s1(n1)

      do i=1,n1
        s2(i)=s1(i)*ss2/ss1
      enddo

      call uglydk(n1,1,1,s1,x1,0.0,0.0,wcx)
      call uglydk(n1,1,1,s1,y1,0.0,0.0,wcy)
      call uglydk(n1,1,1,s1,z1,0.0,0.0,wcz)
      call evaldk(n1,n1,s1,s2,x3,wcx)
      call evaldk(n1,n1,s1,s2,y3,wcy)
      call evaldk(n1,n1,s1,s2,z3,wcz)

      return
      end subroutine wk_interp5

      subroutine wpt_check_plt(ii)
      use m_wake
      use m_dwake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
!     implicit none
      integer i, j, k, n
      integer ii
      real blade_angle
      real veldotpanel(3), wcosa
      real dveldotpanel(3), dwcosa
      real, allocatable :: wvel(:,:,:,:)
      real, allocatable :: dwvel(:,:,:,:)
      CHARACTER*30 :: WPTCHECK, FILENAME, FILENAME2

      IF (.NOT.allocated(wvel)) THEN
        ALLOCATE(wvel(3,ncpw,mrpw,nbld))
        ALLOCATE(dwvel(3,ncpdw,mrpdw,nbld))
      ENDIF

      IF(ii.EQ.1) THEN
        WRITE(FILENAME2, '(I3.3)') icavt
        WPTCHECK='fwa-'//trim(FILENAME2)//'ite.plt'
        open(100+icavt,file=WPTCHECK,status='unknown')
        write(100+icavt,*) 'variables = "x","y","z","u","v","w"
     *,"angle(deg)"'
      ENDIF

      do n = 1,nbld
        blade_angle = -real(n-1)*delk
        write(100+icavt,*)'zone datapacking=block i=', ncpw, ',j='
     *, mrpw
        write(100+icavt,*) 'SOLUTIONTIME= ',ii
        do k = 1, 3
          do i = 1, mrpw
            do j = 1, ncpw
              write(100+icavt,*) wpt(k,j,i,n)
            end do
          end do
        end do
        do k = 1, 3
          do i = 1, mrpw
            do j = 1, ncpw
              wvel(k,j,i,n) = 0.5*(wv(k,j,i,1) + wvi(k,j,i) +
     *                             wv(k,j+1,i,1) + wvi(k,j+1,i))
              if(n.gt.1) then
                if(k.eq.1) then
                  wvel(k,j,i,n) = wvel(k,j,i,1)
                elseif(k.eq.2) then
                  wvel(k,j,i,n) = wvel(k,j,i,1)*cos(blade_angle)
     *                          - wvel(k+1,j,i,1)*sin(blade_angle)
                elseif(k.eq.3) then
                  wvel(k,j,i,n) = wvel(k-1,j,i,1)*sin(blade_angle)
     *                          + wvel(k,j,i,1)*cos(blade_angle)
                endif
              endif
              write(100+icavt,*) wvel(k,j,i,n)
            end do
          end do
        end do
        do i = 1, mrpw
          do j = 1, ncpw
            veldotpanel(1) 
     *                  = dot_product(wvel(:,j,i,n)
     *                               ,wpt(:,j+1,i,n) - wpt(:,j,i,n))
            veldotpanel(2)
     *                  = sqrt(dot_product(wvel(:,j,i,n)
     *                                    ,wvel(:,j,i,n)))
            veldotpanel(3)
     *                  = sqrt(dot_product(wpt(:,j+1,i,n)-wpt(:,j,i,n)
     *                                    ,wpt(:,j+1,i,n)-wpt(:,j,i,n)))
            wcosa = veldotpanel(1)/(veldotpanel(2)*veldotpanel(3))
            if (wcosa .gt. 0.99999) then
              write(100+icavt,*) 0.0 ! In order to avoid NaN value.
            else
              write(100+icavt,*) acos(wcosa)*(180.0/acos(-1.0))
            endif
          enddo
        enddo

        IF(iduct.eq.1) THEN
          write(100+icavt,*)'zone datapacking=block i=', ncpdw, ',j='
     *, mrpdw
          write(100+icavt,*) 'SOLUTIONTIME= ',ii
          do k = 1, 3
            do i = 1, mrpdw
              do j = 1, ncpdw
                write(100+icavt,*) dwpt(k,j,i,n)
              enddo
            enddo
          enddo
          do k = 1, 3
            do i = 1, mrpdw
              do j = 1, ncpdw
                dwvel(k,j,i,n) = 0.5*(dwv(k,j,i,1) + dwvi(k,j,i) +
     *                                dwv(k,j+1,i,1) + dwvi(k,j+1,i))
                if(n.gt.1) then
                  if(k.eq.1) then
                    dwvel(k,j,i,n) = dwvel(k,j,i,1)
                  elseif(k.eq.2) then
                    dwvel(k,j,i,n) = dwvel(k,j,i,1)*cos(blade_angle)
     *                             - dwvel(k+1,j,i,1)*sin(blade_angle)
                  elseif(k.eq.3) then
                    dwvel(k,j,i,n) = dwvel(k-1,j,i,1)*sin(blade_angle)
     *                             + dwvel(k,j,i,1)*cos(blade_angle)
                  endif
                endif
                if(idfwa.ne.1) dwvel(k,j,i,n) = 0.0
                write(100+icavt,*) dwvel(k,j,i,n)
              enddo
            enddo
          enddo
          do i = 1, mrpdw
            do j = 1, ncpdw
              dveldotpanel(1) 
     *                = dot_product(dwvel(:,j,i,n)
     *                             ,dwpt(:,j+1,i,n)-dwpt(:,j,i,n))
              dveldotpanel(2)
     *                = sqrt(dot_product(dwvel(:,j,i,n)
     *                                  ,dwvel(:,j,i,n)))
              dveldotpanel(3)
     *                = sqrt(dot_product(dwpt(:,j+1,i,n)-dwpt(:,j,i,n)
     *                                  ,dwpt(:,j+1,i,n)-dwpt(:,j,i,n)))
              dwcosa = dveldotpanel(1)/(dveldotpanel(2)*dveldotpanel(3))
              if(idfwa.ne.1) dwcosa = 1.0
              if (dwcosa .gt. 0.9999) then
                write(100+icavt,*) 0.0 ! In order to avoid NaN value.
              else
                write(100+icavt,*) acos(dwcosa)*(180.0/acos(-1.0))
              endif
            enddo
          enddo
        ENDIF
      enddo

      IF(ii.EQ.40) THEN
        close(100+icavt)
      ENDIF

      return
      end subroutine wpt_check_plt

      subroutine bldwk_pen757
      use m_wake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      integer n1,n2,i,j
      integer ndte1
      real,allocatable:: xtmp1(:),ytmp1(:),ztmp1(:),rtmp1(:)
      real,allocatable:: xtmp2(:),ytmp2(:),ztmp2(:)
      real,allocatable:: xid1(:),etad1(:),dcx1(:)
      real rout1,rout2,xxx,angtmp,thtmp1
      real,allocatable::xx1(:),yy1(:),zz1(:),rrr(:)
      real,allocatable::stmp2(:)
      real,allocatable::delx(:,:)
      real,allocatable::xx(:),yy(:),zz(:)

      n1=nddat/2+1
      n2=nwpanel+1

      if (.NOT.allocated(xtmp1)) then
        allocate(xtmp1(n2),ytmp1(n2),ztmp1(n2),rtmp1(n2))
        allocate(xtmp2(mrp),ytmp2(mrp),ztmp2(mrp))
        allocate(xid1(n1),etad1(n1),dcx1(4*(n1-1)))
        allocate(xx1(mrp),yy1(mrp),zz1(mrp),rrr(mrp),stmp2(n2))
        allocate(delx(n2,mrp))
        allocate(xx(mrpw),yy(mrpw),zz(mrpw))
      endif

      do j = mrp, mrp-3, -1
        do i=1,n2
          xtmp1(i)=xw(i,j)
          ytmp1(i)=yw(i,j)
          ztmp1(i)=zw(i,j)
          rtmp1(i)=sqrt(ytmp1(i)**2+ztmp1(i)**2)
        enddo

        do i=1,n1
          xid1(i)=xid(n1-i+1)
          etad1(i)=etad(n1-i+1)
        enddo
        call uglydk(n1,1,1,xid1,etad1,0,0,dcx1)

        do i=2,n2
          if (xtmp1(i-1).lt.xid(1).and.xtmp1(i).ge.xid(1)) ndte1=i
        enddo

        do i=2,ndte1-1
          call evaldks(n1,1,xid1,xtmp1(i),rout1,dcx1)
          if (rtmp1(i).gt.rout1) then
            thtmp1 = DANGLE(ztmp1(i),ytmp1(i))
            yw(i,j) = rout1 * cos(thtmp1)
            zw(i,j) = rout1 * sin(thtmp1)
          endif
        enddo
      enddo

C/s S.N.KIM |  In case the cylindrical ductwake is assumed, blade wake
C           |  after duct trailing edge needs to be forced to be 
C           |  located inside the cylinder.
      if(idfwa.ne.1) then
        do i=1,mrp
          rrr(i)=sqrt(yw(ndte1-1,i)**2+zw(ndte1-1,i)**2)
        enddo

        do i=ndte1,n2
          do j=1,mrp
            angtmp=atan2(zw(i,j),yw(i,j))
            yw(i,j)=rrr(j)*cos(angtmp)
            zw(i,j)=rrr(j)*sin(angtmp)
          enddo
        enddo
      endif
C/e S.N.KIM | Aug. 2018

      do i=2,n2
        do j=1,mrp
          xtmp2(j) = xw(i,j)
          ytmp2(j) = yw(i,j)
          ztmp2(j) = zw(i,j)
        enddo
        call wk_interp2(xtmp2,ytmp2,ztmp2,mrp,xx,yy,zz,mrpw)
        do j=1,mrpw
          wpt(1,i,j,1) = xx(j)
          wpt(2,i,j,1) = yy(j)
          wpt(3,i,j,1) = zz(j)
        enddo
      enddo

      call symm_wkp2
      call cyl2cart

      return
      end subroutine bldwk_pen757

      subroutine inflow2fwa
      USE m_wake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION U(3),NH1(3)

      IF(NTSTEP.EQ.0) THEN
        DO IH=1,3
          NH1(IH)=1
        ENDDO
      ELSE
        DO IH=1,3
          NH1(IH) = NHARM(IH)
        ENDDO
      ENDIF

      DO I=1,mrpw
        DO J=1,ncpw
          RCP=SQRT(wpt(2,J,I,1)**2 + wpt(3,J,I,1)**2)
          THP=ATAN2(wpt(3,J,I,1),wpt(2,J,I,1))
          DO IH=1,3
            U(IH)=0.0
            DO JH=1,NH1(IH)
              CALL EVALDKs(NWKCOE,1,XRW,RCP,COEC,XWCUB(1,JH,1,IH))
              CALL EVALDKs(NWKCOE,1,XRW,RCP,COES,XWCUB(1,JH,2,IH))
              T=FLOAT(JH-1)*(THP-TSTEP)
              ST=SIN(T)
              CT1=COS(T)
              U(IH)=U(IH)+COEC*CT1+COES*ST
            ENDDO
          ENDDO
          wvi(1,J,I)=U(1)
          VORW=U(2)
          VOTW=U(3)
          COCP=wpt(2,J,I,1)/RCP
          SICP=wpt(3,J,I,1)/RCP
          WROVS=PI*RCP/ADVCO
          wvi(2,J,I)=VORW*COCP-(WROVS+VOTW)*SICP
          wvi(3,J,I)=VORW*SICP+(WROVS+VOTW)*COCP
        ENDDO
      ENDDO       
 
      RETURN
      END

      subroutine inflow2fwa_duct
      USE m_dwake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION U(3), NH1(3)

      IF(NTSTEP.EQ.0) THEN
        DO IH=1,3
          NH1(IH)=1
        ENDDO
      ELSE
        DO IH=1,3
          NH1(IH)=NHARM(IH)
        ENDDO
      ENDIF

      DO I=1,mrpdw
        DO J=1,ncpdw
          RCP=SQRT(dwpt(2,J,I,1)**2 + dwpt(3,J,I,1)**2)
          THP=ATAN2(dwpt(3,J,I,1),dwpt(2,J,I,1))
          DO IH=1,3
            U(IH)=0.0
            DO JH=1,NH1(IH)
              CALL EVALDKs(NWKCOE,1,XRW,RCP,COEC,XWCUB(1,JH,1,IH))
              CALL EVALDKs(NWKCOE,1,XRW,RCP,COES,XWCUB(1,JH,2,IH))
              T=FLOAT(JH-1)*(THP-TSTEP)
              ST=SIN(T)
              CT1=COS(T)
              U(IH)=U(IH)+COEC*CT1+COES*ST
            ENDDO
          ENDDO
          dwvi(1,J,I)=U(1)
          VORW=U(2)
          VOTW=U(3)
          COCP=dwpt(2,J,I,1)/RCP
          SICP=dwpt(3,J,I,1)/RCP
          WROVS=PI*RCP/ADVCO
          dwvi(2,J,I)=VORW*COCP-(WROVS+VOTW)*SICP
          dwvi(3,J,I)=VORW*SICP+(WROVS+VOTW)*COCP
        ENDDO
      ENDDO

      RETURN
      END

      subroutine bldwk_hub_pen
      use m_wake
      INCLUDE 'PUFCAV.INC' 
      INCLUDE 'PUFCAVB.INC'

      integer i,j
      integer n1,n2
      real,allocatable :: xtmp1(:),ytmp1(:),ztmp1(:),rtmp1(:)
      real,allocatable :: xtmp2(:),ytmp2(:),ztmp2(:)
      real,allocatable :: xx1(:),yy1(:),zz1(:)
      real,allocatable :: stmp2(:)

      n1 = ncpw
      n2 = mrpw

      if (.not.allocated(xtmp1)) then
        allocate(xtmp1(n1),ytmp1(n1),ztmp1(n1),rtmp1(n1))
        allocate(xtmp2(n2),ytmp2(n2),ztmp2(n2))
        allocate(xx1(n2),yy1(n2),zz1(n2))
        allocate(stmp2(n2))
      endif

      do i = 1, ncpw
        xtmp1(i) = wpt(1,i,1,1)
        ytmp1(i) = wpt(2,i,1,1)
        ztmp1(i) = wpt(3,i,1,1)
        rtmp1(i) = sqrt(ytmp1(i)**2 + ztmp1(i)**2)       
      enddo

      do i = 2, ncpw
        CALL EVALDKs(NHIN,1,XHIN,xtmp1(i),RPOD1,XYHINCUB)

        if (RPOD1.GT.rtmp1(i)) then

          xx1 = wpt(1,i,:,1)
          yy1 = wpt(2,i,:,1)
          zz1 = wpt(3,i,:,1)

          rout2 = RPOD1     

          do while (rtmp1(i).lt.rout2)
            rtmp1(i) = rtmp1(i) + 0.001
            call wk_interp3(n2,xx1,yy1,zz1,rtmp1(i),xxx)
            call evaldks(NHIN,1,XHIN,xxx,rout2,XYHINCUB)
          enddo

          call wk_interp4_hub(n2,n2,xx1,yy1,zz1,rtmp1(i),stmp2(i))
          do j = 1, n2
            xx1(j) = wpt(1,i,n2-j+1,1)
            yy1(j) = wpt(2,i,n2-j+1,1)
            zz1(j) = wpt(3,i,n2-j+1,1)
          enddo
          call wk_interp5(n2,xx1,yy1,zz1,stmp2(i),
     &                    xtmp2,ytmp2,ztmp2)
          do j = 1, n2
            xx1(j) = xtmp2(n2-j+1)
            yy1(j) = ytmp2(n2-j+1)
            zz1(j) = ztmp2(n2-j+1)

            wpt(1,i,j,1) = xx1(j)
            wpt(2,i,j,1) = yy1(j)
            wpt(3,i,j,1) = zz1(j)
          enddo

        endif
      enddo
        
      return
      end subroutine bldwk_hub_pen  

      subroutine wk_interp4_hub(n1,ip,x1,y1,z1,r2,s2)
      implicit none
      integer n1,ip
      real x1(n1),y1(n1),z1(n1)
      real xtmp(ip),ytmp(ip),ztmp(ip),rtmp(ip),stmp(ip)
      real x2,y2,z2,r2,s2
      real wcs(4*ip-4)
      real ds,vtmp(3)
      integer i

      stmp(1)=0.0D0
      do i=1,ip
         xtmp(i)=x1(i)
         ytmp(i)=y1(i)
         ztmp(i)=z1(i)
      enddo

      do i=2,ip
        vtmp(1) = xtmp(i)-xtmp(i-1)
        vtmp(2) = ytmp(i)-ytmp(i-1)
        vtmp(3) = ztmp(i)-ztmp(i-1)
        ds = sqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
        stmp(i) = stmp(i-1) + ds
      enddo

      do i=1,ip
        rtmp(i)=sqrt(ytmp(i)**2+ztmp(i)**2)
      enddo

      call uglydk(ip,1,1,rtmp,stmp,0.0,0.0,wcs)
      call evaldks(ip,1,rtmp,r2,s2,wcs)

      s2 = stmp(ip) - s2

      return
      end subroutine wk_interp4_hub
    
C =====================================
      SUBROUTINE PODGWAKE_FWA
C =====================================
      use m_wake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
!     PARAMETER(NWR=MBPZ)
!     COMMON /WAKGEO/ RTBL(11,NWR),XTBL(11),RW(NWR),VAR(NWR),VTR(NWR),
!    *       UASTAR(NWR),UTSTAR(NWR),UAUSTR(NWR),UTUSTR(NWR)
      DIMENSION DRWP(NWPZ)
      

C----- PART 1
      RPOD1 = wr(1,1,1)         
      DO N=2,ncpw
         XWL = wpt(1,N,1,1)     
         CALL EVALDKs(NHIN,1,XHIN,XWL,RPOD2,XYHINCUB)
         DR = RPOD2-RPOD1
         IF(DR.LT.0.0) GO TO 10
         RPOD1 = RPOD2
      ENDDO

 10   XCST  = wpt(1,N-1,1,1)    
      CALL EVALDKs(NHIN,1,XHIN,XCST,RPCST,XYHINCUB)
      DRCST = RPCST-SQRT(wpt(2,N-1,1,1)**2+wpt(3,N-1,1,1)**2)        


C----- PART 2
      DRWP(1)=0.0
      DO 20 N=2,ncpw         
         XWL = wpt(1,N,1,1)     
         IF(XWL.LE.XCST) THEN
            CALL EVALDKs(NHIN,1,XHIN,XWL,RPXL,XYHINCUB)
            RWWL = SQRT(wpt(2,N,1,1)**2 + wpt(3,N,1,1)**2)
            DRWP(N) = RPXL-RWWL 
         ELSE
            DRWP(N) = DRCST
         ENDIF
 20   ENDDO


C----- PART 3
      DO M = 1, mrpw
        DO N = 2, ncpw
          wr(N,M,1) = wr(N,M,1) + DRWP(N) 
          wpt(2,N,M,1) = wr(N,M,1)*COS(wth(N,M,1))
          wpt(3,N,M,1) = wr(N,M,1)*SIN(wth(N,M,1))
        ENDDO
      ENDDO  
      DO K = 2, nbld
        wr(:,:,k) = wr(:,:,1)
        DO M = 1, mrpw
          DO N = 2, ncpw
            wpt(2,N,M,K) = wr(N,M,K)*cos(wth(N,M,K))
            wpt(3,N,M,K) = wr(N,M,K)*sin(wth(N,M,K))
          ENDDO
        ENDDO
      ENDDO


      RETURN
      END

      subroutine bldwk2hub_pen
      use m_wake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      integer n1, n2, n3, i, j
      real rout1, xxx, dk
      real,allocatable :: rtmp1(:,:)
      real,allocatable :: xih1(:), rih1(:), hcx1(:)

      n1 = nhbx + 1
      n2 = ncpw 
      n3 = mrpw
      dk = 0.001 ! Height of the numerical fence.

      if (.not.allocated(rtmp1)) then
        allocate(rtmp1(n2,n3))
        allocate(xih1(n1),rih1(n1),hcx1(4*(n1-1)))
      endif

      do i = 1, n1
        xih1(i) = xh(i,1)
        rih1(i) = sqrt(yh(i,1)**2 + zh(i,1)**2)
      enddo
      call uglydk(n1,1,1,xih1,rih1,0,0,hcx1)

      do i = 2, n2
        if(wpt(1,i,1,1).GT.xh(nhbx+1,1)) go to 888
        do j = 1, n3
          rtmp1(i,j) = sqrt(wpt(2,i,j,1)**2 + wpt(3,i,j,1)**2)
          call evaldks(n1,1,xih1,wpt(1,i,j,1),rout1,hcx1)
          if (rtmp1(i,j).lt.rout1) then
            xxx = dangle(wpt(3,i,j,1),wpt(2,i,j,1))
            wpt(2,i,j,1) = (rout1+dk)*cos(xxx)
            wpt(3,i,j,1) = (rout1+dk)*sin(xxx)
          endif
        enddo
 888    continue
      enddo

      return
      end subroutine bldwk2hub_pen

C ===============================================
      SUBROUTINE DUCTWAKGEO_FIXEDPANEL_FWA
C ===============================================
      USE m_dwake
      INCLUDE 'PUFCAV.INC'
!------- Art Director Added This Line, 03/22/2016 -------!
      DIMENSION XID12(200), ETAD12(200)
      DIMENSION XID22(201), ETAD22(201)
      DIMENSION XXT(202), WWT(202), XWWTCUB(802)
      DIMENSION THTRSW(300,MDUCT+1)
      REAL,ALLOCATABLE:: XD_TMP1(:),YD_TMP1(:),ZD_TMP1(:)
      REAL,ALLOCATABLE:: XD_TMP2(:),YD_TMP2(:),ZD_TMP2(:)
      if(.NOT.allocated(XD_TMP1)) then
        allocate(XD_TMP1(mduct+1),YD_TMP1(mduct+1),ZD_TMP1(mduct+1))
        allocate(XD_TMP2(mrdw+1),YD_TMP2(mrdw+1),ZD_TMP2(mrdw+1))
      endif
!------- Art Director Added This Line, 03/22/2016
      IF(IDGEO .EQ. 0) THEN

        THETA1 = DANGLE(ZD(1,1),YD(1,1))
        THETA2 = DANGLE(ZDW(1,1),YDW(1,1))
        DTHETA = THETA1 - THETA2

        DO M = 1, MDUCTP
          DO N = 1, NDWKP 
            RRR = SQRT(YDW(N,M)**2 + ZDW(N,M)**2)
            THETA11 = DANGLE(ZDW(N,M),YDW(N,M))
            YDW(N,M) = RRR*COS(THETA11 + DTHETA)
            ZDW(N,M) = RRR*SIN(THETA11 + DTHETA)
          ENDDO
          XDW(1,M) = XD(1,M)
          YDW(1,M) = YD(1,M)
          ZDW(1,M) = ZD(1,M)
        ENDDO

        THETA2 = DANGLE(dwpt(3,1,1,1),dwpt(2,1,1,1))
        DTHETA = THETA1 - THETA2

        DO M = 1, mrpdw
          DO N = 1, ncpdw
            RRR = SQRT(dwpt(2,N,M,1)**2 + dwpt(3,N,M,1)**2)
            THETA11 = DANGLE(dwpt(3,N,M,1),dwpt(2,N,M,1))
            dwpt(2,N,M,1) = RRR*COS(THETA11 + DTHETA)
            dwpt(3,N,M,1) = RRR*SIN(THETA11 + DTHETA)
          ENDDO
        ENDDO
c        dwpt(1,1,1,1) = XD(1,1)
c        dwpt(2,1,1,1) = YD(1,1)
c        dwpt(3,1,1,1) = ZD(1,1)
c        dwpt(1,1,mrpdw,1) = XD(1,MDUCTP)
c        dwpt(2,1,mrpdw,1) = YD(1,MDUCTP)
c        dwpt(3,1,mrpdw,1) = ZD(1,MDUCTP)
        DO M = 1, mduct + 1
          XD_TMP1(M) = XD(1,M)
          YD_TMP1(M) = YD(1,M)
          ZD_TMP1(M) = ZD(1,M)
        ENDDO
        CALL wk_interp2(XD_TMP1,YD_TMP1,ZD_TMP1,MDUCT+1,
     *                  XD_TMP2,YD_TMP2,ZD_TMP2,mrdw+1)
        DO M = 1, mrdw+1
          dwpt(1,1,M,1) = XD_TMP2(M)
          dwpt(2,1,M,1) = YD_TMP2(M)
          dwpt(3,1,M,1) = ZD_TMP2(M)
        ENDDO

        CALL symm_dwkp2

      ELSE

        RRR = SQRT(YD(1,1)**2 + ZD(1,1)**2)

        NDWK = 100
        NDWKP = NDWK + 1
        DWL = DWAKEL
        DTN = PI / REAL(NDWK)

        DELTAK = TWOPI / NBLADE

        IF(IDGEO .EQ. 0) THEN
           PPD = XPI(NX)
           PITCH2 = 0.5 * PI - ATAN(PPD/PI)
        ELSEIF(IDGEO .EQ. 1) THEN
           PITCH2 = DUCTPT
        ENDIF

        DO M = 1 , MDUCTP
           XDW(1,M) = XD(1,M)
           YDW(1,M) = YD(1,M)
           ZDW(1,M) = ZD(1,M)

           THETA1 = DANGLE(ZDW(1,M),YDW(1,M))

           DO N = 2 , NDWKP
              XDI = 0.5 * DWL * (1. - COS(DTN*(N-1)))
              XDW(N,M) = XDW(1,M) + XDI

              DTHETA = XDI * TAN(PITCH2) / RRR
              THETA2 = THETA1 + DTHETA

              YDW(N,M) = RRR * COS(THETA2)
              ZDW(N,M) = RRR * SIN(THETA2)
           ENDDO
        ENDDO

      ENDIF

C --- Duct Wake Ends with Disk
      DO M = 1 , MDUCTP
         XDWDK(M,2) = XDW(NDWKP,M)
         YDWDK(M,2) = YDW(NDWKP,M)
         ZDWDK(M,2) = ZDW(NDWKP,M)

         XDWDK(M,1) = XDWDK(M,2)
         YDWDK(M,1) = 0.0
         ZDWDK(M,1) = 0.0
      ENDDO

      RETURN
      END

C************************************************************************
C  S.N.KIM  |  The following three subroutines (cal_bld_vort,
C           |  duct_wake_v, bld_wake_v) are brought here from 'bldindv.f' 
C Aug. 2018 |  to satisfy compatibility with 'm_wake' and 'm_dwake' 
C           |  modules. 
C************************************************************************
      subroutine cal_bld_vort
      USE m_dwake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      integer i,j,l
C/s S.N.KIM | includes for ductwake alignment
      integer imrw_tmp, incw_tmp
      real vtmp(3), ds_int
      real,dimension(mr) :: idx1,idx2,s1,s2,DELP_int
      real,dimension(4*(mr-1)) :: cubspl, cubdelp
C/e S.N.KIM | includes for ductwake alignment
      if (.NOT. allocated(gamonb_s)) then
        allocate(gamonb_s(nc,mr))
        allocate(gamonb_c(nc,mr+1))
C/s S.N.KIM | includes HUB effect | 08-15-2016
        if(ihub.ne.0) then
          allocate(gamonh_s(nhbx,mhbt))
          allocate(gamonh_c(nhbx+1,mhbt))
        endif
C/e S.N.KIM | includes HUB effect | 08-15-2016
        if(iduct.eq.1) then
          allocate(gamond_s(nduct,mduct))
          allocate(gamond_c(nduct,mduct+1))
          allocate(gamonwd_s(ndwk+20,mduct))    !=(ncdw,mduct)
          allocate(gamonwd_c(ndwk+20,mduct+1))  !=(ncdw,mduct+1)
        endif
C/s S.N.KIM | includes ductwake alignment | 08-15-2018
        if(idfwa.eq.1) then
          allocate(gamonw_s(nwpanel+20,mr))   !=(ncpw,mr)
          allocate(gamonw_c(nwpanel+20,mr+1)) !=(ncpw,mr+1) 
        endif
C/e S.N.KIM | includes ductwake alignment | 08-15-2018
      end if

      do i = 1, mr
        do j = 2, nc
          l1 = indexb(j-1,i)
          l2 = indexb(j,i)
          gamonb_s(j,i)=pot(l1)-pot(l2)
        end do
        l1 = indexb(1,i)
        l2 = indexb(nc,i)
        gamonb_s(1,i)=pot(l2)-pot(l1)-DELPonW(i)
      end do

      do i = 2, mr
        do j = 1,nc
          l1 = indexb(j,i-1)
          l2 = indexb(j,i)
          gamonb_c(j,i)=pot(l2)-pot(l1)
        end do
      end do

      do j = 1,nc
        l1 = indexb(j,1)
        l2 = indexb(j,mr)
        gamonb_c(j,1) = pot(l1)
        gamonb_c(j,mr+1) = -pot(l2)
      end do

C/s S.N.KIM | incldues HUB effect | 08-15-2016
      if(ihub.ne.0) then
        do i = 1, nhbx
          do j = 2, mhbt
            l1 = indexh(i,j)
            l2 = indexh(i,j-1)
            gamonh_s(i,j)=pot(l1)-pot(l2)
          enddo
          l1 = indexh(i,1)
          l2 = indexh(i,mhbt)
          gamonh_s(i,1)=pot(l1)-pot(l2)
        enddo
 
        do i = 2, nhbx
          do j = 1, mhbt
            l1 = indexh(i,j)
            l2 = indexh(i-1,j)
            gamonh_c(i,j)=pot(l2)-pot(l1)
          enddo
        enddo
 
        do j = 1, mhbt
          l1 = indexh(1,j)
          l2 = indexh(nhbx,j)
          gamonh_c(1,j) = -pot(l1)
          gamonh_c(nhbx+1,j) = pot(l2) 
        enddo
      endif

c      if(ihub.ne.0) then
c        do i = 1, mhbt
c          do j = 2, nhbx
c            l1 = indexh(j-1,i)
c            l2 = indexh(j,i)
c            gamonh_s(j,i)=pot(l1)-pot(l2)
c          enddo
c          l1 = indexh(1,i)
cc          l2 = indexh(nhbx,i)
c          gamonh_s(1,i)=-pot(l1)
c        enddo
c
c        do i = 2, mhbt
c          do j = 1, nhbx
c            l1 = indexh(j,i-1)
c            l2 = indexh(j,i)
c            gamonh_c(j,i)=pot(l2)-pot(l1)
c          enddo
c        enddo
c        do j = 1, nhbx
c          l1 = indexh(j,1)
c          l2 = indexh(j,mhbt)
c          gamonh_c(j,1) = pot(l1)
c          gamonh_c(j,mhbt+1) = -pot(l2)
c        enddo
c      endif
C/s S.N.KIM | incldues HUB effect | 08-15-2016

      if(iduct.eq.1) then
        do i = 1, mduct
          do j = 2, nduct
            l1 = indexd(j-1,i)
            l2 = indexd(j,i)
            gamond_s(j,i)=pot(l1)-pot(l2)
          end do
          l1 = indexd(1,i)
          l2 = indexd(nduct,i)
          gamond_s(1,i)=pot(l2)-pot(l1)-DELPD(i)
        end do

        do i = 2, mduct
          do j = 1,nduct
            l1 = indexd(j,i-1)
            l2 = indexd(j,i)
            gamond_c(j,i)=pot(l2)-pot(l1)
          end do
        end do

        do j = 1,nduct
          l1 = indexd(j,1)
          l2 = indexd(j,mduct)
          gamond_c(j,1) = pot(l1)
          gamond_c(j,mduct+1) = -pot(l2)
        end do

        do i = 1, mrdw
          do j = 1, ncdw-1
            gamonwd_s(j,i) = 0.0
          end do
          gamonwd_s(ncdw,i) = 0.0 !phidw1(i)
        end do

        do i = 2, mrdw
          do j = 1, ncdw
            gamonwd_c(j,i) = phidw1(i) - phidw1(i-1) !DELPD(i)-DELPD(i-1)
          end do
        end do

        do j = 1, ncdw
          gamonwd_c(j,1) = phidw1(1) !DELPD(1)
          gamonwd_c(j,mrdw+1) = -phidw1(mrdw) !-DELPD(mduct)
        end do
      endif

C/s S.N.KIM | Ductwake alignment | Aug. 2018
      IF(IDFWA.EQ.1) THEN
C--- Interpolate DELP for ductwake.
        DO i = 1, mr
          idx1(i) = real(i-1)/real(mr-1)
        ENDDO
        imrw_tmp = mr/ALIGN_DIV
        incw_tmp = nwpanel + 20
        DO i = 1, imrw_tmp
          idx2(i) = real(i-1)/real(imrw_tmp-1)
        ENDDO

        s1(1) = 0.0d0
        DO i = 2, mr
          vtmp(1) = (xw(1,i) - xw(1,i-1))/2 + (xw(1,i+1) - xw(1,i))/2
          vtmp(2) = (yw(1,i) - yw(1,i-1))/2 + (yw(1,i+1) - yw(1,i))/2
          vtmp(3) = (zw(1,i) - zw(1,i-1))/2 + (zw(1,i+1) - zw(1,i))/2
          ds_int = sqrt(dot_product(vtmp,vtmp))
          s1(i) = s1(i-1) + ds_int
        ENDDO
        CALL uglydk(mr,0,0,idx1,s1,0.0,0.0,cubspl)
        CALL evaldk(mr,imrw_tmp,idx1,idx2,s2,cubspl)
        CALL uglydk(mr,0,0,s1,DELP,0.0,0.0,cubdelp)
        CALL evaldk(mr,imrw_tmp,s1,s2,DELP_int,cubdelp)
C--- Interpolate DELP for ductwake.

        do i = 1, imrw_tmp
          do j = 1, incw_tmp-1
            gamonw_s(j,i) = 0.0
          enddo
          gamonw_s(incw_tmp,i) = 0.0 !DELP_int(i)
        enddo

        do i = 2, imrw_tmp
          do j = 1, incw_tmp
            gamonw_c(j,i) = DELP_int(i) - DELP_int(i-1)
          enddo
        enddo
 
        do j = 1, incw_tmp
          gamonw_c(j,1) = DELP_int(1)
          gamonw_c(j,imrw_tmp+1) = -DELP_int(imrw_tmp)
        enddo
      END IF
C/e S.N.KIM | Ductwake alignment | Aug. 2018

      return
      end subroutine cal_bld_vort

      subroutine duct_wake_v(nipt,x0,y,z,ux,ur,ut,deltab)
      USE m_dwake
      use omp_lib
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3)
      real x1(3),x2(3)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1,inf
      real,dimension(nipt) :: r, th, x0,r0,th0
      real the, the2, deltab, deltab2
      real*8 xppd(3),x1d(3),x2d(3),vtmpd(3),ddelta
      real x_tmp

      r0 = sqrt(y*y+z*z)
      th0 = atan2(z,y)

      inf = 0
!$OMP parallel private(i,m,n,kk,x1,x2,the,the2,vtmp,xpp,vdr,vdt
!$OMP&                 ,x_tmp,deltab2) 
!$OMP&shared(dwpt,delk,deltab,gamonwd_s,gamonwd_c,
!$OMP&       nblade,mrdw,ncdw,nipt)
!$OMP&shared(x,r,th,x0,r0,th0,ux,ur,ut)

C initiate x,r,th so that affinity is satisfied
!$OMP DO
      do i = 1, nipt
        x(i)=x0(i)
        r(i)=r0(i)
        th(i)=th0(i)
      enddo
!$OMP END DO

!$OMP DO
      do i = 1, nipt
        !spanwise vorticies
        do m=1,mrdw !mduct !mrdw
          do n=1,ncdw !ndwk !ncdw
            x1(1)=dwpt(1,n,m,1) !xdw(n,m) !dwpt(1,n,m,1)
            x1(2)=dwpt(2,n,m,1) !ydw(n,m) !dwpt(2,n,m,1)
            x1(3)=dwpt(3,n,m,1) !zdw(n,m) !dwpt(3,n,m,1)

            x2(1)=dwpt(1,n,m+1,1) !xdw(n,m+1) !dwpt(1,n,m+1,1)
            x2(2)=dwpt(2,n,m+1,1) !ydw(n,m+1) !dwpt(2,n,m+1,1)
            x2(3)=dwpt(3,n,m+1,1) !zdw(n,m+1) !dwpt(3,n,m+1,1)

C -- ductwake strength gets weaker with axial distance due to diffustion | S.N.KIM Aug. 2018.
            x_tmp = 0.5*(x1(1) + x2(1))
            if (x_tmp .gt. 0.5) then
              deltab2 = deltab + (x_tmp-0.5)*0.05
            else
              deltab2 = deltab
            endif
C -- ductwake strength gets weaker with axial distance due to diffustion | S.N.KIM Aug. 2018.

            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th(i)+the
              xpp(1) = x(i)
              xpp(2) = r(i)*cos(the2)
              xpp(3) = r(i)*sin(the2)

              call vseg_sd(xpp, x1, x2, deltab2, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)+(gamonwd_s(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)+(gamonwd_s(n,m)*vdr/TWOPI)
              ut(i)=ut(i)+(gamonwd_s(n,m)*vdt/TWOPI)
            end do
          enddo
        end do
        !chordwise vorticies
        do m=1,mrdw+1 !mduct+1 !mrdw+1
          do n=1,ncdw !ndwk !ncdw
            x1(1)=dwpt(1,n,m,1) !xdw(n,m) !dwpt(1,n,m,1)
            x1(2)=dwpt(2,n,m,1) !ydw(n,m) !dwpt(2,n,m,1)
            x1(3)=dwpt(3,n,m,1) !zdw(n,m) !dwpt(3,n,m,1)

            x2(1)=dwpt(1,n+1,m,1) !xdw(n+1,m) !dwpt(1,n+1,m,1)
            x2(2)=dwpt(2,n+1,m,1) !ydw(n+1,m) !dwpt(2,n+1,m,1)
            x2(3)=dwpt(3,n+1,m,1) !zdw(n+1,m) !dwpt(3,n+1,m,1)

C -- ductwake strength gets weaker with axial distance due to diffustion | S.N.KIM Aug. 2018.
            x_tmp = 0.5*(x1(1) + x2(1))
            if (x_tmp .gt. 0.5) then
              deltab2 = deltab + (x_tmp-0.5)*0.05
            else
              deltab2 = deltab
            endif
C -- ductwake strength gets weaker with axial distance due to diffustion | S.N.KIM Aug. 2018.

            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th(i)+the
              xpp(1) = x(i)
              xpp(2) = r(i)*cos(the2)
              xpp(3) = r(i)*sin(the2)

              call vseg_sd(xpp, x1, x2, deltab2, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)+(gamonwd_c(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)+(gamonwd_c(n,m)*vdr/TWOPI)
              ut(i)=ut(i)+(gamonwd_c(n,m)*vdt/TWOPI)
            end do
          enddo
        end do
      end do
!$OMP END DO
!$OMP end parallel
      return
      end subroutine duct_wake_v

C********************************************************************
C/s S.N.KIM | includes the calculation of the induced velocity 
C           | from blade wake to ductwake. | Aug. 2018
C********************************************************************
      subroutine bld_wake_v(nipt,x0,y,z,ux,ur,ut,deltab)
      USE m_wake
      use omp_lib
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3)
      real x1(3),x2(3)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1,inf
      real,dimension(nipt) :: r, th, x0,r0,th0
      real the, the2, deltab, deltab2
      real*8 xppd(3),x1d(3),x2d(3),vtmpd(3),ddelta
      real x_tmp

      r0 = sqrt(y*y+z*z)
      th0 = atan2(z,y)

      inf = 0
!$OMP parallel private(i,m,n,kk,x1,x2,the,the2,vtmp,xpp,vdr,vdt
!$OMP&                 ,x_tmp,deltab2)
!$OMP&shared(wpt,delk,deltab,gamonw_s,gamonw_c,
!$OMP&       nblade,mrw,ncw,nipt)
!$OMP&shared(x,r,th,x0,r0,th0,ux,ur,ut)

C initiate x,r,th so that affinity is satisfied
!$OMP DO
      do i = 1, nipt
        x(i)=x0(i)
        r(i)=r0(i)
        th(i)=th0(i)
      enddo
!$OMP END DO

!$OMP DO
      do i = 1, nipt
        !spanwise vorticies
        do m=1,mrw !mr !mrw
          do n=1,ncw !nwpanel !ncw
            x1(1)=wpt(1,n,m,1) !xw(n,m)  !wpt(1,n,m,1)
            x1(2)=wpt(2,n,m,1) !yw(n,m)  !wpt(2,n,m,1)
            x1(3)=wpt(3,n,m,1) !zw(n,m)  !wpt(3,n,m,1)

            x2(1)=wpt(1,n,m+1,1) !xw(n,m+1)  !wpt(1,n,m+1,1)
            x2(2)=wpt(2,n,m+1,1) !yw(n,m+1)  !wpt(2,n,m+1,1)
            x2(3)=wpt(3,n,m+1,1) !zw(n,m+1)  !wpt(3,n,m+1,1)

C -- blade wake strength gets weaker with axial distance due to diffusion | S.N.KIM Aug. 2018.
            x_tmp = 0.5*(x1(1) + x2(1))
            if (x_tmp .gt. 0.5) then
              deltab2 = deltab + (x_tmp-0.5)*0.05
            else
              deltab2 = deltab
            endif
C -- blade wake strength gets weaker with axial distance due to diffusion | S.N.KIM Aug. 2018.

            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th(i)+the
              xpp(1) = x(i)
              xpp(2) = r(i)*cos(the2)
              xpp(3) = r(i)*sin(the2)

              call vseg_sd(xpp, x1, x2, deltab2, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)-(gamonw_s(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)-(gamonw_s(n,m)*vdr/TWOPI)
              ut(i)=ut(i)-(gamonw_s(n,m)*vdt/TWOPI)
            end do
          enddo
        end do
        !chordwise vorticies
        do m=1,mrw+1 !mrp  !mrw+1
          do n=1,ncw !nwpanel !ncw
            x1(1)=wpt(1,n,m,1) !xw(n,m)  !wpt(1,n,m,1)
            x1(2)=wpt(2,n,m,1) !yw(n,m)  !wpt(2,n,m,1)
            x1(3)=wpt(3,n,m,1) !zw(n,m)  !wpt(3,n,m,1)

            x2(1)=wpt(1,n+1,m,1) !xw(n+1,m)  !wpt(1,n+1,m,1)
            x2(2)=wpt(2,n+1,m,1) !yw(n+1,m)  !wpt(2,n+1,m,1)
            x2(3)=wpt(3,n+1,m,1) !zw(n+1,m)  !wpt(3,n+1,m,1)

C -- blade wake strength gets weaker with axial distance due to diffusion | S.N.KIM Aug. 2018.
            x_tmp = 0.5*(x1(1) + x2(1))
            if (x_tmp .gt. 0.5) then
              deltab2 = deltab + (x_tmp-0.5)*0.05
            else
              deltab2 = deltab
            endif
C -- blade wake strength gets weaker with axial distance due to diffusion | S.N.KIM Aug. 2018.

            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th(i)+the
              xpp(1) = x(i)
              xpp(2) = r(i)*cos(the2)
              xpp(3) = r(i)*sin(the2)

              call vseg_sd(xpp, x1, x2, deltab2, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)-(gamonw_c(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)-(gamonw_c(n,m)*vdr/TWOPI)
              ut(i)=ut(i)-(gamonw_c(n,m)*vdt/TWOPI)
            end do
          enddo
        end do
      end do
!$OMP END DO
!$OMP end parallel
      return
      end subroutine bld_wake_v
C********************************************************************
C/e S.N.KIM | includes the calculation of the induced velocity 
C           | from blade wake to ductwake. | Aug. 2018
C********************************************************************

      subroutine run_alignp_unsteady(icav,imrw)
      use m_wake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      real,allocatable::x(:,:),y(:,:),z(:,:)
      real,allocatable:: xbtmp_o1(:),ybtmp_o1(:),zbtmp_o1(:)
      real,allocatable:: xbtmp_o2(:),ybtmp_o2(:),zbtmp_o2(:)
      real xbtmp1(imrw+1),ybtmp1(imrw+1),zbtmp1(imrw+1)
      real xbtmp2(imrw+1),ybtmp2(imrw+1),zbtmp2(imrw+1)
      real phiwtmp(imrw), phiwtmp2(imrw)
      real vtmp(3), dl
      real beta1,err1
      integer i,nipt,nmiss,icav, imrw, ii, iii
      character*29 :: wkcheck2,filename2
      write(*,*) 'Unsteady Full Wake Alignment Iteration Step =', NTSTEP
      if(.NOT.allocated(x)) then
        allocate(x(imrw+1,2),y(imrw+1,2),z(imrw+1,2))
        allocate(xbtmp_o1(mr+1),ybtmp_o1(mr+1),zbtmp_o1(mr+1))
        allocate(xbtmp_o2(mr+1),ybtmp_o2(mr+1),zbtmp_o2(mr+1))
      end if
      do i = 1, mr+1
        xbtmp_o1(i) = xb(1,i)
        ybtmp_o1(i) = yb(1,i)
        zbtmp_o1(i) = zb(1,i)
        xbtmp_o2(i) = xb(nc,i)
        ybtmp_o2(i) = yb(nc,i)
        zbtmp_o2(i) = zb(nc,i)
      end do

      DO KK = 1, NBLADE
        IREC = NTPOS(KK)
        CALL READ2(46,IREC,TEMP5,NWPANEL*MR)
        DO I = 1, NWPANEL
          DO J = 1, MR
            I35 = INDEXW2(I,J)
            DELPonW(J) = TEMP5(I35)
          ENDDO
          call wk_interp1(xbtmp_o1,ybtmp_o1,zbtmp_o1,delponw,mrp,
     &                    xbtmp1,ybtmp1,zbtmp1,phiwtmp,imrw+1)
          call trans_delponw_uns(phiwtmp,imrw)
          gamw1(kk,I,:) = gamw0
        ENDDO
      ENDDO 

      DO KK = 1, NBLADE
        IREC = NTPOS(KK)
        CALL READ2(46,IREC,TEMP5,NWPANEL*MR) 
        DO I = 1, NWPANEL-1
          DO J = 1, MR
            I35 = INDEXW2(I,J)
            I36 = INDEXW2(I+1,J)
            DELPonW(J) = TEMP5(I35)
            DELPonW2(J) = TEMP5(I36)
          ENDDO
          call wk_interp1(xbtmp_o1,ybtmp_o1,zbtmp_o1,delponw,mrp,
     &                    xbtmp1,ybtmp1,zbtmp1,phiwtmp,imrw+1)
          call wk_interp1(xbtmp_o1,ybtmp_o1,zbtmp_o1,delponw2,mrp,
     &                    xbtmp1,ybtmp1,zbtmp1,phiwtmp2,imrw+1)
          DO K = 1, MRW
             gamw2(kk,I+1,K) = phiwtmp2(K) - phiwtmp(K)
          ENDDO 
        ENDDO
      ENDDO

      nipt = mrpw*ncpw  ! = (mrw+1)*(ncw+1)
                        ! mrw = MR/ALIGN_DIV,  ncw = nwpanel + 20
      call cal_bld_vort
      call indvel4pU(nwpanel,mr)  ! After this, all the 'wv(:,ncpw,mrpw,1)' are set.
      call dump_wakgeo
      call velfrombld(nipt,wpx_dumy,wpy_dumy,wpz_dumy
     &               ,wvxb,wvrb,wvtb,nmiss,ivbmiss,delta_w
     &               ,iduct,ihub)! nipt = mrpw * ncpw
      call merge_vel2_uns
      call align_wkps_unsteady
      call out_wk  ! 'wake.dat' file is generated
      call out_bld  ! nc+1 X mr+1 xb, yb, zb data will be written on 'bld.dat'
      end subroutine run_alignp_unsteady

      subroutine align_wkps_unsteady
      use m_wake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      INTEGER NMZ ! SABUSABU  
      ALLOCATABLE :: GEOX(:),GEOY(:),GEOZ(:) ! SABUSABU

      ALLOCATABLE :: DUMV(:),DUMX(:),DUMY(:),DUMZ(:)
      ALLOCATABLE :: DUMV1(:),DUMX1(:),DUMY1(:),DUMZ1(:)
      ALLOCATABLE :: DUMVX(:),DUMVY(:),DUMVZ(:)
      ALLOCATABLE :: VXYZ(:,:,:)
      ALLOCATABLE :: XW2(:,:),YW2(:,:),ZW2(:,:)
      ALLOCATABLE :: RW2(:,:),TW2(:,:)

      REAL, DIMENSION (3) :: V_MEAN,DEL_SI,U_N,PT_I

      ALLOCATE(DUMV(MRMAX),DUMX(MRMAX),DUMY(MRMAX),DUMZ(MRMAX))
      ALLOCATE(DUMV1(MRMAX),DUMX1(MRMAX),DUMY1(MRMAX),DUMZ1(MRMAX))
      ALLOCATE(DUMVX(4*MRMAX),DUMVY(4*MRMAX),DUMVZ(4*MRMAX))
      ALLOCATE(VXYZ(NWMAX,MRMAXP,3))
      ALLOCATE(XW2(NWMAXP,MRMAXP),YW2(NWMAXP,MRMAXP),ZW2(NWMAXP,MRMAXP))
      ALLOCATE(RW2(NWMAXP,MRMAXP),TW2(NWMAXP,MRMAXP))

      NMZ = MBPZ * NWPZ ! SABUSABU
      ALLOCATE(GEOX(NMZ),GEOY(NMZ),GEOZ(NMZ)) ! SABUSABU

      DTPROP = DELTAT/RAD
      DTF = DTPROP*ADVCO/180.0
      DTI = DTF
      NWPP = NWPANEL + 1 !KKKK NCPW

      IST = NWPANEL - 4
      NITER = 4*ICAVT
      NN1 = 1
      NSWW = IST

      IF(NTSTEP .EQ. 0 .AND. NITER .LT. IST) THEN
        NN1 = NITER - 3
        NSWW = NITER
      ENDIF
   
      DO N = 1, NWPP 
        DO M = 1, MRPW
          XW2(N,M) = wpt(1,N,M,1)
          YW2(N,M) = wpt(2,N,M,1)
          ZW2(N,M) = wpt(3,N,M,1)

          RW2(N,M) = sqrt(wpt(2,N,M,1)**2 + wpt(3,N,M,1)**2)
          TW2(N,M) = atan2(wpt(3,N,M,1),wpt(2,N,M,1))
        ENDDO
      ENDDO

      IF(IHUB .NE. 0) NHBUB = NHBU + NH

C/S.N.KIM - I commented out this part sine this is the old unsteady
Calignment scheme.
      DO M = MRPW, 1, -1
        DO N = NWPP-1, NN1, -1
          XW2(N+1,M) = XW2(N,M) 
     %               + (0.5*wv(1,N,M,1) + 0.5*wv(1,N+1,M,1))*DTI
          YW2(N+1,M) = RW2(N,M)*COS(TW2(N,M) + PI*DTI/ADVCO 
     %               + womg(N,M)*DTI + wvtb(N,M)*DTI/RW2(N,M))
     %               + (0.5*wv(2,N,M,1) + 0.5*wv(2,N+1,M,1))*DTI
          ZW2(N+1,M) = RW2(N,M)*SIN(TW2(N,M) + PI*DTI/ADVCO 
     %               + womg(N,M)*DTI + wvtb(N,M)*DTI/RW2(N,M))
     %               + (0.5*wv(3,N,M,1) + 0.5*wv(3,N+1,M,1))*DTI
        ENDDO
      ENDDO

C/ SABUSABU -------------------------------------------------------
C/ Seungnam Kim sets underrelaxation factor to improve convergence.
      CALL READ2(143,NTPOS(1),GEOX,NSAPW)
      CALL READ2(144,NTPOS(1),GEOY,NSAPW)
      CALL READ2(145,NTPOS(1),GEOZ,NSAPW)
      ICC = MRP
      DO N = 2, NWK
        DO M = 1, MRP
          ICC = ICC + 1
          XW2(N,M) = XW2(N,M)*0.90 + GEOX(ICC)*0.10
          YW2(N,M) = YW2(N,M)*0.90 + GEOY(ICC)*0.10
          ZW2(N,M) = ZW2(N,M)*0.90 + GEOZ(ICC)*0.10
        ENDDO
      ENDDO
C/ SABUSABU -------------------------------------------------------

 1111 DO N = 1, NWPANEL+1 ! KKKK NCPW
        DO M = 1, MRPW
          wpt(1,N,M,1) = XW2(N,M)
          wpt(2,N,M,1) = YW2(N,M)
          wpt(3,N,M,1) = ZW2(N,M)
        ENDDO
      ENDDO

C This subroutine pushes wake panels penetrating hub above the hub
C surface.
      IF (IHUB.NE.0) CALL BLDWK2HUB_PEN
C/S.N.KIM - Wake-Hub Penetration Control

      DEALLOCATE( DUMV,DUMX,DUMY,DUMZ )
      DEALLOCATE( DUMV1,DUMX1,DUMY1,DUMZ1 )
      DEALLOCATE( DUMVX,DUMVY,DUMVZ )
      DEALLOCATE( VXYZ )
      DEALLOCATE( XW2,YW2 )
      DEALLOCATE( ZW2 )

      return
      end subroutine align_wkps_unsteady

      SUBROUTINE INFLOWK3_NO_WROVS
************************************************************************
*     INFLOW3: INFLOW velocities at the wake control points            *
*                                                                      *
*     Date     Comment or Revision                                     *
*     -------- -------------------                                     *
*     10/13/99   HSLEE                                                 *
*     1/9/2017   SNKIM                                                 *
************************************************************************
      USE m_wake
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION U(3),NH1(3)

      IF(NTSTEP.EQ.0) THEN
        DO 10 IH=1,3
          NH1(IH)=1
CCC          NH1(IH)=NHARM(IH)
   10   CONTINUE
      ELSE
        DO 15 IH=1,3
          NH1(IH)=NHARM(IH)
   15   CONTINUE
      END IF
c
      DO I = 1, mrpw
        DO J = 1, ncpw
          RCP=SQRT(wpt(2,J,I,1)**2+wpt(3,J,I,1)**2)  
          THP=ATAN2(wpt(3,J,I,1),wpt(2,J,I,1))
          DO 30 IH=1,3
            U(IH)=0.0
            DO 20 JH=1,NH1(IH)
              CALL EVALDKs(NWKCOE,1,XRW,RCP,COEC,XWCUB(1,JH,1,IH))
              CALL EVALDKs(NWKCOE,1,XRW,RCP,COES,XWCUB(1,JH,2,IH))
              T=FLOAT(JH-1)*(THP-TSTEP)
              ST=SIN(T)
              CT1=COS(T)
              U(IH)=U(IH)+COEC*CT1+COES*ST
   20       CONTINUE
   30     CONTINUE
          wvi(1,J,I)=U(1)
          VORW=U(2)
          VOTW=U(3)
          COCP=wpt(2,J,I,1)/RCP !SQRT(wpt(2,J,I,1)**2+wpt(3,J,I,1)**2)
          SICP=wpt(3,J,I,1)/RCP !SQRT(wpt(2,J,I,1)**2+wpt(3,J,I,1)**2)
c          WROVS=PI*RCP/ADVCO  ! R*Omega component is considered in subroutine 'align_wkps_unsteady'
          womg(J,I)=VOTW/RCP   ! Rotational Angle due to tangential velcity.
          wvi(2,J,I)=VORW*COCP !-VOTW*SICP
          wvi(3,J,I)=VORW*SICP !+VOTW*COCP
        ENDDO
      ENDDO

      RETURN
      END
