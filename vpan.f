      module m_constant
        implicit none
        real*8, parameter :: PI = 3.141592653589793d0
        real*8, parameter :: TWOPI = 2.0d0*PI
        real*8, parameter :: ZERO=0.0d0
        real*8, parameter :: HLF=0.5d0
        real*8, parameter :: HALF=0.5d0
        real*8, parameter :: ONE=1.0d0
        real*8, parameter :: TWO=2.0d0
        real*8, parameter :: THREE=3.0d0
        real*8, parameter :: DEG = PI/180d0
        real*8, parameter :: Tenth = 0.1d0 
        integer,parameter :: dm = 3 !dimension 3d
      end module
! =======================================================================
      module m_constant_s
        implicit none
        real, parameter :: PI = 3.141592653589793
        real, parameter :: TWOPI = 2.00*PI
        real, parameter :: ZERO=0.00
        real, parameter :: HLF=0.50
        real, parameter :: HALF=0.50
        real, parameter :: ONE=1.00
        real, parameter :: TWO=2.00
        real, parameter :: THREE=3.00
        real, parameter :: DEG = PI/1800
        real, parameter :: Tenth = 0.10 
        integer,parameter :: dm = 3 !dimension 3d
      end module

      subroutine svpan(xp, x, delta, v)
      use m_constant
      implicit none
      real, intent(in) ::xp(dm)
      real, intent(in) ::x(dm, 4)
      real, intent(in) :: delta
      real, intent(out):: v(dm)

      real*8 ::dxp(dm)
      real*8 ::dx(dm, 4)
      real*8 ::ddelta
      real*8:: dv(dm)

      integer i
      real*8 vtmp(dm)

      dxp = dble(xp)
      dx  = dble(x)
      ddelta  = dble(delta)

      dv = ZERO
      do i = 1, 3
        call vseg_d(dxp,dx(:,i),dx(:,i+1),ddelta ,vtmp,0)
        dv = dv+ vtmp
      end do
        call vseg_d(dxp,dx(:,4),dx(:,1),ddelta ,vtmp,0)
        dv = dv+ vtmp
      dv = dv/TWOPI

      v = sngl(dv)
      return
      end subroutine svpan
! =======================================================================

! =======================================================================
      subroutine svpan_124(xp, x, delta, v)
      use m_constant
      implicit none
      real, intent(in) ::xp(dm)
      real, intent(in) ::x(dm, 4)
      real, intent(in) :: delta
      real, intent(out):: v(dm)

      real*8 ::dxp(dm)
      real*8 ::dx(dm, 4)
      real*8 ::ddelta
      real*8:: dv(dm)

      integer i
      real*8 vtmp(dm)

      dxp = dble(xp)
      dx  = dble(x)
      ddelta  = dble(delta)

      dv = ZERO
      call vseg_d(dxp,dx(:,1),dx(:,2),ddelta ,vtmp,0)
      dv = dv+ vtmp
      call vseg_d(dxp,dx(:,2),dx(:,3),ddelta ,vtmp,0)
      dv = dv+ vtmp
      call vseg_d(dxp,dx(:,4),dx(:,1),ddelta ,vtmp,0)
      dv = dv+ vtmp
      dv = dv/TWOPI

      v = sngl(dv)
      return
      end subroutine svpan_124
! =======================================================================
! =======================================================================
      subroutine svpan_234(xp, x, delta, v)
      use m_constant
      implicit none
      real, intent(in) ::xp(dm)
      real, intent(in) ::x(dm, 4)
      real, intent(in) :: delta
      real, intent(out):: v(dm)

      real*8 ::dxp(dm)
      real*8 ::dx(dm, 4)
      real*8 ::ddelta
      real*8:: dv(dm)

      integer i
      real*8 vtmp(dm)

      dxp = dble(xp)
      dx  = dble(x)
      ddelta  = dble(delta)

      dv = ZERO
      call vseg_d(dxp,dx(:,2),dx(:,3),ddelta ,vtmp,0)
      dv = dv+ vtmp
      call vseg_d(dxp,dx(:,3),dx(:,4),ddelta ,vtmp,0)
      dv = dv+ vtmp
      call vseg_d(dxp,dx(:,4),dx(:,1),ddelta ,vtmp,0)
      dv = dv+ vtmp
      dv = dv/TWOPI

      v = sngl(dv)
      return
      end subroutine svpan_234
! =======================================================================
! =======================================================================
      subroutine svpan_24(xp, x, delta, v)
      use m_constant
      implicit none
      real, intent(in) ::xp(dm)
      real, intent(in) ::x(dm, 4)
      real, intent(in) :: delta
      real, intent(out):: v(dm)

      real*8 ::dxp(dm)
      real*8 ::dx(dm, 4)
      real*8 ::ddelta
      real*8:: dv(dm)

      integer i
      real*8 vtmp(dm)

      dxp = dble(xp)
      dx  = dble(x)
      ddelta  = dble(delta)

      dv = ZERO
      call vseg_d(dxp,dx(:,2),dx(:,3),ddelta ,vtmp,0)
      dv = dv+ vtmp
      call vseg_d(dxp,dx(:,4),dx(:,1),ddelta ,vtmp,0)
      dv = dv+ vtmp
      dv = dv/TWOPI

      v = sngl(dv)
      return
      end subroutine svpan_24
! =======================================================================

! =======================================================================
      subroutine svpan_1(xp, x, delta, v)
      use m_constant
      implicit none
      real, intent(in) ::xp(dm)
      real, intent(in) ::x(dm, 4)
      real, intent(in) :: delta
      real, intent(out):: v(dm)

      real*8 ::dxp(dm)
      real*8 ::dx(dm, 4)
      real*8 ::ddelta
      real*8:: dv(dm)

      integer i
      real*8 vtmp(dm)

      dxp = dble(xp)
      dx  = dble(x)
      ddelta  = dble(delta)

      dv = ZERO
      call vseg_d(dxp,dx(:,1),dx(:,2),ddelta ,vtmp,0)
      dv = dv+ vtmp
      dv = dv/TWOPI

      v = sngl(dv)
      return
      end subroutine svpan_1
! =======================================================================

! =======================================================================
      subroutine svpan_3(xp, x, delta, v)
      use m_constant
      implicit none
      real, intent(in) ::xp(dm)
      real, intent(in) ::x(dm, 4)
      real, intent(in) :: delta
      real, intent(out):: v(dm)

      real*8 ::dxp(dm)
      real*8 ::dx(dm, 4)
      real*8 ::ddelta
      real*8:: dv(dm)

      integer i
      real*8 vtmp(dm)

      dxp = dble(xp)
      dx  = dble(x)
      ddelta  = dble(delta)

      dv = ZERO
      call vseg_d(dxp,dx(:,3),dx(:,4),ddelta ,vtmp,0)
      dv = dv+ vtmp
      dv = dv/TWOPI

      v = sngl(dv)
      return
      end subroutine svpan_3
! =======================================================================
! =======================================================================
      subroutine dvpan(xp, x, delta, v)
      use m_constant
      implicit none
      real*8, intent(in) ::xp(dm)
      real*8, intent(in) ::x(dm, 4)
      real*8, intent(in) :: delta
!     integer,intent(in) :: inf
      real*8 , intent(out):: v(dm)

      integer i
      real*8 vtmp(dm)

      v = ZERO
      do i = 1, 4
        call vseg_d(xp, x(:,i), x(:,i+1), delta ,vtmp,0)
        v = v+ vtmp
      end do
      v = v/TWOPI

      return
      end subroutine dvpan
! =======================================================================
      subroutine vseg_d(xp, x1, x2, delta, v, inf)
      use m_constant
      implicit none
      real*8 , intent(in) :: xp(dm), x1(dm), x2(dm)  
      real*8,  intent(in) :: delta != 0.01d0
      integer, intent(in) :: inf
      real*8 , intent(out):: v(dm)

      real*8,  parameter  :: radius = 5.0d0
      real*8,  parameter  :: tol =    1.0d-8

      real*8 :: r0(dm), r1(dm), r2(dm), r1xr2(dm), fvec(dm), h(dm)
      real*8 :: lr0,l2r0, lr1, lr2, l2r1xr2, tmp, f

      r0 = x2 - x1
      l2r0= dot_product(r0,r0)
      if (l2r0 .lt. tol) then
        v = ZERO
        return
      endif
      lr0= dsqrt(l2r0)
      fvec =HALF * (x2+x1)-xp
      f=dsqrt(dot_product(fvec,fvec))

      if((f/lr0.ge.radius).and.(inf .ne. 0)) then
        tmp=HALF/f**3
        h(1) =r0(3)*fvec(2)-r0(2)*fvec(3)
        h(2) =r0(1)*fvec(3)-r0(3)*fvec(1)
        h(3) =r0(2)*fvec(1)-r0(1)*fvec(2)
        v=tmp*h
        return
      endif


      r1 = xp - x1
      lr1= dsqrt(dot_product(r1,r1))
      r2 = xp - x2
      lr2= dsqrt(dot_product(r2,r2))
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
      end subroutine vseg_d 

      subroutine vseg_sd(xp, x1, x2, delta, v, inf)
      use m_constant
      implicit none
      real , intent(in) :: xp(dm), x1(dm), x2(dm)  
      real,  intent(in) :: delta != 0.01d0
      integer, intent(in) :: inf
      real , intent(out):: v(dm)

      real,  parameter  :: radius = 5.0
      real*8,  parameter  :: tol =    1.0e-8

      real :: r0(dm), r1(dm), r2(dm), r1xr2(dm), fvec(dm), h(dm)
      real :: lr0,l2r0, lr1, lr2, l2r1xr2, tmp, f

      r0 = x2 - x1
      l2r0= dot_product(r0,r0)
      if (l2r0 .lt. tol) then
        v = ZERO
        return
      endif
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
      end subroutine vseg_sd 

