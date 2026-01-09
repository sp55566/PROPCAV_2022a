!Ye TIAN 04/08/2012
!calculate the blade induced flow field, and convert into (x,r,th)
!coordinate systme as harmonics:  u=sum(A(i,x,r)*cos(i*th),i=1,nh).
!**Notice**
! The subroutines in this file CANNOT handle the rudder case
!****NOTICE 2********************************************
!in order to get the delphi on the wake after IPK, the kutta.f has to be
!modified. after LINE:140, before calling PRSDIF(1), the value of DELP
!has to be stored in a different array.
!It is convenient to add a DELPonW(:) array in the  circ module in the
!propcav.f file

      module m_field_harmonics
      implicit none
      integer fh_nh    !number of harmonics
      integer fh_mr    !number of radial stations
      integer fh_nx    !number of axial stations
      real    fh_amp_factor !=real(nblade)
      real,allocatable:: fh_x0(:),fh_r0(:) !the coordinates of the initial station, dimension(fh_mr)
      real,allocatable:: fh_x(:)           !the axial offset of all the stations, dimension(fh_nx)

      real,allocatable:: fhs_ux_a0(:,:)  !a0 coefficient, dimension(fh_mr,fh_nx), blade source induced
      real,allocatable:: fhs_ux_a(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)
      real,allocatable:: fhs_ux_b(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)
      real,allocatable:: fhs_ur_a0(:,:)  !a0 coefficient, dimension(fh_mr,fh_nx), blade source induced
      real,allocatable:: fhs_ur_a(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)
      real,allocatable:: fhs_ur_b(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)
      real,allocatable:: fhs_ut_a0(:,:)  !a0 coefficient, dimension(fh_mr,fh_nx), blade source induced
      real,allocatable:: fhs_ut_a(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)
      real,allocatable:: fhs_ut_b(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)

      real,allocatable:: fhd_ux_a0(:,:)  !a0 coefficient, dimension(fh_mr,fh_nx), blade source induced
      real,allocatable:: fhd_ux_a(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)
      real,allocatable:: fhd_ux_b(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)
      real,allocatable:: fhd_ur_a0(:,:)  !a0 coefficient, dimension(fh_mr,fh_nx), blade source induced
      real,allocatable:: fhd_ur_a(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)
      real,allocatable:: fhd_ur_b(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)
      real,allocatable:: fhd_ut_a0(:,:)  !a0 coefficient, dimension(fh_mr,fh_nx), blade source induced
      real,allocatable:: fhd_ut_a(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)
      real,allocatable:: fhd_ut_b(:,:,:) !a coefficients, dimension(fh_nh,fh_mr,fh_nx)
      end module m_field_harmonics

      module m_rpan_data
        type t_geom_rpan
          real:: xv1(4)
          real:: yv1(4)
          real:: sq1(15)
          real:: side1(4)
        end type t_geom_rpan
      end module m_rpan_data


      subroutine bld_source_v(nipt,x,y,z,ux,ur,ut)
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3), xpan(3, 4)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1
      real,dimension(nipt) :: r, th
      real the, the2

      r = sqrt(y*y+z*z)
      th = atan2(z,y)

      ux = 0.0
      ur = 0.0
      ut = 0.0
      imr1 = 1
      do i = 1, nipt
        do j = 1 , npanb
          do k = 1 , 4
            xv(k) = xvp(j,k)
            yv(k) = yvp(j,k)
            side(k) = sid(j,k)
          end do
          do k = 1 , 15
            s(k) = ss(j,k)
          enddo

          do kk = 1 , nblade
            the=-delk*real(kk-1)
            the2=th(i)+the
            xpp(1) = x(i)
            xpp(2) = r(i)*cos(the2)
            xpp(3) = r(i)*sin(the2)

            xloc = 0.0
            yloc = 0.0
            zloc = 0.0

            do k = 1 , 3
              xloc = xloc+(xpp(k)-xct(j,k))*dir(j,1,k)
              yloc = yloc+(xpp(k)-xct(j,k))*dir(j,2,k)
              zloc = zloc+(xpp(k)-xct(j,k))*dir(j,3,k)
            enddo

            imr = imr1
            call rpan(xloc,yloc,zloc,chrleps(j),
     %                fs,fd,fsx,fsy,fdx,fdy,fdz,1,imr)

            vsx=fsx*dir(j,1,1)+fsy*dir(j,2,1)-fd*dir(j,3,1)
            vsy=fsx*dir(j,1,2)+fsy*dir(j,2,2)-fd*dir(j,3,2)
            vsz=fsx*dir(j,1,3)+fsy*dir(j,2,3)-fd*dir(j,3,3)

            vsr= vsy*cos(the2)+vsz*sin(the2)
            vst=-vsy*sin(the2)+vsz*cos(the2)

            ux(i)=ux(i)+(dpdn(j)*vsx/4.0/pi)
            ur(i)=ur(i)+(dpdn(j)*vsr/4.0/pi)
            ut(i)=ut(i)+(dpdn(j)*vst/4.0/pi)
          end do
        enddo
      end do

      return
      end subroutine bld_source_v

      subroutine bld_source_v2(nipt,x,y,z,ux,ur,ut)
      use m_rpan_data
      use omp_lib
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      type(t_geom_rpan)::geom_data
      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3), xpan(3, 4)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1
      real,dimension(nipt) :: r, th
      real the, the2, chl


      r = sqrt(y*y+z*z)
      th = atan2(z,y)

      ux = 0.0
      ur = 0.0
      ut = 0.0
      imr1 = 1
!$OMP parallel private(i,j,k,kk,geom_data,the,the2,
!$OMP&                 fs,fd,fsx,fsy,fdx,fdy,fdz,chl,imr,
!$OMP& xloc,yloc,zloc,vsx,vsy,vsz,vsr,vst,xpp,vdr,vdt)
!$OMP&shared(xvp,yvp,sid,ss,x,r,th,chrleps,dir,dpdn,imr1
!$OMP&       ux,ur,ut,nblade,npanb,nipt)

!$OMP DO
      do i = 1, nipt
        do j = 1 , npanb
          do k = 1 , 4
            geom_data%xv1(k) = xvp(j,k)
            geom_data%yv1(k) = yvp(j,k)
            geom_data%side1(k) = sid(j,k)
          end do
          do k = 1 , 15
            geom_data%sq1(k) = ss(j,k)
          enddo
C Yiran Su 20180522 avoided false sharing problem by putting the shared
C                   variable to a private variable. This can greatly
C                   improve the parallel efficiency
          chl=chrleps(j)

          do kk = 1 , nblade
            the=-delk*real(kk-1)
            the2=th(i)+the
            xpp(1) = x(i)
            xpp(2) = r(i)*cos(the2)
            xpp(3) = r(i)*sin(the2)

            xloc = 0.0
            yloc = 0.0
            zloc = 0.0

            do k = 1 , 3
              xloc = xloc+(xpp(k)-xct(j,k))*dir(j,1,k)
              yloc = yloc+(xpp(k)-xct(j,k))*dir(j,2,k)
              zloc = zloc+(xpp(k)-xct(j,k))*dir(j,3,k)
            enddo

            imr = imr1

            call RPAN_PAR(xloc,yloc,zloc,chl,
     %                fs,fd,fsx,fsy,fdx,fdy,fdz,1,imr,geom_data)

            vsx=fsx*dir(j,1,1)+fsy*dir(j,2,1)-fd*dir(j,3,1)
            vsy=fsx*dir(j,1,2)+fsy*dir(j,2,2)-fd*dir(j,3,2)
            vsz=fsx*dir(j,1,3)+fsy*dir(j,2,3)-fd*dir(j,3,3)

            vsr= vsy*cos(the2)+vsz*sin(the2)
            vst=-vsy*sin(the2)+vsz*cos(the2)

            ux(i)=ux(i)+(dpdn(j)*vsx/4.0/pi)
            ur(i)=ur(i)+(dpdn(j)*vsr/4.0/pi)
            ut(i)=ut(i)+(dpdn(j)*vst/4.0/pi)
          end do
        enddo
      end do
!$OMP END DO
!$OMP end parallel
      return
      end subroutine bld_source_v2



      subroutine bld_dipole_v(nipt,x,y,z,ux,ur,ut,deltab)
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3), xpan(3, 4)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1
      real,dimension(nipt) :: r, th
      real the, the2, deltab
      real*8 xppd(3),x1d(3),x2d(3),vtmpd(3),ddelta

      ddelta = dble(deltab)
      r = sqrt(y*y+z*z)
      th = atan2(z,y)
      ux = 0.0
      ur = 0.0
      ut = 0.0
      imr1 = 1
      do i = 1, nipt
        do m=mr,1,-1
          do n=1,nc
            l = indexb(n,m)
            xpan(1,1)=xb(n,m+1)
            xpan(2,1)=yb(n,m+1)
            xpan(3,1)=zb(n,m+1)
            xpan(1,2)=xb(n,m)
            xpan(2,2)=yb(n,m)
            xpan(3,2)=zb(n,m)
            xpan(1,3)=xb(n+1,m)
            xpan(2,3)=yb(n+1,m)
            xpan(3,3)=zb(n+1,m)
            xpan(1,4)=xb(n+1,m+1)
            xpan(2,4)=yb(n+1,m+1)
            xpan(3,4)=zb(n+1,m+1)

            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th(i)+the
              xpp(1) = x(i)
              xpp(2) = r(i)*cos(the2)
              xpp(3) = r(i)*sin(the2)


              call svpan(xpp, xpan, deltab, vtmp)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)-(pot(l)*vtmp(1))
              ur(i)=ur(i)-(pot(l)*vdr)
              ut(i)=ut(i)-(pot(l)*vdt)
!last panel of a section
              if (n.eq.nc) then
                xppd = dble(xpp)
                x1d  = dble(xpan(:,3))
                x2d  = dble(xpan(:,4))
                call vseg_d(xppd,x1d,x2d,ddelta,vtmpd,0)
                vtmp = sngl(vtmpd /TWOPI)
                vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
                vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

                ux(i)=ux(i)+(DELPonW(m)*vtmp(1))
                ur(i)=ur(i)+(DELPonW(m)*vdr)
                ut(i)=ut(i)+(DELPonW(m)*vdt)
              endif
            end do
          enddo
        end do
      end do

      return
      end subroutine bld_dipole_v

      subroutine bld_dipole_v2(nipt,x0,y,z,ux,ur,ut,deltab)
      use omp_lib

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3), xpan(3, 4)
      real x1(3),x2(3)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1,inf
      real,dimension(nipt) :: r, th, x0,r0,th0
      real the, the2, deltab
      real*8 xppd(3),x1d(3),x2d(3),vtmpd(3),ddelta

      r0 = sqrt(y*y+z*z)
      th0 = atan2(z,y)

c      ux = 0.0
c      ur = 0.0
c      ut = 0.0

      inf = 0
!$OMP parallel private(i,m,n,kk,x1,x2,the,the2,vtmp,xpp,vdr,vdt)
!$OMP&shared(xb,yb,zb,delk,deltab,gamonb_s,gamonb_c,
!$OMP&       nblade,mr,nc,nipt)
!$OMP&shared(x,r,th,x0,r0,th0,ux,ur,ut)

C initiate x,r,th so that affinity is satisfied
!$OMP DO
      do i = 1, nipt
        x(i)=x0(i)
        r(i)=r0(i)
        th(i)=th0(i)
      end do
!$OMP END DO

!$OMP DO
      do i = 1, nipt
        !spanwise vorticies
        do m=1,mr
          do n=1,nc
            x1(1)=xb(n,m)
            x1(2)=yb(n,m)
            x1(3)=zb(n,m)

            x2(1)=xb(n,m+1)
            x2(2)=yb(n,m+1)
            x2(3)=zb(n,m+1)
            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th(i)+the
              xpp(1) = x(i)
              xpp(2) = r(i)*cos(the2)
              xpp(3) = r(i)*sin(the2)

              call vseg_sd(xpp, x1, x2, deltab, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)-(gamonb_s(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)-(gamonb_s(n,m)*vdr/TWOPI)
              ut(i)=ut(i)-(gamonb_s(n,m)*vdt/TWOPI)

            end do

          enddo
        end do

        !chordwise vorticies
        do m=1,mr+1
          do n=1,nc
            x1(1)=xb(n,m)
            x1(2)=yb(n,m)
            x1(3)=zb(n,m)

            x2(1)=xb(n+1,m)
            x2(2)=yb(n+1,m)
            x2(3)=zb(n+1,m)
            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th(i)+the
              xpp(1) = x(i)
              xpp(2) = r(i)*cos(the2)
              xpp(3) = r(i)*sin(the2)


              call vseg_sd(xpp, x1, x2, deltab, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)-(gamonb_c(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)-(gamonb_c(n,m)*vdr/TWOPI)
              ut(i)=ut(i)-(gamonb_c(n,m)*vdt/TWOPI)
            end do
          enddo
        end do
      end do
!$OMP END DO
!$OMP end parallel
      return
      end subroutine bld_dipole_v2


C********************************************************
C/s S.N.KIM | includes HUB effect | 08-15-2018
C********************************************************  
      subroutine hub_source_v2(nipt,x,y,z,ux,ur,ut)
      use m_rpan_data
      use omp_lib
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      type(t_geom_rpan)::geom_data
      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3), xpan(3, 4)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1
      real,dimension(nipt) :: r, th
      real the, the2, chl

      r = sqrt(y*y+z*z)
      th = atan2(z,y)

      imr1 = 1
!$OMP parallel private(i,j,k,kk,geom_data,the,the2,
!$OMP&                 fs,fd,fsx,fsy,fdx,fdy,fdz,chl,imr,
!$OMP& xloc,yloc,zloc,vsx,vsy,vsz,vsr,vst,xpp,vdr,vdt)
!$OMP&shared(xvp,yvp,sid,ss,x,r,th,chrleps,dir,dpdn,imr1
!$OMP&       ux,ur,ut,nblade,mhbt,nhbx,nipt)

!$OMP DO
      do i = 1, nipt
        do m=1,mhbt
          do n=1,nhbx
            j=indexh(n,m)
          do k = 1 , 4
            geom_data%xv1(k) = xvp(j,k)
            geom_data%yv1(k) = yvp(j,k)
            geom_data%side1(k) = sid(j,k)
          end do
          do k = 1 , 15
            geom_data%sq1(k) = ss(j,k)
          enddo

          chl = chrleps(j)

          do kk = 1 , nblade
            the=-delk*real(kk-1)
            the2=th(i)+the
            xpp(1) = x(i)
            xpp(2) = r(i)*cos(the2)
            xpp(3) = r(i)*sin(the2)

            xloc = 0.0
            yloc = 0.0
            zloc = 0.0

            do k = 1 , 3
              xloc = xloc+(xpp(k)-xct(j,k))*dir(j,1,k)
              yloc = yloc+(xpp(k)-xct(j,k))*dir(j,2,k)
              zloc = zloc+(xpp(k)-xct(j,k))*dir(j,3,k)
            enddo

            imr = imr1
            call RPAN_PAR(xloc,yloc,zloc,chl,
     %                fs,fd,fsx,fsy,fdx,fdy,fdz,1,imr,geom_data)

            vsx=fsx*dir(j,1,1)+fsy*dir(j,2,1)-fd*dir(j,3,1)
            vsy=fsx*dir(j,1,2)+fsy*dir(j,2,2)-fd*dir(j,3,2)
            vsz=fsx*dir(j,1,3)+fsy*dir(j,2,3)-fd*dir(j,3,3)

            vsr= vsy*cos(the2)+vsz*sin(the2)
            vst=-vsy*sin(the2)+vsz*cos(the2)

            ux(i)=ux(i)+(dpdn(j)*vsx/4.0/pi)
            ur(i)=ur(i)+(dpdn(j)*vsr/4.0/pi)     ! j = indexd(n,m)
            ut(i)=ut(i)+(dpdn(j)*vst/4.0/pi)
          end do
        enddo
      end do
      enddo
!$OMP END DO
!$OMP end parallel
      return
      end subroutine hub_source_v2

      subroutine hub_dipole_v2(nipt,x0,y,z,ux,ur,ut,deltab)
      use omp_lib

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3), xpan(3, 4)
      real x1(3),x2(3)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1,inf
      real,dimension(nipt) :: r, th, x0, r0, th0
      real the, the2, deltab
      real*8 xppd(3),x1d(3),x2d(3),vtmpd(3),ddelta

      r0 = sqrt(y*y+z*z)
      th0 = atan2(z,y)

      inf = 0
!$OMP parallel private(i,m,n,kk,x1,x2,the,the2,vtmp,xpp,vdr,vdt)
!$OMP&shared(xh,yh,zh,delk,deltab,gamonh_s,gamonh_c,
!$OMP&       nblade,mhbt,nhbx,nipt)
!$OMP&shared(x,r,th,x0,r0,th0,ux,ur,ut)

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
        do n=1,nhbx
          do m=1,mhbt
            x1(1)=xh(n+1,m)
            x1(2)=yh(n+1,m)
            x1(3)=zh(n+1,m)

            x2(1)=xh(n,m)
            x2(2)=yh(n,m)
            x2(3)=zh(n,m)
            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th(i)+the
              xpp(1) = x(i)
              xpp(2) = r(i)*cos(the2)
              xpp(3) = r(i)*sin(the2)

              call vseg_sd(xpp, x1, x2, deltab, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)-(gamonh_s(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)-(gamonh_s(n,m)*vdr/TWOPI)
              ut(i)=ut(i)-(gamonh_s(n,m)*vdt/TWOPI)
            end do
          enddo
        end do
        !chordwise vorticies
        do n=1,nhbx+1
         do m=1,mhbt
            x1(1)=xh(n,m)
            x1(2)=yh(n,m)
            x1(3)=zh(n,m)

            x2(1)=xh(n,m+1)
            x2(2)=yh(n,m+1)
            x2(3)=zh(n,m+1)
            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th(i)+the
              xpp(1) = x(i)
              xpp(2) = r(i)*cos(the2)
              xpp(3) = r(i)*sin(the2)


              call vseg_sd(xpp, x1, x2, deltab, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)-(gamonh_c(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)-(gamonh_c(n,m)*vdr/TWOPI)
              ut(i)=ut(i)-(gamonh_c(n,m)*vdt/TWOPI)

            end do
          enddo
        end do
      end do
!$OMP END DO
!$OMP end parallel
      return
      end subroutine hub_dipole_v2
C********************************************************
C/e S.N.KIM | includes HUB effect | 08-15-2018
C********************************************************  

c-----------------------------------------
c  H.F. added following lines.       2015
c  S.N.KIM OpenMP paralled.     Aug. 2018 
c-----------------------------------------
      subroutine duct_source_v(nipt,x,y,z,ux,ur,ut)
      use m_rpan_data
      use omp_lib
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      type(t_geom_rpan)::geom_data
      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3), xpan(3, 4)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1
      real,dimension(nipt) :: r, th
      real the, the2, chl

      r = sqrt(y*y+z*z)
      th = atan2(z,y)

c      ux = 0.0
c      ur = 0.0
c      ut = 0.0

      imr1 = 1
!$OMP parallel private(i,j,k,kk,geom_data,the,the2,
!$OMP&                 fs,fd,fsx,fsy,fdx,fdy,fdz,chl,imr,
!$OMP& xloc,yloc,zloc,vsx,vsy,vsz,vsr,vst,xpp,vdr,vdt)
!$OMP&shared(xvp,yvp,sid,ss,x,r,th,chrleps,dir,dpdn,imr1
!$OMP&       ux,ur,ut,nblade,mduct,nduct,nipt)


!$OMP DO
      do i = 1, nipt
        do m=1,mduct
          do n=1,nduct
            j=indexd(n,m)
          do k = 1 , 4
            geom_data%xv1(k) = xvp(j,k)
            geom_data%yv1(k) = yvp(j,k)
            geom_data%side1(k) = sid(j,k)
          end do
          do k = 1 , 15
            geom_data%sq1(k) = ss(j,k)
          enddo
C Seungnam Kim 20180602 avoided false sharing problem by putting the
C                       shared variable to a private variable. This can
C                       greatly improve the parallel efficiency
C                       - brough this part from blade part.
          chl=chrleps(j)

          do kk = 1 , nblade
            the=-delk*real(kk-1)
            the2=th(i)+the
            xpp(1) = x(i)
            xpp(2) = r(i)*cos(the2)
            xpp(3) = r(i)*sin(the2)

            xloc = 0.0
            yloc = 0.0
            zloc = 0.0

            do k = 1 , 3
              xloc = xloc+(xpp(k)-xct(j,k))*dir(j,1,k)
              yloc = yloc+(xpp(k)-xct(j,k))*dir(j,2,k)
              zloc = zloc+(xpp(k)-xct(j,k))*dir(j,3,k)
            enddo

            imr = imr1
            call RPAN_PAR(xloc,yloc,zloc,chl,
     %                fs,fd,fsx,fsy,fdx,fdy,fdz,1,imr,geom_data)

            vsx=fsx*dir(j,1,1)+fsy*dir(j,2,1)-fd*dir(j,3,1)
            vsy=fsx*dir(j,1,2)+fsy*dir(j,2,2)-fd*dir(j,3,2)
            vsz=fsx*dir(j,1,3)+fsy*dir(j,2,3)-fd*dir(j,3,3)

            vsr= vsy*cos(the2)+vsz*sin(the2)
            vst=-vsy*sin(the2)+vsz*cos(the2)

            ux(i)=ux(i)+(dpdn(j)*vsx/4.0/pi)
            ur(i)=ur(i)+(dpdn(j)*vsr/4.0/pi)     ! j = indexd(n,m)
            ut(i)=ut(i)+(dpdn(j)*vst/4.0/pi)
          end do
        enddo
      end do
      enddo
!$OMP ENDDO
!$OMP end parallel
      return
      end subroutine duct_source_v


      subroutine duct_dipole_v(nipt,x0,y,z,ux,ur,ut,deltab)
      use omp_lib

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3), xpan(3, 4)
      real x1(3),x2(3)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1,inf
      real,dimension(nipt) :: r, th, x0,r0,th0
      real the, the2, deltab
      real*8 xppd(3),x1d(3),x2d(3),vtmpd(3),ddelta

      r0 = sqrt(y*y+z*z)
      th0 = atan2(z,y)

c      ux=0
c      ur=0
c      ut=0

      inf = 0
!$OMP parallel private(i,m,n,kk,x1,x2,the,the2,vtmp,xpp,vdr,vdt
!$OMP&                 ,r1t,r2t,rr0t,r1xr2t,lr0t,hht)
!$OMP&shared(xd,yd,zd,delk,deltab,gamond_s,gamond_c,
!$OMP&       nblade,mduct,nduct,nipt)
!$OMP&shared(x,r,th,x0,r0,th0,ux,ur,ut)

C initiate x,r,th so that affinity is satisfied
!$OMP DO
      do i = 1, nipt
        x(i)=x0(i)
        r(i)=r0(i)
        th(i)=th0(i)
      end do
!$OMP END DO

!$OMP DO
      do i = 1, nipt
        !spanwise vorticies
        do m=1,mduct
          do n=1,nduct
            x1(1)=xd(n,m)
            x1(2)=yd(n,m)
            x1(3)=zd(n,m)

            x2(1)=xd(n,m+1)
            x2(2)=yd(n,m+1)
            x2(3)=zd(n,m+1)
            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th(i)+the
              xpp(1) = x(i)
              xpp(2) = r(i)*cos(the2)
              xpp(3) = r(i)*sin(the2)

              call vseg_sd(xpp, x1, x2, deltab, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)+(gamond_s(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)+(gamond_s(n,m)*vdr/TWOPI)
              ut(i)=ut(i)+(gamond_s(n,m)*vdt/TWOPI)
            end do
          enddo
        end do
        !chordwise vorticies
        do m=1,mduct+1
         do n=1,nduct
            x1(1)=xd(n,m)
            x1(2)=yd(n,m)
            x1(3)=zd(n,m)

            x2(1)=xd(n+1,m)
            x2(2)=yd(n+1,m)
            x2(3)=zd(n+1,m)
            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th(i)+the
              xpp(1) = x(i)
              xpp(2) = r(i)*cos(the2)
              xpp(3) = r(i)*sin(the2)

              call vseg_sd(xpp, x1, x2, deltab, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)+(gamond_c(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)+(gamond_c(n,m)*vdr/TWOPI)
              ut(i)=ut(i)+(gamond_c(n,m)*vdt/TWOPI)

            end do
          enddo
        end do
      end do
!$OMP END DO
!$OMP end parallel 
      return
      end subroutine duct_dipole_v

      subroutine duct_ind
        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC'
        integer n_temp !number of points to evaluate velocities
        real x_temp(8),y_temp(8),z_temp(8)
        real ux_temp(8),ur_temp(8),ut_temp(8)
        integer i
        real deltab
        deltab = 0.01
        n_temp=8

        do i=1,8
          x_temp(i)=0.16
          y_temp(i)=0.10*(i-1)+0.17
          z_temp(i)=0.0
        enddo

        call cal_bld_vort

        call duct_source_v(n_temp,x_temp,y_temp,z_temp,
     &                    ux_temp,ur_temp,ut_temp)

        open(1106,FILE='duct_ind.dat',STATUS='UNKNOWN')
!        write(1106,'(4F8.5)') (y_temp(i),ux_temp(i),
!     &                        ur_temp(i),ut_temp(i),i=1,8)
!        write(1106,*) ' '

        call duct_dipole_v(n_temp,x_temp,y_temp,z_temp,
     &                    ux_temp,ur_temp,ut_temp,deltab)
         call duct_wake_v(n_temp,x_temp,y_temp,z_temp,
     &                   ux_temp,ur_temp,ut_temp,deltab)

        write(1106,'(4F8.5)') (y_temp(i),ux_temp(i),
     &                        ur_temp(i),ut_temp(i),i=1,8)

!        write(1106,'(I4,F8.5)') (i,DELPD(i),i=1,mduct)
!        close(1106)

        return

       end subroutine duct_ind



      subroutine harmonics(nsr, f, th0, nh, a0, a, b)
      implicit none
      real, parameter :: eps=0.00000000001
      integer :: nsr   !number of sectors
      real    :: f(nsr), th0 !f, the periodic function
                             !th0, the phase angle for the first data
                             !point
      integer :: nh    !number of harmonics
      real    :: a0, a(nh), b(nh)

      real    :: c(nh, nsr), s(nh, nsr)
      real    :: t, step
      real    :: PI

      integer :: n, k

      PI = 2.0*acos(0.0)
      step = 2.0*PI/real(NSR)

      do n = 1, nh
        do k = 1, nsr
          t = n*(step*(k-1.0)+th0)
          c(n,k)=cos(t)
          s(n,k)=sin(t)
        end do
      end do

      a0 = 0.0
      a  = 0.0
      b  = 0.0

      do k = 1, nsr
        a0   = a0 + f(k)
      end do
      do n = 1, nh
        do k = 1, nsr
          a(n) = a(n) + f(k) * c(n,k)
          b(n) = b(n) + f(k) * s(n,k)
        end do
      end do
      a0 = a0 / nsr
      a  = a  / nsr * 2.0
      b  = b  / nsr * 2.0
      if (abs(a0).lt.eps) a0=0.0
      do n = 1, nh
        if (abs(a(n)).lt.eps) a(n)=0.0
        if (abs(b(n)).lt.eps) b(n)=0.0
      end do
      return
      end subroutine harmonics

      subroutine rearrage_data(nsr, y, z, th, id)
      implicit none
      real, parameter :: eps=0.00001
      integer :: nsr
      real    :: y(nsr), z(nsr)
      integer :: id(nsr)
      real    :: th(nsr), TWOPI, tmp

      integer :: i, j, tmpi

      TWOPI = 4.0*acos(0.0)


      do i = 1, nsr
        id(i) = i
      end do

      th = atan2(z,y)
      do i = 1, nsr
        if (th(i) .lt. -eps) then
          th(i) = th(i) + TWOPI
        endif
!       write(17,*) 'th', i, th(i)
      end do

!     write(17,*) '========='

!-----Start sorting for incremental angle
      do i = 1, nsr-1
        do j = i + 1, nsr
          if (th(i) .gt. th(j)) then
            !---swap th(i) and th(j)---
            tmp   = th(i)
            th(i) = th(j)
            th(j) = tmp
            !---swap id(i) and id(j)---
            tmpi   = id(i)
            id(i) = id(j)
            id(j) = tmpi
          endif
        end do
      end do

!     write(17,*) 'after sort';
!     do i = 1, nsr
!       write(17,*) id(i),th(i)
!     end do
!-----end of sorting for incremental angle
      return
      end subroutine rearrage_data

      subroutine har2v(nipt,x,y,z,ux,ur,ut,nmiss,miss_id)
      use m_field_harmonics
      implicit none
      integer nipt,nmiss
      real,dimension(nipt)::x,y,z,ux,ur,ut
      real r,x_bnd, alfa_r,beta_r, xdiff_tmp, th
      real alfa_x,beta_x
      real c11,c12,c21,c22
      integer miss_id(nipt)
      integer i, j, k,id_r,id_x
      real,allocatable:: fh_xcoef(:)

      miss_id = 0
      ux = 0.0
      ur = 0.0
      ut = 0.0
      if (.NOT. allocated(fh_xcoef)) then
        allocate(fh_xcoef(5*(fh_mr-1)))
      end if
      do i = 1, nipt
        r = sqrt(y(i)*y(i)+z(i)*z(i))
        if (r.gt.fh_r0(fh_mr)) then
          miss_id(i)=1
          cycle
        end if
        do j = 1, fh_mr-1
          if ((r.ge.fh_r0(j)).and.(r.lt.fh_r0(j+1))) then
            alfa_r=(r-fh_r0(j))/(fh_r0(j+1)-fh_r0(j))
            beta_r=(fh_r0(j+1)-r)/(fh_r0(j+1)-fh_r0(j))
            x_bnd = fh_x0(j)*beta_r+fh_x0(j+1)*alfa_r
            id_r = j
            exit
          end if
        end do
        if (x(i).lt.x_bnd) then
          miss_id(i)=1
          cycle
        end if
        xdiff_tmp = x(i)-x_bnd
        if (xdiff_tmp.ge.fh_x(fh_nx)) then
          ux(i) = 0.0
          ur(i) = 0.0
          ut(i) = 0.0
          cycle
        end if
        do k = 1, fh_nx-1
          if ((xdiff_tmp .ge. fh_x(k))
     &      .and.(xdiff_tmp.le.fh_x(k+1))) then
            alfa_x = (xdiff_tmp-fh_x(k))/(fh_x(k+1)-fh_x(k))
            beta_x = (fh_x(k+1)-xdiff_tmp)/(fh_x(k+1)-fh_x(k))
            id_x = k
            exit
          endif
        end do
        th = atan2(z(i),y(i))*fh_amp_factor
        c11 = beta_r*beta_x
        c12 = alfa_r*beta_x
        c21 = beta_r*alfa_x
        c22 = alfa_r*alfa_x
        ux(i)=c11*fhd_ux_a0(id_r,id_x)
     &       +c12*fhd_ux_a0(id_r+1,id_x)
     &       +c21*fhd_ux_a0(id_r,id_x+1)
     &       +c22*fhd_ux_a0(id_r+1,id_x+1)
        ur(i)=c11*fhd_ur_a0(id_r,id_x)
     &       +c12*fhd_ur_a0(id_r+1,id_x)
     &       +c21*fhd_ur_a0(id_r,id_x+1)
     &       +c22*fhd_ur_a0(id_r+1,id_x+1)
        ut(i)=c11*fhd_ut_a0(id_r,id_x)
     &       +c12*fhd_ut_a0(id_r+1,id_x)
     &       +c21*fhd_ut_a0(id_r,id_x+1)
     &       +c22*fhd_ut_a0(id_r+1,id_x+1)
        do k = 1, fh_nh
          ux(i)=ux(i) + cos(real(k)*th)*
     &       (c11*fhd_ux_a(k,id_r,id_x)
     &       +c12*fhd_ux_a(k,id_r+1,id_x)
     &       +c21*fhd_ux_a(k,id_r,id_x+1)
     &       +c22*fhd_ux_a(k,id_r+1,id_x+1))
          ur(i)=ur(i) + cos(real(k)*th)*
     &       (c11*fhd_ur_a(k,id_r,id_x)
     &       +c12*fhd_ur_a(k,id_r+1,id_x)
     &       +c21*fhd_ur_a(k,id_r,id_x+1)
     &       +c22*fhd_ur_a(k,id_r+1,id_x+1))
          ut(i)=ut(i) + cos(real(k)*th)*
     &       (c11*fhd_ut_a(k,id_r,id_x)
     &       +c12*fhd_ut_a(k,id_r+1,id_x)
     &       +c21*fhd_ut_a(k,id_r,id_x+1)
     &       +c22*fhd_ut_a(k,id_r+1,id_x+1))
        end do
        do k = 1, fh_nh
          ux(i)=ux(i) + sin(real(k)*th)*
     &       (c11*fhd_ux_b(k,id_r,id_x)
     &       +c12*fhd_ux_b(k,id_r+1,id_x)
     &       +c21*fhd_ux_b(k,id_r,id_x+1)
     &       +c22*fhd_ux_b(k,id_r+1,id_x+1))
          ur(i)=ur(i) + sin(real(k)*th)*
     &       (c11*fhd_ur_b(k,id_r,id_x)
     &       +c12*fhd_ur_b(k,id_r+1,id_x)
     &       +c21*fhd_ur_b(k,id_r,id_x+1)
     &       +c22*fhd_ur_b(k,id_r+1,id_x+1))
          ut(i)=ut(i) + sin(real(k)*th)*
     &       (c11*fhd_ut_b(k,id_r,id_x)
     &       +c12*fhd_ut_b(k,id_r+1,id_x)
     &       +c21*fhd_ut_b(k,id_r,id_x+1)
     &       +c22*fhd_ut_b(k,id_r+1,id_x+1))
        end do
      end do
      nmiss = sum(miss_id)
      write(*,*) 'missed:',nmiss
      return
      end subroutine har2v


      subroutine velfrombld(nipt,x,y,z,ux,ur,ut,nmiss,miss_id,deltab
     *                      ,iduct,ihub)
      implicit none
      integer nipt,nmiss
      real,dimension(nipt)::x,y,z,ux,ur,ut
      integer miss_id(nipt)
      integer i,iduct,ihub
      real deltab
      deltab = 0.07
      call bld_source_v2(nipt,x,y,z,ux,ur,ut)
      call bld_dipole_v2(nipt,x,y,z,ux,ur,ut,0.09) !deltab)
C/s S.N.KIM | includes HUB effect | 08-15-2016
      if (ihub.ne.0) then
        call hub_source_v2(nipt,x,y,z,ux,ur,ut)
        call hub_dipole_v2(nipt,x,y,z,ux,ur,ut,0.09) !deltab)
      end if
C/e S.N.KIM | includes HUB effect | 08-15-2016
      if (iduct.eq.1) then
        call duct_source_v(nipt,x,y,z,ux,ur,ut)
        call duct_dipole_v(nipt,x,y,z,ux,ur,ut,0.09) !deltab)
        call duct_wake_v(nipt,x,y,z,ux,ur,ut,deltab)
      end if

      return
      end subroutine velfrombld

      subroutine veltoduct(niptd,x,y,z,ux,ur,ut,nmissd,miss_id,deltab
     *                      ,idopt)
      implicit none
      integer niptd,nmissd
      real,dimension(niptd)::x,y,z,ux,ur,ut
      integer miss_id(niptd)
      integer i,idopt
      real deltab
      deltab = 0.07
      if (idopt.eq.1) then
        ux = 0.0
        ur = 0.0
        ut = 0.0
      endif
      if (idopt.eq.1) go to 216
      call bld_source_v2(niptd,x,y,z,ux,ur,ut)
      call bld_dipole_v2(niptd,x,y,z,ux,ur,ut,0.09) !deltab)
      CALL bld_wake_v(niptd,x,y,z,ux,ur,ut,deltab)
 216  continue
      call duct_source_v(niptd,x,y,z,ux,ur,ut)
      call duct_dipole_v(niptd,x,y,z,ux,ur,ut,0.09) !deltab)
      return
      end subroutine veltoduct

      subroutine bld_source_v_missed(nipt,x,y,z,ux,ur,ut,miss_id)
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      integer miss_id(nipt)
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3), xpan(3, 4)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1,ipt
      real :: r, th
      real the, the2

      imr1 = 1
      do i = 1, nipt
        if (miss_id(i).ne.1) cycle
        r = sqrt(y(i)*y(i)+z(i)*z(i))
        th = atan2(z(i),y(i))
        do j = 1 , npanb
          do k = 1 , 4
            xv(k) = xvp(j,k)
            yv(k) = yvp(j,k)
            side(k) = sid(j,k)
          end do
          do k = 1 , 15
            s(k) = ss(j,k)
          enddo

          do kk = 1 , nblade
            the=-delk*real(kk-1)
            the2=th+the
            xpp(1) = x(i)
            xpp(2) = r*cos(the2)
            xpp(3) = r*sin(the2)

            xloc = 0.0
            yloc = 0.0
            zloc = 0.0

            do k = 1 , 3
              xloc = xloc+(xpp(k)-xct(j,k))*dir(j,1,k)
              yloc = yloc+(xpp(k)-xct(j,k))*dir(j,2,k)
              zloc = zloc+(xpp(k)-xct(j,k))*dir(j,3,k)
            enddo

            imr = imr1
            call rpan(xloc,yloc,zloc,chrleps(j),
     %                fs,fd,fsx,fsy,fdx,fdy,fdz,1,imr)

            vsx=fsx*dir(j,1,1)+fsy*dir(j,2,1)-fd*dir(j,3,1)
            vsy=fsx*dir(j,1,2)+fsy*dir(j,2,2)-fd*dir(j,3,2)
            vsz=fsx*dir(j,1,3)+fsy*dir(j,2,3)-fd*dir(j,3,3)

            vsr= vsy*cos(the2)+vsz*sin(the2)
            vst=-vsy*sin(the2)+vsz*cos(the2)

            ux(i)=ux(i)+(dpdn(j)*vsx/4.0/pi)
            ur(i)=ur(i)+(dpdn(j)*vsr/4.0/pi)
            ut(i)=ut(i)+(dpdn(j)*vst/4.0/pi)
          end do
        enddo
      end do

      return
      end subroutine bld_source_v_missed

      subroutine bld_dipole_v2_missed(nipt,x,y,z,ux,ur,ut,deltab
     &                               ,miss_id)
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      integer nipt !number of points to evaluate velocities
      real,dimension(nipt) :: x,y,z,ux,ur,ut
      integer miss_id(nipt)
      real XPP(3), xptmp(3), vtmp(3), vtmp1(3), xpan(3, 4)
      real x1(3),x2(3)
      real xloc,yloc,zloc
      integer i,n, m, kk, imr, imr1,inf
      real :: r, th
      real the, the2, deltab
      real*8 xppd(3),x1d(3),x2d(3),vtmpd(3),ddelta

      inf = 0
!$OMP parallel private(i,m,n,kk,x1,x2,the,the2,vtmp,xpp,vdr,vdt)
!$OMP&shared(xb,yb,zb,delk,deltab,gamonb_s,gamonb_c,x,r,th,
!$OMP&       ux,ur,ut,nblade,mr,nc,nipt)
!$OMP DO
      do i = 1, nipt
        !spanwise vorticies
        if (miss_id(i).ne.1) cycle
        r = sqrt(y(i)*y(i)+z(i)*z(i))
        th = atan2(z(i),y(i))
        do m=1,mr
          do n=1,nc
            x1(1)=xb(n,m)
            x1(2)=yb(n,m)
            x1(3)=zb(n,m)

            x2(1)=xb(n,m+1)
            x2(2)=yb(n,m+1)
            x2(3)=zb(n,m+1)
            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th+the
              xpp(1) = x(i)
              xpp(2) = r*cos(the2)
              xpp(3) = r*sin(the2)

              call vseg_sd(xpp, x1, x2, deltab, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)-(gamonb_s(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)-(gamonb_s(n,m)*vdr/TWOPI)
              ut(i)=ut(i)-(gamonb_s(n,m)*vdt/TWOPI)
            end do
          enddo
        end do
        !chordwise vorticies
        do m=1,mr+1
          do n=1,nc
            x1(1)=xb(n,m)
            x1(2)=yb(n,m)
            x1(3)=zb(n,m)

            x2(1)=xb(n+1,m)
            x2(2)=yb(n+1,m)
            x2(3)=zb(n+1,m)
            do kk = 1 , nblade
              the=-delk*real(kk-1)
              the2=th+the
              xpp(1) = x(i)
              xpp(2) = r*cos(the2)
              xpp(3) = r*sin(the2)


              call vseg_sd(xpp, x1, x2, deltab, vtmp, inf)

              vdr= vtmp(2)*cos(the2)+vtmp(3)*sin(the2)
              vdt=-vtmp(2)*sin(the2)+vtmp(3)*cos(the2)

              ux(i)=ux(i)-(gamonb_c(n,m)*vtmp(1)/TWOPI)
              ur(i)=ur(i)-(gamonb_c(n,m)*vdr/TWOPI)
              ut(i)=ut(i)-(gamonb_c(n,m)*vdt/TWOPI)
            end do
          enddo
        end do
      end do
!$OMP END DO
!$OMP end parallel

      return
      end subroutine bld_dipole_v2_missed


C--- ------   MIT  --  RPAN FOR HIGH ASPECT RATIO &TOL E-8 ------------C
C                                                                      C
C       COPYRIGHT (C) 1985 MASSACHUSETTS INSTITUTE OF TECHNOLOGY       C
C                                                                      C
C----------------------  RELEASE DATE 9/1/86   ------------------------C
C                                                                      C
C     RPAN SUBROUTINE EVALUATES THE SOURCE POTENTIAL FS AND            C
C     NORMAL DIPOLE POTENTIAL FD  FOR GREEN FUNCTION 1/R, OVER         C
C     A QUADRILATERAL PANEL WITH VERTICES LOCATED AT XV(N),YV(N)       C
C     AND FIELD POINT (X,Y,Z).  IF INPUT PARAMETER ID.NE.0 SUBROUTINE  C
C     ALSO EVALUATES THE DERIVATIVES FSX,FSY,FDX,FDY,FDZ.              C
C                                                                      C
C     IF NOT USING DERIVATIVE OPTION LINES  INDICATED BY COMMENTS      C
C        MAY BE DELETED AND ID PARAMETER REMOVED FROM INPUTS           C
C                                                                      C
C     GEOMETRY OF PANEL IS DEFINED IN COMMON BLOCK INCLUDING VERTEX    C
C     COORDINATES XV(4),YV(4), LENGTH OF SIDES SIDE(4), AND            C
C     MOMENTS OF INERTIA IN THE ARRAY S(15) DEFINED AS FOLLOWS         C
C     (ROW=POWER OF X, COLUMN=POWER OF Y)                              C
C                                                                      C
C          S(1)   S(2)   S(3)   S(4)   S(5)                            C
C          S(6)   S(7)   S(8)   S(9)                                   C
C          S(10)  S(11)  S(12)                                         C
C          S(13)  S(14)                                                C
C          S(15)                                                       C
C                                                                      C
C       E.G. S(1)=AREA, S(3)=2ND Y-MOMENT, S(6)=1ST X-MOMENT           C
C                                                                      C
C     IF TWO VERTICES ARE CLOSER THAN TOL (1E**-8) THEN A              C
C     TRIANGLE IS ASSUMED.                                             C
C                                                                      C
C     FIELD POINT COORDINATES (X,Y,Z) ARE INPUT VIA SUBROUTINE         C
C     DUMMY VARIABLES.  PLANE Z=0 COINCIDES WITH THE PANEL.  IF        C
C     ABS(Z).LT.TOL (10**-8 UNLESS CHANGED IN DATA) THEN VALUE OF      C
C     Z IS SET EQUAL TO TOL TO AVOID COMPUTATIONAL ERRORS IF THE       C
C     FIELD POINT IS ON OR CLOSE TO THE PANEL.  FOR THIS REASON THE    C
C     CONVENTION MUST BE FOLLOWED THAT DOMAIN Z.GT.0 IS FLUID REGION   C
C     ADJACENT TO THIS PANEL, AND PANEL VERTICES ARE NUMBERED FROM     C
C     1 TO 4 IN CLOCKWISE DIRECTION.  SEE HESS AND SMITH PAPERS FOR    C
C     THIS CONVENTION.  SEE MANUSCRIPT "DISTRIBUTIONS OF SOURCES       C
C     AND NORMAL DIPOLES OVER A QUADRILATERAL PANEL" BY NEWMAN FOR     C
C     ADDITIONAL DETAILS.                                              C
C                                                                      C
C     EVALUATION IS BASED ON MULTIPOLE APPROXIMATIONS WITH 2ND         C
C     MOMENTS INCLUDED, FOR R*R GREATER THAN 150*CHRLENS, OTHERWISE    C
C     4TH MOMENTS FOR 40*CHRLENS OR EXACT FORMULATION IF LESS.  THIS   C
C     SHOULD GIVE 6D ACCURACY UNLESS PANEL ASPECT RATIO EXCEEDS 10.    C
C                                                                      C
C                 HART TABLE 5090 IS USED FOR THE ARCTANGENT.          C
C                                                                      C
C     FOLLOWING TIMES APPLY ON VAX 11/750, MILLISECONDS PER RETURN     C
C                                                                      C
C     ID    EXACT    4TH MOMENT    2ND MOMENT                          C
C     0      1.4       0.4           0.2                               C
C     1      2.6       1.0           0.4                               C
C                                                                      C
C  9/1/86 REVISED TO FIX DERIVATIVES NEAR EXTENSION LINES, ALSO        C
C    STREAMLINED MAIN CODE AND INSERTED SPECIAL SECTIONS FOR FIELD     C
C    POINT NEAR PANEL EDGE OR NORMAL TO A VERTEX.  ALSO REVISED SOME   C
C    PARAMETER NOTATIONS TO ACCORD WITH PAPER.  IN REVISED FORM        C
C    FIRST FORM OF (3.9) IS USED NEAR PANEL EDGES, FDZ EVALUATED       C
C    FROM DERIVATIVE OF (3.10), AND (2.7) IS USED INSTEAD OF (2.14)    C
C    EXCEPT NEAR NORMALS TO VERTICES.                                  C
C                                                                      C
C      HSLEE (051500) : Change to Double precision Version.
C
C  IFLAG: input 1 to skip hyperboloid panels' calculations             C
C         input 0 to use hyperboloid panels for near field points      C
C                 THE IFLAG output will be 2 and no velocities         C
C                 calculations.                                        C
C----------------------------------------------------------------------C

      SUBROUTINE RPAN_PAR(Xi,Yi,Zi,CHRLENSi, FSs,FDd,FSXi,FSYi,
     %                  FDXi,FDYi,FDZi,ID,IFLAG, geom_data)
      use m_rpan_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      type(t_geom_rpan)::geom_data
      REAL XI,YI,ZI,CHRLENSI,FSS,FDD,FSXI,FSYI,FDXI,FDYI,FDZI
      DIMENSION XV(4),YV(4),S(15),SIDE(4)
      DIMENSION R(4),RR(4),RI(4),XRI(4),YRI(4),FE(4),
     *          B(5),XMXV(4),YMYV(4),N1(4)

C----------------------------------------------------------------------
C
C  DATA ARRAY B FOR 6D RATIONAL APPROXIMATION OF ARCTANGENT
C  SEE HART, ET AL, COMPUTER APPROXIMATIONS, WILEY, 1968, TABLE 5090
C
C----------------------------------------------------------------------
      DATA B/ 2.4091197D-01, 3.7851122D+00, 5.6770721D+00,
     *        5.6772854D+00, 5.6770747D+00/
C      DATA PI/ 3.1415927D+00 /, PI2/ 1.5707963D+00 /
C      DATA TWOPI/ 6.2831853D+00 /
      DATA ONE/ 1.D+00/, A3/ 3.D+00/, A5/ 5.D+00/
      DATA A7/ 7.D+00/, A9/ 9.D+00/, A11/ 11.D+00/
      DATA A14/ 14.D+00/, A35/ 35.D+00/
      DATA A49/ 49.D+00/, A63/ 63.D+00/, A99/ 99.D+00/
C      DATA ONE10/.1D+00/, ONE6/.1666667D+00/, ONE3/ .3333333D+00/
C      DATA FIVE3/ 1.666667D+00/, SEVEN3/2.333333D+00/, ZERO/ 0.0D0/
      DATA ZERO/ 0.0D0/
      DATA FT3/ 4.666667D+00/, TOL/ 1D-08 /
      DATA N1/ 2, 3, 4, 1 /

      pi = dacos(-1.d0)
      pi2 = 0.5d0*pi
      twopi = 2.d0*pi
      one10 = 0.1d0
      one6 = 1.d0/6.d0
      one3 = 1.d0/3.d0
      five3 = 5.d0/3.d0
      seven3 = 7.d0/3.d0


C --- Change Input variable to double precision

      x = dble(xi)
      y = dble(yi)
      z = dble(zi)
      chrlens = dble(chrlensi)
      do i = 1 , 4
      xv(i) = dble(geom_data%xv1(i))
      yv(i) = dble(geom_data%yv1(i))
      side(i) = dble(geom_data%side1(i))
      enddo
      do i = 1 , 15
        s(i) = dble(geom_data%sq1(i))
      enddo

      XMXC=X-S(6)
      YMYC=Y-S(2)
      XX=XMXC*XMXC
      YY=YMYC*YMYC
      ZZ=Z*Z
      RRC=XX+YY+ZZ
      IF (RRC.LT.100.d0*CHRLENS) THEN
         IF(IFLAG.EQ.1) THEN
            GO TO 11
         ELSE IF (IFLAG.EQ.0) THEN
            IFLAG=2
            RETURN
         END IF
      END IF
C----------------------------------------------------------------------
C
C   TWO-TERM MULTIPOLE EXPANSION INCLUDING SECOND MOMENTS
C
C----------------------------------------------------------------------
      R2=ONE/RRC
      R1=DSQRT(R2)
      R3=R1*R2
      R5=R3*R2
      ZR2=Z*R2
      XY=XMXC*YMYC
      SS1=S(1)*R1
      SS3=-(S(3)+S(10))*R3
      SS5=(XX*S(10)+XY*S(7)+YY*S(3))*R5
      FS=SS1+ONE3*SS3+SS5
      FDSUM=SS1+SS3+A5*SS5
      FD=ZR2*FDSUM
C----------------------------------------------------------------------
C    DELETE NEXT  14  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 8
      RSS3=R2*SS1
      SSX3=-XMXC*RSS3
      SSY3=-YMYC*RSS3
      SSX5=(XMXC*(S(3)+A3*S(10))+YMYC*S(7))*R5
      SSY5=(YMYC*(S(10)+A3*S(3))+XMXC*S(7))*R5
      A5R2=A5*R2
      RSS7=-A5R2*SS5
      SSX7=XMXC*RSS7
      SSY7=YMYC*RSS7
      FSX=SSX3+SSX5+SSX7
      FSY=SSY3+SSY5+SSY7
      FDX=ZR2*(A3*SSX3+A5*SSX5+A7*SSX7)
      FDY=ZR2*(A3*SSY3+A5*SSY5+A7*SSY7)
      ZZR4=ZR2*ZR2
      FDZ=R2*FDSUM-ZZR4*(A3*SS1+A5*SS3+A35*SS5)
  8   IF (RRC.GT.15.d0*CHRLENS) GO TO 99
C----------------------------------------------------------------------
C
C    THIRD AND FOURTH MOMENTS ADDED FOR RRC/AREA BETWEEN 40 AND 150
C
C----------------------------------------------------------------------
      S914=S(9)+S(14)
      S813=S(8)+S(13)
      S411=S(4)+S(11)
      S512=S(5)+S(12)
      S1215=S(12)+S(15)
      R7=R5*R2
      R9=R7*R2
      SS5=(-XMXC*S813-YMYC*S411+ONE10*(S512+S1215))*R5
      SS7=(FIVE3*((XMXC*XX*S(13)+YMYC*YY*S(4))+A3*XY*(XMXC*S(11)
     *+YMYC*S(8)))-XX*S1215-YY*S512-XY*S914)*R7
      SS9=(A7*(ONE6*(XX*XX*S(15)+YY*YY*S(5))+XX*YY*S(12))
     *+SEVEN3*XY*(XX*S(14)+YY*S(9)))*R9
      FS=FS+SS5+SS7+SS9
      FDSUM=A5*SS5+A7*SS7+A9*SS9
      FD=FD+ZR2*FDSUM
C----------------------------------------------------------------------
C    DELETE NEXT  20  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 99
      TXY=XY+XY
      SSX5=-S813*R5
      SSY5=-S411*R5
      RSS7=A5R2*SS5
      SSX7=(A5*(XX*S(13)+TXY*S(11)+YY*S(8))-S1215*(XMXC+XMXC)
     *-YMYC*S914)*R7-XMXC*RSS7
      SSY7=(A5*(YY*S(4)+XX*S(11)+TXY*S(8))-S512*(YMYC+YMYC)
     *-XMXC*S914)*R7-YMYC*RSS7
      RSS9=A7*SS7*R2
      SSX9=(FT3*XMXC*XX*S(15)+A14*XMXC*YY*S(12)+A49*YMYC*(XX*S(14)
     *+ONE3*YY*S(9)))*R9-XMXC*RSS9
      SSY9=(FT3*YMYC*YY*S(5)+A14*YMYC*XX*S(12)+A49*XMXC*(YY*S(9)
     *+ONE3*XX*S(14)))*R9-YMYC*RSS9
      RSS11=A9*SS9*R2
      SSX11=-XMXC*RSS11
      SSY11=-YMYC*RSS11
      FSX=FSX+SSX5+SSX7+SSX9+SSX11
      FSY=FSY+SSY5+SSY7+SSY9+SSY11
      FDX=FDX+ZR2*(A5*SSX5+A7*SSX7+A9*SSX9+A11*SSX11)
      FDY=FDY+ZR2*(A5*SSY5+A7*SSY7+A9*SSY9+A11*SSY11)
      FDZ=FDZ+R2*FDSUM-ZZR4*(A35*SS5+A63*SS7+A99*SS9)
      GO TO 99
C----------------------------------------------------------------------
C
C    NEAR-FIELD SECTION USES EXACT FORMULATION
C      SET Z=TOL IF Z.LT.TOL TO AVOID INDETERMINACY ON PANEL
C      ZVTX IS USED TO DETERMINE PROXIMITY TO VERTEX NORMALS
C      MFLAG=1 IF NEAR VERTEX NORMALS
C
C----------------------------------------------------------------------
  11  FD=ZERO
      FS=ZERO
      ABZ=DABS(Z)
      IF(ABZ.GT.TOL) GO TO 12
        Z=TOL
        ZZ=Z*Z
        ABZ=TOL
  12  ZVTX=1.005D0*ABZ
      MFLAG=0
C----------------------------------------------------------------------
C    DELETE NEXT 5 LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      FSX=ZERO
      FSY=ZERO
      FDX=ZERO
      FDY=ZERO
      FDZ=ZERO
C----------------------------------------------------------------------
C
C    LOOP FOR CORNER FUNCTIONS
C
C----------------------------------------------------------------------
      DO 13 N=1,4
        XMXV(N)=X-XV(N)
        YMYV(N)=Y-YV(N)
        XX=XMXV(N)*XMXV(N)
        YY=YMYV(N)*YMYV(N)
        FE(N)=ZZ+XX
        RR(N)=FE(N)+YY
        R(N)=DSQRT(RR(N))
        IF (R(N).LT.ZVTX) MFLAG=1
        RI(N)=ONE/R(N)
C----------------------------------------------------------------------
C    DELETE NEXT   3   LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 13
        XRI(N)=RI(N)*XMXV(N)
        YRI(N)=RI(N)*YMYV(N)
   13 CONTINUE
C----------------------------------------------------------------------
C
C    LOOP FOR SIDE FUNCTIONS AND SUMS OVER FOUR SIDES
C
C----------------------------------------------------------------------
      DO 33 N=1,4
        IF (SIDE(N).LT.TOL) GO TO 33
        SIDI=ONE/SIDE(N)
        CT=(XV(N1(N))-XV(N))*SIDI
        ST=(YV(N1(N))-YV(N))*SIDI
        V=XMXV(N)*ST-YMYV(N)*CT
        VV=V*V
        RADS=VV+ZZ
        U1=XMXV(N)*CT+YMYV(N)*ST
        U2=XMXV(N1(N))*CT+YMYV(N1(N))*ST
        RSUM=R(N)+R(N1(N))
        FLAG=RI(N)*RI(N1(N))*U1*U2
C----------------------------------------------------------------------
C
C       FLAG=1 ON EXTENSIONS, -1 ON SIDES
C         IN FOLLOWING SUBSECTIONS FS,FSX,FSY,FDZ ARE EVALUATED FROM
C           LAST FORM OF (3.9) IN NORMAL CASE
C         ELSE
C           FIRST FORM OF (3.9) IF NEAR SIDE OF PANEL
C
C----------------------------------------------------------------------
        IF (FLAG.GT.-.99D0) THEN
          RSP=RSUM+SIDE(N)
          RSM=RSUM-SIDE(N)
          FLN=DLOG(RSP/RSM)
        ELSE
          RU1=R(N)+U1
          RU2=R(N1(N))-U2
          RADI=ONE/RADS
          FLN=DLOG(RU1*RU2*RADI)
        ENDIF
        FS=FS+V*FLN
C----------------------------------------------------------------------
C    DELETE FOLLOWING LINES TO NEXT ENDIF TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 14
        IF (FLAG.GT.-.99D0) THEN
          FAC=V*(SIDE(N)+SIDE(N))/(RSP*RSM)
          FSX=FSX+FLN*ST-FAC*(XRI(N)+XRI(N1(N)))
          FSY=FSY-FLN*CT-FAC*(YRI(N)+YRI(N1(N)))
          FDZ=FDZ-FAC*(RI(N)+RI(N1(N)))
        ELSE
          RU1I=ONE/RU1
          RU2I=ONE/RU2
          FA=RU1I-RU2I
          FB=-(V+V)*RADI
          FSX=FSX+FLN*ST+V*(FA*CT+FB*ST+RU1I*XRI(N)+RU2I*XRI(N1(N)))
          FSY=FSY-FLN*CT+V*(FA*ST-FB*CT+RU1I*YRI(N)+RU2I*YRI(N1(N)))
          FDZ=FDZ+FB+V*(RU1I*RI(N)+RU2I*RI(N1(N)))
        ENDIF
C----------------------------------------------------------------------
C
C        IN FOLLOWING SUBSECTIONS FACTORS IN (2.15) ARE EVALUATED FROM
C          (2.7) IN NORMAL CASE
C        ELSE
C          (2.14) IF NEAR NORMAL TO A VERTEX
C
C----------------------------------------------------------------------
   14   IF (MFLAG.EQ.0) THEN
          S1=V*R(N)
          C1=ABZ*U1
          S2=V*R(N1(N))
          C2=ABZ*U2
        ELSE
          FH1=XMXV(N)*YMYV(N)
          FH2=XMXV(N1(N))*YMYV(N1(N))
          S1=FE(N)*ST-FH1*CT
          C1=ABZ*R(N)*CT
          S2=FE(N1(N))*ST-FH2*CT
          C2=ABZ*R(N1(N))*CT
        ENDIF
        S12=S1*C2-S2*C1
        C12=C1*C2+S1*S2
C----------------------------------------------------------------------
C
C    EVALUATE THIRD ARCTANGENT IN (2.15)
C         ANGLE (MODULO PI) BETWEEN -PI/4 AND PI/4
C       ELSE
C         USE INVERSE COTANGENT AND ADD/SUBTRACT PI/2
C
C----------------------------------------------------------------------
        IF (DABS(S12).LE.DABS(C12)) THEN
          U=S12/C12
          IF (C12.LT.ZERO) FD=FD+DSIGN(PI,S12)
        ELSE
          U=-C12/S12
          FD=FD+DSIGN(PI2,S12)
        ENDIF
        UU=U*U
        FD=FD+U*((B(1)*UU+B(2))*UU+B(3))/((UU+B(4))*UU+B(5))
C----------------------------------------------------------------------
C
C    FOLLOWING THREE SECTIONS EVALUATE FDX,FDY FOR
C       FIELD POINT NEAR NORMAL TO VERTEX (MFLAG=1)
C       NORMAL CASE (FLAG.LT.0.99)
C       NEAR EXTENSIONS OF SIDES (FLAG.GE.0.99)
C
C    DELETE ALL FOLLOWING LINES ABOVE LABEL 33 TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 33
        IF (MFLAG.EQ.0) GO TO 20
          FAC=C1/((C1*C1+S1*S1)*RR(N))
        IF(Z.LT.ZERO) FAC=-FAC   ! Fix sign of derivatives, Newman
          FDX=FDX+(RR(N)*V+FH1*U1)*FAC
          FDY=FDY-FE(N)*U1*FAC
          FAC=C2/((C2*C2+S2*S2)*RR(N1(N)))
        IF(Z.LT.ZERO) FAC=-FAC   ! Fix sign of derivatives, Newman
          FDX=FDX-(RR(N1(N))*V+FH2*U2)*FAC
          FDY=FDY+FE(N1(N))*U2*FAC
          GO TO 33
  20    IF (FLAG.LT.0.99D0) THEN
          U1V=U1*V
          FAC=Z/(C1*C1+S1*S1)
          FDX=FDX+(U1V*XRI(N)+R(N)*YMYV(N))*FAC
          FDY=FDY+(U1V*YRI(N)-R(N)*XMXV(N))*FAC
          U2V=U2*V
          FAC=Z/(C2*C2+S2*S2)
          FDX=FDX-(U2V*XRI(N1(N))+R(N1(N))*YMYV(N1(N)))*FAC
          FDY=FDY-(U2V*YRI(N1(N))-R(N1(N))*XMXV(N1(N)))*FAC
        ELSE
          ZS=Z*SIDE(N)
          USUM=U1+U2
          VRADS=V*RADS
          SFAC=VRADS*USUM
          SFS=-SFAC*ZS
          SFA=SFAC*C12
          CFAC=U2*R(N)+U1*R(N1(N))
          SFB=SFAC*CFAC
          CCF=C12*CFAC
          PA=(CCF+CCF)*VRADS-SFA*RSUM-SFB*ZZ*USUM
          PB=CCF*USUM*(VV+VV+RADS)-SFB*(S1+S1)*R(N1(N))
          PC=-SFA*U2-SFB*VV*R(N1(N))
          PD=-SFA*U1-SFB*VV*R(N)
          FAC=ZS/(CCF*CCF+SFS*SFS)
          FDX=FDX-(PA*CT+PB*ST+PC*XRI(N)+PD*XRI(N1(N)))*FAC
          FDY=FDY-(PA*ST-PB*CT+PC*YRI(N)+PD*YRI(N1(N)))*FAC
        ENDIF
  33  CONTINUE
      IF (FD.LT.ZERO) FD=FD+TWOPI
      IF (Z.LT.ZERO) FD=-FD
      FS=FS-Z*FD
C----------------------------------------------------------------------
C    DELETE NEXT  3  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 99
      FSX=FSX-Z*FDX
      FSY=FSY-Z*FDY

99    CONTINUE

      fdd = sngl(fd)
      fss = sngl(fs)
      if(id .ne. 0) then
        fsxi = sngl(fsx)
        fsyi = sngl(fsy)
        fdxi = sngl(fdx)
        fdyi = sngl(fdy)
        fdzi = sngl(fdz)
      endif

      RETURN
      END

