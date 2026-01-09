C***********************************************************************
C
C     DOUBLE PRECISION
C
C***********************************************************************
C-----------------------------------------------------------------------
C
C  Calculate panel information for camber
C
C-----------------------------------------------------------------------
C
C  INPUT :
C  -----
C  NNN          : NUMBER OF PANELS
C  XV(I,J,K) ,  : ARRAY STORING THE K-TH REFERENCE COORDINATE
C                 OF THE J-TH VERTEX OF THE I-TH PANEL, WHERE
C                    I = 1,..,nnn
C                    J = 1,..,4           : THE PANEL VERTICES ARE
C                                           NUMBERED IN THE CLOCKWISE
C                                           SENSE WHEN VIEWED FROM THE
C                                           FLUID DOMAIN.
C                    K = 1,2,3 <=> x,y,z.
C
C
C  OUTPUT :  GEOMETRICAL PANEL DATA AND NORMAL VELOCITIES TRANSFERED
C  ------    VIA THE COMMON BLOCK / GEOMT /.
C
C           IER : ERROR INDEX
C               =  1 : ACCEPTABLE PANEL
C               =  0 : UNACCEPTABLE PANEL
C
C
C-----------------------------------------------------------------------

      module GEOCAMB

        private :: geom_camber

        public  :: UNS_CAM_GEO_INIT, indexCamber

        double precision, allocatable :: xvpc(:,:),yvpc(:,:),ssc(:,:),
     &                      sidc(:,:),chlenc(:),xctc(:,:),dirc(:,:,:)
!     &                    ,gsx(:,:,:),gsw(:,:),vel(:,:),delu(:),delv(:)
!     &                    ,cosphi(:),sinphi(:),ul(:,:),vl(:,:),wl(:,:)

        double precision :: vpp(4,3),xctp(3),dirp(3,3),s(15),side(4),
     &                      gaussx(4,3),gaussw(4),deu,dev,cosph,sinph,
     &                      ulc(3),vlc(3)

      contains

        subroutine UNS_CAM_GEO_INIT
          use M_BLADEM, only : xbm, ybm, zbm
          use M_UNS_GLOBAL_PAR, only : pc_nh, pc_mr
          implicit none
          double precision, dimension(pc_mr*pc_nh,4,3) :: xg
          integer n, m, i, ier
          do n = 1, pc_nh
            do m = 1, pc_mr
              i = indexCamber(n,m)
              xg(i,1,1) = xbm(n,m)
              xg(i,1,2) = ybm(n,m)
              xg(i,1,3) = zbm(n,m)
              xg(i,2,1) = xbm(n,m+1)
              xg(i,2,2) = ybm(n,m+1)
              xg(i,2,3) = zbm(n,m+1)
              xg(i,3,1) = xbm(n+1,m+1)
              xg(i,3,2) = ybm(n+1,m+1)
              xg(i,3,3) = zbm(n+1,m+1)
              xg(i,4,1) = xbm(n+1,m)
              xg(i,4,2) = ybm(n+1,m)
              xg(i,4,3) = zbm(n+1,m)
            end do
          end do
          call geom_camber(pc_mr*pc_nh, xg, ier)
          if (ier == 0) then
            write(*,*) "Unacceptable camber panel"
            stop
          end if
        end subroutine

        subroutine geom_camber(nnn,xv,ier)
          implicit none
          double precision :: xv(nnn,4,3)
          integer :: ier, nnn, i, j, k

          if (allocated(xvpc)) then
            deallocate(xvpc, yvpc, ssc, sidc,chlenc, xctc, dirc)
          end if
          allocate(xvpc(nnn,4), yvpc(nnn,4), ssc(nnn,15), sidc(nnn,4))
          allocate(chlenc(nnn), xctc(nnn,3), dirc(nnn,3,3))

          do i=1,nnn

            call gpanel2(xv(i,:,:), chlenc(i), ier)

            if(ier.eq.0) then
               write(*,*) ' ------ i: ',i
               do j=1,4
                 write(*,'(3f12.6)') (xv(i,j,k), k=1,3)
               end do
               return
             end if

            ssc(i,1:15) = s(1:15)
            xvpc(i,1:4) = vpp(1:4,1)
            yvpc(i,1:4) = vpp(1:4,2)
            sidc(i,1:4) = side(1:4)
            xctc(i,:) = xctp(:)
            dirc(i,:,:) = dirp(:,:)

!          gsw(i,1:4) = gaussw(1:4)
!          delu(i) = deu
!          delv(i) = dev
!          cosphi(i) = cosph
!          sinphi(i) = sinph
!          gsx(i,:,:) = gaussx(:,:)
!          vl(i,:) = vlc(:)
!          ul(i,:) = ulc(:)
C-----------------------------------------------------------------------
C        NORMAL VELOCITIES AT THE PANEL CENTROIDS normal vector outward from the fluid
C-----------------------------------------------------------------------
!          vel(i,1) = dirp(3,1)
!          vel(i,2) = dirp(3,2)
!          vel(i,3) = dirp(3,3)
!          vel(i,4) = xctp(2)*dirp(3,3) - xctp(3)*dirp(3,2)
!          vel(i,5) = xctp(3)*dirp(3,1) - xctp(1)*dirp(3,3)
!          vel(i,6) = xctp(1)*dirp(3,2) - xctp(2)*dirp(3,1)
          end do

        end subroutine

        integer function indexCamber(n, m)
          use M_UNS_GLOBAL_PAR, only : pc_mr
          implicit none
          integer n, m
          indexCamber = (n - 1) * pc_mr + m
        end function

      end module

