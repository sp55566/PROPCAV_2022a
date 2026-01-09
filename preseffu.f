! ----------------------------------------------------------------------------- !
!                       CREATED BY YIRAN SU 02/18/2018                          !
!                        Calculate Effective Pressure                           !
! ----------------------------------------------------------------------------- !

      subroutine PRES_EFF_CALC(EPRES)
        use PC2NSVEL
        include 'PUFCAV.INC'
        include 'PUFCAVB.INC'
        real,dimension(nc,mr) :: EPRES, QEINT
        real,dimension(nc,mr) :: tuex1, tuey1, tuez1
        real,dimension(nc,mr) :: tuex0, tuey0, tuez0

        ! Initialization
        EPRES = 0.0
        QEINT = 0.0
        deltt = 2.0 / NTPREV * ADVCO

        ! Get effective wake velocity
        do j = 1, mr
          do i = 1, nc
            k = indexb(i,j)
            tttt = atan2(xct(k,3), xct(k,2))

            k = IDXREV
            tuex1(i,j) = UEFX(i,j,k)
            tuey1(i,j) = UEFR(i,j,k)*cos(tttt) - UEFT(i,j,k)*sin(tttt)
            tuez1(i,j) = UEFR(i,j,k)*sin(tttt) + UEFT(i,j,k)*cos(tttt)

            k = IDXREV - 1
            if (k.eq.0) k = NTPREV
            tuex0(i,j) = UEFX(i,j,k)
            tuey0(i,j) = UEFR(i,j,k)*cos(tttt) - UEFT(i,j,k)*sin(tttt)
            tuez0(i,j) = UEFR(i,j,k)*sin(tttt) + UEFT(i,j,k)*cos(tttt)
          end do
        end do

        ! Calculate effective pressure
        do j = 1, mr

          ! One side of blade
          i = nh + 1
          tuele = sqrt(tuex1(i,j)**2 + tuey1(i,j)**2 + tuez1(i,j)**2)
          do i = (nh+2), nc
            tuelc = sqrt(tuex1(i,j)**2 + tuey1(i,j)**2 + tuez1(i,j)**2)
            tuext = (tuex1(i,j) - tuex0(i,j)) / deltt
            tueyt = (tuey1(i,j) - tuey0(i,j)) / deltt
            tuezt = (tuez1(i,j) - tuez0(i,j)) / deltt

            k1 = indexb(i-1,j)
            k2 = indexb(  i,j)
            ttx = xct(k2,1) - xct(k1,1)
            tty = xct(k2,2) - xct(k1,2)
            ttz = xct(k2,3) - xct(k1,3)

            QEINT(i,j)= QEINT(i-1,j) + tuext*ttx + tueyt*tty + tuezt*ttz
            EPRES(i,j)= tuele**2 - tuelc**2 - 2.0 * QEINT(i,j)
          end do

          ! The other side of blade
          i = nh
          tuele = sqrt(tuex1(i,j)**2 + tuey1(i,j)**2 + tuez1(i,j)**2)
          do i = (nh-1), 1, -1
            tuelc = sqrt(tuex1(i,j)**2 + tuey1(i,j)**2 + tuez1(i,j)**2)
            tuext = (tuex1(i,j) - tuex0(i,j)) / deltt
            tueyt = (tuey1(i,j) - tuey0(i,j)) / deltt
            tuezt = (tuez1(i,j) - tuez0(i,j)) / deltt

            k1 = indexb(i+1,j)
            k2 = indexb(  i,j)
            ttx = xct(k2,1) - xct(k1,1)
            tty = xct(k2,2) - xct(k1,2)
            ttz = xct(k2,3) - xct(k1,3)

            QEINT(i,j)= QEINT(i+1,j) + tuext*ttx + tueyt*tty + tuezt*ttz
            EPRES(i,j)= tuele**2 - tuelc**2 - 2.0 * QEINT(i,j)
          end do

        end do
      end subroutine
