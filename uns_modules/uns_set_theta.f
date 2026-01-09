! ----------------------------------------------------------------------------- !
!                       CREATED BY YIRAN SU 10/04/2017                          !
!                       Location of the current blades !                        !
!                     Start from the Key blade, 2nd blade, ...                  !
!                   The index value corresponds to the RANS mesh,               !
!                                ,which is opposite to the direction of rotation!
! ----------------------------------------------------------------------------- !
      module M_UNS_THETA
        implicit none
        logical              :: isRotate    ! Does body force field rotate relative to RANS mesh? false means moving mesh in fluent
        double precision     :: rot_rans    ! Ranging from 1 to mesh_nt
        integer, allocatable :: it_rans(:)  ! Ranging from 1 to mesh_nt, not consistent with
        integer, allocatable :: it_prop(:)  ! Ranging from 1 to pc_ntpr


      contains

        subroutine UNS_SET_THETA(nstep)
          use M_UNS_GLOBAL_PAR
          implicit none
          integer nstep, i
          double precision :: delTheta

          delTheta = uns_2pi / mesh_nt  ! This parameter can be set, so that mesh_nt does not have to be the same as PC_NTPR
          isRotate = .false.            ! This parameter can be changed for research purpose; default is .false.


          ! Initialization
          if (nstep .eq. 0) then
            allocate(it_rans(pc_nbld))
            allocate(it_prop(pc_nbld))
            do i = 1, pc_nbld
              it_rans(i) = (i - 1) * mesh_nt / pc_nbld
              if(it_rans(i) < 1) it_rans(i) = it_rans(i) + mesh_nt
              it_prop(i) = 1 + (i - 1) * mesh_nt / pc_nbld
              if(it_prop(i) > mesh_nt) it_prop(i) = it_prop(i)-mesh_nt
            end do
            rot_rans = 0.0d0

          ! Update every time step
          else
            ! Update it_rans and it_prop
            do i = 1, pc_nbld
              if (isRotate) then
                it_rans(i) = it_rans(i) + 1
                if(it_rans(i) > mesh_nt) it_rans(i) = it_rans(i)-mesh_nt
              end if
              it_prop(i) = it_prop(i) + 1
              if(it_prop(i) > mesh_nt) it_prop(i) = it_prop(i) - mesh_nt
            end do
            ! Update rot_rans
            if (.not.isRotate) then
              rot_rans = rot_rans + delTheta;
              if (rot_rans .ge. uns_2pi) rot_rans = rot_rans - uns_2pi
            end if

          end if

          ! Write for Debugging Purpose
          ! write(*,'(8I4)') it_rans(1:4), it_prop(1:4)

        end subroutine

      end module
