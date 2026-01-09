! ----------------------------------------------------------------------------- !
!                       CREATED BY YIRAN SU 10/04/2017                          !
!                   Global parameters for unsteady PC2NS!                       !
! ----------------------------------------------------------------------------- !
      module M_UNS_GLOBAL_PAR
        public  :: UNS_GLOBAL_PAR

        ! Basic parameter in propcav
        integer pc_nbld
        integer pc_ntpr
        integer pc_nc,pc_ncp
        integer pc_mr,pc_mrp
        integer pc_nh,pc_nhp
        integer pc_npn

        ! Parameter for structured mesh
        integer mesh_nx, mesh_nr, mesh_nt

        double precision :: uns_pi, uns_2pi

        ! When Ueff is activated
        integer :: pc2ns_start = 4
        ! Ueff relaxation factor
        double precision :: ueReF = 1.0

        ! Number of PROPCAV inner revolution per time step
        ! Fluent inner iteration number should be (PC_NREV+1)*ITER_PER_REV
        ! 1 blade case pc_nrev = 9, ITER_PER_REV = 2, total fluent inner iteration = 20
        ! 4 blade case pc_nrev = 6, ITER_PER_REV = 3, total fluent inner iteration = 21
        integer :: pc_nrev = 6


        ! Number of seed to interpolate Ut/Uind and the distance between two seeds
        integer :: nSeed = 6
        integer :: useLayer = 6
        double precision :: hSeed = 0.01

      contains

        subroutine UNS_GLOBAL_PAR
          include 'PUFCAV.INC'
          include 'PUFCAVB.INC'
          character (len = 99) :: cTemp

          ! initialize
          pc_nbld = nblade
          pc_ntpr = ntprev
          pc_nc   = nc
          pc_mr   = mr
          pc_nh   = nh
          pc_ncp  = nc + 1
          pc_mrp  = mr + 1
          pc_nhp  = nh + 1
          pc_npn  = npanel

          open(286, file = 'mesh3d.plt', status = 'old')
          read(286, "(1X, A)") cTemp
          read(cTemp, *) mesh_nx, mesh_nr, mesh_nt
          close(286)

          if (mesh_nt.ne.pc_ntpr) write(*,*) 'Error!! Mesh_NT!=NTPREV'

          uns_pi  = acos(0.0d0)*2.0d0
          uns_2pi = uns_pi * 2.0d0

        end subroutine

      end module
