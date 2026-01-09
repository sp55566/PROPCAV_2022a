
      subroutine UNS0_INIT
        use M_UNS_F2M2P_CORR
        use M_UNS_GLOBAL_PAR, only : UNS_GLOBAL_PAR
        use M_UNS_PC2NS_BF, only : UNS_BF_INIT
        use M_UNS_PC2NS_VEL, only : UNS_VEL_INIT
        use M_UNS_THETA, only : UNS_SET_THETA
        use PC2NSVEL, only : ipc2ns
        use GEOCAMB, only : UNS_CAM_GEO_INIT
        implicit none
        ipc2ns = .false.      ! Use input .wak file for inflow
        call UNS_GLOBAL_PAR   ! Initialize global parameters
        call UNS_CAM_GEO_INIT ! Initialize camber panel geometry
        call UNS_READ_MESH    ! Read Mesh
        call UNS_READ_PROPCAV ! Read PROPCAV panels
        call sync_with_fluent ! =1=== SYNC to wait for Fluent cell info
        call UNS_READ_RANS    ! Read RANS Cells
        call UNS_BLD_R2M      ! Create RANS to Mesh Correlation
        call UNS_WRT_R2M      ! Write RANS to Mesh Correlation
        call sync_with_fluent ! =2=== SYNC to tell fluent corr.dat is ready
        call UNS_BLD_P2M      ! Create PROPCAV to Mesh Correlation
        call UNS_BLD_M2P      ! Create Mesh to PROPCAV Correlation
        call UNS_BF_INIT      ! Initialize body force arrays
        call UNS_VEL_INIT     ! Initialize velocity arrays
        call UNS_SET_THETA(0) ! Initialize two theta arrays
      end subroutine


      subroutine UNS1_END_STEADY
        use M_UNS_PC2NS_BF, only : UNS_BF_CALC
        implicit none
        call UNS_BF_CALC(0)   ! Calculate body force and assign it to all locations
      end subroutine


      subroutine UNS2_START_TSTEP(ntstep)
        use M_UNS_GLOBAL_PAR, only : pc2ns_start
        use PC2NSVEL, only : ipc2ns
        use M_UNS_THETA, only : UNS_SET_THETA
        use M_UNS_PC2NS_VEL, only : UNS_UE_EXTRA
        use M_UNS_PC2NS_BF, only : UNS_BF_EXTRA, UNS_WRT_BF
        implicit none
        integer ntstep
        if (ntstep.ge.pc2ns_start) ipc2ns = .true. ! Use Ueff from PC2NS as inflow
        if (ntstep.ne.1) call UNS_SET_THETA(1)  ! Update blade angles
        call UNS_UE_EXTRA                       ! Extrapolate the effective wake at the new location
        call UNS_BF_EXTRA                       ! Extrapolate the body force at the new location
        call UNS_WRT_BF                         ! Write body force ( this is only a initial guess, accuracy does not matter)
        call sync_with_fluent                   ! =3A== SYNC to tell BF/SC field bfout.dat is ready
        call sync_with_fluent                   ! =3B== SYNC to wait for Fluent finishing reading bfout.dat
      end subroutine


      subroutine UNS3_CALC_BF(idxrev)
        use M_UNS_PC2NS_BF, only : UNS_BF_CALC
        implicit none
        integer idxrev
        call UNS_BF_CALC(idxrev)  ! Calculate body force
      end subroutine


      subroutine UNS4_END_REV
        use M_UNS_PC2NS_BF, only : UNS_BF_CALC,UNS_WRT_BF
        use M_UNS_PC2NS_VEL, only : UNS_UT_CALC,UNS_UE_CALC,UNS_UI_CALC,
     &                              WRT_VEL_FIELD
        implicit none
        call UNS_WRT_BF            ! Write body force
        call sync_with_fluent      ! =4=== SYNC to make sure UTOT is ready
        call UNS_UT_CALC           ! Read total flow
        call sync_with_fluent      ! =5=== SYNC to notify reading UTOT is completed and bfout.dat can be overwritten
        call WRT_VEL_FIELD         ! Write velocity output
        call UNS_UI_CALC           ! Calculate Uind inside fluid domain for all the blade
        call UNS_UE_CALC           ! Calculate Effective wake
      end subroutine


C ***************************************************************
C     Synchronize with fluent
C ***************************************************************
      subroutine sync_with_fluent
        implicit none
        logical iex
        integer er
        do
          inquire(file = 'sync.tag', exist = iex)
          if (iex) exit
        end do
        open(unit = 286, file = 'sync.tag', status = 'old', iostat = er)
        if (er.ne.0) write(*,*) 'Error occurred in opening sync.tag!'
        close(unit = 286, status = 'delete', iostat = er)
        if (er.ne.0) write(*,*) 'Error occurred in deleting sync.tag!'
C        write(*,*) 'Sync complete -- Propcav --'
      end subroutine
