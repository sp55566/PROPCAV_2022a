      PROGRAM PROPCAV_UNS
************************************************************************
*                                                                      *
*     PROPCAV_UNS                                                      *
*     --- Time dependent version of PROPCAV                            *
*     --- Make use of ICAVMAX to do the inner iteration                *
*                                                                      *
*     By Yiran Su 09/19/2017                                           *
*                                                                      *
************************************************************************
      use M_UNS_GLOBAL_PAR, only : PC_NREV
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      CHARACTER*29 FNLOG,FNPLT

C-----------------------------------------------------------------------
C     Preset parameters
C-----------------------------------------------------------------------
! - Yiran Su 09/19/2017 --
      IUNS = 1      ! Unsteady mode
      IVISC = 0     ! Boundary layer solver is not supported
      NTREVC = 0    ! Number of cavitating revolutions
      IWET = 1
      NSUB = 5
      NWSUB1 = 4
      TWOPI = 4.0 * DACOS(0.0D0)
C-----------------------------------------------------------------------
C     Release version
C-----------------------------------------------------------------------
      WRITE(*,'(A)') '  '
      WRITE(*,'(A)') '  **********************************************'
      WRITE(*,'(A)') '  *                                            *'
      WRITE(*,'(A)') '  *      PROPCAV RELEASE VERSION 3.3           *'
      WRITE(*,'(A)') '  *      Released    October,  2016            *'
      WRITE(*,'(A)') '  *                                            *'
      WRITE(*,'(A)') '  **********************************************'
      WRITE(*,'(A)') '  '

C-----------------------------------------------------------------------
C     Read and process input data
C-----------------------------------------------------------------------
      WRITE(*,'(A,$)') ' PROPCAV> ENTER DELTAT (in degrees): '
      READ(*,*) DELTAT
      NDLTAT=INT(DELTAT)
      DELTAT=DELTAT*RAD

!s--YE TIAN for m_param
      NSTEP = 360/NDLTAT
      NNDIM = NSTEP
      NDEL  = NDLTAT

C---------What is the cut off for extrapolation? --------------- JY081600
      WRITE(*,*)
      WRITE(*,'(A,$)') ' PROPCAV> ENTER RADII TO CUT THE TIP (RMRTIP): '
      WRITE(*,*)
      READ(*,*) RMRTIP

C --- How many time steps?
      WRITE(*,*) ' PROPCAV> ENTER NUMBER OF REVOLUTIONS FOR FULLY',
     *  ' WETTED COMPUTATION'
      WRITE(*,*) '          (ignored if ISP=1 or ISTEADY=0):'
      READ(*,*) NTREVW

C---------This subroutine reads in the propeller geometry---------------
      CALL PROINP

C--- YE TIAN 07/09/2013----
      NZWSUB=NWZ+NDEL

C --- allocate dynamic memory - Yiran Su 09/19/2017 --
      CALL ALLOC

C ---- Check to see if DELK can be divided by DELTAT.              JY111900
      if(360/ndltat*ndltat .ne. 360) then
        write(*,*) ' ================ WARNING ====================='
        write(*,*) ' --> 360 degree can not divided by input DELTAT'
        write(*,*) ' =============================================='
      endif
      if (ntprev/nblade*nblade .ne. ntprev) then
        write(*,*) ' ================= ERROR! ====================='
        write(*,*) ' --> NTPREV can not divided by input NBLADE    '
        write(*,*) ' =============================================='
        stop
      end if
      if (ITERKT.ne.0) then
        write(*,*)'Warning: To ensure numerical stability, it is highly',
     *            'recommended to set ITERKT = 0 in PROPCAV_UNS runs.'
      end if

C-----------------------------------------------------------------------
C     Generate the geometries of propeller blade (GBLADE), hub (GHUB)
C     and wake (GWAKE)
C-----------------------------------------------------------------------
      CALL CHRLEN(FN,LENCH)
      FNPLT=FN(1:LENCH)//'.plt'
      OPEN(700, FILE=FNPLT, STATUS='UNKNOWN')
      WRITE(700,*) 'VARIABLES = "X", "Y", "Z"'

      CALL PROGEO

C---- Check parameters for dimensioning errors.  Results are printed to ERR.LOG.
      IF(ISCAV.EQ.0) CALL CHECK_PARAM2

C-----------------------------------------------------------------------
C     NTREVW:  Number of revolutions
C     NTPREV:  Number of time steps per revolution
C     NTREV is the total number of revolutions.
C     NTIME is the total number of time steps.
C-----------------------------------------------------------------------
      NTPREV = 360 / NDLTAT !Number of Time Steps in one revolutions (= 360/NDLTAT)
      NTREV=NTREVW
      NTIME=NTREV*NTPREV  ! NTPREV = 360 / NDLTAT

C-----------------------------------------------------------------------
C     Time step  0: steady case
C     Time step  1: beginning the unsteady effect, the propeller is
C     at theta=0.0
C-----------------------------------------------------------------------
C.....NREV:   current number of revolution
C.....IDXREV: number of time step in current revolution

      NREV=1
      WRITE(*,*) '------Start fully-wetted computation---------------'


C --- Run the steady step !!!
      ITSTEP_KEY = 0
      IDXREV_KEY = 0
      ITSTEP = 0
      TSTEP = 0.
      NTSTEP = 0
      IDXREV = 0
      NTSTEP_KEY = 0
      CALL INDPOT(1)
      CALL UNS0_INIT     !!! Initialization of the PC2UNS model
      CALL SIMEQN
      DO M=1,MR
        DPHI(M,IDXREV)=DELP(M)
      END DO
      CALL FORCEV
      CALL WRITEFILEUNS(0,1)
      CALL UNS1_END_STEADY   !!! Calculate body force; write body force file
C --- Finished steady step !!!


C --- Start unsteady calculation
      DO N=1,NTIME  ! NTIME = NTREV * NTPREV
        ICAVMAX = PC_NREV * NBLADE

        NTSTEP_KEY = N
        IF(N.LE.1) THEN
          TSTEP_KEY = 0.
          ITSTEP_KEY = 0
          IDXREV_KEY = N
        ELSE
          TSTEP_KEY = TSTEP_KEY + DELTAT
          ITSTEP_KEY = ITSTEP_KEY + NDLTAT
          IF(ITSTEP_KEY .EQ. 360) THEN
            NREV = NREV + 1
            IDXREV_KEY = 0
            TSTEP_KEY = ZERO
            ITSTEP_KEY = 0
          END IF
          IDXREV_KEY = IDXREV_KEY + 1
        END IF

 1010   FORMAT('+ Time Step=',I4)
        WRITE(*,1010) N

        CALL UNS2_START_TSTEP(N)

C--   Do inner iteration over all blades -------
        DO ICAVT = 1,ICAVMAX
C --- Implement trick - mock NTSTEP, TSTEP, ITSTEP, IDXREV for non-key blades
          ITSTEP = MOD(ITSTEP_KEY + (ICAVT - 1) * 360 / NBLADE, 360)
          TSTEP = MOD(TSTEP_KEY + (ICAVT - 1) * TWOPI / NBLADE, TWOPI)
          NTSTEP = NTSTEP_KEY + (ICAVT - 1) * NTPREV / NBLADE
          IDXREV = MOD(IDXREV_KEY+(ICAVT-1)*NTPREV/NBLADE, NTPREV)
          IF (IDXREV.EQ.0) IDXREV = NTPREV

C --- Modify the matrix for new iteration
          CALL INDPOT(0) ! Do not Calculate I.C. and Read from saved files.

C --- Solve the simultaneous equation and create output files
          CALL SIMEQN
          DO M=1,MR
            DPHI(M,IDXREV)=DELP(M)
          END DO

          CALL UNS3_CALC_BF(IDXREV)    ! Calculate body force
          IF (MOD(ICAVT,NBLADE).EQ.0) THEN
            WRITE(*,'(A5,I2)',ADVANCE='NO') 'Rev #', ICAVT/NBLADE
            CALL UNS4_END_REV
          END IF

C --- Write output file (forces)
          IF(ICAVT.GT.(ICAVMAX-NBLADE)) THEN
            CALL FORCEV
            CALL WRITEFILEUNS(N,ICAVT)
          END IF

C --- Write output file (pressure distribution)
          IF(ICAVT.EQ.(ICAVMAX-NBLADE+1)) THEN
            CALL WRITEPRESSURE(N)
          END IF

        END DO
      END DO

C --- Close all files
      CALL WRITEFILEUNS(-1,0)


C-----------------------------------------------------------------------
C     End of fully-wetted time marching
C-----------------------------------------------------------------------
      WRITE(*,*) '------End fully-wetted computation-----------------'
      END
