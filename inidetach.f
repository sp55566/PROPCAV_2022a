      SUBROUTINE INIDETACH
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      
C-- Read Pressure distribution on Blade at each time step
C   and find the initial cavity detachment points
C

      ADVCO2 = ADVCO**2

      IUNIT = 690

      READ(IUNIT,*)
      READ(IUNIT,*)
      
      DO M = 1 , MR
         READ(IUNIT,*)
         DO N = 1 , NC
            READ(IUNIT,*) DUM1, DUM2
            CPB(N,M) = -DUM2/ADVCO2
         ENDDO
      ENDDO
      
      DO IDXREV = 1 , NTPREV
         DO M = 1 , MR
            READ(IUNIT,*)
            DO N = 1 , NC
               READ(IUNIT,*) DUM1, DUM2
               CPB(N,M) = -DUM2/ADVCO2
            ENDDO
         ENDDO
         CALL DETACH
      ENDDO
      NTSTEP = 1
      
      RETURN
      END

