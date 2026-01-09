       SUBROUTINE CHECK_PARAM

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       dimension ierb(10)

       write(99,*)
       write(99,*) ' --------------------------------------'
       write(99,*) '  (2) Checking ... Blade Part          '
       write(99,*) ' --------------------------------------'
       write(99,*)

       CALL ICLEAR(IERB,10)

       if(nx .gt. nxmax) ierb(1) = 1
       write(99,10) nxmax , nx,ierb(1)

       if(nblade .gt. kz) ierb(2) = 1
       write(99,15) KZ , nblade,ierb(2)

       if(nc .gt. nbz) ierb(3) = 1
       write(99,20) NBZ , nc,ierb(3)

       if(mr .gt. mbz) ierb(4) = 1
       write(99,25) mbz , mr,ierb(4)

       if(isc .eq. 1) then
          if(nsr2.gt.nzsr2) ierb(5)=1
          write(99,30) nzsr2,nsr2,ierb(5)
       endif

 10    format(' Max NX     = ',i3,'    Input NX     = ',i3,
     %      '    Err = ',i1) 
 15    format(' Max NBLADE = ',i3,'    Input NBLADE = ',i3,
     %      '    Err = ',i1) 
 20    format(' Max NC     = ',i3,'    Input NC     = ',i3,
     %      '    Err = ',i1) 
 25    format(' Max MR     = ',i3,'    Input MR     = ',i3,
     %      '    Err = ',i1) 
 30    format(' Max NSR2   = ',i3,'    Input NSR2   = ',i3,
     %      '    Err = ',i1) 

       isum = 0
       do i = 1 , 10
          isum = isum + ierb(i)
       enddo
       if(isum .gt. 0) then
          write(*,*) ' ===================================='
          write(*,*)
          write(*,*) '  You Have ERROR in INPUT DATA FILE! '
          write(*,*) '  Check <ERR.LOG> File to check error'
          write(*,*)
          write(*,*) ' ===================================='
          stop
       endif

       RETURN
       END


C ===================================================================
       SUBROUTINE CHECK_PARAM0
C
C     This subroutine checks the recommended min. values for parameters
C     depending on the special options in PARAM.INC, 
C
C     it only gives you warning if the values are too small 
C     except for the case of the non-specified values.   
C      
C ====================================================================

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

       dimension ierb(20)

       open(99,file='ERR.LOG',status='unknown')

       write(99,*)
       write(99,*) ' --------------------------------------'
       write(99,*) '  (1) Checking ... PARAM.INC           '
       write(99,*) ' --------------------------------------'
       write(99,*)

       CALL ICLEAR(IERB,20)

C --- Global variables : 

C   KZ    : Max. Number of Blade.
C  NXMAX  : Max. Number of Input Radii (Previous it was a fixed value as 15).
C   NBZ   : Max. Number of chordwise panels.
C   MBZ   : Max. Number of Spanwise Panels.
C   NSTEP : Max. Number of Time step per Revolution.
C   NWZ   : Max. Number of streawise panels on wake. 
C 
C   NHDZ  : Max. number of streamwise panels in the downstream hub.  
C   NHUZ  : Max. number of streamwise panels in the upstream hub.  
C   MHBZ  : Max. number of circumferential panels between adjacent blades.
C
C   NZSR  : Max. number of panels in the separated region behind the blade trailing edge
C            (Need only for ISC = 1).
C   MCAVM : Max. Number of circumferential  Panels on the tip vortex.
C  
C Setting for the min. values

       IKZ = 1
       INXMAX = 10
       INBZ = 20
       IF(ICON .EQ. 0) THEN
          IMBZ = 10
       ELSE
          IMBZ = 5
       ENDIF

       IF(IAN .EQ. 2) THEN
          INSTEP = 60
          INWZ = 100
       ELSE
          INSTEP = 30
          INWZ = 60
       ENDIF

       IF(IHUB .EQ. 0) THEN
          INHDZ = 0
          INHUZ = 0
       ELSE
          INHDZ = INWZ
          INHUZ = 10
       ENDIF

       IMHBZ = 5 
       INZSR = 10
       IMCAVM = 8
                 

       IF(KZ .LT. IKZ) IERB(1) = 1
       WRITE(99,10) IKZ, KZ, IERB(1)

       IF(NXMAX .LE. 0) THEN
          IERB(2) = 1
       ELSEIF(NXMAX .GT. 0 .and. NXMAX .LT. INXMAX) THEN
          IERB(2) = -1
       ENDIF
       WRITE(99,20) INXMAX, NXMAX, IERB(2)

       IF(NBZ .LE. 0) THEN
          IERB(3) = 1
       ELSEIF(NBZ .GT. 0 .and. NBZ .LT. INBZ) THEN
          IERB(3) = -1
       ENDIF
       WRITE(99,30) INBZ, NBZ, IERB(3)

       IF(MBZ .LE. 0) THEN
          IERB(4) = 1
       ELSEIF(MBZ .GT. 0 .and. MBZ .LT. IMBZ) THEN
          IF(ICON .EQ. 0) IERB(4) = -1
       ENDIF
       WRITE(99,40) IMBZ, MBZ, IERB(4)

       IF(NSTEP .LE. 0) THEN
          IERB(5) = 1
       ELSEIF(NSTEP .GT. 0 .and. NSTEP .LT. INSTEP) THEN
          IERB(5) = -1
       ENDIF
       WRITE(99,50) INSTEP, NSTEP, IERB(5)

       IF(NWZ .LE. 10) THEN
          IERB(6) = 1
       ELSEIF(NWZ .GT. 10 .and. NWZ .LT. INWZ) THEN
          IERB(6) = -1
       ENDIF
       WRITE(99,60) INWZ, NWZ, IERB(6)

       IF(IHUB .EQ. 0) THEN
          IF(NHDZ .NE. INHDZ) IERB(7) = -1
          IF(NHUZ .NE. INHUZ) IERB(8) = -1
          WRITE(99,70) INHDZ, NHDZ, IERB(7)
          WRITE(99,80) INHUZ, NHUZ, IERB(8)
       ELSE
          IF(NHDZ .LE. 0) THEN
             IERB(7) = 1
          ELSEIF(NHDZ .GT. 0 .AND. NHDZ .LT. INHDZ) THEN
             IERB(7) = -1
          ENDIF
          IF(NHUZ .LE. 0) THEN
             IERB(8) = 1
          ELSEIF(NHUZ .GT. 0 .AND. NHUZ .LT. INHUZ) THEN
             IERB(8) = -1
          ENDIF
          WRITE(99,75) INHDZ, NHDZ, IERB(7)
          WRITE(99,85) INHUZ, NHUZ, IERB(8)
       ENDIF

       IF(MHBZ .LE. 0) THEN
          IERB(9) = 1
       ELSEIF(MHBZ .GT. 0 .AND. MHBZ .LT. IMHBZ) THEN
          IERB(9) = -1
       ENDIF
       WRITE(99,90) IMHBZ, MHBZ, IERB(9)

       IF(ISC .EQ. 0) THEN
          IF(NZSR .LT. 0) THEN
             IERB(10) = 1
          ELSEIF(NZSR .GT. 0) THEN
             IERB(10) = -1
          ENDIF
          WRITE(99,100) 0, NZSR, IERB(10)          
       ELSEIF(ISC .EQ. 1) THEN
          IF(NZSR .LE. 0) THEN
             IERB(10) = 1
          ELSEIF(NZSR .GT. 0 .AND. NZSR .LT. INZSR) THEN
             IERB(10) = -1
          ENDIF 
          WRITE(99,105) INZSR, NZSR, IERB(10)          
       ENDIF

       IF(IAN. NE. 2) THEN
          IF(MCAVM .LT. 0) THEN
             IERB(11) = 1
          ELSEIF(MCAVM .GT. 0) THEN
             IERB(11) = -1
          ENDIF
          WRITE(99,110) 0, MCAVM, IERB(11)          
       ELSE
          IF(MCAVM .LE. 0) THEN
             IERB(11) = 1
          ELSEIF(MCAVM .GT. 0 .AND. MCAVM .LT. IMCAVM) THEN
             IERB(11) = -1
          ENDIF 
          WRITE(99,115) IMCAVM, MCAVM, IERB(11)          
       ENDIF          

 10    format(' Recommended  KZ   >=  ',i3,'    Input KZ     = ',i4,
     %      '    Err = ',i2)

 20    format(' Recommended NXMAX >=  ',i3,'    Input NXMAX  = ',i4,
     %      '    Err = ',i2)

 30    format(' Recommended  NBZ  >=  ',i3,'    Input NBZ    = ',i4,
     %      '    Err = ',i2)

 40    format(' Recommended  MBZ  >=  ',i3,'    Input MBZ    = ',i4,
     %      '    Err = ',i2)

 50    format(' Recommended NSTEP >=  ',i3,'    Input NSTEP  = ',i4,
     %      '    Err = ',i2)

 60    format(' Recommended  NWZ  >=  ',i3,'    Input  NWZ   = ',i4,
     %      '    Err = ',i2,'  * Check part (3) for the correct value')

 70    format(' Recommended NHDZ   =  ',i3,'    Input NHDZ   = ',i4,
     %      '    Err = ',i2,'  * Use Recommended Value ' )

 75    format(' Recommended NHDZ  >=  ',i3,'    Input NHDZ   = ',i4,
     %      '    Err = ',i2)

 80    format(' Recommended NHUZ   =  ',i3,'    Input NHUZ   = ',i4,
     %      '    Err = ',i2,'  * Use Recommended Value ' )

 85    format(' Recommended NHUZ  >=  ',i3,'    Input NHUZ   = ',i4,
     %      '    Err = ',i2)

 90    format(' Recommended MHBZ  >=  ',i3,'    Input MHBZ   = ',i4,
     %      '    Err = ',i2)

 100   format(' Recommended NZSR   =  ',i3,'    Input NZSR   = ',i4,
     %      '    Err = ',i2,'  * Use Recommended Value ' )

 105   format(' Recommended NZSR  >=  ',i3,'    Input NZSR   = ',i4,
     %      '    Err = ',i2)

 110   format(' Recommended MCAVM  =  ',i3,'    Input MCAVM  = ',i4,
     %      '    Err = ',i2,'  * Use Recommended Value ' )

 115   format(' Recommended MCAVM >=  ',i3,'    Input MCAVM  = ',i4,
     %      '    Err = ',i2)
 

       isum1 = 0
       isum2 = 0
       do i = 1 , 11
          IF(ierb(i) .eq. 1) then
             isum1 = isum1 + 1
          ELSEIF(ierb(i) .eq. -1) then
             isum2 = isum2 + 1
          ENDIF
       enddo

       if(isum1 .gt. 0) then
          write(*,*)
          write(*,*)
          write(*,*) '********** ERROR in PARAM.INC File ! ************'
          write(*,*) '*                                               *'
          write(*,*) '*       Check <ERR.LOG> File to find ERROR      *'
          write(*,*) '*                                               *'
          write(*,*) '*************************************************'
          write(*,*)
          stop
       endif

       if(isum2 .gt. 0) then
          write(*,*)
          write(*,*) 
          write(*,*) '********* WARNINGS in PARAM.INC File! ***********'
          write(*,*) '*                                               *'
          write(*,*) '*   Those warnings may not affect solutions,    *'
          write(*,*) '*   however, those need to be changed to save   *'
          write(*,*) '*   array sizes and avoid the possible error.   *'
          write(*,*) '*   Check <ERR.LOG> File to check WARNINGS.     *'
          write(*,*) '*                                               *'
          write(*,*) '*************************************************'
          write(*,*)
       endif


       RETURN
       END




C ======================================================================
      SUBROUTINE CHECK_NWZ
C
C     This routine check the min. number of NWZ for propeller wake
C
C ======================================================================
      use m_WAKGEO
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
!     PARAMETER(NWR=MBPZ)
!s-- YE TIAN 07/07/2013----
!     following common block is moved into the m_WAKGEO module
!     COMMON /WAKGEO/ RTBL(11,NWR),XTBL(11),RW(NWR),VAR(NWR),VTR(NWR),
!    *     UASTAR(NWR),UTSTAR(NWR),UAUSTR(NWR),UTUSTR(NWR)

      integer NWR
!s--- YE TIAN --- 07/09/--- 
! This subroutine is modified in order to get MAXNSW before accessing
! the array xw,yz,zw, which are allocated after knowing MAXNSW.
!e--- YE TIAN --- 07/09/--- 
      real xw1_tmp(mr+1)


      NWR = MBPZ

      if (.NOT.allocated(RTBL)) then
        allocate(RTBL(11,NWR),XTBL(11),RW(NWR),VAR(NWR),VTR(NWR),
     *     UASTAR(NWR),UTSTAR(NWR),UAUSTR(NWR),UTUSTR(NWR))
      end if


C-----------------------------------------------------------------------
C     Parameters of the transition wake
C-----------------------------------------------------------------------

      XFINAL=XULT

      DTPROP=DELTAT/RAD
      DTF=DTPROP*ADVCO/180.0

      MN=MR+1

      DO 10 M=1,MN
         RW(M)=SQRT(YB(1,M)**2+ZB(1,M)**2)
10    CONTINUE

C-----------------------------------------------------------------------
C     Interpolate velocities in the wake
C-----------------------------------------------------------------------

      CALL EVALDK(NX,MN,XR,RW,VAR,VACUB)
      CALL EVALDK(NX,MN,XR,RW,UASTAR,UACUB)
      CALL EVALDK(NX,MN,XR,RW,UAUSTR,UAUCUB)

      DO 60 M=1,MN

!        XW(1,M)=XB(1,M)
         xw1_tmp(m)=XB(1,M)
         
         UAINC=UAUSTR(M)-UASTAR(M)

C-----------------------------------------------------------------------
c     Create the helical wake
C-----------------------------------------------------------------------
         N = 1
         
 1000    N = N  + 1

C --- S.H.CHANG 04/02/2010         
         IF(N .GT. 300) THEN
            WRITE(*,*) 
C            WRITE(*,*) ' ERROR : NWZ requires more than 300 '
            WRITE(*,*) ' WARNING : NWZ requires more than 300 '
            WRITE(*,*)
!Allen Du 01/06/2018 add the option of coupling with cavopt3d 
            open(70401,file="cavopt_parameter.dat")
            write(70401,*)"1"
            close(70401)
!Allen Du 06/16/2019 add the option of coupling with the extension
!scheme
            open(70401,file="ex_check.log")
            close(70401)
            stop
C            STOP
         ENDIF
C --- S.H.CHANG 04/02/2010
         
!        Q=(XW(N-1,M)-XW(1,M))/XFINAL
         Q=(xw1_tmp(m)-XB(1,M))/XFINAL
         
         IF(ICON.EQ.7) THEN
!           XW(N,M)=XW(N-1,M)+(PITCH(M)/ADVCO)*DTF
            xw1_tmp(m)=xw1_tmp(m)+(PITCH(M)/ADVCO)*DTF
         ELSE
!           XW(N,M)=XW(N-1,M)+(VAR(M)+UASTAR(M)+GROW(Q)*UAINC)*DTF
            xw1_tmp(m)=xw1_tmp(m)+(VAR(M)+UASTAR(M)+GROW(Q)*UAINC)*DTF
         END IF
         
         NSW(M)=N
         
!        IF(XW(N,M).GE.XULT) THEN
         IF(xw1_tmp(m).GE.XULT) THEN
            GO TO 50
         ELSE
            GO TO 1000
         ENDIF
         
 50      CONTINUE
         
 60   CONTINUE

      MAXNSW = 0
      DO M = 1, MN
         MAXNSW = MAX(MAXNSW, NSW(M))
      ENDDO

      write(99,*)
      write(99,*) ' --------------------------------------'
      write(99,*) '  (3) Checking ..Minimum value of NWZ  '
      write(99,*) ' --------------------------------------'
      write(99,*)
      
      IERB = 0
      IF(MAXNSW .GT. NWZ) IERB = 1
      WRITE(99,100) MAXNSW, NWZ, IERB   
 100   format(' Required Min. NWZ  =  ',i4,'    Input NWZ  = ',i4,
     %      '    Err = ',i2,'  * Use Required Min. Value ' )       

       if(IERB .EQ. 1) then
          write(*,*)
          write(*,*)
          write(*,*) '*** ERROR : PARAMETER NWZ in PARAM.INC File ! ***'
          write(*,*) '*                                               *'
          write(*,*) '*  Check <ERR.LOG> to find exact value of NWZ   *'
          write(*,*) '*                                               *'
          write(*,*) '*************************************************'
          write(*,*)
          stop
       endif
  
       RETURN
       END


