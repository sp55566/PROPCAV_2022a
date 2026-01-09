C----------------------------------------------------------------------
      SUBROUTINE GWAKETIP
C----------------------------------------------------------------------
C     
C     ORIGINAL VERSION : APRIL 26 1999
C     LATEST REVISION  : OCT. 11 1999 
C          
C     PURPOSE          : GENERATE HELICAL BLADE WAKE GEOMETRY WITH 
C     GEMETRIC PITCH ANGLE   
C     
C     USAGE            : CALL GWAKETIP  
C     
C     ARGUMENTS :
C     INPUT 
C     OUTPUT 
C     
C     REQUIRED  ROUTINES  
C     
C     REVISION        : JANUARY 18,2010
C                     : CHANGED LOCAL ARRAYS TO ALLOCATABLE 
C                     : VV (VIMAL VINAYAN)
C----------------------------------------------------------------------

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

!     PARAMETER(NWR=MBPZ)
CVV   DIMENSION QQ(NWZ),RR1(NWZ),RRAT(NWR)

CVV
      ALLOCATABLE :: QQ(:),RR1(:),RRAT(:)
CVV      
      integer nwr
      nwr = MBPZ

C-----------------------------------------------------------------------
C     PARAMETERS OF THE TRANSITION WAKE
C-----------------------------------------------------------------------
      DTPROP=DELTAT/RAD
      DTF=DTPROP*ADVCO/180.0
      DTRADF=DELTAT
      DTI=DTF
      DTRADI=DTRADF

      MN=MRP
      NWPANEL = NWK - 1
      NWPP = NWPANEL + 1

      RULT = RHULT

CVV   ALLOCATE MEMORY FOR LOCAL ARRAYS
      ALLOCATE(QQ(NWZ),RR1(NWZ),RRAT(NWR))

C--   IF NO HUB IS CONSIDERED, DO NOT ALLOW RADIUS CHANGE OF ULTIIMATE WAKE.
      IF(IHUB .EQ. 0) RULT = RHUB
C-----------------------------------------------------------------------
C     SEPARATED TIP WAKE GEOMETRY 
C-----------------------------------------------------------------------

      DO M = 1 , MN
        XW(1,M) = XB(1,M)
        YW(1,M) = YB(1,M)
        ZW(1,M) = ZB(1,M)
        RRAT(M) = SQRT(YW(1,M)**2+ZW(1,M)**2)
      ENDDO

      DDR = RRAT(MRP) - RRAT(1)
      RR11 = RRAT(1)

      DO M = 1 , MN
        RRAT(M) = (RRAT(M) - RR11)/DDR
      ENDDO
      
!s-- Ye TIAN --- 08/14 debug
!     write(*,*) 'gwaketip.f:75', M, MN
!e-- Ye TIAN --- 08/14 debug
      DO N = 1, NWPP
        IF(N .GT. 1) XW(N,M) = XW(N-1,M) + DTI
        QQ(N) = (XHBFD-XW(N,1))/XHBT
        RR1(N) = RNOSET(QQ(N),RHUB,RHULT)
      ENDDO

      DO M = 1 , MN
        TW=DANGLE(ZW(1,M),YW(1,M))
        DO N = 2 , NWPP
          XW(N,M) = XW(N-1,M) + DTI
          RTW = RR1(N) + (1. - RR1(N))*RRAT(M)
          TW = TW + DTRADI
          YW(N,M) = RTW*COS(TW)
          ZW(N,M) = RTW*SIN(TW)
          NSW(M) = N
        ENDDO
      ENDDO
      
      XUWDK = XW(NWK,MRP)

CVV   DEALLOCATE MEMORY OF LOCAL ARRAYS
      DEALLOCATE(QQ,RR1,RRAT)
 
      RETURN
      END
