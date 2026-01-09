
      SUBROUTINE GWAKE1
************************************************************************
*     GWAKE1: Geometry of the first wake panel (for positioning of the *
*             separated region.                                        *
*                                                                      *
*     Date         Revision/Comment                                    *
*     ----------- ---------------------------------------------------- *
*     JY080701     Copy from subroutine GWAKE.                         *
*                                                                      *
************************************************************************

      use m_BLADEM
      use m_WAKGEO
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

!     PARAMETER(NWR=MBPZ)
!     COMMON /WAKGEO/ RTBL(11,NWR),XTBL(11),RW(NWR),VAR(NWR),VTR(NWR),
!    *     UASTAR(NWR),UTSTAR(NWR),UAUSTR(NWR),UTUSTR(NWR)
      DIMENSION XBM1(MBPZ),YBM1(MBPZ),ZBM1(MBPZ)
!     COMMON/BLADEM/XBM(NBHZP,MBPZ),YBM(NBHZP,MBPZ),ZBM(NBHZP,MBPZ)
      integer nwr
      nwr = MBPZ
!s-- YE TIAN ---- 08/13/2013---
!     write(*,*) 'gwake1.f:26,nwr',nwr
!     stop
      if (.NOT.allocated(RTBL)) then
        allocate(RTBL(11,NWR),XTBL(11),RW(NWR),VAR(NWR),VTR(NWR),
     *     UASTAR(NWR),UTSTAR(NWR),UAUSTR(NWR),UTUSTR(NWR))
      end if
!s-- YE TIAN ---- 08/13/2013---

C-----------------------------------------------------------------------
C     Length of the separated region in terms of degrees.
C-----------------------------------------------------------------------
      DELTAT2=2.*RAD

C-----------------------------------------------------------------------
C     Parameters of the transition wake
C-----------------------------------------------------------------------
      XFINAL=XULT

      DTPROP=DELTAT2/RAD
      DTF=DTPROP*ADVCO/180.0
      DTRADF=DELTAT2

      NN=NCP
      MN=MRP

      NWK=201

C-----------------------------------------------------------------------
C     Define midpoint between upper and lower sides of the blade T.E.
C-----------------------------------------------------------------------
      DO M=1,MN
         IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
            XBM1(M)=HALF*(XB(1,M)+XB(NCP,M))
            YBM1(M)=HALF*(YB(1,M)+YB(NCP,M))
            ZBM1(M)=HALF*(ZB(1,M)+ZB(NCP,M))
         ELSE
            XBM1(M)=XBM(NHP,M)
            YBM1(M)=YBM(NHP,M)
            ZBM1(M)=ZBM(NHP,M)            
         END IF
      END DO

C---------RW is radial distance (excluding x-component) to the panel----
C---------nodes at the trailing edge of the blades CM-------------------
      DO 10 M=1,MN
         RW(M)=SQRT(YBM1(M)**2+ZBM1(M)**2)
10    CONTINUE

C---------RULT is probably radial distance to outermost point on blade--
C---------RHULT is probably rad dist to innermost pt on te of blade CM--

      IF(DCD.EQ.0.0) THEN
         RULT=RW(MN)
         RHULT=RW(1)
      END IF

C-----------------------------------------------------------------------
C     Interpolate velocities in the wake
C-----------------------------------------------------------------------
      CALL EVALDK(NX,MN,XR,RW,VAR,VACUB)
      CALL EVALDK(NX,MN,XR,RW,VTR,VTCUB)
      CALL EVALDK(NX,MN,XR,RW,UASTAR,UACUB)
      CALL EVALDK(NX,MN,XR,RW,UAUSTR,UAUCUB)
      CALL EVALDK(NX,MN,XR,RW,UTSTAR,UTCUB)
      CALL EVALDK(NX,MN,XR,RW,UTUSTR,UTUCUB)

C-----------------------------------------------------------------------
C     Separated tip wake geometry 
C-----------------------------------------------------------------------
      TANP=PITCH(MRP)/PI
      COSP=1.0/SQRT(1.0+TANP*TANP)
      SINP=TANP*COSP
      PHIP=ATAN(TANP)
      BETA_tmp=ATAN(ADVCO*VAR(MRP)/PI)
      BETAB=PHIP
      XDIS=TWO*CHORD(MRP)

      DELTIP=ZERO

C-----------------------------------------------------------------------
C     Compute radii of the transition wake
C-----------------------------------------------------------------------
      CALL RTABLE(RULT,RHULT,DCD,RTBL,XTBL,XULT,RW,MN,XHBTE,XHBT,XHBFD)

C-----------------------------------------------------------------------
C     Generate the transition wake geometry
C-----------------------------------------------------------------------
      DO 60 M=1,MN

C-----------------------------------------------------------------------
C     Initial coordinates of the transition wake from blade grid
C     These are the same coords as the T.E. points on the blade
C-----------------------------------------------------------------------

         XW(1,M)=XBM1(M)
         YW(1,M)=YBM1(M)
         ZW(1,M)=ZBM1(M)

         IF(ICON.NE.5.AND.ICON.NE.8)THEN
            TWOLD=ATAN2(ZW(1,M),YW(1,M))
         ENDIF

         UAINC=UAUSTR(M)-UASTAR(M)
         UTINC=UTUSTR(M)-UTSTAR(M)

C-----------------------------------------------------------------------
c     Create the helical wake
C-----------------------------------------------------------------------
         DO N=2,5

            DT1=DTF
            DTRAD=DTRADF
            
            Q=(XW(N-1,M)-XW(1,M))/XFINAL
            
C.......Use the geometric pitch of the blade for special fan-like
C.......SP-propeller (JY090500)
            IF(ICON.EQ.7) THEN
               XW(N,M)=XW(N-1,M)+(PITCH(M)/ADVCO)*DT1
            ELSE
               XW(N,M)=XW(N-1,M)+(VAR(M)+UASTAR(M)+GROW(Q)*UAINC)*DT1
            END IF
            
            XWTE=XW(N,M)-XW(1,M)
            RWK=RTBL(11,M)
            IF(XWTE.LT.XULT) THEN
               DO 20 J=1,10
                  IF(XWTE.LT.XTBL(J).OR.XWTE.GE.XTBL(J+1)) THEN
                     CONTINUE
                  ELSE
                     RWK=RTBL(J,M)+(RTBL(J+1,M)-RTBL(J,M))
     *                 *(XWTE-XTBL(J))/(XTBL(J+1)-XTBL(J))
                     GO TO 30
                  END IF
 20            CONTINUE
            END IF
 30         CONTINUE

            IF((ICON.EQ.5).OR.(ICON.EQ.6).OR.(ICON.EQ.8))THEN
               ZW(N,M)=ZW(N-1,M)
               YW(N,M)=YW(N-1,M)
            ELSE
               TW=ATAN2(ZW(N-1,M),YW(N-1,M))
               IF(TW.LT.TWOLD) TW=TW+TWOPI
               TW=TW+DTRAD+(VTR(M)+UTSTAR(M)+GROW(Q)*UTINC)*DT1/RWK
               YW(N,M)=RWK*COS(TW)
               ZW(N,M)=RWK*SIN(TW)
            ENDIF

         END DO
         
 60   CONTINUE

      RETURN
C))))))))))))))))))))) End of subroutine GWAKE1 ((((((((((((((((((((((((
      END


