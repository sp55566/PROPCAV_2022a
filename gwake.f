      SUBROUTINE GWAKE
************************************************************************
*     GWAKE: Geometry of the WAKE                                      *
*      --- Generate geometries of the transition wake and the ultimate *
*          wake                                                        *
*                                                                      *
*     Date         Revision/Comment                                    *
*     ----------- ---------------------------------------------------- *
*     12/04/98    I merged gwakesc.f into this routine so all we       *
*                 generate the wake panels for the fully wetted and    *
*                 cavitating geometry in the same manner.              *
*                                                                      *
************************************************************************

      use m_WAKGEO
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      integer NWR
      NWR = MBPZ
!     PARAMETER(NWR=MBPZ)
!     COMMON /WAKGEO/ RTBL(11,NWR),XTBL(11),RW(NWR),VAR(NWR),VTR(NWR),
!    *     UASTAR(NWR),UTSTAR(NWR),UAUSTR(NWR),UTUSTR(NWR)

!s-- YE TIAN ---- 08/13/2013---
!     write(*,*) 'NWR=',NWR
      if (.NOT.allocated(RTBL)) then
        allocate(RTBL(11,NWR),XTBL(11),RW(NWR),VAR(NWR),VTR(NWR),
     *     UASTAR(NWR),UTSTAR(NWR),UAUSTR(NWR),UTUSTR(NWR))
      end if
!s-- YE TIAN ---- 08/13/2013---
C-----------------------------------------------------------------------
C     Parameters of the transition wake
C-----------------------------------------------------------------------
      XFINAL=XULT

      DTPROP=DELTAT/RAD
      DTF=DTPROP*ADVCO/180.0
      DTRADF=DELTAT

      NN=NCP
      MN=MRP

      NWK=NWPZ

C---------RW is radial distance (excluding x-component) to the panel----
C---------nodes at the trailing edge of the blades CM-------------------

      DO 10 M=1,MN
         RW(M)=SQRT(YB(1,M)**2+ZB(1,M)**2)
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
      BETAtmp=ATAN(ADVCO*VAR(MRP)/PI)
      BETAB=PHIP
      XDIS=TWO*CHORD(MRP)

      DELTIP=ZERO
      XTW=XB(1,MRP)
      YTW=YB(1,MRP)
      ZTW=ZB(1,MRP)

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

C.....If the section have finite T.E. thickness, start the wake at......
C.....the midpoint of the T.E. (JY120498)...............................

            IF(M.EQ.MN) THEN
               XW(1,M)=XTW
               YW(1,M)=YTW
               ZW(1,M)=ZTW
            ELSE
               XW(1,M)=XB(1,M)
               YW(1,M)=YB(1,M)
               ZW(1,M)=ZB(1,M)
            END IF

C..........Same for ICON=8 as well (JY110100)
         IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8)THEN
            TWOLD=ATAN2(ZW(1,M),YW(1,M))
         ENDIF

         UAINC=UAUSTR(M)-UASTAR(M)
         UTINC=UTUSTR(M)-UTSTAR(M)

C-----------------------------------------------------------------------
c     Create the helical wake
C-----------------------------------------------------------------------

         DO 40 N=2,NWK

            DT1=DTF
            DTRAD=DTRADF

            Q=(XW(N-1,M)-XW(1,M))/XFINAL

C..........Use the geometric pitch of the blade for special fan-like
C..........SP-propeller (JY090500)
            IF(ICON.EQ.7) THEN
               XW(N,M)=XW(N-1,M)+(PITCH(M)/ADVCO)*DT1
            ELSE
               IF(IDUCT .EQ. 1) THEN
                  XW(N,M)=XW(N-1,M)+(PITCH(M)/ADVCO)*DT1
               ELSE
                  XW(N,M)=XW(N-1,M)+(VAR(M)+UASTAR(M)+GROW(Q)*UAINC)*DT1
               ENDIF
            END IF

            XWTE=XW(N,M)-XW(1,M)
            RWK=RTBL(11,M)

            IF(XWTE.LT.XULT) THEN
               DO 20 J=1,10
                  IF(XWTE.LT.XTBL(J).OR.XWTE.GE.XTBL(J+1)) THEN
                     CONTINUE
                  ELSE
                     RWK=RTBL(J,M)+(RTBL(J+1,M)-RTBL(J,M))
     *                  *(XWTE-XTBL(J))/(XTBL(J+1)-XTBL(J))
                     GO TO 30
                  END IF
20             CONTINUE
            END IF
30          CONTINUE

            IF((ICON.EQ.5).OR.(ICON.EQ.6).OR.(ICON.EQ.8))THEN
               ZW(N,M)=ZW(N-1,M)
               YW(N,M)=YW(N-1,M)
            ELSE
               TW=ATAN2(ZW(N-1,M),YW(N-1,M))
               IF(TW.LT.TWOLD) TW=TW+TWOPI
               IF(IDUCT .EQ. 1) THEN
                  TW = TW + DTRAD
               ELSE
                  TW=TW+DTRAD+(VTR(M)+UTSTAR(M)+GROW(Q)*UTINC)*DT1/RWK
               ENDIF
               YW(N,M)=RWK*COS(TW)
               ZW(N,M)=RWK*SIN(TW)
            ENDIF

            TWOLD=TW
            NSW(M)=N
            IF(XW(N,M).GE.XUWDK) THEN
               GO TO 50
            END IF
40       CONTINUE

         WRITE(*,'(''!!!!! WARNING !!!!! TRANSITION WAKE SHORT'')')
         WRITE(*,'(''      M = '',I5,''  XW = '',F8.3)') M,XW(NSW(M),M)

C........Next four lines added by JY072200
         WRITE(*,'('' NSW(M) = '',I5,'' NWK = '',I5)') NSW(M),NWK
         WRITE(*,'('' Decrease XUWDK in .geo, OR '')')
         WRITE(*,'('' Increase XWK in gwake.f! '')')
         STOP

C........End of the transition wake at this radius
50       CONTINUE

C........Next IF statement added by JY072200
         IF(NSW(M).GT.NWZ) THEN
            WRITE(*,'('' M = '',I5,'' NSW(M) = '',I5,'' NWZ = '',I5)') 
     *           M,NSW(M),NWZ
            WRITE(*,'('' Decrease XUWDK in .geo, OR '')')
            WRITE(*,'('' Increase XWZ --> NSW(M) in PARAM.INC! '')')
            STOP
         END IF

         FK=( XW(NSW(M),M)-XUWDK)/( XW(NSW(M),M)-XW(NSW(M)-1,M) )
         XW(NSW(M),M)=XUWDK

C.......The next definition for TW is only valid when it's not
C.......a foil. (JY110500)
         IF(ICON.NE.5.AND.ICON.NE.6.AND.ICON.NE.8)THEN
            IF(IDUCT .EQ. 1) THEN
               TW = TW - FK*DTRAD
            ELSE
               TW=TW-FK*(DTRAD+(VTR(M)+UTSTAR(M)+GROW(Q)*UTINC)*DT1/RWK)
            ENDIF
            YW(NSW(M),M)=RWK*COS(TW)
            ZW(NSW(M),M)=RWK*SIN(TW)
         ELSE
            ZW(NSW(M),M)=ZW(NSW(M)-1,M)
            YW(NSW(M),M)=YW(NSW(M)-1,M)            
         END IF

C........stop the program and write out the coordinates of the last
C........wake panel if any of the value is zero.
         IF(XW(NSW(M),M).EQ.0.0.AND.YW(NSW(M),M).EQ.0.0.AND.
     *        ZW(NSW(M),M).EQ.0.0) THEN
            WRITE(*,*) 'PROBLEM WITH LAST WAKE PANEL!'
            WRITE(*,*) 'XW,YW,ZW:',XW(NSW(M),M),YW(NSW(M),M),
     *           ZW(NSW(M),M)
            STOP
         END IF
60    CONTINUE

      DO 11 M=1,MR
         IF(NSW(M).EQ.NSW(M+1))THEN
            NWPAN(M)=NSW(M)-1
         ELSE
            NWPAN(M)=MIN(NSW(M),NSW(M+1))
            NWPAN(M)=NWPAN(M)-1
         ENDIF
 11   CONTINUE


         
C-----------------------------------------------------------------------
C     Ultimate wake geometry
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     Here I undo the comment on TANBUW=(XVA(NX)+....) and commented out
C     TANBUW=.3707.                                             JY022798
C-----------------------------------------------------------------------

C....Use the geometric pitch of the blade for special fan-like
C....SP-propeller (JY090500)
      IF(ICON.EQ.7) THEN
         TANBUW=PITCH(MRP)/PI
      ELSE
         TANBUW=(XVA(NX)+XUAU(NX))/((PI*RULT/ADVCO)+XUTU(NX)+XVT(NX))
C      TANBUW=0.3707
      END IF

      CALL ULTWAK

      RETURN
C))))))))))))))))))))) End of subroutine GWAKE (((((((((((((((((((((((((
      END


