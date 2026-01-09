      SUBROUTINE GWAKE_DUCT
************************************************************************
*     GWAKE: Geometry of the WAKE for DUCT PROPELLER                   *
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
!     PARAMETER(NWR=MBPZ)
!     COMMON /WAKGEO/ RTBL(11,NWR),XTBL(11),RW(NWR),VAR(NWR),VTR(NWR),
!    *     UASTAR(NWR),UTSTAR(NWR),UAUSTR(NWR),UTUSTR(NWR)
      DIMENSION XID1(101),ETAD1(101),DCUBICX1(400)
      integer NWR
      NWR=MBPZ

C-----------------------------------------------------------------------
C     Parameters of the transition wake
C-----------------------------------------------------------------------
      XFINAL=XULT

      DTPROP=DELTAT/RAD
C      DTF=0.5*DTPROP*ADVCO/180.0
C      DTRADF=0.5*DELTAT

      DTF=DTPROP*ADVCO/180.0
      DTRADF=DELTAT

      NN=NCP
      MN=MRP

C --- S.H.CHANG 03/02/2010 for checking with Panel Model
      NWK=NWPZ
C      NWK=120
C --- S.H.CHANG 03/02/2010 for checking with Panel Model
C---------RW is radial distance (excluding x-component) to the panel----
C---------nodes at the trailing edge of the blades CM-------------------

      DO 10 M=1,MN
         RW(M)=SQRT(YB(1,M)**2+ZB(1,M)**2)
10    CONTINUE

C --------------------------------------------------------
C     Data for Duct Inner Surface
C --------------------------------------------------------
      NN1 = NDDAT/2+1

      DO N = 1 , NN1
         N1 = NN1-N+1
         XID1(N) = XID(N1)
         ETAD1(N) = ETAD(N1)
      ENDDO

      CALL UGLYDK(NN1,1,1,XID1,ETAD1,ZERO,ZERO,DCUBICX1)

      CALL EVALDKs(NN1,1,XID1,XB(1,MRP),RRR1,DCUBICX1)

      RPROP = SQRT(YB(1,MRP)**2+ZB(1,MRP)**2)

      IF(RRR1 .LE. RPROP) THEN
         RW(MRP) = RRR1
      ENDIF

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
C     Compute radii of the transition wake
C-----------------------------------------------------------------------
      CALL RTABLE_DUCT(RULT,RHULT,RTBL,XTBL,XULT,RW,MN,XID1,ETAD1,NN1)

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

         XW(1,M)=XB(1,M)
         YW(1,M)=YB(1,M)
         ZW(1,M)=ZB(1,M)

         TWOLD=ATAN2(ZW(1,M),YW(1,M))

C --- S.H.CHANG 03/31/2010 CHECKING PANEL METHOD
         UAINC= UAUSTR(M)-UASTAR(M)
         UTINC= UTUSTR(M)-UTSTAR(M)
C         UAINC= 0.0   !FOR MPUF-3A/BEM
C         UTINC= 0.0   !FOR MPUF-3A/BEM
C --- S.H.CHANG 03/31/2010 CHECKING PANEL METHOD

C-----------------------------------------------------------------------
c     Create the helical wake
C-----------------------------------------------------------------------

         DO 40 N=2,NWK

            DT1=DTF
            DTRAD=DTRADF

            Q=(XW(N-1,M)-XW(1,M))/XFINAL ! XFINAL = XULT

            IF(IDALIGN .EQ. 0) THEN
               XW(N,M) = XW(N-1,M) + (PITCH(M)/ADVCO)*DT1
            ELSE
               XW(N,M)=XW(N-1,M)+(VAR(M)+UASTAR(M)+GROW(Q)*UAINC)*DT1
            ENDIF

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

            TW=ATAN2(ZW(N-1,M),YW(N-1,M))
            IF(TW.LT.TWOLD) TW=TW+TWOPI
            IF(IDALIGN .EQ. 0) THEN
               TW = TW + DTRAD
            ELSE
               TW=TW+DTRAD+(VTR(M)+UTSTAR(M)+GROW(Q)*UTINC)*DT1/RWK
            ENDIF

            YW(N,M)=RWK*COS(TW)
            ZW(N,M)=RWK*SIN(TW)

            TWOLD=TW
            NSW(M)=N

C --- S.H.CHANG 03/02/2010 for checking with Panel Model
            IF(XW(N,M).GE.XUWDK) THEN
               GO TO 50
            END IF
C --- S.H.CHANG 03/02/2010 for checking with Panel Model

40       CONTINUE



C --- S.H.CHANG 03/02/2010 for checking with Panel Model
         WRITE(*,'(''!!!!! WARNING !!!!! TRANSITION WAKE SHORT'')')
         WRITE(*,'(''      M = '',I5,''  XW = '',F8.3)') M,XW(NSW(M),M)
C --- S.H.CHANG 03/02/2010 for checking with Panel Model

50       CONTINUE

         IF(NSW(M).GT.NWZ) THEN
            WRITE(*,'('' M = '',I5,'' NSW(M) = '',I5,'' NWZ = '',I5)')
     *           M,NSW(M),NWZ
            WRITE(*,'('' Decrease XUWDK in .geo, OR '')')
            WRITE(*,'('' Increase XWZ --> NSW(M) in PARAM.INC! '')')
            STOP
         END IF

         FK=(XW(NSW(M),M)-XUWDK)/(XW(NSW(M),M)-XW(NSW(M)-1,M))
         XW(NSW(M),M)=XUWDK

C ----- SNKIM cleaned out useless lines and comments, and modified the 'TW' -----
         IF(IDALIGN .EQ. 0) THEN
            TW = TW - FK*DTRAD
         ELSE
            TW=TW-FK*(DTRAD+(VTR(M)+UTSTAR(M)+GROW(Q)*UTINC)*DT1/RWK)
         ENDIF

C         TW=TW-FK*(DTRAD+(VTR(M)+UTSTAR(M)+GROW(Q)*UTINC)*DT1/RWK)

         YW(NSW(M),M)=RWK*COS(TW)
         ZW(NSW(M),M)=RWK*SIN(TW)
C ----- SNKIM cleaned out useless lines and comments, and modified the 'TW' -----

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

      TANBUW=(XVA(NX)+XUAU(NX))/((PI*RULT/ADVCO)+XUTU(NX)+XVT(NX))

      CALL ULTWAK

      RETURN
      END
