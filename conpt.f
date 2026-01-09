      SUBROUTINE CONPT
************************************************************************
*     CONPT: CONtrol PoinTs set-up                                     *
*     Construct blade panels and set up control points                 *
*                                                                      *
*     Date     Comment or Revision                                     *
*     -------- -------------------                                     *
*     CM010598 Changes made by H.S. Lee below for the hydrofoil        *
*              geometry generation (ICON = 5)                          *
*     JY090799 Moved the part that reads the data from *.wak file to   *
*              subroutine readwak.f.                                   *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      CHARACTER*29 :: CPTCHECK, FILENAME
C-----------------------------------------------------------------------
C
C     Construct blade panels,              V
C      clockwise viewing      (N,M+1) 2 *--^----* 3 (N+1,M+1)
C      from fluid domain                |  |   |
C                                     U <--*   |
C                               (N,M) 1 *------* 4 (N+1,M)
C
C-----------------------------------------------------------------------

      DO 20 N=1,NC
        DO 10 M=1,MR
          L=INDEXB(N,M)
          XG(L,1,1)=XB(N,M)
          XG(L,1,2)=YB(N,M)
          XG(L,1,3)=ZB(N,M)
          XG(L,2,1)=XB(N,M+1)
          XG(L,2,2)=YB(N,M+1)
          XG(L,2,3)=ZB(N,M+1)
          XG(L,3,1)=XB(N+1,M+1)
          XG(L,3,2)=YB(N+1,M+1)
          XG(L,3,3)=ZB(N+1,M+1)
          XG(L,4,1)=XB(N+1,M)
          XG(L,4,2)=YB(N+1,M)
          XG(L,4,3)=ZB(N+1,M)

          XG1(L,1,1)=XB(N,M)
          XG1(L,1,2)=YB(N,M)
          XG1(L,1,3)=ZB(N,M)
          XG1(L,2,1)=XB(N,M+1)
          XG1(L,2,2)=YB(N,M+1)
          XG1(L,2,3)=ZB(N,M+1)
          XG1(L,3,1)=XB(N+1,M+1)
          XG1(L,3,2)=YB(N+1,M+1)
          XG1(L,3,3)=ZB(N+1,M+1)
          XG1(L,4,1)=XB(N+1,M)
          XG1(L,4,2)=YB(N+1,M)
          XG1(L,4,3)=ZB(N+1,M)
   10   CONTINUE
   20 CONTINUE

      IF(IHUB.NE.0) THEN
C-----------------------------------------------------------------------
C
C        Construct hub panels,                 V
C         clockwise viewing         (N,M) 2 *--^----* 3 (N+1,M)
C         from fluid domain                 |  |    |
C                                        U <---*    |
C                                 (N,M+1) 1 *-------* 4 (N+1,M+1)
C
C-----------------------------------------------------------------------
        DO 40 N=1,NHBX
          DO 30 M=1,MHBT
            L=INDEXH(N,M)
            XG(L,1,1)=XH(N,M+1)
            XG(L,1,2)=YH(N,M+1)
            XG(L,1,3)=ZH(N,M+1)
            XG(L,2,1)=XH(N,M)
            XG(L,2,2)=YH(N,M)
            XG(L,2,3)=ZH(N,M)
            XG(L,3,1)=XH(N+1,M)
            XG(L,3,2)=YH(N+1,M)
            XG(L,3,3)=ZH(N+1,M)
            XG(L,4,1)=XH(N+1,M+1)
            XG(L,4,2)=YH(N+1,M+1)
            XG(L,4,3)=ZH(N+1,M+1)

            XG1(L,1,1)=XH(N,M+1)
            XG1(L,1,2)=YH(N,M+1)
            XG1(L,1,3)=ZH(N,M+1)
            XG1(L,2,1)=XH(N,M)
            XG1(L,2,2)=YH(N,M)
            XG1(L,2,3)=ZH(N,M)
            XG1(L,3,1)=XH(N+1,M)
            XG1(L,3,2)=YH(N+1,M)
            XG1(L,3,3)=ZH(N+1,M)
            XG1(L,4,1)=XH(N+1,M+1)
            XG1(L,4,2)=YH(N+1,M+1)
            XG1(L,4,3)=ZH(N+1,M+1)
   30     CONTINUE
   40   CONTINUE
      END IF

      IF(IDUCT .NE. 0) THEN
         DO M=1,MDUCT
            DO N=1,NDUCT
               L=INDEXD(N,M)

               XG(L,1,1)=XD(N,M+1)
               XG(L,1,2)=YD(N,M+1)
               XG(L,1,3)=ZD(N,M+1)
               XG(L,2,1)=XD(N,M)
               XG(L,2,2)=YD(N,M)
               XG(L,2,3)=ZD(N,M)
               XG(L,3,1)=XD(N+1,M)
               XG(L,3,2)=YD(N+1,M)
               XG(L,3,3)=ZD(N+1,M)
               XG(L,4,1)=XD(N+1,M+1)
               XG(L,4,2)=YD(N+1,M+1)
               XG(L,4,3)=ZD(N+1,M+1)

C/Seungnam Kim, Shifting the radius of control points on the duct
c               LL=N+NDUCT*(M-1)
c               WRITE(*,*) 'LL=',LL,N,M
               XG1(L,1,1)=XD(N,M+1)
               XG1(L,1,2)=YD(N,M+1)
               XG1(L,1,3)=ZD(N,M+1)
               XG1(L,2,1)=XD(N,M)
               XG1(L,2,2)=YD(N,M)
               XG1(L,2,3)=ZD(N,M)
               XG1(L,3,1)=XD(N+1,M)
               XG1(L,3,2)=YD(N+1,M)
               XG1(L,3,3)=ZD(N+1,M)
               XG1(L,4,1)=XD(N+1,M+1)
               XG1(L,4,2)=YD(N+1,M+1)
               XG1(L,4,3)=ZD(N+1,M+1)

               IF(N .LE. NDAFT) THEN
                 THETA4=DANGLE(XG1(L,4,3),XG1(L,4,2))
                 THETA3=DANGLE(XG1(L,3,3),XG1(L,3,2))
                 RAD1=SQRT(XG1(L,1,2)**2 + XG1(L,1,3)**2)
                 RAD2=SQRT(XG1(L,2,2)**2 + XG1(L,2,3)**2)

                 XG1(L,1,2)=RAD1*COS(THETA4)
                 XG1(L,1,3)=RAD1*SIN(THETA4)
                 XG1(L,2,2)=RAD2*COS(THETA3)
                 XG1(L,2,3)=RAD2*SIN(THETA3)
               ELSEIF(N .GE. NDUCT-NDAFT+1) THEN
                 THETA1=DANGLE(XG1(L,1,3),XG1(L,1,2))
                 THETA2=DANGLE(XG1(L,2,3),XG1(L,2,2))
                 RAD4=SQRT(XG1(L,4,2)**2 + XG1(L,4,3)**2)
                 RAD3=SQRT(XG1(L,3,2)**2 + XG1(L,3,3)**2)

                 XG1(L,4,2)=RAD4*COS(THETA1)
                 XG1(L,4,3)=RAD4*SIN(THETA1)
                 XG1(L,3,2)=RAD3*COS(THETA2)
                 XG1(L,3,3)=RAD3*SIN(THETA2)
               ENDIF
C/Seungnam Kim, Shifting the radius of control points on the duct

            ENDDO
         ENDDO
      END IF

      IF(ITUN .NE. 0) THEN
         DO  N=1,NAXT
            DO  M=1,MTUNEL
               L=INDEXTN(N,M)

               XG(L,1,1)=XTUN(N,M)
               XG(L,1,2)=YTUN(N,M)
               XG(L,1,3)=ZTUN(N,M)
               XG(L,2,1)=XTUN(N,M+1)
               XG(L,2,2)=YTUN(N,M+1)
               XG(L,2,3)=ZTUN(N,M+1)
               XG(L,3,1)=XTUN(N+1,M+1)
               XG(L,3,2)=YTUN(N+1,M+1)
               XG(L,3,3)=ZTUN(N+1,M+1)
               XG(L,4,1)=XTUN(N+1,M)
               XG(L,4,2)=YTUN(N+1,M)
               XG(L,4,3)=ZTUN(N+1,M)

               XG1(L,1,1)=XTUN(N,M)
               XG1(L,1,2)=YTUN(N,M)
               XG1(L,1,3)=ZTUN(N,M)
               XG1(L,2,1)=XTUN(N,M+1)
               XG1(L,2,2)=YTUN(N,M+1)
               XG1(L,2,3)=ZTUN(N,M+1)
               XG1(L,3,1)=XTUN(N+1,M+1)
               XG1(L,3,2)=YTUN(N+1,M+1)
               XG1(L,3,3)=ZTUN(N+1,M+1)
               XG1(L,4,1)=XTUN(N+1,M)
               XG1(L,4,2)=YTUN(N+1,M)
               XG1(L,4,3)=ZTUN(N+1,M)

            ENDDO
         ENDDO
      END IF

C/s S.N.KIM | For new version of PROPCAV released in 2018, tip vortex model 
C             is omitted for periodic unsteady runs.
C******************************************************************************
cC-- Begin Tip HSLEE(10/12/99)
c
c      if(ian .eq. 2) then
c
c          DO N=1,NTHX
c            DO M=1,MCVT
c              L=INDEXT(N,M)
c              XG(L,1,1)=XCH(N+1,M+1)
c              XG(L,1,2)=YCH(N+1,M+1)
c              XG(L,1,3)=ZCH(N+1,M+1)
c              XG(L,2,1)=XCH(N+1,M)
c              XG(L,2,2)=YCH(N+1,M)
c              XG(L,2,3)=ZCH(N+1,M)
c              XG(L,3,1)=XCH(N,M)
c              XG(L,3,2)=YCH(N,M)
c              XG(L,3,3)=ZCH(N,M)
c              XG(L,4,1)=XCH(N,M+1)
c              XG(L,4,2)=YCH(N,M+1)
c              XG(L,4,3)=ZCH(N,M+1)
c
c              XG1(L,1,1)=XCH(N+1,M+1)
c              XG1(L,1,2)=YCH(N+1,M+1)
c              XG1(L,1,3)=ZCH(N+1,M+1)
c              XG1(L,2,1)=XCH(N+1,M)
c              XG1(L,2,2)=YCH(N+1,M)
c              XG1(L,2,3)=ZCH(N+1,M)
c              XG1(L,3,1)=XCH(N,M)
c              XG1(L,3,2)=YCH(N,M)
c              XG1(L,3,3)=ZCH(N,M)
c              XG1(L,4,1)=XCH(N,M+1)
c              XG1(L,4,2)=YCH(N,M+1)
c              XG1(L,4,3)=ZCH(N,M+1)
c            ENDDO
c          ENDDO
c
c          DO N=1,NCVX
c            DO M=1,MCVT
c              L=INDEXC(N,M)
c              XG(L,1,1)=XVC(N+1,M+1)
c              XG(L,1,2)=YVC(N+1,M+1)
c              XG(L,1,3)=ZVC(N+1,M+1)
c              XG(L,2,1)=XVC(N+1,M)
c              XG(L,2,2)=YVC(N+1,M)
c              XG(L,2,3)=ZVC(N+1,M)
c              XG(L,3,1)=XVC(N,M)
c              XG(L,3,2)=YVC(N,M)
c              XG(L,3,3)=ZVC(N,M)
c              XG(L,4,1)=XVC(N,M+1)
c              XG(L,4,2)=YVC(N,M+1)
c              XG(L,4,3)=ZVC(N,M+1)
c
c              XG1(L,1,1)=XVC(N+1,M+1)
c              XG1(L,1,2)=YVC(N+1,M+1)
c              XG1(L,1,3)=ZVC(N+1,M+1)
c              XG1(L,2,1)=XVC(N+1,M)
c              XG1(L,2,2)=YVC(N+1,M)
c              XG1(L,2,3)=ZVC(N+1,M)
c              XG1(L,3,1)=XVC(N,M)
c              XG1(L,3,2)=YVC(N,M)
c              XG1(L,3,3)=ZVC(N,M)
c              XG1(L,4,1)=XVC(N,M+1)
c              XG1(L,4,2)=YVC(N,M+1)
c              XG1(L,4,3)=ZVC(N,M+1)
c            ENDDO
c          ENDDO
c
c      endif
c
c
cC -- End Tip (10/12/99)
C***********************************************************************
C/e S.N.KIM | Aug. 2018.

C-----------------------------------------------------------------------
C     Calculate control point locations and all the geometric
C     characteristics of the panels (local coordinates of the
C     four vertices centered at the centeroid, normal vectors,
C     moments of panels )
C-----------------------------------------------------------------------

C/Seungnam Kim, Shifting the radius of control points on the duct
      IF (IDUCT.EQ.1 .AND. IREPANEL.EQ.1) THEN
        CALL GEOM3D(NPANEL,XG1,CHRLEPS,IER)
        DO M = 1, MDUCT
          DO N = 1, NDUCT
            IF(N .LE. NDAFT .or. N .GE. NDUCT-NDAFT+1) THEN
              L = INDEXD(N,M)
              R111(L) = SQRT(XCT(L,2)**2+XCT(L,3)**2)
            ENDIF
          ENDDO
        ENDDO
      END IF
C/Seungnam Kim, Shifting the radius of control points on the duct

      CALL GEOM3D(NPANEL,XG,CHRLEPS,IER)
      IF(IER.EQ.0) THEN
        WRITE(*,'(A)') ' UNACCEPTABLE PANELS IN CONPT'
        STOP
      END IF

C     Correction of the radius of control points
      IF(IDUCT .EQ. 1 .AND. IDGEO .EQ. 1 .AND. DUCTPT .NE. 0.0) THEN
         DO M = 1 , MDUCT
            DO N = 1 , NDUCT
               L1 = INDEXD(N,M)
               THE1 = ATAN2(XCT(L1,3),XCT(L1,2))
               XCT(L1,2) = RRCP(N) * COS(THE1)
               XCT(L1,3) = RRCP(N) * SIN(THE1)
            ENDDO
         ENDDO
      ENDIF

C/Seungnam Kim, Shifting of the control points on the duct
      IF (IDUCT.EQ.1 .AND. IREPANEL.EQ.1) THEN
        DO M = 1, MDUCT
          DO N = 1, NDUCT
            IF(N .LE. NDAFT .or. N .GE. NDUCT-NDAFT+1) THEN
              L1 = INDEXD(N,M)
              THE1 = DANGLE(XCT(L1,3),XCT(L1,2))
              XCT(L1,2) = R111(L1)*COS(THE1)
              XCT(L1,3) = R111(L1)*SIN(THE1)
            ENDIF
          ENDDO
        ENDDO
        WRITE(*,*)'Control Points are adjusted due to repaneling.'
      END IF
C/Seungnam Kim, Shifting of the control points on the duct

C.....Find control points...............................................
      DO 60 J=1,NPANEL

         IF(ISP .NE. 0 ) CALL ROTATE2(-1,SPANGLE,XCT(J,1),XCT(J,2))

         RCP=SQRT(XCT(J,2)**2+XCT(J,3)**2)
         THP=ATAN2(XCT(J,3),XCT(J,2))
         XCTP(J,1,1)=XCT(J,1)
         XCTP(J,2,1)=XCT(J,2)
         XCTP(J,3,1)=XCT(J,3)
C----------------------------------------------------------------------
C    Next five lines + ENDIF added by HL on 010198 for hydrofoil case
C -- CASE for ICON = 5
C----------------------------------------------------------------------

c      IF(rhub .eq. 0.0 .or. ICON .EQ. 5) THEN
c         XCTP(J,1,2) = XCT(J,1)
c         XCTP(J,2,2) = -XCT(J,2)
c         XCTP(J,3,2) = XCT(J,3)
c      ELSE
C----------------------------------------------------------------------
C    Next five lines + ENDIF added by Shreenaath Natarajan on 040403 for hydrofoil case
C -- CASE for ICON = 5 and IHULL=1
C----------------------------------------------------------------------

      IF(IHULL .EQ. 1 .and. ICON .EQ. 5) THEN
         XCTP(J,1,2) = XCT(J,1)
         XCTP(J,2,2) = 2. -XCT(J,2)
         XCTP(J,3,2) = XCT(J,3)
      ELSE

C.......Calculate the control points on the other blades.................
C   Control points are located at the oposite position of real blade geometry
C   since we use image of control points for the cal. of inf. coeff's.

           IF(NBLADE.GT.1) THEN
              DO 50 KK=2,NBLADE
                 XCTP(J,1,KK)=XCT(J,1)
                 XCTP(J,2,KK)=RCP*COS(THP+DELK*(KK-1))
                 XCTP(J,3,KK)=RCP*SIN(THP+DELK*(KK-1))
 50           CONTINUE
           END IF

           IF(ISP .NE. 0) THEN
              CALL ROTATE2(0,SPANGLE,XCT(J,1),XCT(J,2))
              DO KK = 1 , NBLADE
                 CALL ROTATE2(0,SPANGLE,XCTP(J,1,KK),XCTP(J,2,KK))
              END DO
           END IF

      ENDIF

 60   CONTINUE


C -- Do not write anything here for unsteady case -- Yiran 09/19/2017 --
      IF (IUNS.EQ.1) RETURN

!s-- YE TIAN 08/13/2013 -----
!    The following block is implemented in the subroutien
!    write_fpg in the indpot.f file, which is called just before the end
!    of the program
!
C------------------------------(S.H.CHANG 02/25/2010)-------------------------
C    OUTPUT GEOMETRIES ON THE BLADE AND HUB FOR HULLFPP
!     WRITE(732,*) NBLADE,MR,MRP,NC,NCP,NTPREV,NCTIME,IHUB
!     WRITE(732,*) NRECL
!     WRITE(732,*) ADVCO

!     IF(IHUB.NE.0)  WRITE(732,*) NHBX,MHBT
!     WRITE(732,*) NPANEL

!     DO L = 1,NPANEL
!        DO M = 1,4
!           DO N = 1,3
!              WRITE(732,*) XG(L,M,N)
!           END DO
!        END DO
!     END DO
C------------------------------(S.H.CHANG 02/25/2010)-------------------------
!e-- YE TIAN 08/13/2013 -----

      RETURN
      END
