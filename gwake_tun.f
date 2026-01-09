      SUBROUTINE GWAKE_TUN
C************************************************************************
C     Gwake for Tunnel
C
C     No induced velocities are included in wake geometry 
C************************************************************************

      use m_BLADEM
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
!     COMMON/BLADEM/XBM(NBHZP,MBPZ),YBM(NBHZP,MBPZ),ZBM(NBHZP,MBPZ)


C-----------------------------------------------------------------------
C     Parameters of the transition wake
C-----------------------------------------------------------------------

      XFINAL=XULT
      RRHULT = SQRT( YB(1,1)**2+ZB(1,1)**2 )

      DTPROP=DELTAT/RAD
      DTF=DTPROP*ADVCO/180.0
      DTRADF=DELTAT

      NN=NCP
      MN=MRP

      NWK=NWPZ

      DO 60 M=1,MN

         XW(1,M)=XB(1,M)
         YW(1,M)=YB(1,M)
         ZW(1,M)=ZB(1,M)

         TWOLD=ATAN2(ZW(1,M),YW(1,M))

         DO 40 N=2,NWK

            DT1=DTF
            DTRAD=DTRADF

            XW(N,M) = XW(N-1,M) + (PITCH(M)/ADVCO)*DT1

            XWTE=XW(N,M)-XW(1,M)

            DDR = (1 - RRHULT) / REAL(MR)
            RWK = RRHULT + DDR*REAL(M-1)

            TW=ATAN2(ZW(N-1,M),YW(N-1,M))

            IF(TW.LT.TWOLD) TW=TW+TWOPI
            TW = TW + DTRAD

            YW(N,M)=RWK*COS(TW)
            ZW(N,M)=RWK*SIN(TW)

            TWOLD=TW
            NSW(M)=N

            IF(XW(N,M).GE.XUWDK) THEN
               GO TO 50
            END IF
40       CONTINUE

         WRITE(*,'(''!!!!! WARNING !!!!! TRANSITION WAKE SHORT'')')
         WRITE(*,'(''      M = '',I5,''  XW = '',F8.3)') M,XW(NSW(M),M)

50       CONTINUE

         IF(NSW(M).GT.NWZ) THEN
            WRITE(*,'('' M = '',I5,'' NSW(M) = '',I5,'' NWZ = '',I5)') 
     *           M,NSW(M),NWZ
            WRITE(*,'('' Decrease XUWDK in .geo, OR '')')
            WRITE(*,'('' Increase XWZ --> NSW(M) in PARAM.INC! '')')
            STOP
         END IF

         FK=( XW(NSW(M),M)-XUWDK)/( XW(NSW(M),M)-XW(NSW(M)-1,M) )
         XW(NSW(M),M)=XUWDK

         TW = TW - FK*DTRAD
         
         YW(NSW(M),M)=RWK*COS(TW)
         ZW(NSW(M),M)=RWK*SIN(TW)

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

      TANBUW=XVA(NX)/((PI*RULT/ADVCO)+XVT(NX))

      RHULT = SQRT(YW(NSW(1),1)**2+ZW(NSW(1),1)**2)
      RULT = SQRT(YW(NSW(MRP),MRP)**2+ZW(NSW(MRP),MRP)**2)

      CALL ULTWAK

      RETURN
      END
