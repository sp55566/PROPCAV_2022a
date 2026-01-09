C =====================================
      SUBROUTINE PODGWAKE
C =====================================
      use m_WAKGEO
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
!     PARAMETER(NWR=MBPZ)
!     COMMON /WAKGEO/ RTBL(11,NWR),XTBL(11),RW(NWR),VAR(NWR),VTR(NWR),
!    *       UASTAR(NWR),UTSTAR(NWR),UAUSTR(NWR),UTUSTR(NWR)
      DIMENSION RWK(NWPZ,MBPZ),DRWP(NWPZ)
      integer NWR
      NWR=MBPZ

      MN  = MRP
      NWK = NWPZ
      WRITE(*,*) MRP,NWPZ,NSW(1)

      DO M=1,MN
         RWK(1,M) = SQRT(YW(1,M)**2+ZW(1,M)**2)
      ENDDO
      
C----- PART 1
      RPOD1 = RWK(1,1)
      DO N=2,NSW(1)
         XWL = XW(N,1)
         CALL EVALDKs(NHIN,1,XHIN,XWL,RPOD2,XYHINCUB)
         DR = RPOD2-RPOD1
         IF(DR.LT.0.0) GO TO 10
         RPOD1 = RPOD2
      ENDDO

 10   XCST  = XW(N-1,1)
      CALL EVALDKs(NHIN,1,XHIN,XCST,RPCST,XYHINCUB)
      DRCST = RPCST-SQRT(YW(N-1,1)**2+ZW(N-1,1)**2)
      
C----- PART 2
      DRWP(1)=0.0
      DO 20 N=2,NSW(1)
         XWL = XW(N,1)
         IF(XWL.LE.XCST) THEN
            CALL EVALDKs(NHIN,1,XHIN,XWL,RPXL,XYHINCUB)
            RWWL = SQRT(YW(N,1)**2+ZW(N,1)**2)
            DRWP(N) = RPXL-RWWL 
         ELSE
            DRWP(N) = DRCST
         ENDIF
 20   ENDDO

C----- PART 3 
c----- EVALUATE INDUCED VELOCITIES
C-----------------------------------------------------------------------
C     Interpolate velocities in the wake
C----------------------------------------------------------------------
      XFINAL = XULT
      DTPROP = DELTAT/RAD
      DTF    = DTPROP*ADVCO/180.0
      DTRADF = DELTAT
      DT1    = DTF
      DO 21 M = 1,MN
         DO 22 N=2,NSW(M)
            TWOLD = ATAN2(ZW(1,M),YW(1,M))
C           Hong fixed the bug     11/28/07
            IF(N.LE.NSW(1)) THEN
              RWAKE = SQRT(YW(N,M)**2+ZW(N,M)**2)+DRWP(N)
            ELSE
              RWAKE = SQRT(YW(N,M)**2+ZW(N,M)**2)+DRWP(NSW(1))
            ENDIF
            CALL EVALDKs(NX,1,XR,RWAKE,WVAR,VACUB)
            CALL EVALDKs(NX,1,XR,RWAKE,WVTR,VTCUB)
            CALL EVALDKs(NX,1,XR,RWAKE,WUASTAR,UACUB)
            CALL EVALDKs(NX,1,XR,RWAKE,WUAUSTR,UAUCUB)
            CALL EVALDKs(NX,1,XR,RWAKE,WUTSTAR,UTCUB)
            CALL EVALDKs(NX,1,XR,RWAKE,WUTUSTR,UTUCUB)
            
            WUAINC = WUAUSTR-WUASTAR
            WUTINC = WUTUSTR-WUTSTAR
            Q  = (XW(N-1,M)-XW(1,M))/XFINAL
            TW = ATAN2(ZW(N-1,M),YW(N-1,M))
            IF(TW.LT.TWOLD)TW=TW+TWOPI
            TW = TW+DTRADF+(WVTR+WUTSTAR+GROW(Q)*WUTINC)*DT1/RWAKE
            YW(N,M)=RWAKE*COS(TW)
            ZW(N,M)=RWAKE*SIN(TW)
 22      ENDDO
 21   ENDDO
 

      RETURN
C))))))))))))))))))))) End of subroutine GWAKE (((((((((((((((((((((((((
      END


