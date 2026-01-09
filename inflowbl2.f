      SUBROUTINE INFLOWBL2
************************************************************************
*     INFLOWBL: calculate INFLOW velocities at all wake control points *
*     Created by  Hong SUN          11/05/05                           * 
*                          CAVITATING  CASE                            *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      DIMENSION U(3),NH1(3)
C     *     ,       VOX2(NSCWZ),VOY2(NSCWZ),VOZ2(NSCWZ)
C     *     ,       VOX3(NSCWZ),VOY3(NSCWZ),VOZ3(NSCWZ)

      ALLOCATABLE :: VOX2(:),VOY2(:),VOZ2(:)
      ALLOCATABLE :: VOX3(:),VOY3(:),VOZ3(:)

CVV
      ALLOCATE(VOX2(NSCWZ),VOY2(NSCWZ),VOZ2(NSCWZ))
      ALLOCATE(VOX3(NSCWZ),VOY3(NSCWZ),VOZ3(NSCWZ)) 
CVV 

      IF(NTSTEP.EQ.0.OR.ISTEADY.EQ.0) THEN
        DO 10 IH=1,3
          NH1(IH)=1               
 10     CONTINUE
      ELSE
         DO 15 IH=1,3
            NH1(IH)=NHARM(IH)               
 15      CONTINUE
      END IF


C***************************** On Wake Sub Panels **********************
      DO 40 J=1,NPWAKS

        RCP=SQRT(XCTPWs(J,2,1)**2+XCTPWs(J,3,1)**2)
        THP=ATAN2(XCTPWs(J,3,1),XCTPWs(J,2,1))
C.......Compute inflow velocities based on the wake harmonics...........
        DO 30 IH=1,3
          U(IH)=0.0
          DO 20 JH=1,NH1(IH)
C...........Spline wake harmonics to control radii (assume same.........
C...........number of panels for each blade)............................
            CALL EVALDKs(NWKCOE,1,XRW,RCP,COEC,XWCUB(1,JH,1,IH))
            CALL EVALDKs(NWKCOE,1,XRW,RCP,COES,XWCUB(1,JH,2,IH))
            T=FLOAT(JH-1)*(THP-TSTEP)
            ST=SIN(T)
            CT1=COS(T)
            U(IH)=U(IH)+COEC*CT1+COES*ST
 20       CONTINUE
 30     CONTINUE
        VOX3(J)=U(1)
        VORW=U(2)
        VOTW=U(3)
        COCP=XCTPWs(J,2,1)/RCP
        SICP=XCTPWs(J,3,1)/RCP
        WROVS=PI*RCP/ADVCO
        VOY3(J) = VORW*COCP-(WROVS+VOTW)*SICP
        VOZ3(J) = VORW*SICP+(WROVS+VOTW)*COCP
C---------------------------- Hydrofoil Case ---------------------------
        IF((ICON.EQ.5).OR.(ICON.EQ.6).OR.(ICON.EQ.8))THEN
           VOX3(J) = COS(ALPHA*RAD)
           VOY3(J) = 0.
           VOZ3(J) = SIN(ALPHA*RAD)
        ENDIF
 40   CONTINUE

C           On-coming velocity in local coordinate system
          DO 60 N=1,NTRA
            DO 50 M=1,MR
               L=NTRA*(MR-M)+N
               VXIWs(N,M) = VOX3(L)*DIRWsCP(L,1,1,1)
     *                     +VOY3(L)*DIRWsCP(L,1,2,1)
     *                     +VOZ3(L)*DIRWsCP(L,1,3,1)
               VETAWs(N,M)= VOX3(L)*DIRWsCP(L,2,1,1)
     *                     +VOY3(L)*DIRWsCP(L,2,2,1)
     *                     +VOZ3(L)*DIRWsCP(L,2,3,1)
               VINFSWs(N,M)=VOX3(L)**2+VOY3(L)**2+VOZ3(L)**2
 50         CONTINUE
 60       CONTINUE

C***************************** On Wake Panels **************************
      NWBIG=NWMIN*MR
      DO 140 J=1,NWBIG

        RCP=SQRT(XCTPW(J,2,1)**2+XCTPW(J,3,1)**2)
        THP=ATAN2(XCTPW(J,3,1),XCTPW(J,2,1))
C.......Compute inflow velocities based on the wake harmonics...........
        DO 130 IH=1,3
          U(IH)=0.0
          DO 120 JH=1,NH1(IH)
C...........Spline wake harmonics to control radii (assume same.........
C...........number of panels for each blade)............................
            CALL EVALDKs(NWKCOE,1,XRW,RCP,COEC,XWCUB(1,JH,1,IH))
            CALL EVALDKs(NWKCOE,1,XRW,RCP,COES,XWCUB(1,JH,2,IH))
            T=FLOAT(JH-1)*(THP-TSTEP)
            ST=SIN(T)
            CT1=COS(T)
            U(IH)=U(IH)+COEC*CT1+COES*ST
 120      CONTINUE
 130    CONTINUE
        VOX2(J) = U(1)
        VORW=U(2)
        VOTW=U(3)
        COCP=XCTPW(J,2,1)/RCP
        SICP=XCTPW(J,3,1)/RCP
        WROVS=PI*RCP/ADVCO
        VOY2(J) = VORW*COCP-(WROVS+VOTW)*SICP
        VOZ2(J) = VORW*SICP+(WROVS+VOTW)*COCP
C---------------------------- Hydrofoil Case ---------------------------
        IF((ICON.EQ.5).OR.(ICON.EQ.6).OR.(ICON.EQ.8))THEN
           VOX2(J) = COS(ALPHA*RAD)
           VOY2(J) = 0.
           VOZ2(J) = SIN(ALPHA*RAD)
        ENDIF
 140  CONTINUE

C           On-coming velocity in local coordinate system
           DO 160 N=1,NWMIN
             DO 150 M=1,MR
               L=INDEXW(N,M)               !!!
               VXIW(N,M) =  VOX2(L)*DIRWCP(L,1,1,1)
     *                     +VOY2(L)*DIRWCP(L,1,2,1)
     *                     +VOZ2(L)*DIRWCP(L,1,3,1)
               VETAW(N,M) = VOX2(L)*DIRWCP(L,2,1,1)
     *                     +VOY2(L)*DIRWCP(L,2,2,1)
     *                     +VOZ2(L)*DIRWCP(L,2,3,1)
               VINFSW(N,M)=VOX2(L)**2+VOY2(L)**2+VOZ2(L)**2
 150        CONTINUE
 160      CONTINUE

CVV
      DEALLOCATE(VOX2,VOY2,VOZ2)
      DEALLOCATE(VOX3,VOY3,VOZ3)
CVV

      RETURN
C>>>>>>>>>>>>>>>>>>>>>>End of subroutine INFLOWBL2>>>>>>>>>>>>>>>>>>>>>>
      END
