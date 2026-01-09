      SUBROUTINE INFLOWBL
************************************************************************
*     INFLOWBL: calculate INFLOW velocities at all wake control points *
*     Created by  Hong SUN          11/05/05                           * 
*                             WETTED  CASE                             *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      DIMENSION U(3),NH1(3)
     *     ,       VOX2(NSCWZ),VOY2(NSCWZ),VOZ2(NSCWZ)
     *     ,       VOX3(NSCWZ),VOY3(NSCWZ),VOZ3(NSCWZ)

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
        VOX3(J) = U(1)
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
          DO 60 N=1,NWSUB
            DO 50 M=1,MR
               L=NWSUB*(M-1)+N
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
      DO 140 J=1,NPWAKE

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
               L=IDXWAK(N,M)
               VXIW(N,M) =  VOX2(L)*DIRWCP(L,1,1,1)
     *                     +VOY2(L)*DIRWCP(L,1,2,1)
     *                     +VOZ2(L)*DIRWCP(L,1,3,1)
               VETAW(N,M) = VOX2(L)*DIRWCP(L,2,1,1)
     *                     +VOY2(L)*DIRWCP(L,2,2,1)
     *                     +VOZ2(L)*DIRWCP(L,2,3,1)
               VINFSW(N,M)=VOX2(L)**2+VOY2(L)**2+VOZ2(L)**2
 150        CONTINUE
 160      CONTINUE



      IF(IDUCT.NE.0) THEN
C************************* On Duct Wake Sub Panels ********************
      DO 240 J=1,NPWAKSD

        RCP=SQRT(XCTPDWs(J,2,1)**2+XCTPDWs(J,3,1)**2)
        THP=ATAN2(XCTPDWs(J,3,1),XCTPDWs(J,2,1))
C.......Compute inflow velocities based on the wake harmonics...........
        DO 230 IH=1,3
          U(IH)=0.0
          DO 220 JH=1,NH1(IH)
C...........Spline wake harmonics to control radii (assume same.........
C...........number of panels for each blade)............................
            CALL EVALDKs(NWKCOE,1,XRW,RCP,COEC,XWCUB(1,JH,1,IH))
            CALL EVALDKs(NWKCOE,1,XRW,RCP,COES,XWCUB(1,JH,2,IH))
            T=FLOAT(JH-1)*(THP-TSTEP)
            ST=SIN(T)
            CT1=COS(T)
            U(IH)=U(IH)+COEC*CT1+COES*ST
 220       CONTINUE
 230     CONTINUE
   
        VOX3(J) = U(1)
        VORW=U(2)
        VOTW=U(3)
        COCP=XCTPDWs(J,2,1)/RCP
        SICP=XCTPDWs(J,3,1)/RCP
C        WROVS=PI*RCP/ADVCO
        WROVS= 0.0
        VOY3(J) = VORW*COCP-(WROVS+VOTW)*SICP
        VOZ3(J) = VORW*SICP+(WROVS+VOTW)*COCP
 240   CONTINUE

C           On-coming velocity in local coordinate system
          DO 260 N=1,NWSUB
            DO 250 M=1,MDUCT
               L=NWSUB*(M-1)+N
               VXIDWs(N,M) = VOX3(L)*DIRDWsCP(L,1,1,1)
     *                     + VOY3(L)*DIRDWsCP(L,1,2,1)
     *                     + VOZ3(L)*DIRDWsCP(L,1,3,1)
               VETADWs(N,M)= VOX3(L)*DIRDWsCP(L,2,1,1)
     *                     + VOY3(L)*DIRDWsCP(L,2,2,1)
     *                     + VOZ3(L)*DIRDWsCP(L,2,3,1)
               VINFSDWs(N,M)=VOX3(L)**2+VOY3(L)**2+VOZ3(L)**2
 250         CONTINUE
 260       CONTINUE

C************************** On Duct Wake Panels ***********************
      DO 340 J=1,NPWAKED

        RCP=SQRT(XCTPDW(J,2,1)**2+XCTPDW(J,3,1)**2)
        THP=ATAN2(XCTPDW(J,3,1),XCTPDW(J,2,1))
C.......Compute inflow velocities based on the wake harmonics...........
        DO 330 IH=1,3
          U(IH)=0.0
          DO 320 JH=1,NH1(IH)
C...........Spline wake harmonics to control radii (assume same.........
C...........number of panels for each blade)............................
            CALL EVALDKs(NWKCOE,1,XRW,RCP,COEC,XWCUB(1,JH,1,IH))
            CALL EVALDKs(NWKCOE,1,XRW,RCP,COES,XWCUB(1,JH,2,IH))
            T=FLOAT(JH-1)*(THP-TSTEP)
            ST=SIN(T)
            CT1=COS(T)
            U(IH)=U(IH)+COEC*CT1+COES*ST
 320      CONTINUE
 330    CONTINUE
     
        VOX2(J) = U(1)
        VORW=U(2)
        VOTW=U(3)
        COCP=XCTPDW(J,2,1)/RCP
        SICP=XCTPDW(J,3,1)/RCP
C        WROVS=PI*RCP/ADVCO
        WROVS= 0.0
        VOY2(J) = VORW*COCP-(WROVS+VOTW)*SICP
        VOZ2(J) = VORW*SICP+(WROVS+VOTW)*COCP
 340  CONTINUE

C           On-coming velocity in local coordinate system
           DO 360 N=1,NDWK
             DO 350 M=1,MDUCT
               L=INDEXWD(N,M)
               VXIDW(N,M) =  VOX2(L)*DIRDWCP(L,1,1,1)
     *                     + VOY2(L)*DIRDWCP(L,1,2,1)
     *                     + VOZ2(L)*DIRDWCP(L,1,3,1)
               VETADW(N,M) = VOX2(L)*DIRDWCP(L,2,1,1)
     *                     + VOY2(L)*DIRDWCP(L,2,2,1)
     *                     + VOZ2(L)*DIRDWCP(L,2,3,1)
               VINFSDW(N,M)=VOX2(L)**2+VOY2(L)**2+VOZ2(L)**2
 350        CONTINUE
 360      CONTINUE
      ENDIF


      RETURN
C>>>>>>>>>>>>>>>>>>>>>>End of subroutine INFLOWBL>>>>>>>>>>>>>>>>>>>>>>>
      END
