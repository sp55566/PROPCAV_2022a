      SUBROUTINE VISPOT2(JM)
C **********************************************************************.
C     Calculate the total potential along each WAKE STRIP.             *
C                      Cavitating   Case                               *
C                                                                      *
C            By Hong Sun                June, 2006                     *
C **********************************************************************.
 
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      INCLUDE 'PUFBL.INC'  
      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)
C      DIMENSION TEMPSUB(NPAWZ+MBZ,KZ)
      ALLOCATABLE :: TEMPSUB(:,:)

CVV
      ALLOCATE(TEMPSUB(NPAWZ+MBZ,KZ))
CVV
            
C     First Calculate SINGULARITY geometry info
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
 10     CONTINUE
 20   CONTINUE

      IF(IHUB.NE.0) THEN
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
 30       CONTINUE
 40     CONTINUE
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

            ENDDO
         ENDDO
      END IF

C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
cC.....FOR TIP VORTEX
c      IF(IAN .EQ. 2) THEN
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
c            ENDDO
c          ENDDO
c      ENDIF
C/e S.N.KIM | Aug. 2018.

      CALL GEOM3D(NPANEL,XG,CHRLEPS,IER)
      IF(IER.EQ.0) THEN
        WRITE(*,'(A)') ' UNACCEPTABLE PANELS IN CONPT'
        STOP
      END IF

C     Compute influence coefficients due to the BLADE
C*********************************************************************

      DO 160 M = MR, 1, -1
       DO 150 N = 1, NC
        J = INDEXB(N,M)
        DO 110 K=1,4
          XV(K)=XVP(J,K)
          YV(K)=YVP(J,K)
          SIDE(K)=SID(J,K)
 110    CONTINUE

        DO 120 K=1,15
          S(K)=SS(J,K)
 120    CONTINUE

        XM1(1)=XB(N,M)
        XM1(2)=YB(N,M) 
        XM1(3)=ZB(N,M)
        XM2(1)=XB(N,M+1)
        XM2(2)=YB(N,M+1)
        XM2(3)=ZB(N,M+1)
        XM3(1)=XB(N+1,M+1)
        XM3(2)=YB(N+1,M+1)
        XM3(3)=ZB(N+1,M+1)
        XM4(1)=XB(N+1,M)
        XM4(2)=YB(N+1,M)
        XM4(3)=ZB(N+1,M)   
        IMR0=0 
            
        ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))+ABS(XM2(3)-XM1(3))
        ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))+ABS(XM3(3)-XM2(3))
        ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))+ABS(XM4(3)-XM3(3))
        ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))+ABS(XM1(3)-XM4(3))
        ERRMIN=AMIN1(ER1,ER2,ER3,ER4) 
        IF(ERRMIN.LE.1.0E-6) THEN
          IMR0=1
        END IF
        
        DO 149 KK = 1, NBLADE
   
           IF(ISTEADY.GT.0)THEN      
          IREC = NTPOS(KK)
          CALL READ2(45,IREC,POT,NPANEL)      
          CALL READ2(47,IREC,DPDNC,NPANEL)         
         ENDIF
         
C------- LOOP 130: FOR WAKE SUBPANELS    (Hong Sun, 04/07/2005)     
        DO 130 IN = 1,NTRA
          L = NTRA*(MR-JM)+IN

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,1,K)
             YLOC=YLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,2,K)
             ZLOC=ZLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPWs(L,IXYZ,KK)
            ENDDO
               CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF
                  
          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          END IF

c         Calculate Total potential
          PPOTWs(L) = PPOTWs(L)+ (-POT(J)*FD +DPDNC(J)*FS)/PI/4.
      
 130    CONTINUE

C------- LOOP 140: FOR WAKE PANELS 
        DO 140 IN = 1,NWMIN
          L = INDEXW(IN,JM)

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,1,K)
             YLOC=YLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,2,K)
             ZLOC=ZLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPW(L,IXYZ,KK)
            ENDDO
               CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF
                  
          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          END IF

c         Calculate Total potential
          PPOTW(L) = PPOTW(L)+ (-POT(J)*FD +DPDNC(J)*FS)/PI/4.
 
 140    CONTINUE 

 149   CONTINUE       
 150   CONTINUE
 160  CONTINUE

C     Compute influence coefficients due to the hub
      IF(IHUB.NE.0) THEN
       DO 260 N = 1, NHBX
        DO 250 M = 1, MHBT
         J = INDEXH(N,M)
         DO 210 K=1,4
          XV(K)=XVP(J,K)
          YV(K)=YVP(J,K)
          SIDE(K)=SID(J,K)
 210     CONTINUE

         DO 220 K=1,15
          S(K)=SS(J,K)
 220     CONTINUE

        XM1(1)=XH(N,M+1)
        XM1(2)=YH(N,M+1) 
        XM1(3)=ZH(N,M+1)
        XM2(1)=XH(N,M)
        XM2(2)=YH(N,M)
        XM2(3)=ZH(N,M)
        XM3(1)=XH(N+1,M)
        XM3(2)=YH(N+1,M)
        XM3(3)=ZH(N+1,M)
        XM4(1)=XH(N+1,M+1)
        XM4(2)=YH(N+1,M+1)
        XM4(3)=ZH(N+1,M+1)   
        IMR0=0 
            
        ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))+ABS(XM2(3)-XM1(3))
        ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))+ABS(XM3(3)-XM2(3))
        ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))+ABS(XM4(3)-XM3(3))
        ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))+ABS(XM1(3)-XM4(3))
        ERRMIN=AMIN1(ER1,ER2,ER3,ER4) 
        IF(ERRMIN.LE.1.0E-6) THEN
          IMR0=1
        END IF

C--------LOOP 230: FOR WAKE SUBPANELS  (Hong Sun, 04/07/2005) 
        DO 249 KK = 1, NBLADE
 
        IF(ISTEADY.GT.0)THEN          
          IREC = NTPOS(KK)    
          CALL READ2(45,IREC,POT,NPANEL)      
          CALL READ2(47,IREC,DPDNC,NPANEL)           
        ENDIF      
              
         DO 230 IN = 1,NTRA       
          L = NTRA*(MR-JM)+IN

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,1,K)
             YLOC=YLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,2,K)
             ZLOC=ZLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPWs(L,IXYZ,KK)
            ENDDO
            CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF
                  
          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          END IF 

C         Calculate Total potential
          PPOTWs(L) = PPOTWs(L)+ (-POT(J)*FD +DPDNC(J)*FS)/PI/4.
 
 230     CONTINUE 


C--------LOOP 240: FOR WAKE PANELS   
         DO 240 IN = 1,NWMIN 
          L = INDEXW(IN,JM)

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,1,K)
             YLOC=YLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,2,K)
             ZLOC=ZLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPW(L,IXYZ,KK)
            ENDDO
            CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF
                  
          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          END IF

c         Calculate Total potential
          PPOTW(L) = PPOTW(L)+ (-POT(J)*FD +DPDNC(J)*FS)/PI/4.
 
 240     CONTINUE 

 249    CONTINUE      
 250    CONTINUE
 260   CONTINUE
      ENDIF

C     Compute influence coefficients due to the duct
C*********************************************************************
      IF(IDUCT.NE.0) THEN
       DO 360 N = 1, NDUCT
        DO 350 M = 1, MDUCT
         J = INDEXD(N,M)
         DO 310 K=1,4
          XV(K)=XVP(J,K)
          YV(K)=YVP(J,K)
          SIDE(K)=SID(J,K)
 310     CONTINUE

         DO 320 K=1,15
          S(K)=SS(J,K)
 320     CONTINUE

        XM1(1)=XD(N,M+1)
        XM1(2)=YD(N,M+1) 
        XM1(3)=ZD(N,M+1)
        XM2(1)=XD(N,M)
        XM2(2)=YD(N,M)
        XM2(3)=ZD(N,M)
        XM3(1)=XD(N+1,M)
        XM3(2)=YD(N+1,M)
        XM3(3)=ZD(N+1,M)
        XM4(1)=XD(N+1,M+1)
        XM4(2)=YD(N+1,M+1)
        XM4(3)=ZD(N+1,M+1)   
        IMR0=0 
            
        ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))+ABS(XM2(3)-XM1(3))
        ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))+ABS(XM3(3)-XM2(3))
        ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))+ABS(XM4(3)-XM3(3))
        ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))+ABS(XM1(3)-XM4(3))
        ERRMIN=AMIN1(ER1,ER2,ER3,ER4) 
        IF(ERRMIN.LE.1.0E-6) THEN
          IMR0=1
        END IF

C--------LOOP 330: FOR WAKE SUBPANELS  (Hong Sun, 06/04/2006) 
        DO 349 KK = 1, NBLADE
 
        IF(ISTEADY.GT.0)THEN         
          IREC = NTPOS(KK)          
          CALL READ2(45,IREC,POT,NPANEL)      
          CALL READ2(47,IREC,DPDNC,NPANEL) 
        ENDIF      
              
         DO 330 IN = 1,NTRA
          L = NTRA*(MR-JM)+IN
C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,1,K)
             YLOC=YLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,2,K)
             ZLOC=ZLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPWs(L,IXYZ,KK)
            ENDDO
            CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF
                  
          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          END IF 

C         Calculate Total potential
          PPOTWs(L) = PPOTWs(L)+ (-POT(J)*FD +DPDNC(J)*FS)/PI/4.
 
 330     CONTINUE 


C--------LOOP 340: FOR WAKE PANELS   
         DO 340 IN = 1,NWMIN 
          L = INDEXW(IN,JM)

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,1,K)
             YLOC=YLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,2,K)
             ZLOC=ZLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPW(L,IXYZ,KK)
            ENDDO
            CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF
                  
          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          END IF

c         Calculate Total potential
          PPOTW(L) = PPOTW(L)+ (-POT(J)*FD +DPDNC(J)*FS)/PI/4.
 
 340     CONTINUE 

 349    CONTINUE      
 350    CONTINUE
 360   CONTINUE
      ENDIF

C     Compute influence coefficients due to the TUNNEL
C*********************************************************************
      IF(ITUN.NE.0) THEN
       DO 460 N = 1, NAXT
        DO 450 M = 1, MTUNEL
         J = INDEXTN(N,M)
         DO 410 K=1,4
          XV(K)=XVP(J,K)
          YV(K)=YVP(J,K)
          SIDE(K)=SID(J,K)
 410     CONTINUE

         DO 420 K=1,15
          S(K)=SS(J,K)
 420     CONTINUE

        XM1(1)=XTUN(N,M)
        XM1(2)=YTUN(N,M) 
        XM1(3)=ZTUN(N,M)
        XM2(1)=XTUN(N,M+1)
        XM2(2)=YTUN(N,M+1)
        XM2(3)=ZTUN(N,M+1)
        XM3(1)=XTUN(N+1,M+1)
        XM3(2)=YTUN(N+1,M+1)
        XM3(3)=ZTUN(N+1,M+1)
        XM4(1)=XTUN(N+1,M)
        XM4(2)=YTUN(N+1,M)
        XM4(3)=ZTUN(N+1,M)   
        IMR0=0 
            
        ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))+ABS(XM2(3)-XM1(3))
        ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))+ABS(XM3(3)-XM2(3))
        ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))+ABS(XM4(3)-XM3(3))
        ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))+ABS(XM1(3)-XM4(3))
        ERRMIN=AMIN1(ER1,ER2,ER3,ER4) 
        IF(ERRMIN.LE.1.0E-6) THEN
          IMR0=1
        END IF

C--------LOOP 330: FOR WAKE SUBPANELS  (Hong Sun, 06/04/2006) 
        DO 449 KK = 1, NBLADE
  
        IF(ISTEADY.GT.0)THEN        
          IREC = NTPOS(KK)     
          CALL READ2(45,IREC,POT,NPANEL)      
          CALL READ2(47,IREC,DPDNC,NPANEL)     
        ENDIF      
              
         DO 430 IN = 1,NTRA
          L = NTRA*(MR-JM)+IN
C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,1,K)
             YLOC=YLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,2,K)
             ZLOC=ZLOC+(XCTPWs(L,K,KK)-XCT(J,K))*DIR(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPWs(L,IXYZ,KK)
            ENDDO
            CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF
                  
          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          END IF 

C         Calculate Total potential
          PPOTWs(L) = PPOTWs(L)+ (-POT(J)*FD +DPDNC(J)*FS)/PI/4.
 
 430     CONTINUE 


C--------LOOP 440: FOR WAKE PANELS   
         DO 440 IN = 1,NWMIN 
          L = INDEXW(IN,JM)

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,1,K)
             YLOC=YLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,2,K)
             ZLOC=ZLOC+(XCTPW(L,K,KK)-XCT(J,K))*DIR(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPW(L,IXYZ,KK)
            ENDDO
            CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF
                  
          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          END IF

c         Calculate Total potential
          PPOTW(L) = PPOTW(L)+ (-POT(J)*FD +DPDNC(J)*FS)/PI/4.
 
 440     CONTINUE 

 449    CONTINUE      
 450    CONTINUE
 460   CONTINUE
      ENDIF

C*********************************************************************
C               Wake Sources          Hong Sun 11/02/2006
C*********************************************************************
      NPWAKS=NTRA*MR
      DO 520 N=1,NTRA
        DO 510 M=1,MR
          L=NTRA*(MR-M)+N
          XGW(L,1,1)=XWS(N,M)
          XGW(L,1,2)=YWS(N,M)
          XGW(L,1,3)=ZWS(N,M)
          XGW(L,2,1)=XWS(N,M+1)
          XGW(L,2,2)=YWS(N,M+1)
          XGW(L,2,3)=ZWS(N,M+1)
          XGW(L,3,1)=XWS(N+1,M+1)
          XGW(L,3,2)=YWS(N+1,M+1)
          XGW(L,3,3)=ZWS(N+1,M+1)
          XGW(L,4,1)=XWS(N+1,M)
          XGW(L,4,2)=YWS(N+1,M)
          XGW(L,4,3)=ZWS(N+1,M)
 510     CONTINUE
 520   CONTINUE

      CALL GEO3DW(NPWAKS,XGW,CHRLEWS,IER)

      IF(IER.EQ.0) THEN
         WRITE(*,*) 'UNACCEPTABLE PANEL IN SCWKINF'
         STOP
      END IF
C
C
      DO 560 M = MR, 1, -1
       DO 550 N = 1, NTRA
        J = NTRA*(MR-M)+N
        DO K=1,4
          XV(K)=XVPW(J,K)
          YV(K)=YVPW(J,K)
          SIDE(K)=SIDW(J,K)
        ENDDO

        DO K=1,15
          S(K)=SSW(J,K)
        ENDDO

        XM1(1)=XWS(N,M)
        XM1(2)=YWS(N,M) 
        XM1(3)=ZWS(N,M)
        XM2(1)=XWS(N,M+1)
        XM2(2)=YWS(N,M+1)
        XM2(3)=ZWS(N,M+1)
        XM3(1)=XWS(N+1,M+1)
        XM3(2)=YWS(N+1,M+1)
        XM3(3)=ZWS(N+1,M+1)
        XM4(1)=XWS(N+1,M)
        XM4(2)=YWS(N+1,M)
        XM4(3)=ZWS(N+1,M)   
        IMR0=0 

        ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))+ABS(XM2(3)-XM1(3))
        ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))+ABS(XM3(3)-XM2(3))
        ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))+ABS(XM4(3)-XM3(3))
        ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))+ABS(XM1(3)-XM4(3))
        ERRMIN=AMIN1(ER1,ER2,ER3,ER4) 
        IF(ERRMIN.LE.1.0E-6) THEN
          IMR0=1
        END IF

        DO 549 KK=1,NBLADE
      
C        Calculate the dipole strength on wake subpanels      
         IF(KK.EQ.1)THEN  
       
           IF(N.LE.(4*N1SUB)) THEN 
           JSUB = N/N1SUB+1
           IF(MOD(N,N1SUB).EQ.0) JSUB=JSUB-1 
              JMAC=MR*(JSUB-1)+M
            JMAC1=MR*JSUB+M
                
            SLOPE=(TEMP4(JMAC1)-TEMP4(JMAC))/(XW(JSUB+1,M)-XW(JSUB,M))
C            TEMPSUB(J,KK)=SLOPE*(XCTPWs(J,1,1)-XW(JSUB,M))+TEMP4(JMAC)
            TEMPSUB(J,KK)=SLOPE*(0.5*(XWS(N,M)+XWS(N+1,M))-XW(JSUB,M))
     &                     +TEMP4(JMAC)
C            TEMPSUB(J,KK)=(0.5-(FLOAT(N)-0.5)/FLOAT(4*N1SUB))*TEMP4(JMAC)
        
         ELSE
            JSUB=N-4*N1SUB+4
               JMAC=MR*(JSUB-1)+M
            TEMPSUB(J,KK) = TEMP4(JMAC)                  
            ENDIF
       
c         IF(JM.EQ.8.AND.M.EQ.8)THEN
c             WRITE(139,*)J,XCTPWs(J,1,1),TEMPSUB(J,KK)         !!!!!!!!
c         ENDIF
       
       ELSE
           IREC = NTPOS(KK)      
           CALL READ2(48,IREC,SORW,NPWAKS)      
       ENDIF 
                         
C------- LOOP 530: FOR WAKE SUBPANELS    (Hong Sun, 11/02/2006)     
        DO 530 IN = 1,NTRA
          L = NTRA*(MR-JM)+IN

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPWs(L,K,KK)-XCTW(J,K))*DIRW(J,1,K)
             YLOC=YLOC+(XCTPWs(L,K,KK)-XCTW(J,K))*DIRW(J,2,K)
             ZLOC=ZLOC+(XCTPWs(L,K,KK)-XCTW(J,K))*DIRW(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)      
                
          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPWs(L,IXYZ,KK)
            ENDDO
               CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF

           IF(L.EQ.J .AND. KK.EQ.1) THEN
             IF(IMR.EQ.1) THEN
                FD=TWOPI
             END IF
             IF(FD.LT.0.) THEN
                FD=4.0*PI+FD
             END IF
                     
             IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.6.0) THEN               
                FD=2.0*PI
             END IF
           ELSE       
             IF(ABS(FD).GT.6.28) THEN                       
                FD=0.0
             ENDIF
           END IF
        
c         Calculate Total potential

            PPOTWs(L)=PPOTWs(L)+SORW(J)*FS/PI/4.
          IF(KK.EQ.1) THEN
            PPOTWs(L)=PPOTWs(L)-TEMPSUB(J,KK)*FD/PI/4. 
c             WRITE(138,*)L,J,TEMPSUB(J,KK) 
            ENDIF 
 530    CONTINUE

C------- LOOP 540: FOR WAKE PANELS 
        DO 540 IN = 1,NWMIN
          L = INDEXW(IN,JM)

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPW(L,K,KK)-XCTW(J,K))*DIRW(J,1,K)
             YLOC=YLOC+(XCTPW(L,K,KK)-XCTW(J,K))*DIRW(J,2,K)
             ZLOC=ZLOC+(XCTPW(L,K,KK)-XCTW(J,K))*DIRW(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPW(L,IXYZ,KK)
            ENDDO
               CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF
                  
          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          END IF

c         Calculate Total potential          
 
          PPOTW(L)=PPOTW(L)+SORW(J)*FS/PI/4.         

 540    CONTINUE 
 549    CONTINUE    

 550   CONTINUE
 560  CONTINUE
  
C
C---------------------------------------------------------------------
C
      NWBIG = NWMIN*MR
      DO 620 N=1, NWMIN
        DO 610 M=1,MR
          L=INDEXW(N,M)
          XGW(L,1,1)=XW(N,M)
          XGW(L,1,2)=YW(N,M)
          XGW(L,1,3)=ZW(N,M)
          XGW(L,2,1)=XW(N,M+1)
          XGW(L,2,2)=YW(N,M+1)
          XGW(L,2,3)=ZW(N,M+1)
          XGW(L,3,1)=XW(N+1,M+1)
          XGW(L,3,2)=YW(N+1,M+1)
          XGW(L,3,3)=ZW(N+1,M+1)
          XGW(L,4,1)=XW(N+1,M)
          XGW(L,4,2)=YW(N+1,M)
          XGW(L,4,3)=ZW(N+1,M)
 610     CONTINUE
 620   CONTINUE
        
      CALL GEO3DW(NWBIG,XGW,CHRLEWS,IER)

      IF(IER.EQ.0) THEN
         WRITE(*,*) 'UNACCEPTABLE PANEL IN SCWKINF'
         STOP
      END IF
C
C
      DO 660 M = MR, 1, -1
       DO 650 N = 1, NWMIN
        J = INDEXW(N,M)      
      L1 = MR*(N-1)+M                  !index for wake dipole strength
      
        DO K=1,4
          XV(K)=XVPW(J,K)
          YV(K)=YVPW(J,K)
          SIDE(K)=SIDW(J,K)
        ENDDO

        DO K=1,15
          S(K)=SSW(J,K)
        ENDDO

        XM1(1)=XW(N,M)
        XM1(2)=YW(N,M) 
        XM1(3)=ZW(N,M)
        XM2(1)=XW(N,M+1)
        XM2(2)=YW(N,M+1)
        XM2(3)=ZW(N,M+1)
        XM3(1)=XW(N+1,M+1)
        XM3(2)=YW(N+1,M+1)
        XM3(3)=ZW(N+1,M+1)
        XM4(1)=XW(N+1,M)
        XM4(2)=YW(N+1,M)
        XM4(3)=ZW(N+1,M)   
        IMR0=0 

        ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))+ABS(XM2(3)-XM1(3))
        ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))+ABS(XM3(3)-XM2(3))
        ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))+ABS(XM4(3)-XM3(3))
        ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))+ABS(XM1(3)-XM4(3))
        ERRMIN=AMIN1(ER1,ER2,ER3,ER4) 
        IF(ERRMIN.LE.1.0E-6) THEN
          IMR0=1
        END IF
      
c       IF(JM.EQ.8.AND.M.EQ.8)THEN
c           WRITE(139,*)J,XCTPW(J,1,1),TEMP4(L1)                   !!!!!!!!
c       ENDIF
       
        DO 649 KK=1, NBLADE
      
        IF(KK.GT.1)THEN        
          IREC = NTPOS(KK)     
          CALL READ2(46,IREC,TEMP5,NWBIG) 
         ENDIF
        
                     
C------- LOOP 630: FOR WAKE SUBPANELS    (Hong, 11/02/2006)     
        DO 630 IN = 1,NTRA
          L = NTRA*(MR-JM)+IN

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPWs(L,K,KK)-XCTW(J,K))*DIRW(J,1,K)
             YLOC=YLOC+(XCTPWs(L,K,KK)-XCTW(J,K))*DIRW(J,2,K)
             ZLOC=ZLOC+(XCTPWs(L,K,KK)-XCTW(J,K))*DIRW(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPWs(L,IXYZ,KK)
            ENDDO
               CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF
                  
          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          END IF

c         Calculate Total potential
       IF(KK.EQ.1.AND.N.GT.NSUB) THEN
         PPOTWs(L) = PPOTWs(L)- TEMP4(L1)*FD/PI/4.
         ELSE IF(KK.GT.1) THEN
           PPOTWs(L) = PPOTWs(L)- TEMP5(L1)*FD/PI/4.       
         ENDIF
 630    CONTINUE

C------- LOOP 640: FOR WAKE PANELS 
        DO 640 IN = 1,NWMIN
          L = INDEXW(IN,JM)

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPW(L,K,KK)-XCTW(J,K))*DIRW(J,1,K)
             YLOC=YLOC+(XCTPW(L,K,KK)-XCTW(J,K))*DIRW(J,2,K)
             ZLOC=ZLOC+(XCTPW(L,K,KK)-XCTW(J,K))*DIRW(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPW(L,IXYZ,KK)
            ENDDO
               CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF

           IF(L.EQ.J .AND. KK.EQ.1) THEN
             IF(IMR.EQ.1) THEN
                FD=TWOPI
             END IF
             IF(FD.LT.0.) THEN
                FD=4.0*PI+FD
             END IF
                     
             IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.6.0) THEN               
                FD=2.0*PI
             END IF
           ELSE       
             IF(ABS(FD).GT.6.28) THEN                       
                FD=0.0
             ENDIF
           END IF                  


c         Calculate Total potential          
         IF(KK.EQ.1) THEN
           PPOTW(L) = PPOTW(L)- TEMP4(L1)*FD/PI/4.
         ELSE 
           PPOTW(L) = PPOTW(L)- TEMP5(L1)*FD/PI/4.       
         ENDIF

 640    CONTINUE 
 649    CONTINUE    

 650   CONTINUE
 660  CONTINUE



      IF(IDUCT.NE.0.AND.IDOPT.EQ.0) THEN
C*********************************************************************
C                Duct Wake Dipoles             Hong 11/02/2006
C*********************************************************************
C                DUCT WAKE SUB PANELS
C
c      
c      NPWAKSD=NWSUB*MDUCT
c      DO 720 N=1,NWSUB
c        DO 710 M=1,MDUCT
c          L=NWSUB*(M-1)+N
c          XGW(L,1,1)=XWSD(N,M)
c          XGW(L,1,2)=YWSD(N,M)
c          XGW(L,1,3)=ZWSD(N,M)
c          XGW(L,2,1)=XWSD(N,M+1)
c          XGW(L,2,2)=YWSD(N,M+1)
c          XGW(L,2,3)=ZWSD(N,M+1)
c          XGW(L,3,1)=XWSD(N+1,M+1)
c          XGW(L,3,2)=YWSD(N+1,M+1)
c          XGW(L,3,3)=ZWSD(N+1,M+1)
c          XGW(L,4,1)=XWSD(N+1,M)
c          XGW(L,4,2)=YWSD(N+1,M)
c          XGW(L,4,3)=ZWSD(N+1,M)
c 710     CONTINUE
c 720   CONTINUE
c
c      CALL GEO3DW(NPWAKSD,XGW,CHRLEWS,IER)
c
c      IF(IER.EQ.0) THEN
c         WRITE(*,*) 'UNACCEPTABLE PANEL IN SCWKINF'
c         STOP
c      END IF
C
C
c      DO 760 M = MDUCT, 1, -1
c       DO 750 N = 1, NWSUB
c        J = NWSUB*(M-1)+N
c        DO K=1,4
c          XV(K)=XVPW(J,K)
c          YV(K)=YVPW(J,K)
c          SIDE(K)=SIDW(J,K)
c        ENDDO
c
c        DO K=1,15
c          S(K)=SSW(J,K)
c        ENDDO
c
c        XM1(1)=XWSD(N,M)
c        XM1(2)=YWSD(N,M) 
c        XM1(3)=ZWSD(N,M)
c        XM2(1)=XWSD(N,M+1)
c        XM2(2)=YWSD(N,M+1)
c        XM2(3)=ZWSD(N,M+1)
c        XM3(1)=XWSD(N+1,M+1)
c        XM3(2)=YWSD(N+1,M+1)
c        XM3(3)=ZWSD(N+1,M+1)
c        XM4(1)=XWSD(N+1,M)
c        XM4(2)=YWSD(N+1,M)
c        XM4(3)=ZWSD(N+1,M)   
c        IMR0=0 
c
c        ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))+ABS(XM2(3)-XM1(3))
c        ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))+ABS(XM3(3)-XM2(3))
c        ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))+ABS(XM4(3)-XM3(3))
c        ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))+ABS(XM1(3)-XM4(3))
c        ERRMIN=AMIN1(ER1,ER2,ER3,ER4) 
c        IF(ERRMIN.LE.1.0E-6) THEN
c          IMR0=1
c        END IF
c
c        DO 749 KK=1,NBLADE
c
c        JSUB = N/NWSUB1+1
c        IF(MOD(N,NWSUB1).EQ.0) JSUB=JSUB-1 
c        JMAC=INDEXWD(JSUB,M)
c        JMAC1=INDEXWD(JSUB+1,M)
c        
c        SLOPE=(TEMPD4(JMAC1)-TEMPD4(JMAC)) / (XDW(JSUB+1,M)-XDW(JSUB,M))
c        TEMPSUB(J,KK) = SLOPE*(XCTPDWs(J,1,1)-XDW(JSUB,M))+TEMPD4(JMAC)
c                      
c        IF(KK.GT.1)THEN        
c          IREC = NTPOS(KK)         
c          NWBIG = NDWK*MDUCT  
c          CALL READ2(49,IREC,TEMP5,NWBIG) 
c                 
c          SLOPE=(TEMP5(JMAC1)-TEMP5(JMAC)) / (XDW(JSUB+1,M)-XDW(JSUB,M))
c          TEMPSUB(J,KK) = SLOPE*(XCTPDWs(J,1,1)-XDW(JSUB,M))+TEMP5(JMAC)          
c         ENDIF      
c
c                 
cC------- LOOP 730: FOR WAKE SUBPANELS    (Hong, 11/02/2006)     
c        DO 730 IN = 1,NWSUB
c          L = NWSUB*(JM-1)+IN
c
C.............Transfer control points to local coordinate...............
c          XLOC=ZERO
c          YLOC=ZERO
c          ZLOC=ZERO
c          DO K=1,3
c             XLOC=XLOC+(XCTPWs(L,K,KK)-XCTW(J,K))*DIRW(J,1,K)
c             YLOC=YLOC+(XCTPWs(L,K,KK)-XCTW(J,K))*DIRW(J,2,K)
c             ZLOC=ZLOC+(XCTPWs(L,K,KK)-XCTW(J,K))*DIRW(J,3,K)
c          ENDDO  
c
c          IMR=IMR0
c          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,FSY,
c     *         FDX,FDY,FDZ,0,IMR)
c
c          IF(IMR.EQ.2) THEN
c            DO IXYZ=1,3
c               XMC(IXYZ)=XCTPWs(L,IXYZ,KK)
c            ENDDO
c               CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
c          END IF
c      
c          IF(ABS(FD).GT.6.28) THEN                       
c            FD=0.0
c          ENDIF
c
c         Calculate Total potential
c 
c          PPOTWs(L) = PPOTWs(L)- TEMPSUB(J,KK)*FD/PI/4.
c 
c 730    CONTINUE
c
cC------- LOOP 740: FOR WAKE PANELS 
c        DO 740 IN = 1,NDWK
c          L = INDEXWD(IN,JM)
c
cC.............Transfer control points to local coordinate...............
c          XLOC=ZERO
c          YLOC=ZERO
c          ZLOC=ZERO
c          DO K=1,3
c             XLOC=XLOC+(XCTPW(L,K,KK)-XCTW(J,K))*DIRW(J,1,K)
c             YLOC=YLOC+(XCTPW(L,K,KK)-XCTW(J,K))*DIRW(J,2,K)
c             ZLOC=ZLOC+(XCTPW(L,K,KK)-XCTW(J,K))*DIRW(J,3,K)
c          ENDDO  
c
c          IMR=IMR0
c          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(J),FS,FD,FSX,FSY,
c     *         FDX,FDY,FDZ,0,IMR)
c
c          IF(IMR.EQ.2) THEN
c            DO IXYZ=1,3
c               XMC(IXYZ)=XCTPW(L,IXYZ,KK)
c            ENDDO
c               CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
c          END IF
c                  
c          IF(ABS(FD).GT.6.28) THEN                       
c            FD=0.0
c          END IF
c
cc         Calculate Total potential          
c 
c          PPOTW(L) = PPOTW(L)- TEMPSUB(J,KK)*FD/PI/4. 
c
c 740    CONTINUE 
c 749    CONTINUE    
c
c 750   CONTINUE
c 760  CONTINUE
C
C                  MACRO WAKE PANELS
C
      NWBIGD = NDWK*MDUCT
      DO 820 N=1, NDWK
        DO 810 M=1,MDUCT
          L=INDEXWD(N,M)
          XGWD(L,1,1)=XDW(N,M+1)
          XGWD(L,1,2)=YDW(N,M+1)
          XGWD(L,1,3)=ZDW(N,M+1)
          XGWD(L,2,1)=XDW(N,M)
          XGWD(L,2,2)=YDW(N,M)
          XGWD(L,2,3)=ZDW(N,M)
          XGWD(L,3,1)=XDW(N+1,M)
          XGWD(L,3,2)=YDW(N+1,M)
          XGWD(L,3,3)=ZDW(N+1,M)
          XGWD(L,4,1)=XDW(N+1,M+1)
          XGWD(L,4,2)=YDW(N+1,M+1)
          XGWD(L,4,3)=ZDW(N+1,M+1) 
 810     CONTINUE
 820   CONTINUE
        
      CALL GEO3DWD(NWBIGD,XGWD,CHRLEWSD,IER)

      IF(IER.EQ.0) THEN
         WRITE(*,*) 'UNACCEPTABLE PANEL IN SCWKINF'
         STOP
      END IF
C
C
      DO 860 M = 1, MDUCT 
       DO 850 N = 1, NDWK             
        J = INDEXWD(N,M)
        DO K=1,4
          XV(K)=XVPWD(J,K)
          YV(K)=YVPWD(J,K)
          SIDE(K)=SIDWD(J,K)
        ENDDO

        DO K=1,15
          S(K)=SSWD(J,K)
        ENDDO

        XM1(1)=XDW(N,M+1)
        XM1(2)=YDW(N,M+1) 
        XM1(3)=ZDW(N,M+1)
        XM2(1)=XDW(N,M)
        XM2(2)=YDW(N,M)
        XM2(3)=ZDW(N,M)
        XM3(1)=XDW(N+1,M)
        XM3(2)=YDW(N+1,M)
        XM3(3)=ZDW(N+1,M)
        XM4(1)=XDW(N+1,M+1)
        XM4(2)=YDW(N+1,M+1)
        XM4(3)=ZDW(N+1,M+1)   
        IMR0=0 

        DO 849 KK=1,NBLADE
      
        IF(KK.GT.1)THEN        
          IREC = NTPOS(KK)     
          CALL READ2(49,IREC,TEMP5,NWBIGD) 
         ENDIF
        
C------- LOOP 830: FOR WAKE SUBPANELS    (Hong, 11/02/2006)     
        DO 830 IN = 1,NTRA
          L = NTRA*(MR-JM)+IN 

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPWs(L,K,KK)-XCTWD(J,K))*DIRWD(J,1,K)
             YLOC=YLOC+(XCTPWs(L,K,KK)-XCTWD(J,K))*DIRWD(J,2,K)
             ZLOC=ZLOC+(XCTPWs(L,K,KK)-XCTWD(J,K))*DIRWD(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWSD(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPWs(L,IXYZ,KK)
            ENDDO
               CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF
                  
          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          END IF

c         Calculate Total potential
         IF(KK.EQ.1) THEN
           PPOTWs(L) = PPOTWs(L)- TEMPD4(J)*FD/PI/4.
         ELSE 
           PPOTWs(L) = PPOTWs(L)- TEMP5(J)*FD/PI/4.       
         ENDIF
 830    CONTINUE

C------- LOOP 840: FOR WAKE PANELS 
        DO 840 IN = 1,NWMIN
          L = INDEXW(IN,JM)

C.............Transfer control points to local coordinate...............
          XLOC=ZERO
          YLOC=ZERO
          ZLOC=ZERO
          DO K=1,3
             XLOC=XLOC+(XCTPW(L,K,KK)-XCTWD(J,K))*DIRWD(J,1,K)
             YLOC=YLOC+(XCTPW(L,K,KK)-XCTWD(J,K))*DIRWD(J,2,K)
             ZLOC=ZLOC+(XCTPW(L,K,KK)-XCTWD(J,K))*DIRWD(J,3,K)
          ENDDO  

          IMR=IMR0
          CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWSD(J),FS,FD,FSX,FSY,
     *         FDX,FDY,FDZ,0,IMR)

          IF(IMR.EQ.2) THEN
            DO IXYZ=1,3
               XMC(IXYZ)=XCTPW(L,IXYZ,KK)
            ENDDO
               CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
          END IF

          IF(ABS(FD).GT.6.28) THEN                       
            FD=0.0
          ENDIF
  

c         Calculate Total potential          
         IF(KK.EQ.1) THEN
           PPOTW(L) = PPOTW(L)- TEMPD4(J)*FD/PI/4.
         ELSE 
           PPOTW(L) = PPOTW(L)- TEMP5(J)*FD/PI/4.       
         ENDIF

 840    CONTINUE 
 849    CONTINUE    

 850   CONTINUE
 860  CONTINUE

      ENDIF !IDUCT, IDOPT 

CVV
      DEALLOCATE(TEMPSUB)
CVV

      RETURN
      END  

