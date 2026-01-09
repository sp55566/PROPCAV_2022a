      SUBROUTINE VISPOT(JM)
C **********************************************************************.
C     Calculate the total potential along each WAKE STRIP.             *
C                        Wetted   Case                                 *
C                                                                      *
C            By Hong Sun                April, 2006                    *
C **********************************************************************
 
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
! XM YU 06/2012
      npwakev=nwv*mr
      npwaks=nwsub*mr
      do i=1,npwakev
         ppotw(i)=0.0
      enddo
      do i=1,npwaks
         ppotws(i)=0.0
      enddo 
! XM YU 06/2012
      
      
           
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
C/e S.N.KIM | Aug. 2018

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
          CALL READ2(47,IREC,DPDN,NPANEL)         
         ENDIF
               
C------- LOOP 130: FOR WAKE SUBPANELS    (Hong Sun, 04/07/2005)     
        DO 130 IN = 1,NWSUB
          L = NWSUB*(JM-1)+IN

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
          PPOTWs(L) = PPOTWs(L)+ (-POT(J)*FD +DPDN(J)*FS)/PI/4.
       
 130    CONTINUE

C------- LOOP 140: FOR WAKE PANELS 
C        DO 140 IN = NSUB+1,NWMIN          Hong check May, 2008
        DO 140 IN = 1,NWMIN
c XM YU 4/18/2012
          L = IDXWAK(IN,JM)
c           L=(mr-jm)*nwmin+IN
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
          PPOTW(L) = PPOTW(L)+ (-POT(J)*FD +DPDN(J)*FS)/PI/4.

 140    CONTINUE 

 149   CONTINUE       
 150   CONTINUE
 160  CONTINUE

C     Compute influence coefficients due to the hub
C*********************************************************************
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
          CALL READ2(47,IREC,DPDN,NPANEL)           
        ENDIF      
               
         DO 230 IN = 1,NWSUB       
          L = NWSUB*(JM-1)+IN

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
          PPOTWs(L) = PPOTWs(L)+ (-POT(J)*FD +DPDN(J)*FS)/PI/4.

 230     CONTINUE 


C--------LOOP 240: FOR WAKE PANELS   
C         DO 240 IN = NSUB+1,NWMIN   Hong check May, 2008
         DO 240 IN = 1,NWMIN        
          L = IDXWAK(IN,JM)

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
          PPOTW(L) = PPOTW(L)+ (-POT(J)*FD +DPDN(J)*FS)/PI/4.

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
          CALL READ2(47,IREC,DPDN,NPANEL) 
        ENDIF      
              
         DO 330 IN = 1,NWSUB       
          L = NWSUB*(JM-1)+IN

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
          PPOTWs(L) = PPOTWs(L)+ (-POT(J)*FD +DPDN(J)*FS)/PI/4.

 330     CONTINUE 


C--------LOOP 340: FOR WAKE PANELS   
C         DO 340 IN = NSUB+1,NWMIN     Hong check May, 2008
         DO 340 IN = 1,NWMIN        
          L = IDXWAK(IN,JM)

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
          PPOTW(L) = PPOTW(L)+ (-POT(J)*FD +DPDN(J)*FS)/PI/4.

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

C--------LOOP 430: FOR WAKE SUBPANELS  (Hong Sun, 06/04/2006) 
        DO 449 KK = 1, NBLADE
      
        IF(ISTEADY.GT.0)THEN        
          IREC = NTPOS(KK)     
          CALL READ2(45,IREC,POT,NPANEL)      
          CALL READ2(47,IREC,DPDN,NPANEL)     
        ENDIF      
              
         DO 430 IN = 1,NWSUB       
          L = NWSUB*(JM-1)+IN

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
          PPOTWs(L) = PPOTWs(L)+ (-POT(J)*FD +DPDN(J)*FS)/PI/4.


 430     CONTINUE 


C--------LOOP 440: FOR WAKE PANELS   
C         DO 440 IN = NSUB+1,NWMIN    Hong check May, 2008
         DO 440 IN = 1,NWMIN 
          L = IDXWAK(IN,JM)

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
          PPOTW(L) = PPOTW(L)+ (-POT(J)*FD +DPDN(J)*FS)/PI/4.


 440     CONTINUE 

 449    CONTINUE      
 450    CONTINUE
 460   CONTINUE
      ENDIF


C*********************************************************************
C                Wake Dipoles               Hong 11/02/2006
C*********************************************************************
C                WAKE SUB PANELS
C
      
      NPWAKS=NWSUB*MR
      DO 520 N=1,NWSUB
        DO 510 M=1,MR
          L=NWSUB*(M-1)+N
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
       DO 550 N = 1, NWSUB
        J = NWSUB*(M-1)+N
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
      
        DO 549 KK = 1, 1
      
          IF(KK.EQ.1)THEN        
            JSUB = N/NWSUB1+1
          IF(MOD(N,NWSUB1).EQ.0) JSUB=JSUB-1 

          JMAC=INDEXW2(JSUB,M)
          JMAC1=INDEXW2(JSUB+1,M)
          
          SLOPE=(TEMP4(JMAC1)-TEMP4(JMAC))/(XW(JSUB+1,M)-XW(JSUB,M))
c          TEMPSUB(J,KK) = SLOPE*(XCTPWs(J,1,1)-XW(JSUB,M))+TEMP4(JMAC)
          TEMPSUB(J,KK)=SLOPE*(0.5*(XWS(N,M)+XWS(N+1,M))-XW(JSUB,M))
     &                +TEMP4(JMAC)      
        
c          IF(JM.EQ.8.AND.M.EQ.8.AND.KK.EQ.1)THEN
c              WRITE(139,*)J,XCTPWs(J,1,1),TEMPSUB(J,KK)         
c          ENDIF
        ENDIF
                 
C------- LOOP 530: FOR WAKE SUBPANELS    (Hong, 11/02/2006)     
        DO 530 IN = 1,NWSUB
          L = NWSUB*(JM-1)+IN

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

           IF(L.EQ.J.AND. KK.EQ.1) THEN
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
           
      IF(KK.EQ.1) PPOTWs(L)=PPOTWs(L)-TEMPSUB(J,KK)*FD/PI/4.        
        
 530    CONTINUE
 549    CONTINUE  
 
 550   CONTINUE
 560  CONTINUE
C
C                  MACRO WAKE PANELS
C
      NWBIG = NWMIN*MR
      DO 620 N=1, NWMIN
        DO 610 M=1,MR
c XM YU
          L=IDXWAK(N,M)
c          l=(mr-m)*nwmin+n
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
c          j=(mr-m)*nwmin+n
        J = IDXWAK(N,M)
      L1 = INDEXW2(N,M)
      
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
       
c       IF(JM.EQ.9.AND.M.EQ.9)THEN
c           WRITE(139,*) XCTPW(J,1,1),TEMP4(L1)                   !!!!!!!!
c       ENDIF

        DO 649 KK=1,NBLADE
   
        IF(KK.GT.1)THEN        
          IREC = NTPOS(KK)     
          CALL READ2(46,IREC,TEMP5,NWBIG) 
         ENDIF
        
C------- LOOP 630: FOR WAKE SUBPANELS    (Hong, 11/02/2006)     
        DO 630 IN = 1,NWSUB
          L=NWSUB*(JM-1)+IN

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
          L = IDXWAK(IN,JM)
c           l=nwmin*(mr-jm)+in
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
C             Include Duct Wake Dipoles            Hong 11/02/2006
C*********************************************************************
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
        DO 830 IN = 1,NWSUB
          L=NWSUB*(JM-1)+IN

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
          L = IDXWAK(IN,JM) 

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

      DEALLOCATE(TEMPSUB)
      RETURN
      END

