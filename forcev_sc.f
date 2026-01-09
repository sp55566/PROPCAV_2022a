      SUBROUTINE FORCEV_SC
C***********************************************************************
C     FORCEV: FORCE calculations (Viscous corrections)
C      --- Calculate KT,KQ, and ETA
C
C     FX: 1,2,3: x,y,z components of force on the propeller
C         4,5,6: moments about x,y,z, axese
C         7,8  : KT and KQ
C
C     Date           Revisions or Comments
C     ------         -----------------------
C     JY100901       Copied from forcev_sc.f.  Modified for ISC=1.  
C***********************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      REAL MXLOC1,MYLOC1,MZLOC1
      DIMENSION FXLOC1(MBZ),FYLOC1(MBZ),FZLOC1(MBZ),MXLOC1(MBZ),
     *          MYLOC1(MBZ),MZLOC1(MBZ),RSQ1(MBZ)
      DIMENSION FXD(MBZ),FYD(MBZ),FZD(MBZ),FMXD(MBZ),FMYD(MBZ),FMZD(MBZ)
      DIMENSION BVFX(MBZ),BVFZ(MBZ)

      CALL CLEAR(FBXP,6)
      CALL CLEAR(FHXP,6)
      CALL CLEAR(FBXV,6)
      CALL CLEAR(FHXV,6)

C-----------------------------------------------------------------------
C     Integrate the pressure on the blade
C-----------------------------------------------------------------------
      BFX=0.
      BFY=0.
      BFZ=0.
      BMX=0.
      BMY=0.
      BMZ=0.

      DO 20 M=1,MR
         FX=ZERO
         FY=ZERO
         FZ=ZERO
         TMX=ZERO
         TMY=ZERO
         TMZ=ZERO

         N1=N0(2)
         N2=N0(1)-1

         DO 10 N=N1,N2
            L=INDEXB(N,M)
            CPA=CPBN(N,M)*SS(L,1)

C-----------------------------------------------------------------------
C           Leading edge suction force correction based on Polhamus' 
C             model
C-----------------------------------------------------------------------
            N1=N-NH
            XNML=VEL(L,1)
            YNML=VEL(L,2)
            ZNML=VEL(L,3)
            XRNML=VEL(L,4)
            YRNML=VEL(L,5)
            ZRNML=VEL(L,6)

            FX=FX+CPA*XNML
            FY=FY+CPA*YNML
            FZ=FZ+CPA*ZNML
            TMX=TMX+CPA*XRNML
            TMY=TMY+CPA*YRNML
            TMZ=TMZ+CPA*ZRNML

10       CONTINUE

C.......Add base pressure for ISC=1. (JY071001)
         CPSIGMA=-SIGMA/ADVCO/ADVCO
         CPA=CPSIGMA*SSTE(M)
         XNML=VELTE(M,1)
         YNML=VELTE(M,2)
         ZNML=VELTE(M,3)
         XRNML=VELTE(M,4)
         YRNML=VELTE(M,5)
         ZRNML=VELTE(M,6) 
         FX=FX+CPA*XNML
         FY=FY+CPA*YNML
         FZ=FZ+CPA*ZNML
         TMX=TMX+CPA*XRNML
         TMY=TMY+CPA*YRNML
         TMZ=TMZ+CPA*ZRNML

         BKFX(M)=FX
         BKFY(M)=FY
         BKFZ(M)=FZ
         BKMX(M)=TMX
         BKMY(M)=TMY
         BKMZ(M)=TMZ

20    CONTINUE

C.....Smooth the force -- using MRTIP as the edge for smoothing)
C.....Set ISW=1 to turn of smoothing of forces. (JY080400)
      ISW=0
      IF(ISW.EQ.1) GO TO 1000

      IF(MRTIP.NE.MR) THEN
         MR2P=MRTIP+1
         DO 22 M=1,MRTIP
            DR=RZ(M+1)-RZ(M)
            FXD(M)=BKFX(M)/DR
            FYD(M)=BKFY(M)/DR
            FZD(M)=BKFZ(M)/DR
            FMXD(M)=BKMX(M)/DR
            FMYD(M)=BKMY(M)/DR
            FMZD(M)=BKMZ(M)/DR

            RSQ1(M)=RZPSQ(M)
            FXLOC1(M)=FXD(M)
            FYLOC1(M)=FYD(M)
            FZLOC1(M)=FZD(M)
            MXLOC1(M)=FMXD(M)
            MYLOC1(M)=FMYD(M)
            MZLOC1(M)=FMZD(M)
 22      CONTINUE
         RSQ1(MR2P)=ONE
         FXD(MR2P)=ZERO
         FYD(MR2P)=ZERO
         FZD(MR2P)=ZERO
         FMXD(MR2P)=ZERO
         FMYD(MR2P)=ZERO
         FMZD(MR2P)=ZERO

         DO M=MR2P,MR
            CALL BQUAD(RSQ1,FXD,RZPSQ(M),IM,
     *           COA,COB,COC,YY,MR2P)
            FXLOC1(M)=YY
            CALL BQUAD(RSQ1,FYD,RZPSQ(M),IM,
     *           COA,COB,COC,YY,MR2P)
            FYLOC1(M)=YY
            CALL BQUAD(RSQ1,FZD,RZP(M),IM,
     *           COA,COB,COC,YY,MR2P)
            FZLOC1(M)=YY
            CALL BQUAD(RSQ1,FMXD,RZPSQ(M),IM,
     *           COA,COB,COC,YY,MR2P)
            MXLOC1(M)=YY
            CALL BQUAD(RSQ1,FMYD,RZPSQ(M),IM,
     *           COA,COB,COC,YY,MR2P)
            MYLOC1(M)=YY
            CALL BQUAD(RSQ1,FMZD,RZPSQ(M),IM,
     *           COA,COB,COC,YY,MR2P)
            MZLOC1(M)=YY
         END DO

         DO 26 M=1,MR
            DR=RZ(M+1)-RZ(M)
            BKFX(M)=FXLOC1(M)*DR
            BKFY(M)=FYLOC1(M)*DR
            BKFZ(M)=FZLOC1(M)*DR
            BKMX(M)=MXLOC1(M)*DR
            BKMY(M)=MYLOC1(M)*DR
            BKMZ(M)=MZLOC1(M)*DR
 26      CONTINUE

      END IF

 1000 CONTINUE

      DO 28 M=1,MR
         BFX=BFX+BKFX(M)
         BFY=BFY+BKFY(M)
         BFZ=BFZ+BKFZ(M)
         BMX=BMX+BKMX(M)
         BMY=BMY+BKMY(M)
         BMZ=BMZ+BKMZ(M)
 28   CONTINUE

C-----------------------------------------------------------------------
C     Viscous friction force on the blade
C-----------------------------------------------------------------------
      VBFX=ZERO
      VBFY=ZERO
      VBFZ=ZERO
      VBMX=ZERO
      VBMY=ZERO
      VBMZ=ZERO

      CDB=XCDF
      DO 40 M=1,MR
         FX=ZERO
         FY=ZERO
         FZ=ZERO
         TMX=ZERO
         TMY=ZERO
         TMZ=ZERO

C-----------------------------------------------------------------------
C     Added new variables BVFX & BVFZ so we can calculate lift and
C     Drag coefficients.                                       JY110700
C-----------------------------------------------------------------------
         BVFX(M)=ZERO
         BVFZ(M)=ZERO

C-----------------------------------------------------------------------
C     Viscous correction should be applied to the entire blade for 
C     fully wetted flow.  For cavitating flow, viscous forces should 
C     be integrated over wetted panels only.                    JY052599
C-----------------------------------------------------------------------

C.......Wetted panels on the face side of the blade
         IF(IWET.EQ.1.OR.NOCAV(M,IDXREV,2).EQ.1.OR.
     *        JCV(M,2).EQ.0) THEN
            CALL VISFR(M,N0(2),NH,CDB,FX,FY,FZ,TMX,TMY,TMZ)
         ELSE
            LCAVF=M0(M,2)-JCV(M,2)-NSPP(M,2)-1
            CALL VISFR(M,N0(2),LCAVF,CDB,FX,FY,FZ,TMX,TMY,TMZ)
            M0MF=M0(M,2)
            CALL VISFR(M,M0MF,NH,CDB,FX,FY,FZ,TMX,TMY,TMZ)
            IF(NSPP(M,2).EQ.1) THEN
               FAC=CDB*FRP(M,2)
               ISPF=LCAVF+1
               CALL VISFR(M,ISPF,ISPF,FAC,FX,FY,FZ,TMX,TMY,TMZ)
            END IF
         END IF

C.......Wetted panels on the back side of the blade
         IF(IWET.EQ.1.OR.NOCAV(M,IDXREV,1).EQ.1.OR.
     *        JCV(M,1).EQ.0) THEN
            CALL VISFR(M,NHP,N0(1)-1,CDB,FX,FY,FZ,TMX,TMY,TMZ)
         ELSE
            M0MBM=M0(M,1)-1
            CALL VISFR(M,NHP,M0MBM,CDB,FX,FY,FZ,TMX,TMY,TMZ)
            LCAVB=M0(M,1)+JCV(M,1)+NSPP(M,1)
            CALL VISFR(M,LCAVB,N0(1)-1,CDB,FX,FY,FZ,TMX,TMY,TMZ)
            IF(NSPP(M,1).EQ.1) THEN
               FAC=CDB*FRP(M,1)
               ISPB=LCAVB-1
               CALL VISFR(M,ISPB,ISPB,FAC,FX,FY,FZ,TMX,TMY,TMZ)
            END IF
         END IF
C-----------------------------------------------------------------------

         VBFX=VBFX+FX
         VBFY=VBFY+FY
         VBFZ=VBFZ+FZ
         VBMX=VBMX+TMX
         VBMY=VBMY+TMY
         VBMZ=VBMZ+TMZ

C-----------------------------------------------------------------------
C     Added new variables BVFX & BVFZ so we can calculate lift and
C     Drag coefficients.                                       JY110700
C-----------------------------------------------------------------------
         BVFX(M)=FX
         BVFZ(M)=FZ

40    CONTINUE

C-----------------------------------------------------------------------
C     Integrate the pressure on the hub
C-----------------------------------------------------------------------
      HFX=ZERO
      HFY=ZERO
      HFZ=ZERO
      HMX=ZERO
      HMY=ZERO
      HMZ=ZERO
      DO 60 N=1,NHBX
         FX=ZERO
         FY=ZERO
         FZ=ZERO
         TMX=ZERO
         TMY=ZERO
         TMZ=ZERO
         DO 50 M=1,MHBT
            L=INDEXH(N,M)
            CPA=CPHN(N,M)*SS(L,1)
            FX=FX+CPA*VEL(L,1)
            FY=FY+CPA*VEL(L,2)
            FZ=FZ+CPA*VEL(L,3)
            TMX=TMX+CPA*VEL(L,4)
            TMY=TMY+CPA*VEL(L,5)
            TMZ=TMZ+CPA*VEL(L,6)
50       CONTINUE
         HKFX(N)=FX
         HKFY(N)=FY
         HKFZ(N)=FZ
         HKMX(N)=TMX
         HKMY(N)=TMY
         HKMZ(N)=TMZ
         HFX=HFX+HKFX(N)
         HFY=HFY+HKFY(N)
         HFZ=HFZ+HKFZ(N)
         HMX=HMX+HKMX(N)
         HMY=HMY+HKMY(N)
         HMZ=HMZ+HKMZ(N)
60    CONTINUE

C-----------------------------------------------------------------------
C     Integrate the hub friction force
C-----------------------------------------------------------------------
      VHFX=ZERO
      VHFY=ZERO
      VHFZ=ZERO
      VHMX=ZERO
      VHMY=ZERO
      VHMZ=ZERO
      DO 80 N=1,NHBX
         FX=ZERO
         FY=ZERO
         FZ=ZERO
         TMX=ZERO
         TMY=ZERO
         TMZ=ZERO
         DO 70 M=1,MHBT
            L=INDEXH(N,M)
            CFA=XCDF*VTOTS(L)*SS(L,1)
C-----------------------------------------------------------------------
C     Modified vector of forces to align with that of velocity instead
C     of hub geometry.                                          JY031898
C-----------------------------------------------------------------------
            TXH=UXHTOT(N,M)/SQRT(VTOTS(L))
            TYH=UYHTOT(N,M)/SQRT(VTOTS(L))
            TZH=UZHTOT(N,M)/SQRT(VTOTS(L))
            FX=FX+CFA*TXH
            FY=FY+CFA*TYH
            FZ=FZ+CFA*TZH
            TMX=TMX+CFA*(XCT(L,2)*TZH-XCT(L,3)*TYH)
            TMY=TMY+CFA*(XCT(L,3)*TXH-XCT(L,1)*TZH)
            TMZ=TMX+CFA*(XCT(L,1)*TYH-XCT(L,2)*TXH)
C------End of modification (JY031898)-----------------------------------
70       CONTINUE
         VHFX=VHFX+FX
         VHFY=VHFY+FY
         VHFZ=VHFZ+FZ
         VHMX=VHMX+TMX
         VHMY=VHMY+TMY
         VHMZ=VHMZ+TMZ
80    CONTINUE

C----------------------------------------------------------------------
C     I'm multiplying the forces by 1/8*J^2 and moments by 1/16*J^2 so 
C     the non-dimensionalization is the same as HPUF3A.  This 
C     modification is not needed for TKT and TKQ because it has 
C     already been done above.
C     F/(rho*n^2*D^4)=CPB*SS(*J^2/4)(*1/2) because 
C     Cp(HPUF)=1/2Cp(PROPCAV)
C     Originally, FXV(1) was = bfx+hfx+vbfx+vhfx               JY012198
C----------------------------------------------------------------------


      FBXP(1) = BFX * 0.125*ADVCO**2
      FBXP(2) = BFY * 0.125*ADVCO**2
      FBXP(3) = BFZ * 0.125*ADVCO**2
      FBXP(4) = BMX * 0.0625*ADVCO**2
      FBXP(5) = BMY * 0.0625*ADVCO**2
      FBXP(6) = BMZ * 0.0625*ADVCO**2

      IF(IHUB .NE. 0) THEN
         FHXP(1) = HFX * 0.125*ADVCO**2
         FHXP(2) = HFY * 0.125*ADVCO**2
         FHXP(3) = HFZ * 0.125*ADVCO**2
         FHXP(4) = HMX * 0.0625*ADVCO**2
         FHXP(5) = HMY * 0.0625*ADVCO**2
         FHXP(6) = HMZ * 0.0625*ADVCO**2
      ENDIF


C----------------------------------------------------------------------
C     Eight potential force components                        
C----------------------------------------------------------------------


      FBXV(1) = VBFX * 0.125*ADVCO**2
      FBXV(2) = VBFY * 0.125*ADVCO**2
      FBXV(3) = VBFZ * 0.125*ADVCO**2
      FBXV(4) = VBMX * 0.0625*ADVCO**2
      FBXV(5) = VBMY * 0.0625*ADVCO**2
      FBXV(6) = VBMZ * 0.0625*ADVCO**2

      IF(IHUB .NE. 0) THEN
         FHXV(1) = VHFX * 0.125*ADVCO**2
         FHXV(2) = VHFY * 0.125*ADVCO**2
         FHXV(3) = VHFZ * 0.125*ADVCO**2
         FHXV(4) = VHMX * 0.0625*ADVCO**2
         FHXV(5) = VHMY * 0.0625*ADVCO**2
         FHXV(6) = VHMZ * 0.0625*ADVCO**2
      ENDIF

      DO I = 1 , 6
C -- Blade total force (Potential + Viscous)
         FBXPV(I) = FBXP(I) + FBXV(I)
C -- Hub total force (Potential + Viscous)
         FHXPV(I) = FHXP(I) + FHXV(I)
C ==========Total Potential
         FXP(I)=FBXP(I)+FHXP(I)
C ========= Total Viscous
         FXV(I)=FBXV(I)+FHXV(I)
C ========== Total Forces
         FTOTAL(I) = FXP(I) + FXV(I)
      ENDDO

C----------------------------------------------------------------------
C     Calculate lift and drag coefficient for hydrofoils.      JY110300
C----------------------------------------------------------------------
      IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
         DO M=1,MR
            CHDM=CHORD(M)+CHORD(M+1)
            DELR=RZ(M+1)-RZ(M)
            IF(DELR .EQ. 0.) DELR=1.E-6
            SUMFX=(BKFX(M)+BVFX(M))/CHDM/DELR
            SUMFZ=(BKFZ(M)+BVFZ(M))/CHDM/DELR
         END DO
         CZ=NBLADE*(SUMFZ)
         CX=NBLADE*(SUMFX)
         COSA=COS(ALPHA*RAD)
         SINA=SIN(ALPHA*RAD)
         CLIFT=CZ*COSA-CX*SINA
         CDRAG=CX*COSA+CZ*SINA
      END IF

      RETURN
C))))))))))))))))))))) End of subroutine FORCE (((((((((((((((((((((((((
      END

