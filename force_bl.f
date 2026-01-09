      SUBROUTINE FORCE_BL
C***********************************************************************
C     FORCE_BL: FORCE CALCULATIONS (INCLUDE BOUNDARY LAYER EFFECTS))
C     --- CALCULATE KT,KQ, AND ETA
C     
C     FX: 1,2,3: X,Y,Z COMPONENTS OF FORCE ON THE PROPELLER
C     4,5,6: MOMENTS ABOUT X,Y,Z, AXESE
C     7,8  : KT AND KQ
C     
C     DATE           REVISIONS OR COMMENTS
C     ------         -----------------------
C     06.2006        BY HONG SUN   
C***********************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
C     INCLUDE VISCOUS DRAG COEFFICIENT BY HONG SUN
      INCLUDE 'PUFBL.INC'

      REAL MXLOC1,MYLOC1,MZLOC1
      DIMENSION FXLOC1(MBZ),FYLOC1(MBZ),FZLOC1(MBZ),MXLOC1(MBZ),
     *  MYLOC1(MBZ),MZLOC1(MBZ),RSQ1(MBZ)
      DIMENSION FXD(MBZ),FYD(MBZ),FZD(MBZ),FMXD(MBZ),FMYD(MBZ),FMZD(MBZ)
      DIMENSION BVFX(MBZ),BVFZ(MBZ)
      DIMENSION THKFX(NHMX),THKFY(NHMX),THKFZ(NHMX),
     %  THKMX(NHMX),THKMY(NHMX),THKMZ(NHMX)
      DIMENSION DKFX(NDMAX),DKFY(NDMAX),DKFZ(NDMAX),
     %  DKMX(NDMAX),DKMY(NDMAX),DKMZ(NDMAX)
      DIMENSION CHKFX(NCAVM),CHKFY(NCAVM),CHKFZ(NCAVM),
     %  CHKMX(NCAVM),CHKMY(NCAVM),CHKMZ(NCAVM)
      
      DATA COFLS/1.0/

      CALL CLEAR(FBXP,6)
      CALL CLEAR(FHXP,6)
      CALL CLEAR(FDXP,6)
      CALL CLEAR(FTXP,6)
      CALL CLEAR(FBXV,6)
      CALL CLEAR(FHXV,6)
      CALL CLEAR(FDXV,6)

C-----------------------------------------------------------------------
C     INTEGRATE THE PRESSURE ON THE BLADE
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

        N1=1
        N2=NC

        DO 10 N=N1,N2
          L=INDEXB(N,M)

C-----------------------------------------------------------------------
C     IN THE FORCE CALCULATION, LIMIT THE PRESSURE TO SIGMA IN THE
C     CAVITATING ANALYSIS.  
C-----------------------------------------------------------------------
          DUM=-CPBL(N,M)*ADVCO*ADVCO
          IF(IWET.EQ.0.AND.DUM.GT.SIGMA) THEN
            CPA=-SIGMA/ADVCO/ADVCO*SS(L,1)
          ELSE
            CPA=CPBL(N,M)*SS(L,1)
          END IF

C-----------------------------------------------------------------------
C     LEADING EDGE SUCTION FORCE CORRECTION BASED ON POLHAMUS' 
C     MODEL
C-----------------------------------------------------------------------
          N1=N-NH

          IF(N.GT.NH.AND.VEL(L,3).GT.0.0.AND.ICON.NE.8) THEN
            XNML=COFLS*XCON(N1,M)+(1.0-COFLS)*VEL(L,1)
            YNML=COFLS*YCON(N1,M)+(1.0-COFLS)*VEL(L,2)
            ZNML=COFLS*ZCON(N1,M)+(1.0-COFLS)*VEL(L,3)
            XLEN=SQRT( XNML**2+YNML**2+ZNML**2)
            XNML=XNML/XLEN
            YNML=YNML/XLEN
            ZNML=ZNML/XLEN
            XRNML=XCTP(L,2,1)*ZNML-XCTP(L,3,1)*YNML
            YRNML=XCTP(L,3,1)*XNML-XCTP(L,1,1)*ZNML
            ZRNML=XCTP(L,1,1)*YNML-XCTP(L,2,1)*XNML
            XRLEN=SQRT( XRNML**2+YRNML**2+ZRNML**2)
          ELSE
            XNML=VEL(L,1)
            YNML=VEL(L,2)
            ZNML=VEL(L,3)
            XRNML=VEL(L,4)
            YRNML=VEL(L,5)
            ZRNML=VEL(L,6)
          END IF

          FX=FX+CPA*XNML
          FY=FY+CPA*YNML
          FZ=FZ+CPA*ZNML
          TMX=TMX+CPA*XRNML
          TMY=TMY+CPA*YRNML
          TMZ=TMZ+CPA*ZRNML

 10     CONTINUE

        BKFX(M)=FX
        BKFY(M)=FY
        BKFZ(M)=FZ
        BKMX(M)=TMX
        BKMY(M)=TMY
        BKMZ(M)=TMZ

 20   CONTINUE

C.....SMOOTH THE FORCE -- USING MRTIP AS THE EDGE FOR SMOOTHING)
C.....SET ISW=1 TO TURN OF SMOOTHING OF FORCES. (JY080400)
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
 22     CONTINUE
        RSQ1(MR2P)=ONE
        FXD(MR2P)=ZERO
        FYD(MR2P)=ZERO
        FZD(MR2P)=ZERO
        FMXD(MR2P)=ZERO
        FMYD(MR2P)=ZERO
        FMZD(MR2P)=ZERO

        DO M=MR2P,MR
          CALL BQUAD(RSQ1,FXD,RZPSQ(M),IM,
     *      COA,COB,COC,YY,MR2P)
          FXLOC1(M)=YY
          CALL BQUAD(RSQ1,FYD,RZPSQ(M),IM,
     *      COA,COB,COC,YY,MR2P)
          FYLOC1(M)=YY
          CALL BQUAD(RSQ1,FZD,RZP(M),IM,
     *      COA,COB,COC,YY,MR2P)
          FZLOC1(M)=YY
          CALL BQUAD(RSQ1,FMXD,RZPSQ(M),IM,
     *      COA,COB,COC,YY,MR2P)
          MXLOC1(M)=YY
          CALL BQUAD(RSQ1,FMYD,RZPSQ(M),IM,
     *      COA,COB,COC,YY,MR2P)
          MYLOC1(M)=YY
          CALL BQUAD(RSQ1,FMZD,RZPSQ(M),IM,
     *      COA,COB,COC,YY,MR2P)
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
 26     CONTINUE

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
C     VISCOUS FRICTION FORCE ON THE BLADE
C-----------------------------------------------------------------------
      VBFX=ZERO
      VBFY=ZERO
      VBFZ=ZERO
      VBMX=ZERO
      VBMY=ZERO
      VBMZ=ZERO

      CDB=XCDF
      DO 40 M=1,MR
C........INCLUDE SECTION DRAG COEFFICIENT  BY HONG SUN        
        IF(IVISC.EQ.1) THEN                    
          CDB=CDV(M)
        ENDIF
        FX=ZERO
        FY=ZERO
        FZ=ZERO
        TMX=ZERO
        TMY=ZERO
        TMZ=ZERO
        
C-----------------------------------------------------------------------
C     ADDED NEW VARIABLES BVFX & BVFZ SO WE CAN CALCULATE LIFT AND
C     DRAG COEFFICIENTS.                                       JY110700
C-----------------------------------------------------------------------
        BVFX(M)=ZERO
        BVFZ(M)=ZERO

C-----------------------------------------------------------------------
C     VISCOUS CORRECTION SHOULD BE APPLIED TO THE ENTIRE BLADE FOR 
C     FULLY WETTED FLOW.  FOR CAVITATING FLOW, VISCOUS FORCES SHOULD 
C     BE INTEGRATED OVER WETTED PANELS ONLY.                    
C-----------------------------------------------------------------------

C.......WETTED PANELS ON THE FACE SIDE OF THE BLADE
        IF(IWET.EQ.1.OR.IFACE.EQ.0.OR.NOCAV(M,IDXREV,2).EQ.1.OR.
     *    JCV(M,2).EQ.0) THEN
c XM YU 03/2012
C Generate a new subroutine to calculate the forces on the propeller
c with local friction coefficient cf from XFOIL
          CALL VISFRBL(M,1,NH,FX,FY,FZ,TMX,TMY,TMZ)          
c          CALL VISFRBL(M,1,NH,CDB,FX,FY,FZ,TMX,TMY,TMZ)          
c XM YU 03/2012 end

        ELSE
          LCAVF=M0(M,2)-JCV(M,2)-NSPP(M,2)-1
          CALL VISFR(M,1,LCAVF,CDB,FX,FY,FZ,TMX,TMY,TMZ)          
          M0MF=M0(M,2)
          CALL VISFR(M,M0MF,NH,CDB,FX,FY,FZ,TMX,TMY,TMZ)          
          IF(NSPP(M,2).EQ.1) THEN
            FAC=CDB*FRP(M,2)
            ISPF=LCAVF+1
            CALL VISFR(M,ISPF,ISPF,FAC,FX,FY,FZ,TMX,TMY,TMZ)             
          END IF
        END IF

C.......WETTED PANELS ON THE BACK SIDE OF THE BLADE
        IF(IWET.EQ.1.OR.IFACE.EQ.1.OR.NOCAV(M,IDXREV,1).EQ.1.OR.
     *    JCV(M,1).EQ.0) THEN
c XM YU 03/2012
C Generate a new subroutine to calculate the forces on the propeller
c with local friction coefficient cf from XFOIL
          CALL VISFRBL(M,NHP,NC,FX,FY,FZ,TMX,TMY,TMZ)
c          CALL VISFRBL(M,1,NH,CDB,FX,FY,FZ,TMX,TMY,TMZ)          
c XM YU 03/2012 end

        ELSE
          M0MBM=M0(M,1)-1
          CALL VISFR(M,NHP,M0MBM,CDB,FX,FY,FZ,TMX,TMY,TMZ)
          LCAVB=M0(M,1)+JCV(M,1)+NSPP(M,1)
          CALL VISFR(M,LCAVB,NC,CDB,FX,FY,FZ,TMX,TMY,TMZ)
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

c        write(*,*) 'vis',vbfx

C-----------------------------------------------------------------------
C     ADDED NEW VARIABLES BVFX & BVFZ SO WE CAN CALCULATE LIFT AND
C     DRAG COEFFICIENTS.                                       JY110700
C-----------------------------------------------------------------------
        BVFX(M)=FX
        BVFZ(M)=FZ

 40   CONTINUE

      IF(IHUB .NE. 0) THEN
C-----------------------------------------------------------------------
C     INTEGRATE THE PRESSURE ON THE HUB
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
 60     CONTINUE
        
C-----------------------------------------------------------------------
C     INTEGRATE THE HUB FRICTION FORCE
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
C     MODIFIED VECTOR OF FORCES TO ALIGN WITH THAT OF VELOCITY INSTEAD
C     OF HUB GEOMETRY.                                          JY031898
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
C------END OF MODIFICATION (JY031898)-----------------------------------
 70       CONTINUE
          VHFX=VHFX+FX
          VHFY=VHFY+FY
          VHFZ=VHFZ+FZ
          VHMX=VHMX+TMX
          VHMY=VHMY+TMY
          VHMZ=VHMZ+TMZ
 80     CONTINUE
      ENDIF

      IF(IDUCT .NE. 0) THEN
C-----------------------------------------------------------------------
C     INTEGRATE THE PRESSURE ON THE DUCT
C-----------------------------------------------------------------------
        DFX=ZERO
        DFY=ZERO
        DFZ=ZERO
        DMX=ZERO
        DMY=ZERO
        DMZ=ZERO
        DO N=1,NDUCT
          FX=ZERO
          FY=ZERO
          FZ=ZERO
          TMX=ZERO
          TMY=ZERO
          TMZ=ZERO
          DO M=1,MDUCT
            L=INDEXD(N,M)
            CPA=CPDN(N,M)*SS(L,1)
            FX=FX+CPA*VEL(L,1)
            FY=FY+CPA*VEL(L,2)
            FZ=FZ+CPA*VEL(L,3)
            TMX=TMX+CPA*VEL(L,4)
            TMY=TMY+CPA*VEL(L,5)
            TMZ=TMZ+CPA*VEL(L,6)
          ENDDO
          DKFX(N)=FX
          DKFY(N)=FY
          DKFZ(N)=FZ
          DKMX(N)=TMX
          DKMY(N)=TMY
          DKMZ(N)=TMZ
          DFX=DFX+DKFX(N)
          DFY=DFY+DKFY(N)
          DFZ=DFZ+DKFZ(N)
          DMX=DMX+DKMX(N)
          DMY=DMY+DKMY(N)
          DMZ=DMZ+DKMZ(N)
        ENDDO
        
C-----------------------------------------------------------------------
C     INTEGRATE THE DUCT FRICTION FORCE
C-----------------------------------------------------------------------
        VDFX=ZERO
        VDFY=ZERO
        VDFZ=ZERO
        VDMX=ZERO
        VDMY=ZERO
        VDMZ=ZERO
        DO N=1,NDUCT
          FX=ZERO
          FY=ZERO
          FZ=ZERO
          TMX=ZERO
          TMY=ZERO
          TMZ=ZERO
          DO M=1,MDUCT
            L=INDEXD(N,M)
            CFA=XCDF*VTOTS(L)*SS(L,1) 
            TXH=UXDTOT(N,M)/SQRT(VTOTS(L))
            TYH=UYDTOT(N,M)/SQRT(VTOTS(L))
            TZH=UZDTOT(N,M)/SQRT(VTOTS(L))
C     INCLUDE SECTION DRAG COEFFICIENT  BY HONG SUN
C     IF(IVISC.EQ.1) THEN
C     CFA=CDVD(M)*VTOTS_V(L)*SS(L,1)  
C     TXH=UXDTOT_V(N,M)/SQRT(VTOTS_V(L))
C     TYH=UYDTOT_V(N,M)/SQRT(VTOTS_V(L))
C     TZH=UZDTOT_V(N,M)/SQRT(VTOTS_V(L))             
C     ENDIF                   
            FX=FX+CFA*TXH
            FY=FY+CFA*TYH
            FZ=FZ+CFA*TZH
            TMX=TMX+CFA*(XCT(L,2)*TZH-XCT(L,3)*TYH)
            TMY=TMY+CFA*(XCT(L,3)*TXH-XCT(L,1)*TZH)
            TMZ=TMX+CFA*(XCT(L,1)*TYH-XCT(L,2)*TXH)
          ENDDO
          VDFX=VDFX+FX
          VDFY=VDFY+FY
          VDFZ=VDFZ+FZ
          VDMX=VDMX+TMX
          VDMY=VDMY+TMY
          VDMZ=VDMZ+TMZ
        ENDDO
      ENDIF

C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
c      IF(IAN .EQ. 2) THEN
cC-----------------------------------------------------------------------
cC     INTEGRATE THE PRESSURE ON THE TIP
cC-----------------------------------------------------------------------
c        THFX=ZERO
c        THFY=ZERO
c        THFZ=ZERO
c        THMX=ZERO
c        THMY=ZERO
c        THMZ=ZERO
c        DO N=1,NTHX
c          FX=ZERO
c          FY=ZERO
c          FZ=ZERO
c          TMX=ZERO
c          TMY=ZERO
c          TMZ=ZERO
c          DO M=1,MCVT
c            L=INDEXT(N,M)
c            CPA=CPTHN(N,M)*SS(L,1)
c            FX=FX+CPA*VEL(L,1)
c            FY=FY+CPA*VEL(L,2)
c            FZ=FZ+CPA*VEL(L,3)
c            TMX=TMX+CPA*VEL(L,4)
c            TMY=TMY+CPA*VEL(L,5)
c            TMZ=TMZ+CPA*VEL(L,6)
c          ENDDO
c          THKFX(N)=FX
c          THKFY(N)=FY
c          THKFZ(N)=FZ
c          THKMX(N)=TMX
c          THKMY(N)=TMY
c          THKMZ(N)=TMZ
c          THFX=THFX+THKFX(N)
c          THFY=THFY+THKFY(N)
c          THFZ=THFZ+THKFZ(N)
c          THMX=THMX+THKMX(N)
c          THMY=THMY+THKMY(N)
c          THMZ=THMZ+THKMZ(N)
c        ENDDO
c        
cC-----------------------------------------------------------------------
cC     INTEGRATE THE PRESSURE ON THE TIP
cC-----------------------------------------------------------------------
c        CHFX=ZERO
c        CHFY=ZERO
c        CHFZ=ZERO
c        CHMX=ZERO
c        CHMY=ZERO
c        CHMZ=ZERO
c        DO N=1,NCVX
c          FX=ZERO
c          FY=ZERO
c          FZ=ZERO
c          TMX=ZERO
c          TMY=ZERO
c          TMZ=ZERO
c          DO M=1,MCVT
c            L=INDEXC(N,M)
c            CPA=CPCN(N,M)*SS(L,1)
c            FX=FX+CPA*VEL(L,1)
c            FY=FY+CPA*VEL(L,2)
c            FZ=FZ+CPA*VEL(L,3)
c            TMX=TMX+CPA*VEL(L,4)
c            TMY=TMY+CPA*VEL(L,5)
c            TMZ=TMZ+CPA*VEL(L,6)
c          ENDDO
c          CHKFX(N)=FX
c          CHKFY(N)=FY
c          CHKFZ(N)=FZ
c          CHKMX(N)=TMX
c          CHKMY(N)=TMY
c          CHKMZ(N)=TMZ
c          CHFX=CHFX+CHKFX(N)
c          CHFY=CHFY+CHKFY(N)
c          CHFZ=CHFZ+CHKFZ(N)
c          CHMX=CHMX+CHKMX(N)
c          CHMY=CHMY+CHKMY(N)
c          CHMZ=CHMZ+CHKMZ(N)
c        ENDDO
c
cC-----------------------------------------------------------------------
cC     INTEGRATE THE TIP FRICTION FORCE
cC-----------------------------------------------------------------------
c        VTHFX=ZERO
c        VTHFY=ZERO
c        VTHFZ=ZERO
c        VTHMX=ZERO
c        VTHMY=ZERO
c        VTHMZ=ZERO
c        DO N=1,NTHX
c          FX=ZERO
c          FY=ZERO
c          FZ=ZERO
c          TMX=ZERO
c          TMY=ZERO
c          TMZ=ZERO
c          DO M=1,MCVT
c            L=INDEXT(N,M)
c            CFA=XCDF*VTOTS(L)*SS(L,1)
c            TXH=UXTHTOT(N,M)/SQRT(VTOTS(L))
c            TYH=UYTHTOT(N,M)/SQRT(VTOTS(L))
c            TZH=UZTHTOT(N,M)/SQRT(VTOTS(L))
c            FX=FX+CFA*TXH
c            FY=FY+CFA*TYH
c            FZ=FZ+CFA*TZH
c            TMX=TMX+CFA*(XCT(L,2)*TZH-XCT(L,3)*TYH)
c            TMY=TMY+CFA*(XCT(L,3)*TXH-XCT(L,1)*TZH)
c            TMZ=TMX+CFA*(XCT(L,1)*TYH-XCT(L,2)*TXH)
c          ENDDO
c          
c          VTHFX=VTHFX+FX
c          VTHFY=VTHFY+FY
c          VTHFZ=VTHFZ+FZ
c          VTHMX=VTHMX+TMX
c          VTHMY=VTHMY+TMY
c          VTHMZ=VTHMZ+TMZ
c        ENDDO
c
cC-----------------------------------------------------------------------
cC     INTEGRATE THE TIP FRICTION FORCE
cC-----------------------------------------------------------------------
c        VCHFX=ZERO
c        VCHFY=ZERO
c        VCHFZ=ZERO
c        VCHMX=ZERO
c        VCHMY=ZERO
c        VCHMZ=ZERO
c        DO N=1,NCVX
c          FX=ZERO
c          FY=ZERO
c          FZ=ZERO
c          TMX=ZERO
c          TMY=ZERO
c          TMZ=ZERO
c          DO M=1,MCVT
c            L=INDEXC(N,M)
c            CFA=XCDF*VTOTS(L)*SS(L,1)
c            TXH=UXCHTOT(N,M)/SQRT(VTOTS(L))
c            TYH=UYCHTOT(N,M)/SQRT(VTOTS(L))
c            TZH=UZCHTOT(N,M)/SQRT(VTOTS(L))
c            FX=FX+CFA*TXH
c            FY=FY+CFA*TYH
c            FZ=FZ+CFA*TZH
c            TMX=TMX+CFA*(XCT(L,2)*TZH-XCT(L,3)*TYH)
c            TMY=TMY+CFA*(XCT(L,3)*TXH-XCT(L,1)*TZH)
c            TMZ=TMX+CFA*(XCT(L,1)*TYH-XCT(L,2)*TXH)
c          ENDDO
c          
c          VCHFX=VCHFX+FX
c          VCHFY=VCHFY+FY
c          VCHFZ=VCHFZ+FZ
c          VCHMX=VCHMX+TMX
c          VCHMY=VCHMY+TMY
c          VCHMZ=VCHMZ+TMZ
c        ENDDO         
c      ENDIF
C/e S.N.KIM | Aug. 2018.

C----------------------------------------------------------------------
C     I'M MULTIPLYING THE FORCES BY 1/8*J^2 AND MOMENTS BY 1/16*J^2 SO 
C     THE NON-DIMENSIONALIZATION IS THE SAME AS HPUF3A.  THIS 
C     MODIFICATION IS NOT NEEDED FOR TKT AND TKQ BECAUSE IT HAS 
C     ALREADY BEEN DONE ABOVE.
C     F/(RHO*N^2*D^4)=CPB*SS(*J^2/4)(*1/2) BECAUSE 
C     CP(HPUF)=1/2CP(PROPCAV)
C     ORIGINALLY, FXV(1) WAS = BFX+HFX+VBFX+VHFX               JY012198
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C     EIGHT POTENTIAL FORCE COMPONENTS                        
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

      IF(IDUCT .NE. 0) THEN
        FDXP(1)=DFX*.125*ADVCO**2
        FDXP(2)=DFY*.125*ADVCO**2
        FDXP(3)=DFZ*.125*ADVCO**2
        FDXP(4)=DMX*.0625*ADVCO**2
        FDXP(5)=DMY*.0625*ADVCO**2
        FDXP(6)=DMZ*.0625*ADVCO**2
C     WRITE(*,*) 'DFX TOTAL = ', DFX*NBLADE
      ENDIF

C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
c      IF(IAN .EQ. 2) THEN
c        FTXP(1) = (THFX+CHFX)*.125*ADVCO**2
c        FTXP(2) = (THFY+CHFY)*.125*ADVCO**2
c        FTXP(3) = (THFZ+CHFZ)*.125*ADVCO**2
c        FTXP(4) = (THMX+CHMX)*.0625*ADVCO**2
c        FTXP(5) = (THMY+CHMY)*.0625*ADVCO**2
c        FTXP(6) = (THMZ+CHMZ)*.0625*ADVCO**2
c      ENDIF
C/e S.N.KIM | Aug. 2018.


C----------------------------------------------------------------------
C     EIGHT VISCOUS FORCE COMPONENTS                        
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

      IF(IDUCT .NE. 0) THEN
        FDXV(1) = VDFX*.125*ADVCO**2
        FDXV(2) = VDFY*.125*ADVCO**2
        FDXV(3) = VDFZ*.125*ADVCO**2
        FDXV(4) = VDMX*.0625*ADVCO**2
        FDXV(5) = VDMY*.0625*ADVCO**2
        FDXV(6) = VDMZ*.0625*ADVCO**2
      ENDIF

      DO I = 1 , 6
C     -- BLADE TOTAL FORCE (POTENTIAL + VISCOUS)
        FBXPV(I) = FBXP(I) + FBXV(I)
C     -- HUB TOTAL FORCE (POTENTIAL + VISCOUS)
        FHXPV(I) = FHXP(I) + FHXV(I)
C     -- DUCT TOTAL FORCE (POTENTIAL + VISCOUS)
        FDXPV(I) = FDXP(I) + FDXV(I)
C     ==========TOTAL POTENTIAL
        FXP(I)=FBXP(I)+FHXP(I)+FDXP(I)+FTXP(I)
C     ========= TOTAL VISCOUS
        FXV(I)=FBXV(I)+FHXV(I)+FDXV(I)
C     ========== TOTAL FORCES
        FTOTAL(I) = FXP(I) + FXV(I)
      ENDDO

      IF(IDUCT .NE. 0) THEN
C     WRITE(*,*) ' TOTAL FORCE FROM PRESSURE INTEGRATION ON DUCT'
C     WRITE(*,*) ' POTENTIAL FORCES, POTENTIAL + FRICTION FORCES'
C     WRITE(*,*) 'CTD  = ', DFX*NBLADE, (DFX+VDFX)*NBLADE
C     WRITE(*,*) 'KXD = ',FDXP(1)*NBLADE, FDXPV(1)*NBLADE
C     WRITE(*,*) 'KYD = ',FDXP(2)*NBLADE, FDXPV(2)*NBLADE
C     WRITE(*,*) 'KZD = ',FDXP(3)*NBLADE, FDXPV(3)*NBLADE
C     WRITE(*,*) 'MXD = ',FDXP(4)*NBLADE, FDXPV(4)*NBLADE
C     WRITE(*,*) 'MYD = ',FDXP(5)*NBLADE, FDXPV(5)*NBLADE
C     WRITE(*,*) 'MZD = ',FDXP(6)*NBLADE, FDXPV(6)*NBLADE
C     WRITE(*,*) 'CTD = ', (DFX+VDFX)*NBLADE
C     WRITE(*,*) 'KT  = ', FDXPV(1)*NBLADE
C     WRITE(*,*) 'KQ  = ', FDXPV(4)*NBLADE
        FKTB = NBLADE*FBXPV(1) 
        FKTD = NBLADE*FDXPV(1)
        FKQB = NBLADE*FBXPV(4)
C     WRITE(*,*) 'KT_TOTAL, KT_DUCT, KQ_TOTAL'
C     WRITE(*,*) FKTB+FKTD, FKTD, 10*FKQB
      ENDIF

C----------------------------------------------------------------------
C     CALCULATE LIFT AND DRAG COEFFICIENT FOR HYDROFOILS.      JY110300
C----------------------------------------------------------------------
      IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
        SUMAREA=0.
        SUMFX=0.
        SUMFZ=0.
        DO M=1,MR
          CHDM=CHORD(M)+CHORD(M+1)
          DELR=RZ(M+1)-RZ(M)
          IF(DELR .EQ. 0.) DELR=1.E-6
          SUMFX=SUMFX+(BKFX(M)+BVFX(M))
          SUMFZ=SUMFZ+(BKFZ(M)+BVFZ(M))
          SUMAREA=SUMAREA+CHDM*DELR
        END DO
        CZ=NBLADE*(SUMFZ/SUMAREA)
        CX=NBLADE*(SUMFX/SUMAREA)
        COSA=COS(ALPHA*RAD)
        SINA=SIN(ALPHA*RAD)
        CLIFT=CZ*COSA-CX*SINA
        CDRAG=CX*COSA+CZ*SINA
      END IF

      RETURN
      END

