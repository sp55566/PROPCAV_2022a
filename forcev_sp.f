      SUBROUTINE FORCEV_SP
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
C     JY113099       Copied from FORCEV.F.  Modfied subroutine so the
C                    viscous correction is applied over the wetted,
C                    submerged panels only.
C     JY071700       Modified subroutine to allow for detachment for 
C                    ISP=1.
C     JY111200       Removed Polhamus' correction as well as smoothing
C                    of the forces for ISP=1.
C
C***********************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

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

         IF(ISC.EQ.0) THEN
            N1=1
            N2=NC
         ELSE
            N1=N0(2)
            N2=N0(1)-1
         END IF

         DO 10 N=N1,N2
            L=INDEXB(N,M)
            CPA=CPB(N,M)*SS(L,1)

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

         BKFX(M)=FX
         BKFY(M)=FY
         BKFZ(M)=FZ
         BKMX(M)=TMX
         BKMY(M)=TMY
         BKMZ(M)=TMZ
20    CONTINUE

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
C     Apply viscous correction over wetted, submerged panels only.
C                                                               JY120399
C-----------------------------------------------------------------------
         IF(ICB(M,2,IDXREV).EQ.1) THEN
            IF(ISC.EQ.0) THEN
               N1=IW(1,M,IDXREV)
            ELSE
               N1=MAX(N0(2),IW(1,M,IDXREV))
            END IF
               
            N2=IW(2,M,IDXREV)
            IF(N2.GE.N1) CALL VISFR(M,N1,N2,CDB,FX,FY,FZ,TMX,TMY,TMZ)
         END IF

         IF(ICB(M,1,IDXREV).EQ.1) THEN
            N1=IC(1,M,IDXREV)

            IF(JCV(M,1).EQ.0) THEN
               IF(ISC.EQ.0) THEN
                  N2=IC(2,M,IDXREV)
               ELSE
                  N2=MIN(IC(2,M,IDXREV),N0(1)-1)
               END IF
            ELSE
               N2=M0(M,1)-1
            END IF
            IF(N2.GE.N1) CALL VISFR(M,N1,N2,CDB,FX,FY,FZ,TMX,TMY,TMZ)
         END IF
C-----------------------------------------------------------------------

         VBFX=VBFX+FX
         VBFY=VBFY+FY
         VBFZ=VBFZ+FZ
         VBMX=VBMX+TMX
         VBMY=VBMY+TMY
         VBMZ=VBMZ+TMZ
40    CONTINUE

C----------------------------------------------------------------------
C     Set hub forces equal to zero for ISP=1
C----------------------------------------------------------------------
      HFX=ZERO
      HFY=ZERO
      HFZ=ZERO
      HMX=ZERO
      HMY=ZERO
      HMZ=ZERO

      VHFX=ZERO
      VHFY=ZERO
      VHFZ=ZERO
      VHMX=ZERO
      VHMY=ZERO
      VHMZ=ZERO

C----------------------------------------------------------------------
C     I'm multiplying the forces by 1/8*J^2 and moments by 1/16*J^2 so 
C     the non-dimensionalization is the same as HPUF3A.  This 
C     modification is not needed for TKT and TKQ because it has 
C     already been done above.
C     F/(rho*n^2*D^4)=CPB*SS(*J^2/4)(*1/2) because 
C     Cp(HPUF)=1/2Cp(PROPCAV)
C     Originally, FXV(1) was = bfx+hfx+vbfx+vhfx               JY012198
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C     Eight potential force components                        
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
C     Eight viscous force components                        
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

      RETURN
C))))))))))))))))))))) End of subroutine FORCE (((((((((((((((((((((((((
      END
