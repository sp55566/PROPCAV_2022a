C ===================================
      SUBROUTINE PREKUTTA
C ===================================
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
C      DIMENSION WPRE(NTZ)
      ALLOCATABLE :: WPRE(:)

C      COMMON/BKUTTA/POTEMP(NTZ),WWK(NTZ,MBZ)

C.....Note: If ITMAX>100, You MUST set ITMX=ITMAX+1 in BLIC2.F(JY051800)

      ITMAX=100
      NORD=NPANEL
      NBLOCK=MR
      IF(IHUB .NE. 0) NBLOCK = NBLOCK + NHBX
      IF(IDUCT .NE. 0) NBLOCK = NBLOCK + MDUCT
      IF(ITUN .NE. 0) NBLOCK = NBLOCK + NAXT
c      IF(IAN .EQ. 2) NBLOCK = NBLOCK + NTHX + NCVX

      DO M=1,MR
         NPERB(M)=NC
      ENDDO
      
      NNN = MR
      
      IF(IHUB .NE. 0) THEN
         DO N=1,NHBX
            NPERB(NNN+N)=MHBT
         ENDDO
         NNN = NNN + NHBX
      ENDIF
      
      IF(IDUCT .NE. 0) THEN
         DO N=1,MDUCT
            NPERB(NNN+N)=NDUCT
         ENDDO
         NNN = NNN + MDUCT
      ENDIF
      
      IF(ITUN .NE. 0) THEN
         DO N=1,NAXT
            NPERB(NNN+N)=MTUNEL
         ENDDO
         NNN = NNN + NAXT
      ENDIF

C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.      
c      IF(IAN.EQ.2) THEN
c         DO N=1,NTHX
c            NPERB(NNN+N)=MCVT
c         END DO
c         
c         NNN=NNN+NTHX
c         DO N=1,NCVX
c            NPERB(NNN+N)=MCVT
c         END DO
c         NNN = NNN + NCVX
c      END IF
C/e S.N.KIM | Aug. 2018.
      
C.....Tolerance of the matrix solution is set to be 0.000005
      TOL=5.0E-06

CVV
      ALLOCATE(WPRE(NTZ))
CVV      
      
      IF(ITERKT.NE.0.AND.NTSTEP.EQ.1) THEN
         DO M=1,MR
            DO  I=1,NORD
               WWK(I,M)=0.0
               IF(ISTEADY.EQ.0) THEN
                  WPRE(I)=-W(I,M)
               ELSE
                  WPRE(I)=-HALF*W(I,M)-WSUBIF(I,M)
               END IF
            ENDDO 

            CALL CBLIC1(WWK(1,M),WPRE,TOL,NORD,NBLOCK,NPERB,ITMAX,AA,NTZ)
         ENDDO
      END IF
        
      IF(ITERKT.NE.0) THEN
         DO I = 1, NPANEL
            POTEMP(I) = POT(I)
         ENDDO
         CALL KUTTA(ITERKT)
      END IF


CVV
      DEALLOCATE(WPRE) 
CVV

      RETURN
      END
