      SUBROUTINE GEO3DWo(NPANEL,XV,CHRLEWS,KK,IER,MZ,KZ1)
C***********************************************************************
C
C     DOUBLE PRECISION
C
C***********************************************************************
C-----------------------------------------------------------------------
C     SAME AS GEOM3D
C-----------------------------------------------------------------------
      USE GEOMTWO
      USE UVECWo
      INCLUDE 'PARAM.INC'
!     PARAMETER(MZ=NSCWZ)
!     PARAMETER(KZ1=KZ)
      DOUBLE PRECISION VPP,XCTP,DIRP,S,SIDE,GAUSSX,GAUSSW,DEU,DEV,
     *                 COSPH,SINPH,CHRLENS,DR1,DR2,DR3,VEL4S,VEL5S,
     *                 VEL6S,ULC,VLC
      DOUBLE PRECISION V(4,3)
c      COMMON /GEOMTWo/ XVP(MZ,4,KZ1),YVP(MZ,4,KZ1),XCT(MZ,3,KZ1),
c     *                 DIR(MZ,3,3,KZ1),SS(MZ,15,KZ1),SID(MZ,4,KZ1),
c     *               GSX(MZ,4,3,KZ1),GSW(MZ,4,KZ1),VEL(MZ,6,KZ1),
c     *               DELU(MZ,KZ1),DELV(MZ,KZ1),COSPHI(MZ,KZ1),
c     *               SINPHI(MZ,KZ1)
      COMMON /PANEL/ VPP(4,3),XCTP(3),DIRP(3,3),S(15),SIDE(4),
     *               GAUSSX(4,3),GAUSSW(4),DEU,DEV,COSPH,SINPH
      COMMON/VSPAN/ULC(3),VLC(3)
!     COMMON/UVECWo/ULWo(MZ,3,KZ1)
      DIMENSION XV(MZ,4,3),CHRLEWS(MZ,KZ1)

      if (.NOT.allocated(ULWo)) then
        allocate(ULWo(MZ,3,KZ))
      end if
C-----------------------------------------------------------------------
C     LOOP OVER PANELS
C-----------------------------------------------------------------------

      DO 200 I=1,NPANEL

         DO 20 J=1,4

            DO 10 K=1,3
               V(J,K)=DBLE(XV(I,J,K))
 10         CONTINUE

 20      CONTINUE

C-----------------------------------------------------------------------
C        CALL TO GPANEL
C-----------------------------------------------------------------------

         CALL GPANEL(V,CHRLENS,IER)
         IF(IER.EQ.0) GO TO 99
         CHRLEWS(I,KK) = SNGL(CHRLENS)

         DO 30 K=1,15
            SSWO(I,K,KK) = SNGL(S(K))
 30      CONTINUE

         DO 40 LG=1,4
            XVPWO(I,LG,KK) = SNGL(VPP(LG,1))
            YVPWO(I,LG,KK) = SNGL(VPP(LG,2))
            GSWWO(I,LG,KK) = SNGL(GAUSSW(LG))
            SIDWO(I,LG,KK) = SNGL(SIDE(LG))
 40      CONTINUE

         DELUWO(I,KK) = SNGL(DEU)
         DELVWO(I,KK) = SNGL(DEV)
         COSPHIWO(I,KK) = SNGL(COSPH)
         SINPHIWO(I,KK) = SNGL(SINPH)

C-----------------------------------------------------------------------
C        LOOP OVER x-y-z DIMENSIONS
C-----------------------------------------------------------------------

         DO 70 K=1,3
            XCTWO(I,K,KK) = SNGL(XCTP(K))

            DO 50 K1=1,3
               DIRWO(I,K,K1,KK) = SNGL(DIRP(K,K1))
 50         CONTINUE
            
            DO 60 M1=1,4
               GSXWO(I,M1,K,KK) = SNGL(GAUSSX(M1,K))
 60         CONTINUE
            ULWo(I,K,KK)=SNGL(ULC(K))

 70      CONTINUE

C-----------------------------------------------------------------------
C        NORMAL VELOCITIES AT THE PANEL CENTROIDS
C-----------------------------------------------------------------------

         DR1=DIRP(3,1)
         DR2=DIRP(3,2)
         DR3=DIRP(3,3)

C........NORMAL VECTOR OUTWARD FROM THE FLUID
         VELWO(I,1,KK) = SNGL(DR1)
         VELWO(I,2,KK) = SNGL(DR2)
         VELWO(I,3,KK) = SNGL(DR3)
         VEL4S=XCTP(2)*DR3-XCTP(3)*DR2
         VEL5S=XCTP(3)*DR1-XCTP(1)*DR3
         VEL6S=XCTP(1)*DR2-XCTP(2)*DR1
         VELWO(I,4,KK) = SNGL(VEL4S)
         VELWO(I,5,KK) = SNGL(VEL5S)
         VELWO(I,6,KK) = SNGL(VEL6S)

200   CONTINUE


99    RETURN
      END
