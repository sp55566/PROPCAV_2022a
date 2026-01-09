      SUBROUTINE GEO3DW(NPANEL,XV,CHRLEWS,IER)
C***********************************************************************
C
C     DOUBLE PRECISION
C
C***********************************************************************
C-----------------------------------------------------------------------
C     SAME AS GEOM3D
C-----------------------------------------------------------------------
      USE GEOMTW
      USE UVECW
      INCLUDE 'PARAM.INC'
!     PARAMETER(MZ=NSCWZ)
      DOUBLE PRECISION VPP,XCTP,DIRP,S,SIDE,GAUSSX,GAUSSW,DEU,DEV,
     *                 COSPH,SINPH,CHRLENS,DR1,DR2,DR3,VEL4S,VEL5S,
     *                 VEL6S,ULC,VLC
      DOUBLE PRECISION V(4,3)
c      COMMON /GEOMTW/ XVP(MZ,4),YVP(MZ,4),XCT(MZ,3),DIR(MZ,3,3),
c     *               SS(MZ,15),SID(MZ,4),GSX(MZ,4,3),GSW(MZ,4),
c     *               VEL(MZ,6),DELU(MZ),DELV(MZ),COSPHI(MZ),
c     *               SINPHI(MZ)
      COMMON /PANEL/ VPP(4,3),XCTP(3),DIRP(3,3),S(15),SIDE(4),
     *               GAUSSX(4,3),GAUSSW(4),DEU,DEV,COSPH,SINPH
      COMMON/VSPAN/ULC(3),VLC(3)
C      COMMON/UVECW/ULW(MZ,3)
!     DIMENSION XV(MZ,4,3),CHRLEWS(MZ)
      DIMENSION XV(NSCWZ,4,3),CHRLEWS(NSCWZ)

      INTEGER MZ=NSCWZ

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
         CHRLEWS(I) = SNGL(CHRLENS)

         DO 30 K=1,15
            SSW(I,K) = SNGL(S(K))
 30      CONTINUE

         DO 40 LG=1,4
            XVPW(I,LG) = SNGL(VPP(LG,1))
            YVPW(I,LG) = SNGL(VPP(LG,2))
            GSWW(I,LG) = SNGL(GAUSSW(LG))
            SIDW(I,LG) = SNGL(SIDE(LG))
 40      CONTINUE

         DELUW(I) = SNGL(DEU)
         DELVW(I) = SNGL(DEV)
         COSPHIW(I) = SNGL(COSPH)
         SINPHIW(I) = SNGL(SINPH)

C-----------------------------------------------------------------------
C        LOOP OVER x-y-z DIMENSIONS
C-----------------------------------------------------------------------

         DO 70 K=1,3
            XCTW(I,K) = SNGL(XCTP(K))

            DO 50 K1=1,3
               DIRW(I,K,K1) = SNGL(DIRP(K,K1))
 50         CONTINUE
            
            DO 60 M1=1,4
               GSXW(I,M1,K) = SNGL(GAUSSX(M1,K))
 60         CONTINUE
            ULW(I,K)=SNGL(ULC(K))

 70      CONTINUE

C-----------------------------------------------------------------------
C        NORMAL VELOCITIES AT THE PANEL CENTROIDS
C-----------------------------------------------------------------------

         DR1=DIRP(3,1)
         DR2=DIRP(3,2)
         DR3=DIRP(3,3)

C........NORMAL VECTOR OUTWARD FROM THE FLUID
         VELW(I,1) = SNGL(DR1)
         VELW(I,2) = SNGL(DR2)
         VELW(I,3) = SNGL(DR3)
         VEL4S=XCTP(2)*DR3-XCTP(3)*DR2
         VEL5S=XCTP(3)*DR1-XCTP(1)*DR3
         VEL6S=XCTP(1)*DR2-XCTP(2)*DR1
         VELW(I,4) = SNGL(VEL4S)
         VELW(I,5) = SNGL(VEL5S)
         VELW(I,6) = SNGL(VEL6S)

200   CONTINUE


99    RETURN
      END
