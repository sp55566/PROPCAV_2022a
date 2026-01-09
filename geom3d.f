      SUBROUTINE GEOM3D(NPANEL,XV,CHRLEPS,IER)
C***********************************************************************
C
C     DOUBLE PRECISION
C
C***********************************************************************
C-----------------------------------------------------------------------
C
C  SUBROUTINE GEOM3D GENERARTES THE PANEL GOMETRICAL DATA
C
C-----------------------------------------------------------------------
C
C  INPUT :
C  -----
C  NPANEL       : NUMBER OF PANELS 
C  XV(I,J,K) ,  : ARRAY STORING THE K-TH REFERENCE COORDINATE
C                 OF THE J-TH VERTEX OF THE I-TH PANEL, WHERE
C                    I = 1,..,NPANEL
C                    J = 1,..,4           : THE PANEL VERTICES ARE
C                                           NUMBERED IN THE CLOCKWISE
C                                           SENSE WHEN VIEWED FROM THE
C                                           FLUID DOMAIN.
C                    K = 1,2,3 <=> x,y,z.
C
C
C  OUTPUT :  GEOMETRICAL PANEL DATA AND NORMAL VELOCITIES TRANSFERED
C  ------    VIA THE COMMON BLOCK / GEOMT /.
C
C           IER : ERROR INDEX
C               =  1 : ACCEPTABLE PANEL
C               =  0 : UNACCEPTABLE PANEL
C
C    Date of last revision           Revision
C    ---------------------         -----------
C    1-30-91 Neal Fine        -Put UL,VL,WL in common /PANEL/
C                             -Transfered UL,.. to ULL,... and put
C                              these arrays in common /PCAV/
C    2-26-91 Neal Fine        -changed nmz to 400
C    5-08-91 Neal Fine        -changed nmz to 1600
C    5-14-91 Neal Fine        -changed nmz to 4800
C    8-29-91 Neal Fine        -put VLC in common /VSPAN/
C
C-----------------------------------------------------------------------
      USE GEOMT
      USE UVECB
      INCLUDE 'PARAM.INC'
!     PARAMETER(NMZ=NPANZ)
      DOUBLE PRECISION VPP,XCTP,DIRP,S,SIDE,GAUSSX,GAUSSW,DEU,DEV,
     *                 COSPH,SINPH,CHRLENS,DR1,DR2,DR3,VEL4S,VEL5S,
     *                 VEL6S,VLC,ULC
      DOUBLE PRECISION V(4,3)

c      COMMON /GEOMT/ XVP(NMZ,4),YVP(NMZ,4),XCT(NMZ,3),DIR(NMZ,3,3),
c     *               SS(NMZ,15),SID(NMZ,4),GSX(NMZ,4,3),GSW(NMZ,4),
c     *               VEL(NMZ,6),DELU(NMZ),DELV(NMZ),COSPHI(NMZ),
c     *               SINPHI(NMZ)
      COMMON /PANEL/ VPP(4,3),XCTP(3),DIRP(3,3),S(15),SIDE(4),
     *               GAUSSX(4,3),GAUSSW(4),DEU,DEV,COSPH,SINPH
      COMMON/VSPAN/ULC(3),VLC(3)
c      COMMON/UVECB/ULL(NMZ,3),VLL(NMZ,3),WLL(NMZ,3)
!     DIMENSION XV(NMZ,4,3),CHRLEPS(NMZ)
      DIMENSION XV(NPANZ,4,3),CHRLEPS(NPANZ)
      integer NMZ
      NMZ = NPANZ



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
         IF(IER.EQ.0) THEN
            WRITE(*,*) ' ------ I: ',I
            DO 25 J=1,4
               WRITE(*,'(3F12.6)') (XV(I,J,K),K=1,3)
 25         CONTINUE
            WRITE(*,'(A)') ' '
            GO TO 99
         END IF
         CHRLEPS(I) = SNGL(CHRLENS)
         DO 30 K=1,15
            SS(I,K) = SNGL(S(K))
 30      CONTINUE
         DO 40 LG=1,4
            XVP(I,LG) = SNGL(VPP(LG,1))
            YVP(I,LG) = SNGL(VPP(LG,2))
            GSW(I,LG) = SNGL(GAUSSW(LG))
            SID(I,LG) = SNGL(SIDE(LG))
 40      CONTINUE
         DELU(I) = SNGL(DEU)
         DELV(I) = SNGL(DEV)
         COSPHI(I) = SNGL(COSPH)
         SINPHI(I) = SNGL(SINPH)

C-----------------------------------------------------------------------
C        LOOP OVER x-y-z DIMENSIONS
C        modified 1-30-91 to include unit vectors ULL,VLL....Neal Fine
C-----------------------------------------------------------------------

         DO 70 K=1,3
            XCT(I,K) = SNGL(XCTP(K))
            DO 50 K1=1,3
               DIR(I,K,K1) = SNGL(DIRP(K,K1))
 50         CONTINUE
            DO 60 M1=1,4
                GSX(I,M1,K) = SNGL(GAUSSX(M1,K))
 60          CONTINUE
            VL(I,K)=SNGL(VLC(K))
            UL(I,K)=SNGL(ULC(K))
 70      CONTINUE

C-----------------------------------------------------------------------
C        NORMAL VELOCITIES AT THE PANEL CENTROIDS
C-----------------------------------------------------------------------

         DR1=DIRP(3,1)
         DR2=DIRP(3,2)
         DR3=DIRP(3,3)

C........NORMAL VECTOR OUTWARD FROM THE FLUID
         VEL(I,1) = SNGL(DR1)
         VEL(I,2) = SNGL(DR2)
         VEL(I,3) = SNGL(DR3)
         VEL4S=XCTP(2)*DR3-XCTP(3)*DR2
         VEL5S=XCTP(3)*DR1-XCTP(1)*DR3
         VEL6S=XCTP(1)*DR2-XCTP(2)*DR1
         VEL(I,4) = SNGL(VEL4S)
         VEL(I,5) = SNGL(VEL5S)
         VEL(I,6) = SNGL(VEL6S)

 200  CONTINUE

 99   RETURN
      END
