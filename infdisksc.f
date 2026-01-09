      SUBROUTINE INFDISKSC
************************************************************************
*   INFDISKSC: INFluence coefficients due to the wake disk and the hub *
*              disk at the control points in the supercavity wake      *
*      --- Calculate the influence coefficients due to a source disk   *
*          represents the ultimate wake and a dipole disk represents   *
*          the far-upstream hub                                        *
*      moved to INFDISKSC from INFDISK on 08-14-92 Neal Fine           *
*                                                                      *
*      JY010600    Modified subroutine to also calculate the influence *
*                  of hub dipole disk to supercavitating panels.       *
*      JY111900    The hub disk is not needed for IHUB=2.              *
*      JY052401    Modified subroutine to accomodate new hub options.  *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVC.INC'

C-----------------------------------------------------------------------
C     prepare the geometry of a source disk replacing the ultimate wake
C
C     Clockwise viewing              (N,2) 2 *-----* 3 (N+1,2)
C       from upstream                        |     |
C                                            |     |
C                                    (N,1) 1 *-----* 4 (N+1,1)
C
C-----------------------------------------------------------------------
      IDDUCT2 = NPANB + NPANH
      IDDUCT3 = NPANB + NPANH + NPAND

C/s S.N.KIM | If FWA is used, ultimate wake disk is ignored.
C           | Same for tip loaded propellers (IAN=5).
C           | Tip loaded propellers use FWA, hence do not
C           | have the ultimate wake disk.
      IF (IAN.EQ.6.OR.IREADW.EQ.1) GO TO 1004
C/e S.N.KIM | Aug. 2018.

      DO 100 N=1,NUWDK
         XGW(N,1,1)=XWDK(N,1)
         XGW(N,1,2)=YWDK(N,1)
         XGW(N,1,3)=ZWDK(N,1)
         XGW(N,2,1)=XWDK(N,2)
         XGW(N,2,2)=YWDK(N,2)
         XGW(N,2,3)=ZWDK(N,2)
         XGW(N,3,1)=XWDK(N+1,2)
         XGW(N,3,2)=YWDK(N+1,2)
         XGW(N,3,3)=ZWDK(N+1,2)
         XGW(N,4,1)=XWDK(N+1,1)
         XGW(N,4,2)=YWDK(N+1,1)
         XGW(N,4,3)=ZWDK(N+1,1)
100   CONTINUE
      CALL GEO3DW(NUWDK,XGW,CHRLEWS,IER)
      IF(IER.EQ.0) THEN
         WRITE(*,'(A)') ' UNACCEPTABLE PANELS IN INFDISKSC!'
         STOP
      END IF
C-----------------------------------------------------------------------
C     Calculate the influence coefficients due to the wake disk
C       WDK1(I): induced potential at I due to the wake disk
C-----------------------------------------------------------------------
      DO 110 I=1,NPWAKS
         WDK1(I)=0.
110   CONTINUE
      DO 170 L=1,NUWDK
C........Transfer data to common block /GEOM/
         DO 120 K=1,4
            XV(K)=XVPW(L,K)
            YV(K)=YVPW(L,K)
            SIDE(K)=SIDW(L,K)
120      CONTINUE
         DO 130 K=1,15
            S(K)=SSW(L,K)
130      CONTINUE
C........Calculate the influence coefficients
         DO 160 I=1,NPWAKS
            DO 150 KK=1,NBLADE
C..............Transfer control points to the local coordinate
               XLOC=0.
               YLOC=0.
               ZLOC=0.
               DO 140 K=1,3
                  XLOC=XLOC+(XCPW(I,K,KK)-XCTW(L,K))*DIRW(L,1,K)
                  YLOC=YLOC+(XCPW(I,K,KK)-XCTW(L,K))*DIRW(L,2,K)
                  ZLOC=ZLOC+(XCPW(I,K,KK)-XCTW(L,K))*DIRW(L,3,K)
140             CONTINUE
C..............Compute induced potentials
               IMR=1
               CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(L),FS,FD,
     *                   FSX,FSY,FDX,FDY,FDZ,0,IMR)
               IF(IDUCT .EQ. 1 .AND. IDOPT .EQ. 1 
     %              .AND. I .GT. IDDUCT2 .AND. I .LE. IDDUCT2) FS = 0.0
               WDK1(I)=WDK1(I)+FS
150         CONTINUE
160       CONTINUE
170    CONTINUE

 1004  CONTINUE

C-----------------------------------------------------------------------
C     Calculate the influence of the hub dipole disk on the 
C     supercavitating panels.                                   JY010600
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     Dipole disk representing the far upstream hub
C-----------------------------------------------------------------------

       DO 210 I=1,NPWAKS
          CHDK1(I)=0.
 210   CONTINUE

C....The hub disk is not needed for IHUB=2 or IHUB=4 b/c it's close 
C....far upstream (JY052401)
       IF(IHUB.EQ.0.OR.IHUB.EQ.2.OR.IHUB.EQ.4
     *             .OR.IHUB.EQ.5.OR.IHUB.EQ.6) GO TO 1111

C-----------------------------------------------------------------------
C     prepare the geometry of a dipole disk replacing the far upstream
C       hub
C     Clockwise viewing              (N,2) 2 *-----* 3 (N+1,2)
C       from upstrean                        |     |
C                                            |     |
C                                    (N,1) 1 *-----* 4 (N+1,1)
C
C-----------------------------------------------------------------------
       DO 200 N=1,NHBDK
          XGW(N,1,1)=XHDK(N,1)
          XGW(N,1,2)=YHDK(N,1)
          XGW(N,1,3)=ZHDK(N,1)
          XGW(N,2,1)=XHDK(N,2)
          XGW(N,2,2)=YHDK(N,2)
          XGW(N,2,3)=ZHDK(N,2)
          XGW(N,3,1)=XHDK(N+1,2)
          XGW(N,3,2)=YHDK(N+1,2)
          XGW(N,3,3)=ZHDK(N+1,2)
          XGW(N,4,1)=XHDK(N+1,1)
          XGW(N,4,2)=YHDK(N+1,1)
          XGW(N,4,3)=ZHDK(N+1,1)
 200   CONTINUE
       CALL GEO3DW(NHBDK,XGW,CHRLEWS,IER)
       IF(IER.EQ.0) THEN
          WRITE(*,'(A)') ' UNACCEPTABLE PANELS IN INFDISKSC!'
          STOP
       END IF

C-----------------------------------------------------------------------
C     Calculate the influence coefficients due to the hub disk
C       CHDK1(I): induced potentials at supercavitating panel I due to 
C                 the hub disk
C-----------------------------------------------------------------------

       DO 270 N=1,NHBDK
C
C........Transfer data to common block /GEOM/
          DO 220 K=1,4
             XV(K)=XVPW(N,K)
             YV(K)=YVPW(N,K)
             SIDE(K)=SIDW(N,K)
 220      CONTINUE
          DO 230 K=1,15
             S(K)=SSW(N,K)
 230      CONTINUE
          DO 260 I=1,NPWAKS
             DO 250 KK=1,NBLADE
C
C..............Transfer control points to local coordinate
                XLOC=0.
                YLOC=0.
                ZLOC=0.
                DO 240 K=1,3
                   XLOC=XLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,1,K)
                   YLOC=YLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,2,K)
                   ZLOC=ZLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,3,K)
 240            CONTINUE
C
C..............Compute induced potentials
                IMR=1
                CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(N),FS,FD,
     1               FSX,FSY,FDX,FDY,FDZ,0,IMR)
                IF(IDUCT .EQ. 1 .AND. IDOPT .EQ. 1 
     %             .AND. I .GT. IDDUCT2 .AND. I .LE. IDDUCT2) FD = 0.0
                CHDK1(I)=CHDK1(I)+FD
 250         CONTINUE
 260      CONTINUE
 270   CONTINUE

 1111  CONTINUE

C-----------------------------------------------------------------------
C     Dipole disk representing the far downstream hub
C-----------------------------------------------------------------------
      DO 310 I=1,NPWAKS
         CHDKDT1(I)=0.
310   CONTINUE

C....The hub disk is not needed for IHUB=1 or IHUB=2 b/c it's close 
C....far downstream (JY052401)
      IF(IHUB.EQ.1.OR.IHUB.EQ.2.OR.IHUB.EQ.6.OR.IHUB.EQ.7) GO TO 2111

C-----------------------------------------------------------------------
C     prepare the geometry of a dipole disk replacing the far downstream
C       hub
C     Clockwise viewing              (N,2) 3 *-----* 2 (N+1,2)
C       from upstrean                        |     |
C                                            |     |
C                                    (N,1) 4 *-----* 1 (N+1,1)
C
C-----------------------------------------------------------------------
      DO 300 N=1,NHBDK
         XGW(N,1,1)=XHDKDT(N+1,1)
         XGW(N,1,2)=YHDKDT(N+1,1)
         XGW(N,1,3)=ZHDKDT(N+1,1)
         XGW(N,2,1)=XHDKDT(N+1,2)
         XGW(N,2,2)=YHDKDT(N+1,2)
         XGW(N,2,3)=ZHDKDT(N+1,2)
         XGW(N,3,1)=XHDKDT(N,2)
         XGW(N,3,2)=YHDKDT(N,2)
         XGW(N,3,3)=ZHDKDT(N,2)
         XGW(N,4,1)=XHDKDT(N,1)
         XGW(N,4,2)=YHDKDT(N,1)
         XGW(N,4,3)=ZHDKDT(N,1)
300   CONTINUE
      CALL GEO3DW(NHBDK,XGW,CHRLEWS,IER)
      IF(IER.EQ.0) THEN
         WRITE(*,'(A)') ' UNACCEPTABLE PANELS IN INFDISK!'
         STOP
      END IF

C-----------------------------------------------------------------------
C     Calculate the influence coefficients due to the hub disk
C       CHDK(I): induced potentials at I due to the hub disk
C-----------------------------------------------------------------------
      DO 370 N=1,NHBDK
C
C........Transfer data to common block /GEOM/
         DO 320 K=1,4
            XV(K)=XVPW(N,K)
            YV(K)=YVPW(N,K)
            SIDE(K)=SIDW(N,K)
320      CONTINUE
         DO 330 K=1,15
            S(K)=SSW(N,K)
330      CONTINUE
          DO 360 I=1,NPWAKS
             DO 350 KK=1,NBLADE
C
C..............Transfer control points to local coordinate
                XLOC=0.
                YLOC=0.
                ZLOC=0.
                DO 340 K=1,3
                   XLOC=XLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,1,K)
                   YLOC=YLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,2,K)
                   ZLOC=ZLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,3,K)
 340            CONTINUE
C
C..............Compute induced potentials
                IMR=1
                CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(N),FS,FD,
     1               FSX,FSY,FDX,FDY,FDZ,0,IMR)
                IF(IDUCT .EQ. 1 .AND. IDOPT .EQ. 1 
     %             .AND. I .GT. IDDUCT2 .AND. I .LE. IDDUCT2) FD = 0.0
                CHDKDT1(I)=CHDKDT1(I)+FD
 350         CONTINUE
 360      CONTINUE
 370   CONTINUE

 2111  CONTINUE

       RETURN
C))))))))))))))))))))) End of subroutine INFDISKSC (((((((((((((((((((((
       END                                                             
  
