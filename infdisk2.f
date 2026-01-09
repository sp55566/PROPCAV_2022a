      SUBROUTINE INFDISK2
      INCLUDE 'PUFCAV.INC'

C-----------------------------------------------------------------------
C     prepare the geometry of a dipole disk replacing the far upstream
C       hub
C     Clockwise viewing              (N,2) 2 *-----* 3 (N+1,2)
C       from upstrean                        |     |
C                                            |     |
C                                    (N,1) 1 *-----* 4 (N+1,1)
C
C-----------------------------------------------------------------------
      IDDUCT2 = NPANB + NPANH
      IDDUCT3 = NPANB + NPANH + NPAND

      DO 210 I=1,NPANEL
         CHDK(I)=0.
210   CONTINUE

      IF(IHUB.EQ.2.OR.IHUB.EQ.4.OR.IHUB.EQ.5.OR.IHUB.EQ.6) GO TO 1111
      IF(IHUB.EQ.1.AND.ITUN.EQ.1) GOTO 1111

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
         WRITE(*,'(A)') ' UNACCEPTABLE PANELS IN INFDISK!'
         STOP
      END IF

C-----------------------------------------------------------------------
C     Calculate the influence coefficients due to the hub disk
C       CHDK(I): induced potentials at I due to the hub disk
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
         DO 260 I=1,NPANEL
            DO 250 KK=1,NBLADE
C
C..............Transfer control points to local coordinate
               XLOC=0.
               YLOC=0.
               ZLOC=0.
               DO 240 K=1,3
                  XLOC=XLOC+(XCTP(I,K,KK)-XCTW(N,K))*DIRW(N,1,K)
                  YLOC=YLOC+(XCTP(I,K,KK)-XCTW(N,K))*DIRW(N,2,K)
                  ZLOC=ZLOC+(XCTP(I,K,KK)-XCTW(N,K))*DIRW(N,3,K)
240            CONTINUE
C
C..............Compute induced potentials
               IMR=1
               CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(N),FS,FD,
     1                   FSX,FSY,FDX,FDY,FDZ,0,IMR)
               IF(IDUCT .EQ. 1 .AND. IDOPT .EQ. 1 
     %             .AND. I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) FD = 0.0
 
               CHDK(I)=CHDK(I)+FD
250         CONTINUE
260      CONTINUE
270   CONTINUE

 1111 CONTINUE

C-----------------------------------------------------------------------
C     Dipole disk representing the far downstream hub
C-----------------------------------------------------------------------
      DO 310 I=1,NPANEL
         CHDKDT(I)=0.
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
         DO 360 I=1,NPANEL
            DO 350 KK=1,NBLADE
C
C..............Transfer control points to local coordinate
               XLOC=0.
               YLOC=0.
               ZLOC=0.
               DO 340 K=1,3
                  XLOC=XLOC+(XCTP(I,K,KK)-XCTW(N,K))*DIRW(N,1,K)
                  YLOC=YLOC+(XCTP(I,K,KK)-XCTW(N,K))*DIRW(N,2,K)
                  ZLOC=ZLOC+(XCTP(I,K,KK)-XCTW(N,K))*DIRW(N,3,K)
340            CONTINUE
C
C..............Compute induced potentials
               IMR=1
               CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(N),FS,FD,
     1                   FSX,FSY,FDX,FDY,FDZ,0,IMR)
               IF(IDUCT .EQ. 1 .AND. IDOPT .EQ. 1 
     %             .AND. I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) FD = 0.0
               CHDKDT(I)=CHDKDT(I)+FD
350         CONTINUE
360      CONTINUE
370   CONTINUE

 2111 CONTINUE

      RETURN
      END                                                               
