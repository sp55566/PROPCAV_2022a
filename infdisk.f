      SUBROUTINE INFDISK
************************************************************************
*     INFDISK: INFluence coefficients due to the wake disk and the hub
*              disk
*      --- Calculate the influence coefficients due to a source disk 
*          represents the ultimate wake and a dipole disk represents 
*          the far-upstream hub
*
*      JY111900    The hub disk is not needed for IHUB=2.              
*      JY052401    Modified subroutine to accomodate new hub options.
*                                                                      
************************************************************************

      INCLUDE 'PUFCAV.INC'
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
         WRITE(*,'(A)') ' UNACCEPTABLE PANELS IN INFDISK!'
         STOP
      END IF

C-----------------------------------------------------------------------
C     Calculate the influence coefficients due to the wake disk
C       WDK(I): induced potential at I due to the wake disk
C-----------------------------------------------------------------------
      DO 110 I=1,NPANEL
         WDK(I)=0.
110   CONTINUE

C/s S.N.KIM | When wake is read-in, source disk is ignored to be
C           | consistent with cavity runs and fully aligned wake. 
      IF(IREADW.EQ.1) GO TO 1004
C/e S.N.KIM | Nov. 2018.
      DO 170 L=1,NUWDK
C
C........Transfer data to common block /GEOM/
         DO 120 K=1,4
            XV(K)=XVPW(L,K)
            YV(K)=YVPW(L,K)
            SIDE(K)=SIDW(L,K)
120      CONTINUE
         DO 130 K=1,15
            S(K)=SSW(L,K)
130      CONTINUE
C
C........Calculate the influence coefficients
         DO 160 I=1,NPANEL
            DO 150 KK=1,NBLADE
C
C..............Transfer control points to the local coordinate
               XLOC=0.
               YLOC=0.
               ZLOC=0.
               DO 140 K=1,3
                  XLOC=XLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,1,K)
                  YLOC=YLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,2,K)
                  ZLOC=ZLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,3,K)
140             CONTINUE
C
C..............Compute induced potentials
               IMR=1
               CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(L),FS,FD,
     1                   FSX,FSY,FDX,FDY,FDZ,0,IMR)

               IF(IDUCT .EQ. 1 .AND. IDOPT .EQ. 1 
     %           .AND. I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                     FS = 0.0
               ENDIF
 
               WDK(I)=WDK(I)+FS
150         CONTINUE
160       CONTINUE
170    CONTINUE

1004   CONTINUE

cC------------------------------(S.H.CHANG 02/25/2010)-------------------------  
cC    OUTPUT GEOMETRIES ON THE BLADE AND HUB FOR HULLFPP 
c       WRITE(733,*) NUWDK
c       DO L = 1, NUWDK
c          DO M = 1,4
c             DO N = 1,3
c                WRITE(733,*) XGW(L,M,N)
c             END DO
c          END DO
c       END DO
cC------------------------------(S.H.CHANG 02/25/2010)------------------------- 

       IF(IDUCT .NE. 0) THEN
          DO N=1,MDUCT
             XGW(N,1,1)=XDWDK(N,1)
             XGW(N,1,2)=YDWDK(N,1)
             XGW(N,1,3)=ZDWDK(N,1)
             XGW(N,2,1)=XDWDK(N,2)
             XGW(N,2,2)=YDWDK(N,2)
             XGW(N,2,3)=ZDWDK(N,2)
             XGW(N,3,1)=XDWDK(N+1,2)
             XGW(N,3,2)=YDWDK(N+1,2)
             XGW(N,3,3)=ZDWDK(N+1,2)
             XGW(N,4,1)=XDWDK(N+1,1)
             XGW(N,4,2)=YDWDK(N+1,1)
             XGW(N,4,3)=ZDWDK(N+1,1)
          ENDDO

          CALL GEO3DW(MDUCT,XGW,CHRLEWS,IER)

          IF(IER.EQ.0) THEN
             WRITE(*,'(A)') ' UNACCEPTABLE PANELS IN INFDISK!'
             STOP
          END IF

          DO I=1,NPANEL
             WDKD(I)=0.0
          ENDDO

          DO L=1,MDUCT
             DO K=1,4
                XV(K)=XVPW(L,K)
                YV(K)=YVPW(L,K)
                SIDE(K)=SIDW(L,K)
             ENDDO
             DO K=1,15
                S(K)=SSW(L,K)
             ENDDO

             DO I=1,NPANEL
                DO KK=1,NBLADE
                   XLOC=0.
                   YLOC=0.
                   ZLOC=0.
                   DO K=1,3
                      XLOC=XLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,1,K)
                      YLOC=YLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,2,K)
                      ZLOC=ZLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,3,K)
                   ENDDO
                   IMR=1
                   CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(L),FS,FD,
     %                  FSX,FSY,FDX,FDY,FDZ,0,IMR)

                   IF(I .LE. IDDUCT2 .OR. I .GT. IDDUCT3) THEN
                      IF(IDOPT .EQ. 1) THEN
                         FS = 0.0
                      ENDIF
                   ENDIF

                   WDKD(I)=WDKD(I)+FS
                ENDDO
             ENDDO
          ENDDO
       ENDIF

C-----------------------------------------------------------------------
C     Dipole disk representing the far upstream hub
C-----------------------------------------------------------------------
      DO 210 I=1,NPANEL
         CHDK(I)=0.
210   CONTINUE

      IF(IHUB.EQ.2 .OR. IHUB.EQ.4 .OR. IHUB.EQ.5 .OR. IHUB.EQ.6) GOTO 1111
      IF(IHUB.EQ.1 .AND. ITUN.EQ.1) GOTO 1111

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
     %           .AND. I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) FD = 0.0

               CHDK(I)=CHDK(I)+FD
250         CONTINUE
260      CONTINUE
270   CONTINUE

 1111 CONTINUE

cC------------------------------(S.H.CHANG 02/25/2010)------------------------- 
cC    OUTPUT GEOMETRIES ON THE BLADE AND HUB FOR HULLFPP 
c      WRITE(734,*) NHBDK
c      DO L = 1,NHBDK
c         DO M = 1,4
c           DO N = 1,3
c              WRITE(734,*) XGW(L,M,N)
c           END DO
c         END DO
c      END DO
cC------------------------------(S.H.CHANG 02/25/2010)------------------------- 

C-----------------------------------------------------------------------
C     Dipole disk representing the far downstream hub
C-----------------------------------------------------------------------
      DO 310 I=1,NPANEL
         CHDKDT(I)=0.
310   CONTINUE

      IF(IHUB.EQ.1 .OR.IHUB.EQ.2 .OR.IHUB.EQ.6 .OR.IHUB.EQ.7) GOTO 2111

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

C..............Transfer control points to local coordinate
               XLOC=0.
               YLOC=0.
               ZLOC=0.
               DO 340 K=1,3
                  XLOC=XLOC+(XCTP(I,K,KK)-XCTW(N,K))*DIRW(N,1,K)
                  YLOC=YLOC+(XCTP(I,K,KK)-XCTW(N,K))*DIRW(N,2,K)
                  ZLOC=ZLOC+(XCTP(I,K,KK)-XCTW(N,K))*DIRW(N,3,K)
340            CONTINUE

C..............Compute induced potentials
               IMR=1
               CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(N),FS,FD,
     1                   FSX,FSY,FDX,FDY,FDZ,0,IMR)

               IF(IDUCT .EQ. 1 .AND. IDOPT .EQ. 1 
     %           .AND. I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) FD= 0.0

               CHDKDT(I)=CHDKDT(I)+FD
350         CONTINUE
360      CONTINUE
370   CONTINUE

 2111 CONTINUE

cC------------------------------(S.H.CHANG 02/25/2010)-------------------------
cC    OUTPUT GEOMETRIES ON THE WAKE DISK AND HUB DISK FOR HULLFPP    
c      WRITE(735,*) NHBDK
c      DO L = 1,NHBDK 
c         DO M = 1,4
c            DO N = 1,3
c               WRITE(735,*) XGW(L,M,N)
c            END DO
c         END DO
c      END DO 
cC------------------------------(S.H.CHANG 02/25/2010)-------------------------

      RETURN
C))))))))))))))))))))) End of subroutine INFDISK (((((((((((((((((((((((
      END                                                               
