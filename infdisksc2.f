      SUBROUTINE INFDISKSC2
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVC.INC'

C-----------------------------------------------------------------------
C     Dipole disk representing the far upstream hub
C-----------------------------------------------------------------------

       DO 210 I=1,NPWAKS
          CHDK1(I)=0.
 210   CONTINUE

       IF(IHUB.EQ.0.OR.IHUB.EQ.2.OR.IHUB.EQ.4
     *             .OR.IHUB.EQ.5.OR.IHUB.EQ.6) GO TO 1111 ! People sometimes forget to put this crucial safety if statement... | S.N.KIM, Aug. 2018.

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
       DO 270 N=1,NHBDK
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
                XLOC=0.
                YLOC=0.
                ZLOC=0.
                DO 240 K=1,3
                   XLOC=XLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,1,K)
                   YLOC=YLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,2,K)
                   ZLOC=ZLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,3,K)
 240            CONTINUE
                IMR=1
                CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(N),FS,FD,
     1               FSX,FSY,FDX,FDY,FDZ,0,IMR)
                CHDK1(I)=CHDK1(I)+FD
 250         CONTINUE
 260      CONTINUE
 270   CONTINUE

 1111 CONTINUE

      DO 310 I=1,NPWAKS
         CHDKDT1(I)=0.
310   CONTINUE

      IF(IHUB.EQ.1.OR.IHUB.EQ.2.OR.IHUB.EQ.6.OR.IHUB.EQ.7) GO TO 2111 !S.N.KIM | Aug. 2018.

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
      DO 370 N=1,NHBDK
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

                XLOC=0.
                YLOC=0.
                ZLOC=0.
                DO 340 K=1,3
                   XLOC=XLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,1,K)
                   YLOC=YLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,2,K)
                   ZLOC=ZLOC+(XCPW(I,K,KK)-XCTW(N,K))*DIRW(N,3,K)
 340            CONTINUE

                IMR=1
                CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(N),FS,FD,
     1               FSX,FSY,FDX,FDY,FDZ,0,IMR)
                CHDKDT1(I)=CHDKDT1(I)+FD
 350         CONTINUE
 360      CONTINUE
 370   CONTINUE
 2111  CONTINUE

       RETURN
       END                                                             
  
