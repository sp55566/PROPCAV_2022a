c-----------------------------------------------------------------------
c       ccc
ccc    ccccccc
ccc  cccccccccccc
cccc ccccccccccc               +++++++
ccccccccccc                +++++++++++++++                         CC C
cccccccccccc             ++++++++++++++++++++                 CC CCCCCC
cc ccc cccc             ++++++++++++++++++++++             CCCCCCCCCCCC
**  *  * *    __________________________________________  CCCCCCCCCCC**
**  *  *      I- Massachusetts Institute of Technology      CC* CCC*CCC
**  **        ------------------------------------------    CCCCC *CCC*
*  **           III    III    III    III    III    III        CCCCC CCC
* **            III    III    III    III    III    III            * **
**              III    III    III    III    III    III             * **
*               III    III    III    III    III    III               **
*               III    III    III    III    III    III               **
*               III    III    III    III    III    III               **
*             ==========================================             **
*            ============================================            **
**          ==============================================          ***
*
C---------------------------- MIT PUF-3A5 ------------------------------
C
C           COPYRIGHT (C) MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C                        VERSION 1.0   DEC. 1986
C                        VERSION 1.01  JULY 1990
C                        VERSION 1.2   JUNE 1993
C                        VERSION 1.3   FEBR 1994 (HPUF-3A, v. 1.0)
C---------------------------------------------------------------------
C                  ANALYSIS OF CAVITATING PROPELLERS
C                       IN NONUNIFORM INFLOW
C
C                  PROGRAMMER CHANG-SUP LEE (MAY 1979)
C                    REVISED BY SPYROS KINNAS (JUNE 1983)
C                               CHARLES CORRADO (JUNE 1986)
C                               KATY BUSHMAN (DEC 1986)
C
C                 INCLUDE HUB IMAGE EFFECTS - S.KINNAS (JAN. 1994)
C
C-----------------------------------------------------------------------
C
C ======================================================================
      SUBROUTINE VORSEG(XP,YP,ZP,X1,Y1,Z1,X2,Y2,Z2,VX,VY,VZ,SX,SY,SZ,K)
C ======================================================================
      COMMON/HIMAGE/JHUB,JDUCT
      COMMON/HIMAGE2/RRHUB
      COMMON/DIMAGE/DGAPSIZE,ENDFOIL

      CALL VSEG(XP,YP,ZP,X1,Y1,Z1,X2,Y2,Z2,VX,VY,VZ,SX,SY,SZ,K)

      IF(JHUB .EQ. 1) THEN
         RIMG = RRHUB
         IF(Y1*Y1+Z1*Z1 .NE. 0.0 .AND. Y2*Y2+Z2*Z2 .NE. 0.0) THEN
            CALL VIMSEG(XP,YP,ZP,X1,Y1,Z1,X2,Y2,Z2,RIMG,
     %           VXI,VYI,VZI,SXI,SYI,SZI,K)
            VX=VX+VXI
            VY=VY+VYI
            VZ=VZ+VZI
            SX=SX+SXI
            SY=SY+SYI
            SZ=SZ+SZI 
         ENDIF
      ENDIF

      IF(JDUCT .NE. 0)THEN
         RIMG = 1. + DGAPSIZE
C         IF(Y1*Y1+Z1*Z1 .NE. 0.0 .AND. Y2*Y2+Z2*Z2 .NE. 0.0) THEN
         IF(X1 .LE. ENDFOIL .AND. X2 .LE. ENDFOIL) THEN
            CALL VIMSEG(XP,YP,ZP,X1,Y1,Z1,X2,Y2,Z2,RIMG,
     %           VXI,VYI,VZI,SXI,SYI,SZI,K)
            VX=VX+VXI
            VY=VY+VYI
            VZ=VZ+VZI
            SX=SX+SXI
            SY=SY+SYI
            SZ=SZ+SZI 
         ENDIF
      ENDIF  

      RETURN
      END 


C ======================================================================
      SUBROUTINE VORSGN(XP,YP,ZP,X1,Y1,Z1,X2,Y2,Z2,XN,YN,ZN,VN,SN,K)
C ======================================================================
      COMMON/HIMAGE/JHUB,JDUCT
      COMMON/HIMAGE2/RRHUB
      COMMON/DIMAGE/DGAPSIZE,ENDFOIL

      CALL VSGN(XP,YP,ZP,X1,Y1,Z1,X2,Y2,Z2,XN,YN,ZN,VN,SN,K)

      IF(JHUB .EQ. 1) THEN
         RIMG = RRHUB
         IF(Y1*Y1+Z1*Z1 .NE. 0.0 .AND. Y2*Y2+Z2*Z2 .NE. 0.0) THEN
            CALL VIMSGN(XP,YP,ZP,X1,Y1,Z1,X2,Y2,Z2,RIMG,
     %           XN,YN,ZN,VNI,SNI,K)
            VN=VN+VNI
            SN=SN+SNI
         ENDIF
      ENDIF

      IF(JDUCT .NE. 0)THEN
         RIMG = 1. + DGAPSIZE
C         IF(Y1*Y1+Z1*Z1 .NE. 0.0 .AND. Y2*Y2+Z2*Z2 .NE. 0.0) THEN
         IF(X1 .LE. ENDFOIL .AND. X2 .LE. ENDFOIL) THEN
            CALL VIMSGN(XP,YP,ZP,X1,Y1,Z1,X2,Y2,Z2,RIMG,
     %           XN,YN,ZN,VNI,SNI,K)
            VN=VN+VNI
            SN=SN+SNI
         ENDIF
      ENDIF  

      RETURN
      END 


C ======================================================================
      SUBROUTINE VSEG(XP,YP,ZP,X1,Y1,Z1,X2,Y2,Z2,VX,VY,VZ,SX,SY,SZ,K)
C ======================================================================
      DOUBLE PRECISION E,V22
      DATA ZERO,TENTH,HALF,ONE,TWO,RADIUS/0.0,0.1,0.5,1.0,2.0,5.0/
C-----VORTEX(A),MEAN DISTANCE(F),AND NORMAL(H) VECTORS------------------
      AX=X2-X1
      AY=Y2-Y1
      AZ=Z2-Z1
      A=SQRT(AX**2+AY**2+AZ**2)
      FX=HALF*(X2+X1)-XP
      FY=HALF*(Y2+Y1)-YP
      FZ=HALF*(Z2+Z1)-ZP
      F=SQRT(FX**2+FY**2+FZ**2)
      HX=AZ*FY-AY*FZ
      HY=AX*FZ-AZ*FX
      HZ=AY*FX-AX*FY
      IF(F/A.LT.RADIUS) GO TO 1
C-----FAR FIELD APPROXIMATIONS------------------------------------------
      T=HALF/F**3
      VX=T*HX
      VY=T*HY
      VZ=T*HZ
      IF(K.EQ.0) GO TO 2
      T=-T*A
      SX=T*FX
      SY=T*FY
      SZ=T*FZ
      GO TO 2
C-----NEAR FIELD EXACT SOLUTION-----------------------------------------
 1    B=SQRT((X2-XP)**2+(Y2-YP)**2+(Z2-ZP)**2)
      C=SQRT((X1-XP)**2+(Y1-YP)**2+(Z1-ZP)**2)
      A2=HALF*A
      Q=F/A2
      IF(Q.GT.TENTH) GO TO 3
C-----SOLUTION CLOSE TO VORTEX TO AVIOD ROUNDOFF------------------------
      D=SQRT(HX**2+HY**2+HZ**2)/A
      E=(A**2+C**2-B**2)/(TWO*A)
      V22=((A-E)/B+E/C)/(TWO*D)
      V2=V22
      GO TO 4
C-----SOLUTION BEYOND ROUNDOFF LIMIT------------------------------------
 3    S=HALF*(A2+B+F)
      ARGMNT=AMIN1(SQRT(ABS(((S-A2)*(S-B)/(A2*B)))),1.0)
      BETA=TWO*ASIN(ARGMNT)
      S=HALF*(A2+C+F)
      ARGMNT=AMIN1(SQRT(ABS(((S-A2)*(S-C)/(A2*C)))),1.0)
      ALPHA=TWO*ASIN(ARGMNT)
      CA=COS(ALPHA)
      D=C*SIN(ALPHA)
      IF(D.GT.1.0E-5) GO TO 5
      V2=ZERO
      D=ONE
      GO TO 4
 5    V2=(CA+COS(BETA))/(TWO*D)
 4    T=V2/(A*D)
      VX=T*HX
      VY=T*HY
      VZ=T*HZ
      IF(K.EQ.0) GO TO 2
      T=V2/D
      W=(ONE/B-ONE/C)/(TWO*A)
      IF(Q.LT.TENTH) W=ZERO
      EPS=C*CA/A
      IF(Q.LT.TENTH) EPS=HALF
      DX=X1+EPS*AX-XP
      DY=Y1+EPS*AY-YP
      DZ=Z1+EPS*AZ-ZP
      SX=W*AX-T*DX
      SY=W*AY-T*DY
      SZ=W*AZ-T*DZ
 2    RETURN
      END


C ===================================================================
      SUBROUTINE VSGN(XP,YP,ZP,X1,Y1,Z1,X2,Y2,Z2,XN,YN,ZN,VN,SN,K)
C ===================================================================
      DOUBLE PRECISION E,V22
c-----s.a.kinnas---4/26/2000------------------------------------------
c     change TENTH to RNEAR and RADIUS to RFAR
c      DATA ZERO,TENTH,HALF,ONE,TWO,RADIUS/0.0,0.1,0.1666667,0.5,
c     1  1.0,2.0,5.0/
      DATA ZERO,RNEAR,HALF,ONE,TWO,RFAR/0.0,0.0,0.5,1.0,2.0,50.0/
c------------------------4/26/2000------------------------------------
C-----VORTEX(A),MEAN DISTANCE(F),AND NORMAL(H) VECTORS------------------
      AX=X2-X1
      AY=Y2-Y1
      AZ=Z2-Z1
      A=SQRT(AX**2+AY**2+AZ**2)
      FX=HALF*(X2+X1)-XP
      FY=HALF*(Y2+Y1)-YP
      FZ=HALF*(Z2+Z1)-ZP
      F=SQRT(FX**2+FY**2+FZ**2)
      HX=AZ*FY-AY*FZ
      HY=AX*FZ-AZ*FX
      HZ=AY*FX-AX*FY
      IF(F/A.LT.RFAR) GO TO 1
C-----FAR FIELD APPROXIMATIONS------------------------------------------
      T=HALF/F**3
      VX=T*HX
      VY=T*HY
      VZ=T*HZ
      VN=XN*VX+YN*VY+ZN*VZ
      IF(K.EQ.0) GO TO 2
      T=-T*A
      SX=T*FX
      SY=T*FY
      SZ=T*FZ
      SN=XN*SX+YN*SY+ZN*SZ
      GO TO 2
C-----NEAR FIELD EXACT SOLUTION-----------------------------------------
 1    B=SQRT((X2-XP)**2+(Y2-YP)**2+(Z2-ZP)**2)
      C=SQRT((X1-XP)**2+(Y1-YP)**2+(Z1-ZP)**2)
      A2=HALF*A
      Q=F/A2
      IF(Q.GT.RNEAR) GO TO 3
C-----SOLUTION CLOSE TO VORTEX TO AVIOD ROUNDOFF------------------------
      D=SQRT(HX**2+HY**2+HZ**2)/A
      E=(A**2+C**2-B**2)/(TWO*A)
      V22=((A-E)/B+E/C)/(TWO*D)
      V2=V22
      GO TO 4
C-----SOLUTION BEYOND ROUNDOFF LIMIT------------------------------------
 3    S=HALF*(A2+B+F)
      ARGMNT=AMIN1(SQRT(ABS(((S-A2)*(S-B)/(A2*B)))),1.0)
      BETA=TWO*ASIN(ARGMNT)
      S=HALF*(A2+C+F)
      ARGMNT=AMIN1(SQRT(ABS(((S-A2)*(S-C)/(A2*C)))),1.0)
      ALPHA=TWO*ASIN(ARGMNT)
      CA=COS(ALPHA)
      D=C*SIN(ALPHA)
      IF(D.GT.1.0E-5) GO TO 5
      V2=ZERO
      D=ONE
      GO TO 4
 5    V2=(CA+COS(BETA))/(TWO*D)
 4    T=V2/(A*D)
      VX=T*HX
      VY=T*HY
      VZ=T*HZ
      VN=XN*VX+YN*VY+ZN*VZ
      IF(K.EQ.0) GO TO 2
      T=V2/D
      W=(ONE/B-ONE/C)/(TWO*A)
      IF(Q.LT.RNEAR) W=ZERO
      EPS=C*CA/A
      IF(Q.LT.RNEAR) EPS=HALF
      DX=X1+EPS*AX-XP
      DY=Y1+EPS*AY-YP
      DZ=Z1+EPS*AZ-ZP
      SX=W*AX-T*DX
      SY=W*AY-T*DY
      SZ=W*AZ-T*DZ
      SN=XN*SX+YN*SY+ZN*SZ
 2    RETURN
      END


C ================================================================
      SUBROUTINE SSEG(XP,YP,ZP,X1,Y1,Z1,X2,Y2,Z2,SX,SY,SZ,PHI,K)
C ================================================================
      A1=X2-X1
      A2=Y2-Y1
      A3=Z2-Z1
      A=SQRT(A1**2+A2**2+A3**2)
      B=SQRT((X2-XP)**2+(Y2-YP)**2+(Z2-ZP)**2)
      C=SQRT((X1-XP)**2+(Y1-YP)**2+(Z1-ZP)**2)
      V1=1.0/B-1.0/C
      E=(A**2+C**2-B**2)/(2.0*A)
      BUG=C**2-E**2
      IF(BUG.LE.0.0) GO TO 2
      D=SQRT(BUG)
      IF(E.GE.0.0.AND.E.LE.A) GO TO 3
      TOR=ABS(E)*0.03
      IF(D.GE.TOR) GO TO 3
 2    V2=ABS(A*(A-2.0*E)/(4.0*E**2*(A-E)**2))
      GO TO 1
 3    V2=((A-E)/B+E/C)/(2.0*D**2)
 1    EPS=E/A
      XQ=X1+EPS*A1
      YQ=Y1+EPS*A2
      ZQ=Z1+EPS*A3
      SX=V1*A1/A/2.0-V2*(XQ-XP)
      SY=V1*A2/A/2.0-V2*(YQ-YP)
      SZ=V1*A3/A/2.0-V2*(ZQ-ZP)
      IF(K.EQ.0) GO TO 99
      PHI=-0.07957747*ALOG((A-E+B)/(-E+C))
 99   RETURN
      END

C
C =======================================================================
      SUBROUTINE VIMSEG(XP,YP,ZP,XX1,YY1,ZZ1,XX2,YY2,ZZ2,RIMAGE,
     %                  VX,VY,VZ,SX,SY,SZ,K)
C =======================================================================
      DOUBLE PRECISION E,V22
      DATA ZERO,TENTH,HALF,ONE,TWO,RADIUS/0.0,0.1,0.5,1.0,2.0,5.0/

      IF(RIMAGE .EQ. 0) RETURN

      X2=XX1
      FF1=RIMAGE**2/(YY1**2+ZZ1**2)
      Y2=YY1*FF1
      Z2=ZZ1*FF1
C
      X1=XX2
      FF2=RIMAGE**2/(YY2**2+ZZ2**2)
      Y1=YY2*FF2
      Z1=ZZ2*FF2
C
C------1/10/94 ---------------------------------------------
      AX=X2-X1
      AY=Y2-Y1
      AZ=Z2-Z1
      A=SQRT(AX**2+AY**2+AZ**2)
      FX=HALF*(X2+X1)-XP
      FY=HALF*(Y2+Y1)-YP
      FZ=HALF*(Z2+Z1)-ZP
      F=SQRT(FX**2+FY**2+FZ**2)
      HX=AZ*FY-AY*FZ
      HY=AX*FZ-AZ*FX
      HZ=AY*FX-AX*FY
      IF(F/A.LT.RADIUS) GO TO 1
C-----FAR FIELD APPROXIMATIONS------------------------------------------
      T=HALF/F**3
      VX=T*HX
      VY=T*HY
      VZ=T*HZ
      IF(K.EQ.0) GO TO 2
      T=-T*A
      SX=T*FX
      SY=T*FY
      SZ=T*FZ
      GO TO 2
C-----NEAR FIELD EXACT SOLUTION-----------------------------------------
 1    B=SQRT((X2-XP)**2+(Y2-YP)**2+(Z2-ZP)**2)
      C=SQRT((X1-XP)**2+(Y1-YP)**2+(Z1-ZP)**2)
      A2=HALF*A
      Q=F/A2
      IF(Q.GT.TENTH) GO TO 3
C-----SOLUTION CLOSE TO VORTEX TO AVIOD ROUNDOFF------------------------
      D=SQRT(HX**2+HY**2+HZ**2)/A
      E=(A**2+C**2-B**2)/(TWO*A)
      V22=((A-E)/B+E/C)/(TWO*D)
      V2=V22
      GO TO 4
C-----SOLUTION BEYOND ROUNDOFF LIMIT------------------------------------
 3    S=HALF*(A2+B+F)
      ARGMNT=AMIN1(SQRT(ABS(((S-A2)*(S-B)/(A2*B)))),1.0)
      BETA=TWO*ASIN(ARGMNT)
      S=HALF*(A2+C+F)
      ARGMNT=AMIN1(SQRT(ABS(((S-A2)*(S-C)/(A2*C)))),1.0)
      ALPHA=TWO*ASIN(ARGMNT)
      CA=COS(ALPHA)
      D=C*SIN(ALPHA)
      IF(D.GT.1.0E-5) GO TO 5
      V2=ZERO
      D=ONE
      GO TO 4
 5    V2=(CA+COS(BETA))/(TWO*D)
 4    T=V2/(A*D)
      VX=T*HX
      VY=T*HY
      VZ=T*HZ
      IF(K.EQ.0) GO TO 2
      T=V2/D
      W=(ONE/B-ONE/C)/(TWO*A)
      IF(Q.LT.TENTH) W=ZERO
      EPS=C*CA/A
      IF(Q.LT.TENTH) EPS=HALF
      DX=X1+EPS*AX-XP
      DY=Y1+EPS*AY-YP
      DZ=Z1+EPS*AZ-ZP
      SX=W*AX-T*DX
      SY=W*AY-T*DY
      SZ=W*AZ-T*DZ
 2    RETURN
      END


C =======================================================================
      SUBROUTINE VIMSGN(XP,YP,ZP,XX1,YY1,ZZ1,XX2,YY2,ZZ2,RIMAGE,
     %                  XN,YN,ZN,VN,SN,K)
C =======================================================================
C
      DOUBLE PRECISION E,V22
      DATA ZERO,TENTH,HALF,ONE,TWO,RADIUS/0.0,0.1,0.5,1.0,2.0,5.0/

      IF(RIMAGE .EQ. 0) RETURN

      X2=XX1
      FF1=RIMAGE**2/(YY1**2+ZZ1**2)
      Y2=YY1*FF1
      Z2=ZZ1*FF1
C
      X1=XX2
      FF2=RIMAGE**2/(YY2**2+ZZ2**2)
      Y1=YY2*FF2
      Z1=ZZ2*FF2
C      
C------1/10/94 ---------------------------------------------
      AX=X2-X1
      AY=Y2-Y1
      AZ=Z2-Z1
      A=SQRT(AX**2+AY**2+AZ**2)
      FX=HALF*(X2+X1)-XP
      FY=HALF*(Y2+Y1)-YP
      FZ=HALF*(Z2+Z1)-ZP
      F=SQRT(FX**2+FY**2+FZ**2)
      HX=AZ*FY-AY*FZ
      HY=AX*FZ-AZ*FX
      HZ=AY*FX-AX*FY
      IF(F/A.LT.RADIUS) GO TO 1
C-----FAR FIELD APPROXIMATIONS------------------------------------------
      T=HALF/F**3
      VX=T*HX
      VY=T*HY
      VZ=T*HZ
      VN=XN*VX+YN*VY+ZN*VZ
      IF(K.EQ.0) GO TO 2
      T=-T*A
      SX=T*FX
      SY=T*FY
      SZ=T*FZ
      SN=XN*SX+YN*SY+ZN*SZ
      GO TO 2
C-----NEAR FIELD EXACT SOLUTION-----------------------------------------
 1    B=SQRT((X2-XP)**2+(Y2-YP)**2+(Z2-ZP)**2)
      C=SQRT((X1-XP)**2+(Y1-YP)**2+(Z1-ZP)**2)
      A2=HALF*A
      Q=F/A2
      IF(Q.GT.TENTH) GO TO 3
C-----SOLUTION CLOSE TO VORTEX TO AVIOD ROUNDOFF------------------------
      D=SQRT(HX**2+HY**2+HZ**2)/A
      E=(A**2+C**2-B**2)/(TWO*A)
      V22=((A-E)/B+E/C)/(TWO*D)
      V2=V22
      GO TO 4
C-----SOLUTION BEYOND ROUNDOFF LIMIT------------------------------------
 3    S=HALF*(A2+B+F)
      ARGMNT=AMIN1(SQRT(ABS(((S-A2)*(S-B)/(A2*B)))),1.0)
      BETA=TWO*ASIN(ARGMNT)
      S=HALF*(A2+C+F)
      ARGMNT=AMIN1(SQRT(ABS(((S-A2)*(S-C)/(A2*C)))),1.0)
      ALPHA=TWO*ASIN(ARGMNT)
      CA=COS(ALPHA)
      D=C*SIN(ALPHA)
      IF(D.GT.1.0E-5) GO TO 5
      V2=ZERO
      D=ONE
      GO TO 4
 5    V2=(CA+COS(BETA))/(TWO*D)
 4    T=V2/(A*D)
      VX=T*HX
      VY=T*HY
      VZ=T*HZ
      VN=XN*VX+YN*VY+ZN*VZ
      IF(K.EQ.0) GO TO 2
      T=V2/D
      W=(ONE/B-ONE/C)/(TWO*A)
      IF(Q.LT.TENTH) W=ZERO
      EPS=C*CA/A
      IF(Q.LT.TENTH) EPS=HALF
      DX=X1+EPS*AX-XP
      DY=Y1+EPS*AY-YP
      DZ=Z1+EPS*AZ-ZP
      SX=W*AX-T*DX
      SY=W*AY-T*DY
      SZ=W*AZ-T*DZ
      SN=XN*SX+YN*SY+ZN*SZ
 2    RETURN
      END
