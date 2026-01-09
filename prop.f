C----------------------------------------------------------------------------
      SUBROUTINE PROPEXT(XPIN,YPIN,ZPIN,NP,ISPC,XPOUT,YPOUT,ZPOUT,
     %                   XIN,YIN,XYCUB,NH,IOPT)
C----------------------------------------------------------------------------
      DIMENSION XIN(*),YIN(*),XYCUB(*)
      DIMENSION XPIN(*),YPIN(*),ZPIN(*)
      DIMENSION XPOUT(*),YPOUT(*),ZPOUT(*)
      DIMENSION XCOEF(4000),YCOEF(4000),ZCOEF(4000)
C      DIMENSION SS1(1000),SS2(1000)
      DIMENSION ps(1000)
      DIMENSION XDUM(1000),YDUM(1000),ZDUM(1000)

      CALL D3INTERPOL(XPIN,YPIN,ZPIN,NP,PS,XCOEF,YCOEF,ZCOEF)
      CALL INTERSECT(XPIN,YPIN,ZPIN,NP,XCOEF,YCOEF,ZCOEF,PS,
     %               XSH,YSH,ZSH,XIN,YIN,XYCUB,NH,IOPT)
                  
      CALL D3XTRAPOL(XPIN,YPIN,ZPIN,NP,XSH,YSH,ZSH,XDUM,YDUM,
     %               ZDUM,NN2,IOPT)
      CALL TPBMOD(XDUM,YDUM,ZDUM,NN2,ISPC,XPOUT,YPOUT,ZPOUT,NP)
      
      RETURN
      END
      
C----------------------------------------------------------------------------
      SUBROUTINE INTERSECT(XPIN,YPIN,ZPIN,NP,XCOEF,YCOEF,ZCOEF,PS,
     %                     XOUT,YOUT,ZOUT,XHIN,YHIN,XYHINCUB,NH,IOPT)
C
C     xpin,ypin,zpin      : x,y,z coordinates of blade (radial direction at given 
C                           chordwise node number)to be extrapolated to find the 
C                           intersection points with duct or hub 
C     np                  : No. of points of xpin, ypin, and zpin
C     xcoef, ycoef, zcoef : spline coeficients calculated from UGLYDK 
C                           based on PS and (xpin, ypin, zpin)
C     PS                  : Arclength calculated from xpin,ypin,zpin 
C                           with dimension of NP
C     xout,yout,zout      : coordinates of intersection points with 
C                           duct or hub surface.
C     IOPT                : 0 for HUB, 1 for DUCT                      
C---------------------------------------------------------------------------
      DIMENSION XHIN(*),YHIN(*),XYHINCUB(*)     
      DIMENSION XPIN(*),YPIN(*),ZPIN(*),PS(*)
      DIMENSION XCOEF(*),YCOEF(*),ZCOEF(*)

      EPS = 1.E-6
    
      IF(IOPT .EQ. 1) THEN
         DO I = NP, 1, -1
            CALL EVALDKs(NH,1,XHIN,XPIN(I),YY,XYHINCUB)
            RR = SQRT(YPIN(I)**2 + ZPIN(I)**2)
            
            DIST = RR - YY

            IF(DIST .LT. 0.0) THEN
               IF(I .EQ. NP) THEN
                  IN = NP - 1
               ELSE
                  IN = I
               ENDIF
               I2N = IN + (NP-1)
               I3N = I2N + (NP-1)
               I4N = I3N + (NP-1)
               
               T1 = PS(IN)
               T2 = PS(IN+1)
            
               a1=xcoef(in)
               a2=xcoef(i2n)
               a3=xcoef(i3n)
               a4=xcoef(i4n)
            
               b1=ycoef(in)
               b2=ycoef(i2n)
               b3=ycoef(i3n)
               b4=ycoef(i4n)
            
               c1=zcoef(in)
               c2=zcoef(i2n)
               c3=zcoef(i3n)
               c4=zcoef(i4n)

               GO TO 1000
            ELSEIF(DIST .EQ. 0.0) THEN
               XOUT = XPIN(I)
               YOUT = YPIN(I)
               ZOUT = ZPIN(I)
               
               RETURN
            ENDIF
         ENDDO
         WRITE(*,*) ' Cannot find intersection point !! '
         STOP

C     -- HUB
         
      ELSE

         DO I = 1, NP
            CALL EVALDKs(NH,1,XHIN,XPIN(I),YY,XYHINCUB)
            RR = SQRT(YPIN(I)**2 + ZPIN(I)**2)
            
            DIST = RR - YY 
            
C     -- If Dist > 0,  Vertex point is located ABOVE the HUB surface, 
C     -- If Dist < 0, Find the location where Dist > 0.

            IF(DIST .GT. 0.0) THEN
               IN = I
               I2N = IN + (NP-1)
               I3N = I2N + (NP-1)
               I4N = I3N + (NP-1)
               
               T1 = PS(IN)
               T2 = PS(IN+1)
               
               a1=xcoef(in)
               a2=xcoef(i2n)
               a3=xcoef(i3n)
               a4=xcoef(i4n)
               
               b1=ycoef(in)
               b2=ycoef(i2n)
               b3=ycoef(i3n)
               b4=ycoef(i4n)
               
               c1=zcoef(in)
               c2=zcoef(i2n)
               c3=zcoef(i3n)
               c4=zcoef(i4n)

               GO TO 1000
            ELSEIF(DIST .EQ. 0.0) THEN
               XOUT = XPIN(I)
               YOUT = YPIN(I)
               ZOUT = ZPIN(I)
               
               RETURN
            ENDIF
         ENDDO
         WRITE(*,*) ' Cannot find intersection point !! '
         STOP
      
      ENDIF

 1000 CONTINUE

      IT = 1

 1500 DELT = IT *(T2 - T1)

      IF(IOPT .EQ. 1) THEN
         T = DELT
      ELSE
         T = T2 - DELT
      ENDIF

      X = A1*T**3 + A2*T*T + A3*T + A4
      Y = B1*T**3 + B2*T*T + B3*T + B4
      Z = C1*T**3 + C2*T*T + C3*T + C4
      
      CALL EVALDKs(NH,1,XHIN,X,YY,XYHINCUB)
      IF(IOPT .EQ. 1) THEN
         DIST = YY - SQRT(Y*Y + Z*Z)
      ELSE
         DIST = SQRT(Y*Y + Z*Z) - YY
      ENDIF

      IF(DIST .LE. 0.0) THEN
         GO TO 2000
      ELSE
         IT = IT + 1
         GO TO 1500
      ENDIF
 
 2000 CONTINUE

      I = 0

      IF(IOPT .EQ. 1) THEN
         A = 0.0
         B = A + REAL(IT)*(T2-T1)
      ELSE
         A = T2 - REAL(IT-1)*(T2-T1)
         B = T2 - IT*(T2 - T1)
      ENDIF

      XB = A1*B**3 + A2*B*B + A3*B + A4
      YB = B1*B**3 + B2*B*B + B3*B + B4
      ZB = C1*B**3 + C2*B*B + C3*B + C4

      XA = A1*A**3 + A2*A*A + A3*A + A4
      YA = B1*A**3 + B2*A*A + B3*A + B4
      ZA = C1*A**3 + C2*A*A + C3*A + C4
      
      CALL EVALDKs(NH,1,XHIN,XB,YY,XYHINCUB)

      IF(IOPT .EQ. 1) THEN
         DIST2 = YY - SQRT(YB*YB + ZB*ZB)
      ELSE
         DIST2 = SQRT(YB*YB + ZB*ZB) - YY
      ENDIF

      CALL EVALDKs(NH,1,XHIN,XA,YY,XYHINCUB)

      IF(IOPT .EQ. 1) THEN
         DIST1 = YY - SQRT(YA*YA + ZA*ZA)
      ELSE
         DIST1 = SQRT(YA*YA + ZA*ZA) - YY
      ENDIF

 1100 I = I + 1

      T = 0.5 * (A+B)

      X = A1*T**3 + A2*T*T + A3*T + A4
      Y = B1*T**3 + B2*T*T + B3*T + B4
      Z = C1*T**3 + C2*T*T + C3*T + C4

      CALL EVALDKs(NH,1,XHIN,X,YY,XYHINCUB)
    

      IF(IOPT .EQ. 1) THEN
         DIST = YY - SQRT(Y*Y+Z*Z) 
      ELSE
         DIST = SQRT(Y*Y+Z*Z) - YY
      ENDIF

      IF(ABS(DIST) .LE. 1.E-6) GO TO 1200
 
      IF(DIST1*DIST .LT. 0.0) THEN
         B = T
         DIST2 = DIST
         GO TO 1100
      ELSEIF(DIST1*DIST .GT. 0.0) THEN
         A = T
         DIST1 = DIST
         Go TO 1100
      ENDIF

 1200 XOUT = X
      YOUT = Y
      ZOUT = Z

      RETURN
      END
      
C--------------------------------------------------------------------------
      SUBROUTINE D3INTERPOL(XIN,YIN,ZIN,NIN,SS1,XCOEF,YCOEF,ZCOEF)
C
C     Calculate Spline coeff. for 3-D geometry
C
C     (SS1 .vs. XIN,YIN,ZIN) of No. of data points NIN
C     XCOEF,YCOEF,ZCOEF : spline coeff. of xin,yin,zin
C
C--------------------------------------------------------------------------
        DIMENSION XIN(NIN),YIN(NIN),ZIN(NIN),SS1(NIN)
        DIMENSION XCOEF(4*(NIN-1)),YCOEF(4*(NIN-1)),ZCOEF(4*(NIN-1))

c---input format
        
        DO I = 1, NIN
           SS1(I) = SQRT(YIN(I)**2+ZIN(I)**2)
        ENDDO

        CALL UGLYDK(NIN,1,1,SS1,XIN,0.0,0.0,XCOEF)
        CALL UGLYDK(NIN,1,1,SS1,YIN,0.0,0.0,YCOEF)
        CALL UGLYDK(NIN,1,1,SS1,ZIN,0.0,0.0,ZCOEF)

        RETURN
        END
C----------------------------------------------------------------------------
      SUBROUTINE D3XTRAPOL(XPIN,YPIN,ZPIN,M,X,Y,Z,
     %                          XOUT,YOUT,ZOUT,NN,IOPT)
C----------------------------------------------------------------------------
C----gets prop.geo,duct or hub with boolean true and false resp.
C-----pts of Xscn with hub and duct....
      DIMENSION XPIN(*),YPIN(*),ZPIN(*)
      DIMENSION XOUT(*),YOUT(*),ZOUT(*)


C -- Re-arrange input blade data to include the intersection point.
      IF(IOPT .EQ. 0) THEN
         K = 1
         XOUT(1) = X
         YOUT(1) = Y
         ZOUT(1) = Z
         DO I = 1 , M
            RR = SQRT(YPIN(I)**2+ZPIN(I)**2)
            RR2 = SQRT(Y*Y + Z*Z)
            IF(RR .GT. RR2) THEN
               K = K + 1
               XOUT(K) = XPIN(I)
               YOUT(K) = YPIN(I)
               ZOUT(K) = ZPIN(I)
            ENDIF
         ENDDO         
         NN = K

         GO TO 100

      ELSEIF(IOPT .EQ. 1) THEN
         K = 0
         DO I = 1 , M
            RR = SQRT(YPIN(I)**2+ZPIN(I)**2)
            RR2 = SQRT(Y*Y + Z*Z)
            IF(RR .LE. RR2) THEN
               K = K + 1
               XOUT(K) = XPIN(I)
               YOUT(K) = YPIN(I)
               ZOUT(K) = ZPIN(I)
               NN = K
               IF(RR .EQ. RR2) THEN
                  GO TO 100
               ELSE
                  IF(I .EQ. M) THEN
                     K = K + 1
                     XOUT(K) = X
                     YOUT(K) = Y
                     ZOUT(K) = Z
                     NN = K
                     GO TO 100
                  ENDIF
               ENDIF
            ELSE
               K = K + 1
               XOUT(K) = X
               YOUT(K) = Y
               ZOUT(K) = Z
               NN = K
               GO TO 100
            ENDIF
         ENDDO

      ENDIF

 100  CONTINUE

      RETURN
      END

C-------------------------------------------------------------------
      SUBROUTINE TPBMOD(XPIN,YPIN,ZPIN,NIN1,ISPC,XOUT,YOUT,ZOUT,MP)
      
      DIMENSION XPIN(*),YPIN(*),ZPIN(*)
      DIMENSION XOUT(*),YOUT(*),ZOUT(*)
      DIMENSION SR(1000),SRM(1000)
      DIMENSION XCOEF(4000),YCOEF(4000),ZCOEF(4000)
C------------------------------------------------------------------
C     CALCULATING RADIUS
      PI = ACOS(-1.)

      DO M = 1,NIN1      
        SR(M) = SQRT(YPIN(M)*YPIN(M)+ZPIN(M)*ZPIN(M))
      END DO
     
C     FITTING CURVES
      CALL UGLYDK(NIN1,1,1,SR,XPIN,0.0,0.0,XCOEF)
      CALL UGLYDK(NIN1,1,1,SR,YPIN,0.0,0.0,YCOEF)
      CALL UGLYDK(NIN1,1,1,SR,ZPIN,0.0,0.0,ZCOEF)

      IF(ISPC.EQ.0) THEN
         DTH=PI/REAL(MP-1)
         DO 5 M=1,MP
            SRM(M)=SR(1)+0.5*(SR(NIN1)-SR(1))*(1.0-COS(DTH*(M-1)))
5        CONTINUE
      ELSE IF(ISPC.EQ.1) THEN 
         DELR=(SR(NIN1)-SR(1))/(MP-1)
         DO 10 M=1,MP
            SRM(M)=SR(1)+DELR*(M-1)
10       CONTINUE
      ELSE IF(ISPC.EQ.2) THEN
         DTH=0.5*PI/REAL(MP-1)
         DO 20 M=1,MP
            SRM(M)=SR(1)+(SR(NIN1)-SR(1))*SIN(DTH*(M-1))
20       CONTINUE
      END IF
C--------------------------------------------------------------
C     INTERPOLATING 
      CALL EVALDK(NIN1,MP,SR,SRM,XOUT,XCOEF)
      CALL EVALDK(NIN1,MP,SR,SRM,YOUT,YCOEF)
      CALL EVALDK(NIN1,MP,SR,SRM,ZOUT,ZCOEF)

      RETURN 
      END
      
      
