C =============================================================
      SUBROUTINE DUCTGEO1
C
C     For IDGEO=0  Generate duct panel using propeller Pitch      
C =============================================================
      INCLUDE 'PUFCAV.INC'

      DIMENSION XBB1(200),YBB1(200),ZBB1(200)
      DIMENSION XBB2(200),YBB2(200),ZBB2(200)
      DIMENSION XID1(201),ETAD1(201),DCUBICX1(800)
      DIMENSION XID2(201),ETAD2(201),DCUBICX2(800)
      DIMENSION SSDF1(201),DCUBICS11(800),DCUBICS12(800)
      DIMENSION SSDF2(201),DCUBICS21(800),DCUBICS22(800)
      DIMENSION SSDX1(201),SSDR1(201),SSDT1(201),DCUBICS13(800)
      DIMENSION SSDX2(201),SSDR2(201),SSDT2(201),DCUBICS23(800)
      DIMENSION TBANGL(100),TBANGL2(100),TBANGT(100)
      DIMENSION XBK(200),YBK(200),ZBK(200)
      DIMENSION XTT(200),RTT(200),WTT(200),XRTTCUB(800),XWTTCUB(800)
      DIMENSION XXXXR(200)
      DIMENSION RDTRS(NDUCT+1,MDUCT+1),THTRS(NDUCT+1,MDUCT+1)

      NN1 = NDDAT/2 + 1  ! NN! = 100
      
      DO N = 1 , NN1
         N1 = NN1 - N + 1
         N2 = NN1 + N - 1
         XID1(N) = XID(N1)
         ETAD1(N) = ETAD(N1)
         
         XID2(N) = XID(N2)
         ETAD2(N) = ETAD(N2)
      ENDDO

      SSDF1(1) = 0.0
      SSDF2(1) = 0.0

      DO N = 2 , NN1
         P1 = SQRT( (XID1(N)-XID1(N-1))**2 
     %             +(ETAD1(N)-ETAD1(N-1))**2)
         SSDF1(N) = SSDF1(N-1) + P1

         P2 = SQRT( (XID2(N)-XID2(N-1))**2 
     %             +(ETAD2(N)-ETAD2(N-1))**2)
         SSDF2(N) = SSDF2(N-1) + P2

      ENDDO
     
      CALL UGLYDK(NN1,1,1,XID1,ETAD1,ZERO,ZERO,DCUBICX1)

      CALL UGLYDK(NN1,1,1,SSDF1,XID1,ZERO,ZERO,DCUBICS11)
      CALL UGLYDK(NN1,1,1,SSDF1,ETAD1,ZERO,ZERO,DCUBICS12)

      CALL UGLYDK(NN1,1,1,SSDF2,XID2,ZERO,ZERO,DCUBICS21)
      CALL UGLYDK(NN1,1,1,SSDF2,ETAD2,ZERO,ZERO,DCUBICS22)

      CALL UGLYDK(NN1,1,1,XID2,ETAD2,ZERO,ZERO,DCUBICX2)

      PPD = XPI(NX)
     
      DO N = 1 , NHP    ! NHP = NH + 1, NH = NC/2
         I1 = NHP - N + 1
         I2 = NHP + N - 1

         IF(DUCTGAP .EQ. 0.0) THEN
            XBB1(N) = XB(I1,MRP)
            YBB1(N) = YB(I1,MRP)
            ZBB1(N) = ZB(I1,MRP)
            
            XBB2(N) = XB(I2,MRP)
            YBB2(N) = YB(I2,MRP)
            ZBB2(N) = ZB(I2,MRP)
         ELSE
            XBB1(N) = 0.5*(XB(I1,MRP)+XB(I2,MRP))
            YBB1(N) = 0.5*(YB(I1,MRP)+YB(I2,MRP))
            ZBB1(N) = 0.5*(ZB(I1,MRP)+ZB(I2,MRP))

            XBB2(N) = XBB1(N)
            YBB2(N) = YBB1(N)
            ZBB2(N) = ZBB1(N)
         ENDIF
         XBK(N) = 0.5*(XB(I1,MRP)+XB(I2,MRP))
         YBK(N) = 0.5*(YB(I1,MRP)+YB(I2,MRP))
         ZBK(N) = 0.5*(ZB(I1,MRP)+ZB(I2,MRP))
      ENDDO

C -- Check No. of Panel on After part of duct ..

      DO N = 1 , NSW(MRP)
         IF(XID1(NN1) .GT. XW(N,MRP) .AND. 
     %          XID1(NN1) .LE. XW(N+1,MRP)) THEN
            NDAFTT =N 
         ENDIF
      ENDDO

      NDAF = NDAFTT + 1

      DO N = 1 , NDAF
         XTT(N) = XW(N,MRP)
         RTT(N) = SQRT(YW(N,MRP)**2+ZW(N,MRP)**2)
         WTT(N) = DANGLE(ZW(N,MRP),YW(N,MRP))
         IF(WTT(N) .LT. 0.0) WTT(N) = WTT(N) + TWOPI
      ENDDO

      DO N = 2 , NDAF
         IF(WTT(N-1) .GT. WTT(N)) THEN
            WTT(N) = WTT(N) + TWOPI
         ENDIF
      ENDDO

      CALL UGLYDK(NDAF,1,1,XTT,RTT,ZERO,ZERO,XRTTCUB)    
      CALL UGLYDK(NDAF,1,1,XTT,WTT,ZERO,ZERO,XWTTCUB)    

      NMID = NH  ! NH = NC/2
      NMIDP = NMID+1
      PITCH2 = 0.5 * PI - ATAN(PPD/PI)

      NDUCTH = NDUCT/2
!      NDAFT = (NDUCTH-NMID)/2
!      NDFWD = NDUCTH - NDAFT - NMID

!****** PR, Art Director : Seungnam, Kim****** 07/03/2016
      NDAFT = NDAFTT
      NDFWD = NDUCTH - NDAFT - NMID
!****** PR, Art Director : Seungnam, Kim****** 07/03/2016

      NDUCTP = NDUCT + 1
      NDUCTHP = NDUCTH + 1

      MDUCTP = MDUCT + 1

      DELTAK = TWOPI / NBLADE

C -- Mid Part (Backside to Faceside on blade)

      I11 = NDAFT + NMIDP 

C --new tmpi
      I22 = NDUCTHP + NDFWD
C --new tmpe

      DO N = 1, NMIDP

         THETAF = DANGLE(ZBB1(N), YBB1(N))
         THETAB = DANGLE(ZBB2(N), YBB2(N))

         IF(THETAB .LT. THETAF) THETAB = THETAB + 2.*PI

         DELTAX = (XBB1(N) - XBB2(N))/REAL(MDUCT)

         IF(DUCTGAP .NE. 0.0) THEN
            DTHETA = DELTAK / REAL(MDUCT)
         ELSE
            DDD = THETAB-THETAF
C           Hong added IF statement to avoid wrong DTHETA 03/06/08
            IF(DDD.GT.DELTAK) DDD=0.0
            DTHETA = (DELTAK-DDD) / REAL(MDUCT)
         ENDIF

         DTHETA2 = DELTAK / REAL(MDUCT)

C --new tmpi
         THETA1 = DANGLE(ZBK(N), YBK(N))
         CALL EVALDKs(NN1,1,XID2,XBK(N),RRR2,DCUBICX2)       
C --new tmpe

         DO M  = 1, MDUCTP

            XX1 = XBB2(N) + DELTAX*REAL(M-1)
         
            CALL EVALDKs(NN1,1,XID1,XX1,RRR1,DCUBICX1)

            THETA = THETAB + DTHETA*REAL(M-1)

            XD(I11-N+1,M) = XX1            
            YD(I11-N+1,M) = RRR1 * COS(THETA)
            ZD(I11-N+1,M) = RRR1 * SIN(THETA)

C --new tmpi
            THETA = THETA1 + DTHETA2*REAL(M-1)
            XD(I22+N-1,M) = XBK(N)          
            YD(I22+N-1,M) = RRR2 * COS(THETA)
            ZD(I22+N-1,M) = RRR2 * SIN(THETA)    
     
            IF(N .EQ. 1) THEN
               TBANGT(M) = THETA
               IF(THETA .GT. 2.*PI) TBANGT(M) = TBANGT(M) - 2.*PI
            ELSEIF(N .EQ. NMIDP) THEN
               TBANGL2(M) = THETA
               IF(THETA .GT. 2.*PI) TBANGL2(M) = TBANGL2(M) - 2.*PI
            ENDIF
C --new tmpe

         ENDDO
      ENDDO            

C -- Forward Part
      
C -- Face Side
 
      DTN = PI / REAL(NDFWD)

      I11 = NDAFT + NMIDP
 
      DO M = 1 , MDUCTP

         DO N = 2 , NN1
            IF(XD(I11,M) .GT. XID1(N-1) 
     %           .AND. XD(I11,M) .LE. XID1(N)) THEN
               NFF = N-1
               CALL EVALDKs(NN1,1,XID1,XD(I11,M),RRR1,DCUBICX1)
               P1 = SQRT( (XID1(NFF)-XD(I11,M))**2 
     %              +(ETAD1(NFF)-RRR1)**2)
               SSSS = SSDF1(NFF) + P1
            ENDIF
         ENDDO

         THETAF = DANGLE(ZD(I11,M),YD(I11,M))
         DELX1 = SSSS

         DO N = 1 , NDFWD
            XDI1 = 0.5*DELX1 * (1. - COS(DTN*N))

            XXXX = SSSS - XDI1

            CALL EVALDKs(NN1,1,SSDF1,XXXX,XD(I11+N,M),DCUBICS11)
            CALL EVALDKs(NN1,1,SSDF1,XXXX,RRR1,DCUBICS12)

            DTHETA1 = XDI1 * TAN(PITCH2)
            THETA2 = THETAF - DTHETA1
            
            YD(I11+N,M) = RRR1 * COS(THETA2)
            ZD(I11+N,M) = RRR1 * SIN(THETA2)

            IF(N .EQ. NDFWD) THEN
               TBANGL(M) = THETA2
               IF(THETA2 .GT. 2.*PI) TBANGL(M) = TBANGL(M) - 2.*PI
            ENDIF
         ENDDO
      ENDDO

C - Back Side

C -- New tmpi
  
      I22 = NDUCTHP + NDFWD

      DTN = PI / REAL(NDFWD)

      DO M = 1 , MDUCTP
         THETAL = TBANGL(M)      
         THETAT = TBANGT(M)         

         IF(THETAL .LT. 0.0) THETAL = THETAL + 2. * PI
         IF(THETAT .LT. 0.0) THETAT = THETAT + 2. * PI

         IF(THETAT .LT. THETAL) THETAT = THETAT + 2. * PI

         DO N = 2 , NN1
            IF(XD(I22,M) .GT. XID2(N-1) 
     %           .AND. XD(I22,M) .LE. XID2(N)) THEN
               NAA = N-1
               CALL EVALDKs(NN1,1,XID2,XD(I22,M),RRR2,DCUBICX2)
               P1 = SQRT( (XID2(NAA)-XD(I22,M))**2 
     %              +(ETAD2(NAA)-RRR2)**2)
               SSSS = SSDF2(NAA) + P1
            ENDIF
         ENDDO
         
C         DELX = XD(I22,M) - XD(NDUCTHP,M)
         DELX = SSSS
         DTN2 = (THETAT - THETAL)/REAL(NDFWD)
C         XDI = DELX / REAL(NDFWD)

         THPP  = (THETAT-THETAL) / DELX

         DO N = 1 , NDFWD-1
            XXXX = 0.5*DELX * (1. - COS(DTN*N))
            
            CALL EVALDKs(NN1,1,SSDF2,XXXX,XD(NDUCTHP+N,M),DCUBICS21)
            CALL EVALDKs(NN1,1,SSDF2,XXXX,RRR2,DCUBICS22)

            DTHE = XXXX * THPP
            THETA2 = THETAL + DTHE

            YD(NDUCTHP+N,M) = RRR2 * COS(THETA2)
            ZD(NDUCTHP+N,M) = RRR2 * SIN(THETA2)

          ENDDO
      ENDDO

C -- New tmpe

C ------------------------ S.N.KIM added. ----------------------------
C This part is added for ducted propellers with round blade tip. 
C The old version could not deal with the panel alignment between
C the short chord length around the blade tip and duct inner side.
C --------------------------------------------------------------------
C Front Face Part Renovation
      IF(XCHD(NX) .LE. 0.1) THEN

C        DTN = PI / REAL(NDFWD + NMID) !FULL COSINE
        DTN = HALFPI / REAL(NDFWD + NMID) !HALF COSINE
  
        I11 = NDAFT + 1
   
        DO M = 1 , MDUCTP
  
           DO N = 2 , NN1
              IF(XD(I11,M) .GT. XID1(N-1) 
     %             .AND. XD(I11,M) .LE. XID1(N)) THEN
                 NFF = N-1
                 CALL EVALDKs(NN1,1,XID1,XD(I11,M),RRR1,DCUBICX1)
                 P1 = SQRT( (XID1(NFF)-XD(I11,M))**2 
     %                +(ETAD1(NFF)-RRR1)**2)
                 SSSS = SSDF1(NFF) + P1
              ENDIF
           ENDDO

           THETAF = DANGLE(ZD(I11,M),YD(I11,M))
           DELX1 = SSSS

           DO N = 1 , NDFWD + NMID
C              XDI1 = 0.5 * DELX1 * (1. - COS(DTN*N)) !FULL COSINE
              XDI1 = DELX1 * SIN(DTN*N) !HALF COSINE
  
              XXXX = SSSS - XDI1

              CALL EVALDKs(NN1,1,SSDF1,XXXX,XD(I11+N,M),DCUBICS11)
              CALL EVALDKs(NN1,1,SSDF1,XXXX,RRR1,DCUBICS12)

              DTHETA1 = XDI1 * TAN(PITCH2)
              THETA2 = THETAF - DTHETA1
              
              YD(I11+N,M) = RRR1 * COS(THETA2)
              ZD(I11+N,M) = RRR1 * SIN(THETA2)
  
              IF(N .EQ. NDFWD + NMID) THEN
                 TBANGL(M) = THETA2
                 IF(THETA2 .GT. 2.*PI) TBANGL(M) = TBANGL(M) - 2.*PI
              ENDIF
           ENDDO
        ENDDO

C Front Back Part Renovation
        I22 = NDUCTHP + NDFWD + NMID
  
C        DTN = PI / REAL(NDFWD + NMID) ! FULL COSINE
        DTN = HALFPI / REAL(NDFWD + NMID) ! HALF COSINE
  
        DO M = 1 , MDUCTP
           THETAL = TBANGL(M)      
           THETAT = TBANGL2(M)         
  
           IF(THETAL .LT. 0.0) THETAL = THETAL + 2. * PI
           IF(THETAT .LT. 0.0) THETAT = THETAT + 2. * PI
  
           IF(THETAT .LT. THETAL) THETAT = THETAT + 2. * PI
  
           DO N = 2 , NN1
              IF(XD(I22,M) .GT. XID2(N-1) 
     %             .AND. XD(I22,M) .LE. XID2(N)) THEN
                 NAA = N-1
                 CALL EVALDKs(NN1,1,XID2,XD(I22,M),RRR2,DCUBICX2)
                 P1 = SQRT( (XID2(NAA)-XD(I22,M))**2 
     %                +(ETAD2(NAA)-RRR2)**2)
                 SSSS = SSDF2(NAA) + P1
              ENDIF
           ENDDO
         
C           DELX = XD(I22,M) - XD(NDUCTHP,M)
           DELX = SSSS
           DTN2 = (THETAT - THETAL) / REAL(NDFWD + NMID)
C           XDI = DELX / REAL(NDFWD)
  
           THPP  = (THETAT - THETAL) / DELX
  
           DO N = 1 , NDFWD + NMID - 1
C              XXXX = 0.5 * DELX * (1. - COS(DTN*N)) ! FULL COSINE
              XXXX = DELX * (1. - COS(DTN*N)) !HALF COSINE
  
              CALL EVALDKs(NN1,1,SSDF2,XXXX,XD(NDUCTHP+N,M),DCUBICS21)
C              XD(NDUCTHP+N,M) = XD(NDUCTHP-N,M)
              CALL EVALDKs(NN1,1,SSDF2,XXXX,RRR2,DCUBICS22)

              DTHE = XXXX * THPP
              THETA2 = THETAL + DTHE
  
              YD(NDUCTHP+N,M) = RRR2 * COS(THETA2)
              ZD(NDUCTHP+N,M) = RRR2 * SIN(THETA2)
  
            ENDDO
        ENDDO

      ENDIF
C ----- S.N.KIM added. -----

C -- Face Aft part
!********* DUCTGEO1.f Under The Modification for PR ********!
!********* Art Director : Seungnam Kim *********************!

      DELTAD = DELTAK / REAL(MDUCT)

      I11 = NDAFT+1  ! = 31
! KSN made Half Cosine
      DTN = HALFPI / REAL(NDAFT)

      DO M = 1 , MDUCTP

         DELX = XID1(NN1) - XD(I11,M)

         DO N = 1 , NDAFT

            XDI = DELX * (1.-COS(DTN*N)) ! Full Cosine Spacing.
! KSN sets Const Spacing
!             XDI = (DELX / NDAFT)*N   
! KSN sets Const Spacing
!*********** PR, Seungnam, Kim, 07/03/2016 **********            
            XD(I11-N,M) = XD(I11,M) + XDI ! Spacing Distribution.
            XD(I11-N,M) = XW(N+1,MRP)
            if (N.eq.NDAFT) XD(I11-N,M) = XD(I11,M) + XDI 
!*********** PR, Seungnam, Kim, 07/03/2016 **********
            CALL EVALDKs(NN1,1,XID1,XD(I11-N,M),RRR1,DCUBICX1)
            CALL EVALDKs(NDAF,1,XTT,XD(I11-N,M),RWW,XRTTCUB)  ! XTT(N:1,NDAF) = XW(N,MRP)
            CALL EVALDKs(NDAF,1,XTT,XD(I11-N,M),THETAT,XWTTCUB)
              
            YD(I11-N,M) = RRR1 * COS(THETAT+DELTAD*(M-1))
            ZD(I11-N,M) = RRR1 * SIN(THETAT+DELTAD*(M-1)) ! DELTAD = DELTAK / REAL(MDUCT)

            IF(N .EQ. NDAFT) THEN
               TBANGT(M) = THETAT+DELTAD*(M-1)
               IF(TBANGT(M) .GT. TWOPI) TBANGT(M) = TBANGT(M) - TWOPI
            ENDIF

         ENDDO

      ENDDO

! ------- STICK, PR, Seungnam, Kim, 24/03/2016, STICK ------- !
!      I11 = NDAFT + 1
c      IF (DUCTGAP .EQ. 0.0) THEN
c       DO Y = 1, NDAFTT-1
c        YW(Y+1,MRP) = YD(I11-Y,1)
c        ZW(Y+1,MRP) = ZD(I11-Y,1)
c       ENDDO
c      ENDIF
! ------- STICK, PR, Seungnam, Kim, 24/03/2016, STICK ------- ! 

C -- New tmpi

      I33 = NDUCTHP + NDFWD + NMID

      DO M = 1 , MDUCTP
         THETAL = TBANGL2(M)      
         THETAT = TBANGT(M)         
         
         IF(THETAL .LT. 0.0) THETAL = THETAL + 2. * PI
         IF(THETAT .LT. 0.0) THETAT = THETAT + 2. * PI
         
         IF(THETAT .LT. THETAL) THETAT = THETAT + 2. * PI

         DELX = XD(1,M) - XD(I33,M)
         
         DTN2 = (THETAT - THETAL)/REAL(NDAFT)

         THPP  = (THETAT-THETAL) / DELX

         DO N = 1 , NDAFT

! KSN sets Constant Spacing
            XDI = 0.5* DELX * (1.-COS(DTN*N))
!             XDI = (DELX / NDAFT)*N
! KSN sets Constant Spacing
! ********** PR, Seungnam, Kim, 07/03/2016 ********** 
            XD(I33+N,M) = XD(I33,M) + XDI 
            XD(I33+N,M) = XW(N+1,MRP)
            if (N.eq.NDAFT) XD(I33+N,M) = XD(I33,M) + XDI

            CALL EVALDKs(NN1,1,XID2,XD(I33+N,M),RRR2,DCUBICX2)
        
            XDI = XD(I33+N,M) - XD(I33,M) 
! ********** PR, Seungnam, Kim, 07/03/2016 **********
            DTHE = XDI * THPP
            THETA2 = THETAL + DTHE            

            YD(I33+N,M) = RRR2 * COS(THETA2)
            ZD(I33+N,M) = RRR2 * SIN(THETA2)

          ENDDO

          XD(NDUCTP,M) = XD(1,M)
          YD(NDUCTP,M) = YD(1,M)
          ZD(NDUCTP,M) = ZD(1,M)
      ENDDO

C      IF(XCHD(NX) .LE. 0.1) THEN
C
C         DO M = 1 , MDUCTP
C            SSDF1(1) = 0.0
C            SSDX1(1) = XD(NDUCTHP,M)
C            SSDR1(1) = SQRT( YD(NDUCTHP,M)**2+ZD(NDUCTHP,M)**2 )
C            SSDT1(1) = DANGLE( ZD(NDUCTHP,M),YD(NDUCTHP,M))
C
C            SSDF2(1) = 0.0
C            SSDX2(1) = XD(NDUCTHP,M)
C            SSDR2(1) = SQRT( YD(NDUCTHP,M)**2+ZD(NDUCTHP,M)**2 )
C            SSDT2(1) = DANGLE( ZD(NDUCTHP,M),YD(NDUCTHP,M))
C
C            DO I = 2, NDUCTHP
C               I1 = NDUCTHP - I + 1
C               I2 = NDUCTHP - I + 2
C               
C               SSDF1(I) = SSDF1(I-1) +
C     %              SQRT( (XD(I1,M) - XD(I2,M))**2 + 
C     %                    (YD(I1,M) - YD(I2,M))**2 + 
C     %                    (ZD(I1,M) - ZD(I2,M))**2 ) 
C               SSDX1(I) = XD(I1,M)
C               SSDR1(I) = SQRT( YD(I1,M)**2+ZD(I1,M)**2)
C               SSDT1(I) = DANGLE(ZD(I1,M),YD(I1,M))
C
C               I1 = NDUCTHP + I - 2
C               I2 = NDUCTHP + I - 1
C               
C               SSDF2(I) = SSDF2(I-1) +
C     %              SQRT( (XD(I1,M) - XD(I2,M))**2 + 
C     %                    (YD(I1,M) - YD(I2,M))**2 + 
C     %                    (ZD(I1,M) - ZD(I2,M))**2 ) 
C
C               SSDX2(I) = XD(I2,M)
C               SSDR2(I) = SQRT( YD(I2,M)**2+ZD(I2,M)**2)
C               SSDT2(I) = DANGLE(ZD(I2,M),YD(I2,M))
C            ENDDO
C
C            DO I = 1, NDUCTH
C               IF(SSDT1(I+1) .LT. SSDT1(I)) SSDT1(I+1)=SSDT1(I+1)+TWOPI
C               IF(SSDT2(I+1) .LT. SSDT2(I)) SSDT2(I+1)=SSDT2(I+1)+TWOPI
C            ENDDO
C
C            CALL UGLYDK(NDUCTHP,1,1,SSDF1,SSDX1,ZERO,ZERO,DCUBICS11)
C            CALL UGLYDK(NDUCTHP,1,1,SSDF1,SSDR1,ZERO,ZERO,DCUBICS12)            
C            CALL UGLYDK(NDUCTHP,1,1,SSDF1,SSDT1,ZERO,ZERO,DCUBICS13)            
C
C            CALL UGLYDK(NDUCTHP,1,1,SSDF2,SSDX2,ZERO,ZERO,DCUBICS21)
C            CALL UGLYDK(NDUCTHP,1,1,SSDF2,SSDR2,ZERO,ZERO,DCUBICS22)            
C            CALL UGLYDK(NDUCTHP,1,1,SSDF2,SSDT2,ZERO,ZERO,DCUBICS23)  
C 
C            DTN = PI / REAL( DNFWD + NMID )
C
C            DELX1 = SSDF1( NDFWD + NMID )
C            DELX2 = SSDF2( NDFWD + NMID )
C
C            DO N = 1, NDFWD + NMID
C               XXX1 = 0.5*DELX1 * (1. - COS(DTN*(N-1)))
C               XXX2 = 0.5*DELX2 * (1. - COS(DTN*(N-1)))
C               
C               CALL EVALDKs(NDUCTHP,1,SSDF1,XXX1,XXXX1,DCUBICS11)
C               CALL EVALDKs(NDUCTHP,1,SSDF1,XXX1,RRRR1,DCUBICS12)
C               CALL EVALDKs(NDUCTHP,1,SSDF1,XXX1,TTTT1,DCUBICS13)
C
C               CALL EVALDKs(NDUCTHP,1,SSDF2,XXX2,XXXX2,DCUBICS21)
C               CALL EVALDKs(NDUCTHP,1,SSDF2,XXX2,RRRR2,DCUBICS22)
C               CALL EVALDKs(NDUCTHP,1,SSDF2,XXX2,TTTT2,DCUBICS23)
C
C               XD(NDUCTHP-N+1,M) = XXXX1
C               YD(NDUCTHP-N+1,M) = RRRR1*COS(TTTT1)
C               ZD(NDUCTHP-N+1,M) = RRRR1*SIN(TTTT1)
C
C               XD(NDUCTHP+N-1,M) = XXXX2
C               YD(NDUCTHP+N-1,M) = RRRR2*COS(TTTT2)
C               ZD(NDUCTHP+N-1,M) = RRRR2*SIN(TTTT2)
C            ENDDO
C         ENDDO
C         
C      ENDIF

C ----- KSN, Rotate Duct Geometry -----
c      DO I = 1, nduct + 1
c        DO J = 1, mduct + 1
c          THTRS(I,J) = DANGLE(ZD(I,J),YD(I,J))
c          RDTRS(I,J) = SQRT(ZD(I,J)**2 + YD(I,J)**2)
c
c          YD(I,J) = RDTRS(I,J)*COS(THTRS(I,J)-(DELK/2))
c          ZD(I,J) = RDTRS(I,J)*SIN(THTRS(I,J)-(DELK/2))
c        ENDDO
c      ENDDO
C ----- KSN, Rotate Duct Geometry ----- 

      RETURN
      END       


C =============================================================
      SUBROUTINE DUCTGEO11
C
C     For IDGEO=0  Generate duct panel using propeller Pitch      
C =============================================================
      INCLUDE 'PUFCAV.INC'

      DIMENSION XBB1(200),YBB1(200),ZBB1(200)
      DIMENSION XBB2(200),YBB2(200),ZBB2(200)
      DIMENSION XID1(201),ETAD1(201),DCUBICX1(800)
      DIMENSION XID2(201),ETAD2(201),DCUBICX2(800)
      DIMENSION SSDF1(201),DCUBICS11(800),DCUBICS12(800)
      DIMENSION SSDF2(201),DCUBICS21(800),DCUBICS22(800)
      DIMENSION SSDX1(201),SSDR1(201),SSDT1(201),DCUBICS13(800)
      DIMENSION SSDX2(201),SSDR2(201),SSDT2(201),DCUBICS23(800)
      DIMENSION TBANGL(100),TBANGL2(100),TBANGT(100)
      DIMENSION XBK(200),YBK(200),ZBK(200)
      DIMENSION XTT(200),RTT(200),WTT(200),XRTTCUB(800),XWTTCUB(800)

      NN1 = NDDAT/2 + 1

      DO N = 1 , NN1
         N1 = NN1 - N + 1
         N2 = NN1 + N - 1
         XID1(N) = XID(N1)
         ETAD1(N) = ETAD(N1)
         
         XID2(N) = XID(N2)
         ETAD2(N) = ETAD(N2)
      ENDDO

      SSDF1(1) = 0.0
      SSDF2(1) = 0.0

      DO N = 2 , NN1
         P1 = SQRT( (XID1(N)-XID1(N-1))**2 
     %             +(ETAD1(N)-ETAD1(N-1))**2)
         SSDF1(N) = SSDF1(N-1) + P1

         P2 = SQRT( (XID2(N)-XID2(N-1))**2 
     %             +(ETAD2(N)-ETAD2(N-1))**2)
         SSDF2(N) = SSDF2(N-1) + P2

      ENDDO
     
      CALL UGLYDK(NN1,1,1,XID1,ETAD1,ZERO,ZERO,DCUBICX1)

      CALL UGLYDK(NN1,1,1,SSDF1,XID1,ZERO,ZERO,DCUBICS11)
      CALL UGLYDK(NN1,1,1,SSDF1,ETAD1,ZERO,ZERO,DCUBICS12)

      CALL UGLYDK(NN1,1,1,SSDF2,XID2,ZERO,ZERO,DCUBICS21)
      CALL UGLYDK(NN1,1,1,SSDF2,ETAD2,ZERO,ZERO,DCUBICS22)

      CALL UGLYDK(NN1,1,1,XID2,ETAD2,ZERO,ZERO,DCUBICX2)

      PPD = XPI(NX)
     
      DO N = 1 , NHP
         I1 = NHP - N + 1
         I2 = NHP + N - 1

         IF(DUCTGAP .EQ. 0.0) THEN
            XBB1(N) = XB(I1,MRP)
            YBB1(N) = YB(I1,MRP)
            ZBB1(N) = ZB(I1,MRP)
            
            XBB2(N) = XB(I2,MRP)
            YBB2(N) = YB(I2,MRP)
            ZBB2(N) = ZB(I2,MRP)
         ELSE
            XBB1(N) = 0.5*(XB(I1,MRP)+XB(I2,MRP))
            YBB1(N) = 0.5*(YB(I1,MRP)+YB(I2,MRP))
            ZBB1(N) = 0.5*(ZB(I1,MRP)+ZB(I2,MRP))

            XBB2(N) = XBB1(N)
            YBB2(N) = YBB1(N)
            ZBB2(N) = ZBB1(N)
         ENDIF
         XBK(N) = 0.5*(XB(I1,MRP)+XB(I2,MRP))
         YBK(N) = 0.5*(YB(I1,MRP)+YB(I2,MRP))
         ZBK(N) = 0.5*(ZB(I1,MRP)+ZB(I2,MRP))
      ENDDO

C -- Check No. of Panel on After part of duct ..

      DO N = 1 , NSW(MRP)
         IF(XID1(NN1) .GT. XW(N,MRP) .AND. 
     %          XID1(NN1) .LE. XW(N+1,MRP)) THEN
            NDAFTT =N 
         ENDIF
      ENDDO

      NDAF = NDAFTT + 1

      DO N = 1 , NDAF
         XTT(N) = XW(N,MRP)
         RTT(N) = SQRT(YW(N,MRP)**2+ZW(N,MRP)**2)
         WTT(N) = DANGLE(ZW(N,MRP),YW(N,MRP))
         IF(WTT(N) .LT. 0.0) WTT(N) = WTT(N) + TWOPI
      ENDDO

      DO N = 2 , NDAF
         IF(WTT(N-1) .GT. WTT(N)) THEN
            WTT(N) = WTT(N) + TWOPI
         ENDIF
      ENDDO

      CALL UGLYDK(NDAF,1,1,XTT,RTT,ZERO,ZERO,XRTTCUB)    
      CALL UGLYDK(NDAF,1,1,XTT,WTT,ZERO,ZERO,XWTTCUB)    

      NMID = NH
      NMIDP = NMID+1
      PITCH2 = 0.5 * PI - ATAN(PPD/PI)

      NDUCTH = NDUCT/2
      NDAFT = (NDUCTH-NMID)/2
      NDFWD = NDUCTH - NDAFT - NMID
 
C      write(*,*) ndfwd, nmid, ndaft, nduct

      NDUCTP = NDUCT + 1
      NDUCTHP = NDUCTH + 1

      MDUCTP = MDUCT + 1

      DELTAK = TWOPI / NBLADE

C -- Mid Part (Backside to Faceside on blade)

      I11 = NDAFT + NMIDP 

C --new tmpi
      I22 = NDUCTHP + NDFWD
C --new tmpe

      DO N = 1, NMIDP

         THETAF = DANGLE(ZBB1(N), YBB1(N))
         THETAB = DANGLE(ZBB2(N), YBB2(N))

         IF(THETAB .LT. THETAF) THETAB = THETAB + 2.*PI

         DELTAX = (XBB1(N) - XBB2(N))/REAL(MDUCT)

         IF(DUCTGAP .NE. 0.0) THEN
            DTHETA = DELTAK / REAL(MDUCT)
         ELSE
            DDD = THETAB-THETAF
C           Hong added IF statement to avoid wrong DTHETA 03/06/08
            IF(DDD.GT.DELTAK) DDD=0.0
            DTHETA = (DELTAK-DDD) / REAL(MDUCT)
         ENDIF

         DTHETA2 = DELTAK / REAL(MDUCT)

C --new tmpi
         THETA1 = DANGLE(ZBK(N), YBK(N))
         CALL EVALDKs(NN1,1,XID2,XBK(N),RRR2,DCUBICX2)       
C --new tmpe

         DO M  = 1, MDUCTP

            XX1 = XBB2(N) + DELTAX*REAL(M-1)
         
            CALL EVALDKs(NN1,1,XID1,XX1,RRR1,DCUBICX1)

            THETA = THETAB + DTHETA*REAL(M-1)

            XD(I11-N+1,M) = XX1            
            YD(I11-N+1,M) = RRR1 * COS(THETA)
            ZD(I11-N+1,M) = RRR1 * SIN(THETA)

C --new tmpi
            THETA = THETA1 + DTHETA2*REAL(M-1)
            XD(I22+N-1,M) = XBK(N)          
            YD(I22+N-1,M) = RRR2 * COS(THETA)
            ZD(I22+N-1,M) = RRR2 * SIN(THETA)    
     
            IF(N .EQ. 1) THEN
               TBANGT(M) = THETA
               IF(THETA .GT. 2.*PI) TBANGT(M) = TBANGT(M) - 2.*PI
            ELSEIF(N .EQ. NMIDP) THEN
               TBANGL2(M) = THETA
               IF(THETA .GT. 2.*PI) TBANGL2(M) = TBANGL2(M) - 2.*PI
            ENDIF
C --new tmpe

         ENDDO
      ENDDO            

C -- Forward Part
     
C -- Face Side
 
      DTN = PI / REAL(NDFWD)

      I11 = NDAFT + NMIDP
 
      DO M = 1 , MDUCTP

         DO N = 2 , NN1
            IF(XD(I11,M) .GT. XID1(N-1) 
     %           .AND. XD(I11,M) .LE. XID1(N)) THEN
               NFF = N-1
               CALL EVALDKs(NN1,1,XID1,XD(I11,M),RRR1,DCUBICX1)
               P1 = SQRT( (XID1(NFF)-XD(I11,M))**2 
     %              +(ETAD1(NFF)-RRR1)**2)
               SSSS = SSDF1(NFF) + P1
            ENDIF
         ENDDO

         THETAF = DANGLE(ZD(I11,M),YD(I11,M))
         DELX1 = SSSS

         DO N = 1 , NDFWD
            XDI1 = 0.5*DELX1 * (1. - COS(DTN*N))

            XXXX = SSSS - XDI1

            CALL EVALDKs(NN1,1,SSDF1,XXXX,XD(I11+N,M),DCUBICS11)
            CALL EVALDKs(NN1,1,SSDF1,XXXX,RRR1,DCUBICS12)

            DTHETA1 = XDI1 * TAN(PITCH2)
            THETA2 = THETAF - DTHETA1
            
            YD(I11+N,M) = RRR1 * COS(THETA2)
            ZD(I11+N,M) = RRR1 * SIN(THETA2)

            IF(N .EQ. NDFWD) THEN
               TBANGL(M) = THETA2
               IF(THETA2 .GT. 2.*PI) TBANGL(M) = TBANGL(M) - 2.*PI
            ENDIF
         ENDDO
      ENDDO

C - Back Side

C -- New tmpi
  
      I22 = NDUCTHP + NDFWD

      DTN = PI / REAL(NDFWD)

      DO M = 1 , MDUCTP
         THETAL = TBANGL(M)      
         THETAT = TBANGT(M)         

         IF(THETAL .LT. 0.0) THETAL = THETAL + 2. * PI
         IF(THETAT .LT. 0.0) THETAT = THETAT + 2. * PI

         IF(THETAT .LT. THETAL) THETAT = THETAT + 2. * PI

         DO N = 2 , NN1
            IF(XD(I22,M) .GT. XID2(N-1) 
     %           .AND. XD(I22,M) .LE. XID2(N)) THEN
               NAA = N-1
               CALL EVALDKs(NN1,1,XID2,XD(I22,M),RRR2,DCUBICX2)
               P1 = SQRT( (XID2(NAA)-XD(I22,M))**2 
     %              +(ETAD2(NAA)-RRR2)**2)
               SSSS = SSDF2(NAA) + P1
            ENDIF
         ENDDO
         
C         DELX = XD(I22,M) - XD(NDUCTHP,M)
         DELX = SSSS
         DTN2 = (THETAT - THETAL)/REAL(NDFWD)
C         XDI = DELX / REAL(NDFWD)

         THPP  = (THETAT-THETAL) / DELX

         DO N = 1 , NDFWD-1
            XXXX = 0.5*DELX * (1. - COS(DTN*N))
            
C            XD(NDUCTHP+N,M) = XD(NDUCTHP,M) + XDI*REAL(N) 
C            XD(NDUCTHP+N,M) = XD(NDUCTHP,M) + XDI

            CALL EVALDKs(NN1,1,SSDF2,XXXX,XD(NDUCTHP+N,M),DCUBICS21)
            CALL EVALDKs(NN1,1,SSDF2,XXXX,RRR2,DCUBICS22)

C            CALL EVALDK(NN1,1,XID2,XD(NDUCTHP+N,M),RRR2,DCUBICX2)

C            DTHE = XDI * THPP
            DTHE = XXXX * THPP
            THETA2 = THETAL + DTHE

            YD(NDUCTHP+N,M) = RRR2 * COS(THETA2)
            ZD(NDUCTHP+N,M) = RRR2 * SIN(THETA2)

          ENDDO
      ENDDO

C -- New tmpe

C -- Face Aft part

      DELTAD = DELTAK / REAL(MDUCT)

      I11 = NDAFT+1

      DTN = PI / REAL(NDAFT)

      DO M = 1 , MDUCTP

         DELX = XID1(NN1) - XD(I11,M)

         DO N = 1 , NDAFT

            XDI = 0.5*DELX * (1. - COS(DTN*N))

            XD(I11-N,M) = XD(I11,M) + XDI

C            write(*,*) m , n, xd(i11-n,m)
            CALL EVALDKs(NN1,1,XID1,XD(I11-N,M),RRR1,DCUBICX1)
            CALL EVALDKs(NDAF,1,XTT,XD(I11-N,M),RWW,XRTTCUB)
            CALL EVALDKs(NDAF,1,XTT,XD(I11-N,M),THETAT,XWTTCUB)
              
            YD(I11-N,M) = RRR1 * COS(THETAT+DELTAD*(M-1))
            ZD(I11-N,M) = RRR1 * SIN(THETAT+DELTAD*(M-1))

            IF(N .EQ. NDAFT) THEN
               TBANGT(M) = THETAT+DELTAD*(M-1)
               IF(TBANGT(M) .GT. TWOPI) TBANGT(M) = TBANGT(M) - TWOPI
            ENDIF

         ENDDO

      ENDDO

C -- New tmpi

      I33 = NDUCTHP + NDFWD + NMID

      DO M = 1 , MDUCTP
         THETAL = TBANGL2(M)      
         THETAT = TBANGT(M)         
         
         IF(THETAL .LT. 0.0) THETAL = THETAL + 2. * PI
         IF(THETAT .LT. 0.0) THETAT = THETAT + 2. * PI
         
         IF(THETAT .LT. THETAL) THETAT = THETAT + 2. * PI

         DELX = XD(1,M) - XD(I33,M)
         
         DTN2 = (THETAT - THETAL)/REAL(NDAFT)

         THPP  = (THETAT-THETAL) / DELX

         DO N = 1 , NDAFT

            XDI = 0.5*DELX * (1. - COS(DTN*N))
 
            XD(I33+N,M) = XD(I33,M) + XDI 

            CALL EVALDKs(NN1,1,XID2,XD(I33+N,M),RRR2,DCUBICX2)

            DTHE = XDI * THPP
            THETA2 = THETAL + DTHE            

            YD(I33+N,M) = RRR2 * COS(THETA2)
            ZD(I33+N,M) = RRR2 * SIN(THETA2)

          ENDDO

          XD(NDUCTP,M) = XD(1,M)
          YD(NDUCTP,M) = YD(1,M)
          ZD(NDUCTP,M) = ZD(1,M)
      ENDDO


      IF(XCHD(NX) .LE. 0.1) THEN

         DO M = 1 , MDUCTP
            SSDF1(1) = 0.0
            SSDX1(1) = XD(NDUCTHP,M)
            SSDR1(1) = SQRT( YD(NDUCTHP,M)**2+ZD(NDUCTHP,M)**2 )
            SSDT1(1) = DANGLE( ZD(NDUCTHP,M),YD(NDUCTHP,M))

            SSDF2(1) = 0.0
            SSDX2(1) = XD(NDUCTHP,M)
            SSDR2(1) = SQRT( YD(NDUCTHP,M)**2+ZD(NDUCTHP,M)**2 )
            SSDT2(1) = DANGLE( ZD(NDUCTHP,M),YD(NDUCTHP,M))

            DO I = 2, NDUCTHP
               I1 = NDUCTHP - I + 1
               I2 = NDUCTHP - I + 2
               
               SSDF1(I) = SSDF1(I-1) +
     %              SQRT( (XD(I1,M) - XD(I2,M))**2 + 
     %                    (YD(I1,M) - YD(I2,M))**2 + 
     %                    (ZD(I1,M) - ZD(I2,M))**2 ) 
               SSDX1(I) = XD(I1,M)
               SSDR1(I) = SQRT( YD(I1,M)**2+ZD(I1,M)**2)
               SSDT1(I) = DANGLE(ZD(I1,M),YD(I1,M))

               I1 = NDUCTHP + I - 2
               I2 = NDUCTHP + I - 1
               
               SSDF2(I) = SSDF2(I-1) +
     %              SQRT( (XD(I1,M) - XD(I2,M))**2 + 
     %                    (YD(I1,M) - YD(I2,M))**2 + 
     %                    (ZD(I1,M) - ZD(I2,M))**2 ) 

               SSDX2(I) = XD(I2,M)
               SSDR2(I) = SQRT( YD(I2,M)**2+ZD(I2,M)**2)
               SSDT2(I) = DANGLE(ZD(I2,M),YD(I2,M))
            ENDDO

            DO I = 1, NDUCTH
               IF(SSDT1(I+1) .LT. SSDT1(I)) SSDT1(I+1)=SSDT1(I+1)+TWOPI
               IF(SSDT2(I+1) .LT. SSDT2(I)) SSDT2(I+1)=SSDT2(I+1)+TWOPI
            ENDDO

            CALL UGLYDK(NDUCTHP,1,1,SSDF1,SSDX1,ZERO,ZERO,DCUBICS11)
            CALL UGLYDK(NDUCTHP,1,1,SSDF1,SSDR1,ZERO,ZERO,DCUBICS12)            
            CALL UGLYDK(NDUCTHP,1,1,SSDF1,SSDT1,ZERO,ZERO,DCUBICS13)            

            CALL UGLYDK(NDUCTHP,1,1,SSDF2,SSDX2,ZERO,ZERO,DCUBICS21)
            CALL UGLYDK(NDUCTHP,1,1,SSDF2,SSDR2,ZERO,ZERO,DCUBICS22)            
            CALL UGLYDK(NDUCTHP,1,1,SSDF2,SSDT2,ZERO,ZERO,DCUBICS23)  

            DTN = PI / REAL(NDUCTH)

            DELX1 = SSDF1(NDUCTHP)
            DELX2 = SSDF2(NDUCTHP)

            DO N = 1, NDUCTHP
               XXX1 = 0.5*DELX1 * (1. - COS(DTN*(N-1)))
               XXX2 = 0.5*DELX2 * (1. - COS(DTN*(N-1)))
               
               CALL EVALDKs(NDUCTHP,1,SSDF1,XXX1,XXXX1,DCUBICS11)
               CALL EVALDKs(NDUCTHP,1,SSDF1,XXX1,RRRR1,DCUBICS12)
               CALL EVALDKs(NDUCTHP,1,SSDF1,XXX1,TTTT1,DCUBICS13)

               CALL EVALDKs(NDUCTHP,1,SSDF2,XXX2,XXXX2,DCUBICS21)
               CALL EVALDKs(NDUCTHP,1,SSDF2,XXX2,RRRR2,DCUBICS22)
               CALL EVALDKs(NDUCTHP,1,SSDF2,XXX2,TTTT2,DCUBICS23)

               XD(NDUCTHP-N+1,M) = XXXX1
               YD(NDUCTHP-N+1,M) = RRRR1*COS(TTTT1)
               ZD(NDUCTHP-N+1,M) = RRRR1*SIN(TTTT1)

               XD(NDUCTHP+N-1,M) = XXXX2
               YD(NDUCTHP+N-1,M) = RRRR2*COS(TTTT2)
               ZD(NDUCTHP+N-1,M) = RRRR2*SIN(TTTT2)
            ENDDO
         ENDDO
         
      ENDIF


      RETURN
      END

