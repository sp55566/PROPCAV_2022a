C ===============================================
      SUBROUTINE DUCTWAKGEO
C ===============================================
      INCLUDE 'PUFCAV.INC'
!------- Art Director Added This Line, 03/22/2016 -------!
      DIMENSION XID12(200), ETAD12(200)
      DIMENSION XID22(201), ETAD22(201)
      DIMENSION XXT(202), WWT(202), XWWTCUB(802) 
!------- Art Director Added This Line, 03/22/2016 
      IF(IDGEO .EQ. 0) THEN

       RRR = SQRT(YD(1,1)**2 + ZD(1,1)**2)

       NDWK = 100
C      NDWKPP = NDWK + 1
       DWL = DWAKEL          ! DWAKEL = Duct Wake Length ( Written From Input File)
       DTN = PI / REAL(NDWK)

       DELTAK = TWOPI / NBLADE
! ------- Check Number of Panel on After part of duct ------- !
       NN12 = NDDAT/2 + 1
      
       DO N = 1, NN12
         N12 = NN12 - N + 1
         N22 = NN12 + M - 1
         XID12(N) = XID(N12)
         ETAD12(N) = ETAD(N12)
         XID22(N) = XID(N22)
         ETAD22(N) = ETAD(N22)
       ENDDO      

       DO N = 1, NSW(MRP)
        IF(XID12(NN12) .GT. XW(N,MRP) .AND.
     *  XID12(NN12) .LE. XW(N+1,MRP)) THEN
        NDAFTT2 = N
        ENDIF
       ENDDO
       NDAF2 = NDAFTT2 + 1
C ------- End of Checking Number of Panels on duct ------- ! 
C      IF(IDGEO .EQ. 0) THEN
C         PPD = XPI(NX)             ! Wake After Duct is aligned with Blade
C         PITCH2 = 0.5 * PI - ATAN(PPD/PI)
C      ELSEIF(IDGEO .EQ. 1) THEN
C         PITCH2 = DUCTPT
C      ENDIF
      
!------- Repanelling Project, Art Director : Seungnam, Kim, 03/22/2016 -------!
!------- I changed pitch of ductwake according to that of after part of duct -!
        DO N = 1, NSW(MRP)
          XXT(N) = XW(N,MRP)
          WWT(N) = DANGLE(ZW(N,MRP), YW(N,MRP))
          IF(WWT(N) .LT. 0.0) WWT(N) = WWT(N) + TWOPI
        ENDDO

        DO N = 2, NSW(MRP)
          IF(WWT(N-1) .GT. WWT(N)) THEN
            WWT(N) = WWT(N) + TWOPI
          ENDIF
        ENDDO
 
        CALL UGLYDK(NSW(MRP),1,1,XXT,WWT,ZERO,ZERO,XWWTCUB)
C      NN12 = NDDAT/2 +1
C  
C      DO N = 1, NN12
C         N12 = NN12 - N +1
C         N22 = NN12 + N -1
C         XID12(N) = XID(N12)
C         ETAD12(N) = ETAD(N12)
C         XID22(N) = XID(N22)
C         ETAD22(N) = ETAD(N22)
C      ENDDO     
C      
C      DO N = 1, NSW(MRP)
C         IF(XID12(NN12) .GT. XW(N,MRP) .AND.
C     * XID12(NN12) .LE. XW(N+1,MRP)) THEN
C          NDAFTT2 = N
C         ENDIF
C      ENDDO
C      NDAF2 = NDAFTT2 + 1
      
C      THETAD1 = DANGLE(ZW(NDAF2+1,MRP),YW(NDAF2+1,MRP))
C      THETAD2 = DANGLE(ZW(NDAF2,MRP),YW(NDAF2,MRP))
C      DELTATHETAD = ABS(THETAD1 - THETAD2)
C      DELTAXD = XW(NDAF2+1,MRP) - XW(NDAF2,MRP)
C      PITCH2 = ATAN(RRR * DELTATHETAD / DELTAXD)
!------- Repanelling Project, Art Director : Seungnam, Kim, 03/22/2016 -------!
       DO M = 1 , MDUCTP
          XDW(1,M) = XD(1,M)
          YDW(1,M) = YD(1,M)
          ZDW(1,M) = ZD(1,M)

          THETA1 = DANGLE(ZDW(1,M),YDW(1,M))
          DELXX = XW(NSW(MRP)-1,MRP) - XW(NSW(MRP)-2,MRP)
          DELXXEND = XDW(1,M) + DWL - XW(NSW(MRP)-1,MRP)
          NPMORE = PI / (2*ASIN(DELXX / DELXXEND))
          NDWKP = NSW(MRP) - NDAFTT2 + NPMORE
          NDWK = NDWKP - 1


         
          DO N = 2 , NDWKP   ! NDWKP = NSW(MRP) - NDAFTT2 + NPMORE
           
             IF(N - 1 .LE. (NSW(MRP) - NDAFTT2 - 1)) THEN
               XDW(N,M) = XW(NDAF2 + N -2, MRP) 
             ELSE
               J = N - (NSW(MRP) - NDAFTT2) 
               XDI = DELXXEND * SIN((HALFPI/NPMORE)*J)
               XDW(N,M) = XW(NSW(MRP)-1,MRP) + XDI             
             ENDIF

             IF(XDW(N,M) .LE. XW(NSW(MRP)-1,MRP)) THEN
               CALL EVALDKs(NSW(MRP),1,XXT,XDW(N,M),WKPITCH,XWWTCUB)
               DTHETA = WKPITCH
C              DTHETA = XDI * TAN(PITCH2) / RRR  ! l = r * theta 
               THETA11 = THETA1 - DANGLE (ZDW(1,1),YDW(1,1))
               THETA2  = THETA11 + DTHETA
               YDW(N,M) = RRR * COS(THETA2)
               ZDW(N,M) = RRR * SIN(THETA2)    
              
             ELSE      
              
               DO P = 1, N - 1 ! I spent whole three days to find a bug right here.
                  IF(XW(NSW(MRP)-1,MRP) .GT. XDW(P,M) .AND.
     *               XW(NSW(MRP)-1,MRP) .LE. XDW(P+1,M)) THEN
                   NWAFTT = P
                  ENDIF
               ENDDO

               THETAD1 = DANGLE(ZDW(NWAFTT,1) , YDW(NWAFTT,1))
               THETAD2 = DANGLE(ZDW(NWAFTT-1,1) , YDW(NWAFTT-1,1))
               DELTATHETAD = THETAD1 - THETAD2
               IF(DELTATHETAD.LT.0.0) DELTATHETAD=(2.00*PI)+DELTATHETAD
               DELTAXD = XDW(NWAFTT,M) - XDW(NWAFTT-1,M)
               PITCH2 = DANGLE(RRR * DELTATHETAD , DELTAXD)       
               XDII = XDI
               DTHETA = XDII * TAN(PITCH2) / RRR
               THETA3 = THETA2 + DTHETA
               YDW(N,M) = RRR * COS(THETA3)
               ZDW(N,M) = RRR * SIN(THETA3) 
 
            ENDIF
          ENDDO
       ENDDO

! ------- Stick, Panel Artist : Seungnam Kim, 24/03/2016 ------- !
c      IF (DUCTGAP .EQ. 0.0) THEN
c       DO Y = 1, NSW(MRP) - NDAFTT2
c        YW(Y+NDAFTT2,MRP) = YDW(Y+1,1)
c        ZW(Y+NDAFTT2,MRP) = ZDW(Y+1,1)
c       ENDDO
c      ENDIF
! ------- Stick, Panel Artist : Seungnam Kim, 24/03/2016 ------- !

!----Up to here, I have finished modification on the duct wake. 26/03/2016----!
!------- Repanelling Project, Art Director : Seungnam, Kim, 23/03/2016 -------!
     
      ELSE 
 
       RRR = SQRT(YD(1,1)**2 + ZD(1,1)**2)

       NDWK = 150
       NDWKP = NDWK + 1
       DWL = DWAKEL
       DTN = PI / REAL(NDWK)

       DELTAK = TWOPI / NBLADE

       IF(IDGEO .EQ. 0) THEN
          PPD = XPI(NX)
          PITCH2 = 0.5 * PI - ATAN(PPD/PI)
       ELSEIF(IDGEO .EQ. 1) THEN
          PITCH2 = DUCTPT
       ENDIF

       DO M = 1 , MDUCTP
          XDW(1,M) = XD(1,M)
          YDW(1,M) = YD(1,M)
          ZDW(1,M) = ZD(1,M)

          THETA1 = DANGLE(ZDW(1,M),YDW(1,M))

          DO N = 2 , NDWKP

             XDI = 0.5 * DWL * (1. - COS(DTN*(N-1)))
             XDW(N,M) = XDW(1,M) + XDI

             DTHETA = XDI * TAN(PITCH2) / RRR
             THETA2 = THETA1 + DTHETA

             YDW(N,M) = RRR * COS(THETA2)
             ZDW(N,M) = RRR * SIN(THETA2)
          ENDDO
       ENDDO

      ENDIF

C --- Duct Wake Ends with Disk

      DO M = 1 , MDUCTP
         XDWDK(M,2) = XDW(NDWKP,M)
         YDWDK(M,2) = YDW(NDWKP,M)
         ZDWDK(M,2) = ZDW(NDWKP,M)

         XDWDK(M,1) = XDWDK(M,2) 
         YDWDK(M,1) = 0.0
         ZDWDK(M,1) = 0.0
      ENDDO

      RETURN
      END

C ===============================================
      SUBROUTINE DUCTWAKGEO2
C ===============================================
      INCLUDE 'PUFCAV.INC'


      RRR = SQRT(YD(1,1)**2 + ZD(1,1)**2)

      NDWK = 100
      NDWKP = NDWK + 1
      DWL = DWAKEL
      DTN = PI / REAL(NDWK) !FULL COSINE
c      DTN = HALFPI / REAL(NDWK) !HALF COSINE

      DELTAK = TWOPI / NBLADE

      IF(IDGEO .EQ. 0) THEN
         PPD = XPI(NX)
         PITCH2 = 0.5 * PI - ATAN(PPD/PI)
      ELSEIF(IDGEO .EQ. 1) THEN
         PITCH2 = DUCTPT
      ENDIF

      DO M = 1 , MDUCTP
         XDW(1,M) = XD(1,M)
         YDW(1,M) = YD(1,M)
         ZDW(1,M) = ZD(1,M)

         THETA1 = DANGLE(ZDW(1,M),YDW(1,M))

         DO N = 2 , NDWKP

            XDI = 0.5 * DWL * (1. - COS(DTN*(N-1))) !FULL COSINE
c            XDI = DWL * SIN(DTN*(N-1)) !HALF COSINE
            XDW(N,M) = XDW(1,M) + XDI

            DTHETA = XDI * TAN(PITCH2) / RRR
            THETA2 = THETA1 + DTHETA

            YDW(N,M) = RRR * COS(THETA2)
            ZDW(N,M) = RRR * SIN(THETA2)
         ENDDO
      ENDDO

C --- Duct Wake Ends with Disk

      DO M = 1 , MDUCTP
         XDWDK(M,2) = XDW(NDWKP,M)
         YDWDK(M,2) = YDW(NDWKP,M)
         ZDWDK(M,2) = ZDW(NDWKP,M)

         XDWDK(M,1) = XDWDK(M,2) 
         YDWDK(M,1) = 0.0
         ZDWDK(M,1) = 0.0
      ENDDO

      RETURN
      END

C ===============================================
      SUBROUTINE DUCTWAKGEO_FIXEDPANEL
C ===============================================
      INCLUDE 'PUFCAV.INC'
!------- Art Director Added This Line, 03/22/2016 -------!
      DIMENSION XID12(200), ETAD12(200)
      DIMENSION XID22(201), ETAD22(201)
      DIMENSION XXT(202), WWT(202), XWWTCUB(802) 
      DIMENSION THTRSW(300,MDUCT+1)

!------- Art Director Added This Line, 03/22/2016 
      IF(IDGEO .EQ. 0) THEN

       RRR = SQRT(YD(1,1)**2 + ZD(1,1)**2)
! ------- Check Number of Panel on After part of duct ------- !
       NN12 = NDDAT/2 + 1
      
       DO N = 1, NN12
         N12 = NN12 - N + 1
         N22 = NN12 + M - 1
         XID12(N) = XID(N12)
         ETAD12(N) = ETAD(N12)
         XID22(N) = XID(N22)
         ETAD22(N) = ETAD(N22)
       ENDDO      

       DO N = 1, NSW(MRP)
        IF(XID12(NN12) .GT. XW(N,MRP) .AND.
     *  XID12(NN12) .LE. XW(N+1,MRP)) THEN
        NDAFTT2 = N
        ENDIF
       ENDDO
       NDAF2 = NDAFTT2 + 1
!------- Repanelling Project, Art Director : Seungnam, Kim, 03/22/2016 -------!
!------- I changed pitch of ductwake according to that of after part of duct -!
       DO N = 1, NSW(MRP)
         XXT(N) = XW(N,MRP)
         WWT(N) = DANGLE(ZW(N,MRP), YW(N,MRP))
         IF(WWT(N) .LT. 0.0) WWT(N) = WWT(N) + TWOPI
       ENDDO

       DO N = 2, NSW(MRP)
         IF(WWT(N-1) .GT. WWT(N)) THEN
           WWT(N) = WWT(N) + TWOPI
         ENDIF
       ENDDO

C ----- KSN added IF 
C       IF(ICAVT .LT. 1) THEN 
       CALL UGLYDK(NSW(MRP),1,1,XXT,WWT,ZERO,ZERO,XWWTCUB)
C       ENDIF
C ----- KSN added

!------- Repanelling Project, Art Director : Seungnam, Kim, 03/22/2016 -------!
       
       DO M = 1 , MDUCTP
          XDW(1,M) = XD(1,M)
          YDW(1,M) = YD(1,M)
          ZDW(1,M) = ZD(1,M)

          THETA1 = DANGLE(ZDW(1,M),YDW(1,M))
          NDWK = NWAKEP
          NDWKP = NDWK + 1
         
          DO N = 2 , NDWKP   ! NDWKP = NWAKEP + 1
           
             IF(N - 1 .LE. NSW(MRP) - NDAFTT2) THEN  ! 0 -> -1
               XDW(N,M) = XW(NDAF2 + N - 2, MRP) ! 2 -> 1 
             ELSE
C               J = N - (NSW(MRP) - NDAFTT2) 
C               XDI = DELXXEND * SIN((HALFPI/NPMORE)*J)
               XDI = XW(NSW(MRP),MRP) - XW(NSW(MRP)-1,MRP)
               XDW(N,M) = XDW(N-1,M) + XDI
             ENDIF

             IF(XDW(N,M) .LE. XW(NSW(MRP),MRP)) THEN
               CALL EVALDKs(NSW(MRP),1,XXT,XDW(N,M),WKPITCH,XWWTCUB)
               DTHETA = WKPITCH
               IF(DTHETA.LT.0.0) DTHETA = 2.00 * PI + DTHETA
C              DTHETA = XDI * TAN(PITCH2) / RRR  ! l = r * theta 
               THETA11 = THETA1 - DANGLE (ZDW(1,1),YDW(1,1))
               THETA2  = THETA11 + DTHETA
               YDW(N,M) = RRR * COS(THETA2)
               ZDW(N,M) = RRR * SIN(THETA2)    
               
             ELSE      
              
               DO P = 1, N - 1 ! I spent whole three days to find a bug right here. - S.N.KIM
                  IF(XW(NSW(MRP),MRP) .GT. XDW(P,M) .AND.
     *               XW(NSW(MRP),MRP) .LE. XDW(P+1,M)) THEN
                   NWAFTT = P
                  ENDIF
               ENDDO

               THETAD1 = DANGLE(ZDW(NWAFTT+1,1) , YDW(NWAFTT+1,1))
               THETAD2 = DANGLE(ZDW(NWAFTT,1) , YDW(NWAFTT,1))
               DELTATHETAD = THETAD1 - THETAD2
               IF(DELTATHETAD.LT.0.0) DELTATHETAD=(2.00*PI)+DELTATHETAD
               DELTAXD = XDW(NWAFTT+1,M) - XDW(NWAFTT,M)
               PITCH2 = DANGLE(RRR * DELTATHETAD , DELTAXD)       
C               XDII = XDI
               XDII = XDW(N,M) - XW(NSW(MRP),MRP)
               DTHETA = XDII * TAN(PITCH2) / RRR
               THETA3 = THETA2 + DTHETA
               YDW(N,M) = RRR * COS(THETA3)
               ZDW(N,M) = RRR * SIN(THETA3) 
 
             ENDIF
          ENDDO
       ENDDO

! ------- Stick, Panel Artist : Seungnam Kim, 24/03/2016 ------- !
c      IF (DUCTGAP .EQ. 0.0) THEN
c       DO Y = 1, NSW(MRP) - NDAFTT2
c        YW(Y+NDAFTT2,MRP) = YDW(Y+1,1)
c        ZW(Y+NDAFTT2,MRP) = ZDW(Y+1,1)
c       ENDDO
c      ENDIF
! ------- Stick, Panel Artist : Seungnam Kim, 24/03/2016 ------- !

!----Up to here, I have finished modification on the duct wake. 26/03/2016----!
!------- Repanelling Project, Art Director : Seungnam, Kim, 23/03/2016 -------!
     
      ELSE 
 
       RRR = SQRT(YD(1,1)**2 + ZD(1,1)**2)

       NDWK = 100
       NDWKP = NDWK + 1
       DWL = DWAKEL
       DTN = PI / REAL(NDWK)

       DELTAK = TWOPI / NBLADE

       IF(IDGEO .EQ. 0) THEN
          PPD = XPI(NX)
          PITCH2 = 0.5 * PI - ATAN(PPD/PI)
       ELSEIF(IDGEO .EQ. 1) THEN
          PITCH2 = DUCTPT
       ENDIF

       DO M = 1 , MDUCTP
          XDW(1,M) = XD(1,M)
          YDW(1,M) = YD(1,M)
          ZDW(1,M) = ZD(1,M)

          THETA1 = DANGLE(ZDW(1,M),YDW(1,M))

          DO N = 2 , NDWKP

             XDI = 0.5 * DWL * (1. - COS(DTN*(N-1)))
             XDW(N,M) = XDW(1,M) + XDI

             DTHETA = XDI * TAN(PITCH2) / RRR
             THETA2 = THETA1 + DTHETA

             YDW(N,M) = RRR * COS(THETA2)
             ZDW(N,M) = RRR * SIN(THETA2)
          ENDDO
       ENDDO

      ENDIF

C ----- KSN, Rotate Duct Wake Geometry ----- 
cc      DO I = 2, NDWKP
cc        DO J = 1, MDUCT+1
cc          THTRSW(I,J) = DANGLE(ZDW(I,J),YDW(I,J))
C          RDTRSW(I,J) = SQRT(ZDW(I,J)**2 + YDW(I,J)**2)

cc          YDW(I,J) = RRR*COS(THTRSW(I,J)-(DELK/2))
cc          ZDW(I,J) = RRR*SIN(THTRSW(I,J)-(DELK/2))
cc        ENDDO
cc      ENDDO
C ----- KSN, Rotate Duct Wake Geometry -----

! ------- Stick, Panel Artist : Seungnam Kim, 24/03/2016 ------- !
c      IF (DUCTGAP .EQ. 0.0) THEN
c       DO Y = 1, NSW(MRP) - NDAFTT2
c        XW(Y+NDAFTT2,MRP) = XDW(Y+1,(MDUCT/2)+1)
c        YW(Y+NDAFTT2,MRP) = YDW(Y+1,(MDUCT/2)+1)
c        ZW(Y+NDAFTT2,MRP) = ZDW(Y+1,(MDUCT/2)+1)
c       ENDDO
c      ENDIF
! ------- Stick, Panel Artist : Seungnam Kim, 24/03/2016 ------- !

C --- Duct Wake Ends with Disk
      DO M = 1 , MDUCTP
         XDWDK(M,2) = XDW(NDWKP,M)
         YDWDK(M,2) = YDW(NDWKP,M)
         ZDWDK(M,2) = ZDW(NDWKP,M)

         XDWDK(M,1) = XDWDK(M,2) 
         YDWDK(M,1) = 0.0
         ZDWDK(M,1) = 0.0
      ENDDO


      RETURN
      END

